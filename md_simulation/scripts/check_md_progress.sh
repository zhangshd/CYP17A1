#!/bin/bash
# ============================================================================
# Integrated MD Simulation Progress Monitor
# 检查MD模拟进度和配体稳定性
#
# Usage:
#   ./check_md_progress.sh <ligand_id>           # Single check
#   ./check_md_progress.sh <ligand_id> --watch   # Real-time monitoring (30s)
#   ./check_md_progress.sh <ligand_id> -w        # Same as --watch
#
# Author: zhangshd
# Date: 2024-12-25
# ============================================================================

# ============================================================================
# Configuration
# ============================================================================

LIGAND_ID=${1:-AER601}
WATCH_MODE=0

# Parse arguments
for arg in "$@"; do
    case $arg in
        --watch|-w)
            WATCH_MODE=1
            shift
            ;;
    esac
done

SYSTEM_DIR="/home/zhangsd/repos/CYP17A1/md_simulation/complexes/systems_gpu/${LIGAND_ID}"
AMBER_BIN="/opt/share/Amber/amber22/bin"

# Check if system directory exists
if [ ! -d "$SYSTEM_DIR" ]; then
    echo "ERROR: System directory not found: $SYSTEM_DIR"
    exit 1
fi

cd "$SYSTEM_DIR" || exit 1

# ============================================================================
# Helper Functions
# ============================================================================

print_header() {
    echo "════════════════════════════════════════════════════════════════════"
    echo "  MD Simulation Progress Monitor"
    echo "  System: ${LIGAND_ID}"
    echo "  Time: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "════════════════════════════════════════════════════════════════════"
}

check_job_status() {
    echo -e "\n【1. Job Status】"
    JOB_INFO=$(squeue -u $USER -n ${LIGAND_ID}_md -h -o "%i %T %M %l" 2>/dev/null)
    if [ -n "$JOB_INFO" ]; then
        JOB_ID=$(echo "$JOB_INFO" | awk '{print $1}')
        JOB_STATE=$(echo "$JOB_INFO" | awk '{print $2}')
        JOB_TIME=$(echo "$JOB_INFO" | awk '{print $3}')
        JOB_LIMIT=$(echo "$JOB_INFO" | awk '{print $4}')
        echo "  Job ID: $JOB_ID"
        echo "  State: $JOB_STATE"
        echo "  Runtime: $JOB_TIME / $JOB_LIMIT"
        return 0
    else
        echo "  Status: Not running or completed"
        return 1
    fi
}

check_md_stage() {
    echo -e "\n【2. MD Stage】"
    
    LOG_FILE=$(ls -t md_*.log 2>/dev/null | head -1)
    if [ -z "$LOG_FILE" ]; then
        echo "  ⚠️  No log file found"
        return 1
    fi
    
    echo "  Log: $LOG_FILE"
    
    # Find last completed stage
    LAST_COMPLETED=$(grep "Stage.*completed successfully" "$LOG_FILE" | tail -1 | sed 's/.*Stage //' | sed 's/ completed.*//')
    if [ -n "$LAST_COMPLETED" ]; then
        echo "  Last completed: $LAST_COMPLETED"
    fi
    
    # Find currently running stage
    CURRENT_STAGE=$(grep "Starting" "$LOG_FILE" | tail -1 | sed 's/Starting //' | sed 's/\.\.\.//')
    if [ -n "$CURRENT_STAGE" ]; then
        echo "  Currently running: $CURRENT_STAGE"
    fi
    
    # Check for errors
    if grep -qi "error\|fatal\|vlimit exceeded\|illegal memory" "$LOG_FILE" 2>/dev/null; then
        echo "  ❌ Errors detected:"
        grep -iE "error:|fatal|vlimit exceeded" "$LOG_FILE" | tail -3 | sed 's/^/    /'
        return 1
    else
        echo "  ✅ No errors"
    fi
    
    return 0
}

find_latest_trajectory() {
    # Find the most recent trajectory file based on log
    LOG_FILE=$(ls -t md_*.log 2>/dev/null | head -1)
    
    # Check in reverse order of MD stages
    if [ -f "prod.nc" ] && grep -q "Starting production" "$LOG_FILE" 2>/dev/null; then
        echo "prod.nc|prod|Production"
        return 0
    fi
    
    # Check eq2 stages (10-1)
    for i in {10..1}; do
        STAGE_FILE="eq2_$(printf "%02d" $i).nc"
        if [ -f "$STAGE_FILE" ]; then
            echo "$STAGE_FILE|eq2_$(printf "%02d" $i)|Equilibration 2.$i"
            return 0
        fi
    done
    
    if [ -f "eq1.nc" ]; then
        echo "eq1.nc|eq1|Equilibration 1"
        return 0
    fi
    
    # Check for heating (old systems had heat1/2/3, new systems have heat)
    if [ -f "heat.nc" ]; then
        echo "heat.nc|heat|Heating"
        return 0
    fi
    
    for i in {3..1}; do
        HEAT_FILE="heat${i}.nc"
        if [ -f "$HEAT_FILE" ]; then
            echo "$HEAT_FILE|heat${i}|Heating Stage $i"
            return 0
        fi
    done
    
    echo "||No trajectory found"
    return 1
}

analyze_fe_ligand_distance() {
    local TRAJ_FILE=$1
    local STAGE_NAME=$2
    
    echo -e "\n【3. Fe-Ligand Distance Analysis】"
    echo "  Trajectory: $TRAJ_FILE ($STAGE_NAME)"
    
    # Create cpptraj input for Fe-ligand center-of-mass distance
    cat > tmp_check_dist.in << EOF
parm complex_solv.prmtop
trajin $TRAJ_FILE
distance fe_lig_com :475@FE :476 com out tmp_fe_lig_dist.dat
go
quit
EOF
    
    $AMBER_BIN/cpptraj -i tmp_check_dist.in > /dev/null 2>&1
    
    if [ -f "tmp_fe_lig_dist.dat" ]; then
        # Calculate statistics
        STATS=$(awk 'NR>1 {
            sum+=$2; n++; 
            if(min=="" || $2<min) min=$2; 
            if($2>max) max=$2;
            if($2<=6.0) stable++;
            if($2>10.0) broken++;
        } END {
            printf "%.2f %.2f %.2f %d %d %d", sum/n, min, max, n, stable, broken
        }' tmp_fe_lig_dist.dat)
        
        AVG=$(echo $STATS | awk '{print $1}')
        MIN=$(echo $STATS | awk '{print $2}')
        MAX=$(echo $STATS | awk '{print $3}')
        FRAMES=$(echo $STATS | awk '{print $4}')
        STABLE=$(echo $STATS | awk '{print $5}')
        BROKEN=$(echo $STATS | awk '{print $6}')
        
        echo "  Total frames: $FRAMES"
        echo "  Average distance: $AVG Å"
        echo "  Range: $MIN - $MAX Å"
        
        if [ $FRAMES -gt 0 ]; then
            STABLE_PCT=$(echo "scale=1; $STABLE*100/$FRAMES" | bc)
            BROKEN_PCT=$(echo "scale=1; $BROKEN*100/$FRAMES" | bc)
            echo "  Stable frames (≤6Å): $STABLE (${STABLE_PCT}%)"
            echo "  Dissociated frames (>10Å): $BROKEN (${BROKEN_PCT}%)"
        fi
        
        # Stability assessment
        if (( $(echo "$AVG <= 6.0" | bc -l) )); then
            echo "  ✅ Ligand is STABLE"
            STABILITY="STABLE"
        elif (( $(echo "$AVG <= 10.0" | bc -l) )); then
            echo "  ⚠️  Ligand is DRIFTING"
            STABILITY="DRIFTING"
        else
            echo "  ❌ Ligand DISSOCIATED"
            STABILITY="DISSOCIATED"
        fi
        
        # Show first and last 3 frames
        echo -e "\n  First 3 frames:"
        head -4 tmp_fe_lig_dist.dat | tail -3 | awk '{printf "    Frame %4d: %.2f Å\n", $1, $2}'
        echo "  Last 3 frames:"
        tail -3 tmp_fe_lig_dist.dat | awk '{printf "    Frame %4d: %.2f Å\n", $1, $2}'
        
        rm -f tmp_check_dist.in tmp_fe_lig_dist.dat
    else
        echo "  ❌ Failed to analyze distance"
        STABILITY="UNKNOWN"
    fi
}

analyze_fe_n_bond() {
    local TRAJ_FILE=$1
    local STAGE_NAME=$2
    
    # Check if N22 atom exists (Type II inhibitor)
    if ! grep -q " N22 " complex_dry.pdb 2>/dev/null; then
        return 0
    fi
    
    echo -e "\n【4. Fe-N22 Coordination Bond】"
    echo "  (Type II inhibitor detected)"
    
    # Create cpptraj input for Fe-N22 distance
    cat > tmp_check_fen.in << EOF
parm complex_solv.prmtop
trajin $TRAJ_FILE
distance fe_n22 :475@FE :476@N22 out tmp_fe_n22.dat
go
quit
EOF
    
    $AMBER_BIN/cpptraj -i tmp_check_fen.in > /dev/null 2>&1
    
    if [ -f "tmp_fe_n22.dat" ]; then
        # Calculate statistics
        STATS=$(awk 'NR>1 {
            sum+=$2; n++; 
            if(min=="" || $2<min) min=$2; 
            if($2>max) max=$2;
            if($2<=3.0) stable++;
            if($2>5.0) broken++;
        } END {
            printf "%.3f %.3f %.3f %d %d %d", sum/n, min, max, n, stable, broken
        }' tmp_fe_n22.dat)
        
        AVG=$(echo $STATS | awk '{print $1}')
        MIN=$(echo $STATS | awk '{print $2}')
        MAX=$(echo $STATS | awk '{print $3}')
        FRAMES=$(echo $STATS | awk '{print $4}')
        STABLE=$(echo $STATS | awk '{print $5}')
        BROKEN=$(echo $STATS | awk '{print $6}')
        
        echo "  Total frames: $FRAMES"
        echo "  Average bond length: $AVG Å (ideal: ~2.1 Å)"
        echo "  Range: $MIN - $MAX Å"
        
        if [ $FRAMES -gt 0 ]; then
            STABLE_PCT=$(echo "scale=1; $STABLE*100/$FRAMES" | bc)
            BROKEN_PCT=$(echo "scale=1; $BROKEN*100/$FRAMES" | bc)
            echo "  Coordinated frames (≤3Å): $STABLE (${STABLE_PCT}%)"
            echo "  Broken frames (>5Å): $BROKEN (${BROKEN_PCT}%)"
        fi
        
        # Bond stability assessment
        if (( $(echo "$AVG <= 2.5" | bc -l) )); then
            echo "  ✅ Coordination bond is EXCELLENT"
        elif (( $(echo "$AVG <= 3.5" | bc -l) )); then
            echo "  ✅ Coordination bond is STABLE"
        elif (( $(echo "$AVG <= 5.0" | bc -l) )); then
            echo "  ⚠️  Coordination bond is STRETCHED"
        else
            echo "  ❌ Coordination bond is BROKEN"
        fi
        
        # Show first and last 3 frames
        echo -e "\n  First 3 frames:"
        head -4 tmp_fe_n22.dat | tail -3 | awk '{printf "    Frame %4d: %.3f Å\n", $1, $2}'
        echo "  Last 3 frames:"
        tail -3 tmp_fe_n22.dat | awk '{printf "    Frame %4d: %.3f Å\n", $1, $2}'
        
        rm -f tmp_check_fen.in tmp_fe_n22.dat
    else
        echo "  ⚠️  Failed to analyze Fe-N22 bond"
    fi
}

check_file_sizes() {
    echo -e "\n【5. File Sizes】"
    
    if [ -f "prod.nc" ]; then
        SIZE=$(du -h prod.nc | awk '{print $1}')
        echo "  prod.nc: $SIZE"
    fi
    
    TOTAL_SIZE=$(du -sh . 2>/dev/null | awk '{print $1}')
    echo "  Total system size: $TOTAL_SIZE"
}

print_summary() {
    echo -e "\n════════════════════════════════════════════════════════════════════"
    echo "  SUMMARY"
    echo "  ──────────────────────────────────────────────────────────────────"
    
    if [ "$STABILITY" == "STABLE" ]; then
        echo "  ✅ System is stable and running normally"
        echo "  📊 Ligand remains in binding pocket"
    elif [ "$STABILITY" == "DRIFTING" ]; then
        echo "  ⚠️  Ligand shows some drift but not dissociated"
        echo "  📊 Continue monitoring"
    elif [ "$STABILITY" == "DISSOCIATED" ]; then
        echo "  ❌ WARNING: Ligand has dissociated from binding site"
        echo "  📊 Simulation may need to be restarted with different parameters"
    else
        echo "  ℹ️  Status unknown - check trajectory files"
    fi
    
    echo "════════════════════════════════════════════════════════════════════"
}

# ============================================================================
# Main Check Function
# ============================================================================

run_check() {
    if [ $WATCH_MODE -eq 1 ]; then
        clear
    fi
    
    print_header
    
    # 1. Check job status
    check_job_status
    JOB_RUNNING=$?
    
    # 2. Check MD stage
    check_md_stage
    
    # 3. Find latest trajectory
    TRAJ_INFO=$(find_latest_trajectory)
    if [ $? -ne 0 ]; then
        echo -e "\n⚠️  No trajectory files found yet"
        print_summary
        return 1
    fi
    
    TRAJ_FILE=$(echo "$TRAJ_INFO" | cut -d'|' -f1)
    STAGE_ID=$(echo "$TRAJ_INFO" | cut -d'|' -f2)
    STAGE_DESC=$(echo "$TRAJ_INFO" | cut -d'|' -f3)
    
    if [ -z "$TRAJ_FILE" ]; then
        echo -e "\n⚠️  No trajectory files available for analysis"
        print_summary
        return 1
    fi
    
    # 4. Analyze Fe-ligand distance
    analyze_fe_ligand_distance "$TRAJ_FILE" "$STAGE_DESC"
    
    # 5. Analyze Fe-N bond if applicable
    analyze_fe_n_bond "$TRAJ_FILE" "$STAGE_DESC"
    
    # 6. Check file sizes
    check_file_sizes
    
    # 7. Print summary
    print_summary
    
    return 0
}

# ============================================================================
# Main Execution
# ============================================================================

if [ $WATCH_MODE -eq 1 ]; then
    echo "Starting real-time monitoring mode (Ctrl+C to exit)..."
    sleep 2
    
    while true; do
        run_check
        
        echo -e "\n💤 Next update in 30 seconds..."
        sleep 30
    done
else
    run_check
fi

