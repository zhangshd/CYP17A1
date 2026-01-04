#!/bin/bash
# ============================================================================
# Batch Check MD Simulation Status
# 批量检查多个系统的MD模拟状态
#
# Usage:
#   ./batch_check_status.sh                    # Check all systems
#   ./batch_check_status.sh GRAS_*             # Check specific pattern
#   ./batch_check_status.sh --problems-only    # Only show systems with issues
#
# Author: zhangshd
# Date: 2024-12-25
# ============================================================================

SYSTEMS_DIR="/home/zhangsd/repos/CYP17A1/md_simulation/complexes/systems_round2"
CHECK_SCRIPT="./check_md_progress.sh"
PROBLEMS_ONLY=0

# Parse arguments
if [ "$1" == "--problems-only" ] || [ "$1" == "-p" ]; then
    PROBLEMS_ONLY=1
    shift
fi

# Determine systems to check
if [ $# -gt 0 ]; then
    # User specified pattern
    PATTERN="$1"
else
    # Check all systems
    PATTERN="*"
fi

cd $(dirname "$0")

echo "════════════════════════════════════════════════════════════════════"
echo "  Batch MD Status Check"
echo "  Pattern: $PATTERN"
echo "  Time: $(date '+%Y-%m-%d %H:%M:%S')"
if [ $PROBLEMS_ONLY -eq 1 ]; then
    echo "  Mode: Problems only"
fi
echo "════════════════════════════════════════════════════════════════════"
echo ""

# Counters
TOTAL=0
RUNNING=0
STABLE=0
PROBLEMS=0
UNKNOWN=0

for system_path in ${SYSTEMS_DIR}/${PATTERN}; do
    if [ ! -d "$system_path" ]; then
        continue
    fi
    
    system=$(basename "$system_path")
    
    # Skip backup/test directories
    if [[ "$system" == *"_old"* ]] || [[ "$system" == *"_backup"* ]] || [[ "$system" == *"_test"* ]]; then
        continue
    fi
    
    TOTAL=$((TOTAL + 1))
    
    # Run check
    OUTPUT=$(${CHECK_SCRIPT} "$system" 2>&1)
    
    # Extract key information
    JOB_STATUS=$(echo "$OUTPUT" | grep -A1 "【1. Job Status】" | grep "State:" | awk '{print $2}')
    STAGE=$(echo "$OUTPUT" | grep "Currently running:" | sed 's/.*Currently running: //')
    ERRORS=$(echo "$OUTPUT" | grep -E "❌|ERROR" | head -1)
    STABILITY=$(echo "$OUTPUT" | grep -E "✅ Ligand is STABLE|⚠️  Ligand is DRIFTING|❌ Ligand DISSOCIATED" | head -1)
    
    # Determine status
    HAS_PROBLEM=0
    
    if [ -n "$ERRORS" ]; then
        STATUS="❌ ERROR"
        HAS_PROBLEM=1
        PROBLEMS=$((PROBLEMS + 1))
    elif echo "$STABILITY" | grep -q "DISSOCIATED"; then
        STATUS="❌ DISSOCIATED"
        HAS_PROBLEM=1
        PROBLEMS=$((PROBLEMS + 1))
    elif echo "$STABILITY" | grep -q "DRIFTING"; then
        STATUS="⚠️  DRIFTING"
        HAS_PROBLEM=1
        PROBLEMS=$((PROBLEMS + 1))
    elif echo "$STABILITY" | grep -q "STABLE"; then
        STATUS="✅ STABLE"
        STABLE=$((STABLE + 1))
    else
        STATUS="❓ UNKNOWN"
        UNKNOWN=$((UNKNOWN + 1))
    fi
    
    if [ -n "$JOB_STATUS" ]; then
        RUNNING=$((RUNNING + 1))
    fi
    
    # Print result
    if [ $PROBLEMS_ONLY -eq 0 ] || [ $HAS_PROBLEM -eq 1 ]; then
        printf "%-20s %s" "$system" "$STATUS"
        
        if [ -n "$JOB_STATUS" ]; then
            printf "  [%s]" "$JOB_STATUS"
        fi
        
        if [ -n "$STAGE" ]; then
            printf "  %s" "$STAGE"
        fi
        
        echo ""
        
        # Show error details if present
        if [ $HAS_PROBLEM -eq 1 ] && [ -n "$ERRORS" ]; then
            echo "$ERRORS" | head -2 | sed 's/^/    /'
        fi
    fi
done

echo ""
echo "════════════════════════════════════════════════════════════════════"
echo "  SUMMARY"
echo "  ──────────────────────────────────────────────────────────────────"
echo "  Total systems checked: $TOTAL"
echo "  Currently running: $RUNNING"
echo "  Stable: $STABLE"
echo "  With problems: $PROBLEMS"
echo "  Unknown status: $UNKNOWN"
echo "════════════════════════════════════════════════════════════════════"

