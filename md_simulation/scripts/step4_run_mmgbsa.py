"""
Step 4: MM/GBSA Binding Free Energy Analysis
Calculate binding free energy using MMPBSA.py from AmberTools.
Author: zhangshd
Date: 2024-12-17
"""

import os
import subprocess
from pathlib import Path
from typing import Optional, Tuple

from config import (
    AMBER_HOME, AMBER_SOURCE_CMD, SYSTEMS_DIR, SLURM_CONFIG,
    RECEPTOR_MASK, LIGAND_MASK, STRIP_MASK, MMGBSA_PARAMS,
    ensure_directories
)


def create_strip_trajectory_script(output_dir: Path) -> Path:
    """Create cpptraj script to strip water/ions from trajectory."""
    
    script_content = f"""# cpptraj script: Strip water and ions from trajectory
# Input: prod.nc (solvated trajectory)
# Output: prod_dry.nc (stripped trajectory for MM/GBSA)

parm complex_solv.prmtop
trajin prod.nc

# Strip water and ions
strip {STRIP_MASK}

# Write stripped trajectory
trajout prod_dry.nc netcdf

go
quit
"""
    
    script_file = output_dir / "strip_traj.cpptraj"
    script_file.write_text(script_content)
    return script_file


def create_topology_prep_script(output_dir: Path) -> Path:
    """Create cpptraj script to generate receptor and ligand topologies.
    
    Note: Box removal is done via parmed in the SLURM script (more reliable for GB).
    """
    
    script_content = f"""# cpptraj script: Prepare topologies for MM/GBSA
# Creates: complex_dry.prmtop, receptor.prmtop, ligand.prmtop

# === Complex (stripped, no water/ions) ===
parm complex_solv.prmtop [solv]
parmstrip {STRIP_MASK} parmindex 0
parmwrite out complex_dry.prmtop parmindex 0

# === Receptor only (protein + heme + Fe) ===
parm complex_solv.prmtop [rec]
parmstrip {STRIP_MASK} parmindex 1
parmstrip {LIGAND_MASK} parmindex 1
parmwrite out receptor.prmtop parmindex 1

# === Ligand only ===
parm complex_solv.prmtop [lig]
parmstrip {STRIP_MASK} parmindex 2
parmstrip {RECEPTOR_MASK} parmindex 2
parmwrite out ligand.prmtop parmindex 2

go
quit
"""
    
    script_file = output_dir / "prep_topologies.cpptraj"
    script_file.write_text(script_content)
    return script_file


def create_mmgbsa_input(output_dir: Path) -> Path:
    """Create MM/GBSA input file for MMPBSA.py."""
    
    params = MMGBSA_PARAMS
    
    input_content = f"""# MM/GBSA input file for MMPBSA.py
&general
  startframe={params['startframe']},
  endframe={params['endframe']},
  interval={params['interval']},
  verbose=2,
  keep_files=0,
/
&gb
  igb={params['igb']},
  saltcon={params['saltcon']},
/
"""
    
    input_file = output_dir / "mmgbsa.in"
    input_file.write_text(input_content)
    return input_file


def create_mmgbsa_slurm_script(output_dir: Path, ligand_id: str) -> Path:
    """Create SLURM script for MM/GBSA analysis - matching successful cases."""
    
    abs_output_dir = output_dir.resolve()
    
    script = f"""#!/bin/bash
#SBATCH --job-name={ligand_id}_mmgbsa
#SBATCH --partition=C9654
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --output=mmgbsa_%j.log

# ============================================================================
# MM/GBSA Binding Free Energy Analysis for {ligand_id}
# ============================================================================

echo "================================================================"
echo "MM/GBSA Analysis for {ligand_id}"
echo "Started at: $(date)"
echo "Running on: $(hostname)"
echo "================================================================"

# Load environment
source {AMBER_HOME}/amber.sh
module load fftw/fftw-3.3.8

cd {abs_output_dir}

# Check prerequisites
if [ ! -f "complex_solv.prmtop" ]; then
    echo "ERROR: complex_solv.prmtop not found!"
    exit 1
fi

if [ ! -f "prod.nc" ]; then
    echo "ERROR: prod.nc trajectory not found!"
    exit 1
fi

# Step 1: Prepare topologies
echo ""
echo "Step 1: Preparing topologies..."
echo "================================"

cpptraj -i prep_topologies.cpptraj

# Verify topologies were created
for f in complex_dry.prmtop receptor.prmtop ligand.prmtop; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Failed to create $f"
        exit 1
    fi
    echo "  Created: $f"
done

# Step 1b: Remove periodic box from topologies (required for GB)
echo ""
echo "Step 1b: Removing periodic box info..."
echo "======================================="

/opt/share/miniconda3/bin/python3 << 'PYEOF'
import parmed as pmd
for prmtop in ["complex_dry.prmtop", "receptor.prmtop", "ligand.prmtop"]:
    parm = pmd.load_file(prmtop)
    parm.box = None
    parm.save(prmtop, overwrite=True)
    print("  Fixed:", prmtop)
PYEOF

# Step 2: Strip trajectory
echo ""
echo "Step 2: Stripping trajectory..."
echo "================================"

cpptraj -i strip_traj.cpptraj

if [ ! -f "prod_dry.nc" ]; then
    echo "ERROR: Failed to create stripped trajectory"
    exit 1
fi
echo "  Created: prod_dry.nc"

# Step 3: Run MMPBSA.py
echo ""
echo "Step 3: Running MM/GBSA calculation..."
echo "======================================"

# MMPBSA.py serial version (more reliable than MPI)
MMPBSA.py -O \\
    -i mmgbsa.in \\
    -o FINAL_RESULTS_MMGBSA.dat \\
    -eo ENERGY_DECOMP.csv \\
    -cp complex_dry.prmtop \\
    -rp receptor.prmtop \\
    -lp ligand.prmtop \\
    -y prod_dry.nc

# Check if MMPBSA completed
if [ ! -f "FINAL_RESULTS_MMGBSA.dat" ]; then
    echo "ERROR: MMPBSA.py failed to produce output!"
    exit 1
fi

# Step 4: Extract key results
echo ""
echo "Step 4: Results Summary"
echo "======================="

echo ""
echo "=== Binding Free Energy (kcal/mol) ==="
grep -A 15 "DELTA TOTAL" FINAL_RESULTS_MMGBSA.dat | head -20

echo ""
echo "Completed at: $(date)"
echo "Full results: {abs_output_dir}/FINAL_RESULTS_MMGBSA.dat"

# Clean up large intermediate files (optional)
# rm -f prod_dry.nc _MMPBSA_*

echo ""
echo "MM/GBSA analysis completed successfully!"
"""
    
    script_path = output_dir / "run_mmgbsa.sh"
    script_path.write_text(script)
    script_path.chmod(0o755)
    
    return script_path


def setup_mmgbsa_analysis(ligand_id: str, force: bool = False) -> bool:
    """
    Setup MM/GBSA analysis files.
    
    Args:
        ligand_id: Unique identifier for the ligand
        force: If True, overwrite existing files
        
    Returns:
        True if successful
    """
    system_dir = SYSTEMS_DIR / ligand_id
    
    # Check production trajectory exists
    prod_nc = system_dir / "prod.nc"
    if not prod_nc.exists():
        print(f"  ERROR: Production trajectory not found for {ligand_id}")
        print(f"  Run MD simulation first")
        return False
    
    # Check if already done
    mmgbsa_script = system_dir / "run_mmgbsa.sh"
    if not force and mmgbsa_script.exists():
        print(f"  SKIP: MM/GBSA files already exist for {ligand_id}")
        return True
    
    # Create input files
    print(f"  Creating cpptraj scripts...")
    create_strip_trajectory_script(system_dir)
    create_topology_prep_script(system_dir)
    
    print(f"  Creating MMPBSA input...")
    create_mmgbsa_input(system_dir)
    
    print(f"  Creating SLURM script...")
    create_mmgbsa_slurm_script(system_dir, ligand_id)
    
    print(f"  SUCCESS: MM/GBSA analysis files created")
    return True


def parse_mmgbsa_result(result_file: Path) -> Optional[Tuple[float, float, float]]:
    """
    Parse FINAL_RESULTS_MMGBSA.dat to extract binding energy.
    
    Returns:
        Tuple of (delta_G, std_dev, std_err) or None if not found
    """
    if not result_file.exists():
        return None
    
    content = result_file.read_text()
    
    for line in content.split('\n'):
        if 'DELTA TOTAL' in line:
            parts = line.split()
            if len(parts) >= 4:
                delta_g = float(parts[2])
                std_dev = float(parts[3])
                std_err = float(parts[4]) if len(parts) > 4 else 0.0
                return delta_g, std_dev, std_err
    
    return None


def main():
    """Standalone test."""
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python step4_run_mmgbsa.py <ligand_id> [--force]")
        sys.exit(1)
    
    ligand_id = sys.argv[1]
    force = "--force" in sys.argv or "-f" in sys.argv
    
    success = setup_mmgbsa_analysis(ligand_id, force=force)
    
    if success:
        print(f"\nTo submit MM/GBSA job:")
        print(f"  cd {SYSTEMS_DIR / ligand_id}")
        print(f"  sbatch run_mmgbsa.sh")
    else:
        print("\nMM/GBSA setup failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()
