#!/usr/bin/env python3
"""
Step 3: Setup MD Simulation
Generate AMBER input files and SLURM scripts for MD simulation.

This version replicates the successful MD protocol from systems/ligand_* cases.
Key features:
- NO distance restraint (DISANG)
- Single-stage heating (0→300K in 50ps)
- NO nvt_long stage
- eq2 stages without positional restraints
- 50 ns production

Author: zhangshd
Date: 2024-12-17
Updated: 2024-12-25 (reverted to successful protocol)
"""

import os
from pathlib import Path
from typing import Optional

from config import (
    SYSTEMS_DIR, SLURM_CONFIG, MD_PROTOCOL,
    AMBER_HOME, LIGAND_FE_COORDINATION, FE_COORDINATION_RESTRAINT,
    FE_RESIDUE, LIGAND_RESIDUE, ensure_directories
)

try:
    import parmed as pmd
except ImportError:
    pmd = None


def get_atom_number(system_dir: Path, residue_num: int, atom_name: str) -> Optional[int]:
    """
    Get absolute atom number from PDB file.
    
    Args:
        system_dir: System directory containing complex_dry.pdb
        residue_num: Residue number (1-indexed)
        atom_name: Atom name (e.g., "FE", "N22")
        
    Returns:
        Absolute atom number (1-indexed), or None if not found
    """
    pdb_file = system_dir / "complex_dry.pdb"
    
    if not pdb_file.exists():
        print(f"Warning: PDB file not found: {pdb_file}")
        return None
    
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    # PDB format: columns 7-11 for atom number, 13-16 for atom name, 23-26 for residue number
                    atom_num = int(line[6:11].strip())
                    atom_name_pdb = line[12:16].strip()
                    res_num = int(line[22:26].strip())
                    
                    if res_num == residue_num and atom_name_pdb == atom_name:
                        return atom_num
        
        print(f"Warning: Atom {atom_name} in residue {residue_num} not found in PDB")
        return None
        
    except Exception as e:
        print(f"Warning: Failed to read PDB file: {e}")
        return None


def create_distance_restraint_file(
    ligand_id: str,
    output_dir: Path,
    stage: str
) -> Optional[Path]:
    """
    Create distance restraint file for Type II inhibitors (Fe-N coordination).
    
    Args:
        ligand_id: Ligand identifier
        output_dir: Output directory for restraint file
        stage: MD stage (min1, min2, heat, eq1, eq2, prod)
        
    Returns:
        Path to restraint file, or None if not applicable
    """
    # Check if this is a Type II inhibitor
    if ligand_id not in LIGAND_FE_COORDINATION:
        return None
    
    # Get restraint force constant for this stage
    if stage not in FE_COORDINATION_RESTRAINT:
        return None
    
    force_const = FE_COORDINATION_RESTRAINT[stage]
    
    # No restraint file needed if force constant is 0
    if force_const == 0.0:
        return None
    
    # Get coordination parameters
    coord_params = LIGAND_FE_COORDINATION[ligand_id]
    coord_atom = coord_params["atom"]
    target_dist = coord_params["target_dist"]
    tolerance = coord_params.get("tolerance", 0.3)
    
    # Get absolute atom numbers from PDB file
    fe_atom_num = get_atom_number(output_dir, FE_RESIDUE, "FE")
    ligand_atom_num = get_atom_number(output_dir, LIGAND_RESIDUE, coord_atom)
    
    if fe_atom_num is None or ligand_atom_num is None:
        print(f"Warning: Could not get atom numbers for Fe-{coord_atom} restraint")
        return None
    
    # Calculate flat-bottom restraint bounds
    r2 = target_dist - tolerance
    r3 = target_dist + tolerance
    
    # Create restraint file
    restraint_file = output_dir / f"fe_coord_restraint_{stage}.RST"
    
    content = f"""# Fe-Ligand coordination restraint for Type II inhibitor {ligand_id}
# Stage: {stage}, Force constant: {force_const} kcal/mol/Å²
# Target distance: {target_dist} Å, Tolerance: ±{tolerance} Å
# Flat-bottom potential: no penalty between {r2:.2f}-{r3:.2f} Å
# Fe atom: {fe_atom_num} (residue {FE_RESIDUE})
# Ligand atom: {ligand_atom_num} ({coord_atom} in residue {LIGAND_RESIDUE})

 &rst
  iat={fe_atom_num},{ligand_atom_num},
  r1=0.0, r2={r2:.2f}, r3={r3:.2f}, r4=999.0,
  rk2={force_const:.1f}, rk3={force_const:.1f},
 &end
"""
    
    restraint_file.write_text(content)
    return restraint_file


def create_minimization_inputs(output_dir: Path, ligand_id: str) -> None:
    """Create minimization input files (2 stages) - matching successful cases."""
    
    # Check if Type II inhibitor needs Fe-coordination restraint
    restraint_file_min1 = create_distance_restraint_file(ligand_id, output_dir, "min1")
    restraint_file_min2 = create_distance_restraint_file(ligand_id, output_dir, "min2")
    
    # Min1 input
    min1 = """# Stage 1: Minimization with restraints on heavy atoms
Minimization of complex with restraints
 &cntrl
  imin=1,                   ! Minimization
  maxcyc=10000,             ! Maximum cycles
  ncyc=5000,                ! Steepest descent for first 5000 cycles
  ntb=1,                    ! Constant volume
  ntr=1,                    ! Positional restraints
  restraint_wt=10.0,        ! 10 kcal/mol-A^2 restraint weight
  restraintmask='!@H=',     ! Restrain all heavy atoms
  cut=10.0,                 ! Nonbonded cutoff
  ntpr=500,                 ! Print frequency
  ntxo=2,                   ! NetCDF restart file
"""
    if restraint_file_min1:
        min1 += f"  nmropt=1,                 ! NMR restraints (distance restraint)\n"
    min1 += " /\n"
    
    if restraint_file_min1:
        min1 += f""" &wt type='END' /
DISANG={restraint_file_min1.name}
"""
    
    # Min2 input
    min2 = """# Stage 2: Minimization without restraints
Minimization of complex without restraints
 &cntrl
  imin=1,                   ! Minimization
  maxcyc=10000,             ! Maximum cycles
  ncyc=5000,                ! Steepest descent for first 5000 cycles
  ntb=1,                    ! Constant volume
  ntr=0,                    ! No restraints
  cut=10.0,                 ! Nonbonded cutoff
  ntpr=500,                 ! Print frequency
  ntxo=2,                   ! NetCDF restart file
"""
    if restraint_file_min2:
        min2 += f"  nmropt=1,                 ! NMR restraints (distance restraint)\n"
    min2 += " /\n"
    
    if restraint_file_min2:
        min2 += f""" &wt type='END' /
DISANG={restraint_file_min2.name}
"""
    
    (output_dir / "min1.in").write_text(min1)
    (output_dir / "min2.in").write_text(min2)


def create_heating_input(output_dir: Path, ligand_id: str) -> None:
    """Create single-stage heating input (0K→300K in 50ps) - matching successful cases."""
    
    # Check if Type II inhibitor needs Fe-coordination restraint
    restraint_file = create_distance_restraint_file(ligand_id, output_dir, "heat")
    
    heat = """# Heating: NVT from 0K to 300K in 50ps
Heating of complex
 &cntrl
  imin=0,                   ! MD (not minimization)
  irest=0,                  ! New simulation
  ntx=1,                    ! Read coordinates only
  nstlim=25000,             ! 50 ps (25000 steps * 2 fs)
  dt=0.002,                 ! 2 fs timestep
  ntc=2,                    ! SHAKE on hydrogen bonds
  ntf=2,                    ! No force evaluation on constrained bonds
  ntt=3,                    ! Langevin thermostat
  gamma_ln=2.0,             ! Collision frequency
  tempi=0.0,                ! Initial temperature
  temp0=300.0,              ! Target temperature
  ntb=1,                    ! Constant volume
  ntp=0,                    ! No pressure control
  ntr=1,                    ! Positional restraints
  restraint_wt=5.0,         ! 5 kcal/mol-A^2 restraint weight
  restraintmask='@CA,C,N,O',  ! Restrain backbone
  cut=10.0,                 ! Nonbonded cutoff
  ntpr=500,                 ! Print frequency
  ntwx=500,                 ! Trajectory frequency
  ntwr=500,                 ! Restart frequency
  ntxo=2,                   ! NetCDF restart file
  ioutfm=1,                 ! NetCDF trajectory format
  nmropt=1,                 ! NMR restraints (for temperature ramping)
 /
 &wt
  type='TEMP0',
  istep1=0, istep2=25000,
  value1=0.0, value2=300.0,
 /
 &wt type='END' /
"""
    
    # Add distance restraint for Type II inhibitors
    if restraint_file:
        heat += f"DISANG={restraint_file.name}\n"
    
    (output_dir / "heat.in").write_text(heat)


def create_equilibration_inputs(output_dir: Path, ligand_id: str, num_stages: int = 10) -> None:
    """Create equilibration input files (NPT) - matching successful cases."""
    
    # Check if Type II inhibitor needs Fe-coordination restraint
    restraint_file_eq1 = create_distance_restraint_file(ligand_id, output_dir, "eq1")
    restraint_file_eq2 = create_distance_restraint_file(ligand_id, output_dir, "eq2")
    
    # Stage 1: With backbone restraints (100ps)
    eq1 = """# Equilibration Stage 1: NPT with restraints (100ps)
NPT equilibration with backbone restraints
 &cntrl
  imin=0,                   ! MD
  irest=1,                  ! Restart
  ntx=5,                    ! Read coordinates and velocities
  nstlim=50000,             ! 100 ps
  dt=0.002,                 ! 2 fs timestep
  ntc=2,                    ! SHAKE on hydrogen bonds
  ntf=2,                    ! No force evaluation on constrained bonds
  ntt=3,                    ! Langevin thermostat
  gamma_ln=2.0,             ! Collision frequency
  temp0=300.0,              ! Target temperature
  ntb=2,                    ! Constant pressure
  ntp=1,                    ! Isotropic pressure scaling
  barostat=2,               ! Monte Carlo barostat
  pres0=1.0,                ! Reference pressure (bar)
  taup=2.0,                 ! Pressure relaxation time
  ntr=1,                    ! Positional restraints
  restraint_wt=2.0,         ! 2 kcal/mol-A^2 restraint weight
  restraintmask='@CA,C,N,O',  ! Restrain backbone
  cut=10.0,                 ! Nonbonded cutoff
  ntpr=500,                 ! Print frequency
  ntwx=500,                 ! Trajectory frequency
  ntwr=500,                 ! Restart frequency
  ntxo=2,                   ! NetCDF restart file
  ioutfm=1,                 ! NetCDF trajectory format
"""
    if restraint_file_eq1:
        eq1 += "  nmropt=1,                 ! NMR restraints (distance restraint)\n"
    eq1 += " /\n"
    
    if restraint_file_eq1:
        eq1 += f" &wt type='END' /\nDISANG={restraint_file_eq1.name}\n"
    
    (output_dir / "eq1.in").write_text(eq1)
    
    # Stage 2: Without restraints (10 x 50ps sub-stages)
    # All stages have ntr=0 (no positional restraints)
    for i in range(1, num_stages + 1):
        eq2_stage = f"""# Equilibration Stage 2.{i}: NPT without restraints (50ps)
NPT equilibration without restraints - sub-stage {i}/{num_stages}
 &cntrl
  imin=0,                   ! MD
  irest=1,                  ! Restart
  ntx=5,                    ! Read coordinates and velocities
  nstlim=25000,             ! 50 ps per sub-stage
  dt=0.002,                 ! 2 fs timestep
  ntc=2,                    ! SHAKE on hydrogen bonds
  ntf=2,                    ! No force evaluation on constrained bonds
  ntt=3,                    ! Langevin thermostat
  gamma_ln=2.0,             ! Collision frequency
  temp0=300.0,              ! Target temperature
  ntb=2,                    ! Constant pressure
  ntp=1,                    ! Isotropic pressure scaling
  barostat=2,               ! Monte Carlo barostat
  pres0=1.0,                ! Reference pressure (bar)
  taup=2.0,                 ! Pressure relaxation time
  ntr=0,                    ! No restraints
  cut=10.0,                 ! Nonbonded cutoff
  ntpr=500,                 ! Print frequency
  ntwx=500,                 ! Trajectory frequency
  ntwr=500,                 ! Restart frequency
  ntxo=2,                   ! NetCDF restart file
  ioutfm=1,                 ! NetCDF trajectory format
"""
        if restraint_file_eq2:
            eq2_stage += "  nmropt=1,                 ! NMR restraints (distance restraint)\n"
        eq2_stage += " /\n"
        
        if restraint_file_eq2:
            eq2_stage += f" &wt type='END' /\nDISANG={restraint_file_eq2.name}\n"
        
        (output_dir / f"eq2_{i:02d}.in").write_text(eq2_stage)


def create_production_input(output_dir: Path) -> None:
    """Create production MD input file (50 ns) - matching successful cases."""
    
    prod = """# Production MD: NPT (50 ns)
Production NPT simulation
 &cntrl
  imin=0,                   ! MD
  irest=1,                  ! Restart
  ntx=5,                    ! Read coordinates and velocities
  nstlim=25000000,          ! 50 ns
  dt=0.002,                 ! 2 fs timestep
  ntc=2,                    ! SHAKE on hydrogen bonds
  ntf=2,                    ! No force evaluation on constrained bonds
  ntt=3,                    ! Langevin thermostat
  gamma_ln=2.0,             ! Collision frequency
  temp0=300.0,              ! Target temperature
  ntb=2,                    ! Constant pressure
  ntp=1,                    ! Isotropic pressure scaling
  barostat=2,               ! Monte Carlo barostat
  pres0=1.0,                ! Reference pressure (bar)
  taup=2.0,                 ! Pressure relaxation time
  cut=10.0,                 ! Nonbonded cutoff
  ntpr=5000,                ! Print every 10ps
  ntwx=5000,                ! Save trajectory every 10ps
  ntwr=5000,                ! Restart every 10ps
  ntxo=2,                   ! NetCDF restart file
  ioutfm=1,                 ! NetCDF trajectory format
 /
"""
    
    (output_dir / "prod.in").write_text(prod)


def create_slurm_script(output_dir: Path, ligand_id: str) -> Path:
    """Create SLURM submission script for MD simulation - matching successful cases."""
    
    num_eq2 = MD_PROTOCOL["eq2"]["num_stages"]
    
    script = f"""#!/bin/bash
#SBATCH --job-name={ligand_id}_md
#SBATCH --partition={SLURM_CONFIG["partition"]}
#SBATCH --nodes={SLURM_CONFIG["nodes"]}
#SBATCH --ntasks={SLURM_CONFIG["ntasks"]}
#SBATCH --cpus-per-task={SLURM_CONFIG["cpus_per_task"]}
#SBATCH --gres={SLURM_CONFIG["gres"]}
#SBATCH --output=md_%j.log

# Exit on error (will be handled by check_error function for specific commands)
# set -e  # Disabled to allow custom error handling

# Load CUDA and AMBER environment
module load cuda/11.8 fftw/fftw-3.3.8
source {AMBER_HOME}/amber.sh

# Define variables
PRMTOP="complex_solv.prmtop"
INPCRD="complex_solv.inpcrd"
EXE="pmemd.cuda"

# Error checking function
check_error() {{
    local stage=$1
    if [ $? -ne 0 ]; then
        echo "ERROR: Stage $stage failed!"
        exit 1
    fi
    
    # Check for common AMBER errors in output
    if [ -f "${{stage}}.out" ]; then
        if grep -q "ERROR" "${{stage}}.out" 2>/dev/null; then
            echo "ERROR: AMBER reported error in ${{stage}}.out"
            grep "ERROR" "${{stage}}.out"
            exit 1
        fi
        if grep -q "Calculation halted" "${{stage}}.out" 2>/dev/null; then
            echo "ERROR: Calculation halted in ${{stage}}.out"
            tail -20 "${{stage}}.out"
            exit 1
        fi
    fi
    
    echo "Stage $stage completed successfully."
}}

echo "=========================================="
echo "Starting MD simulation for {ligand_id}"
echo "Date: $(date)"
echo "=========================================="

# Stage 1: Minimization with restraints
echo "Starting minimization stage 1..."
$EXE -O -i min1.in -o min1.out -p $PRMTOP -c $INPCRD -r min1.rst7 -ref $INPCRD
check_error "min1"

# Stage 2: Minimization without restraints
echo "Starting minimization stage 2..."
$EXE -O -i min2.in -o min2.out -p $PRMTOP -c min1.rst7 -r min2.rst7
check_error "min2"

# Stage 3: Heating (NVT)
echo "Starting heating..."
$EXE -O -i heat.in -o heat.out -p $PRMTOP -c min2.rst7 -r heat.rst7 -x heat.nc -ref min2.rst7
check_error "heat"

# Stage 4: Equilibration with restraints (NPT)
echo "Starting equilibration stage 1..."
$EXE -O -i eq1.in -o eq1.out -p $PRMTOP -c heat.rst7 -r eq1.rst7 -x eq1.nc -ref heat.rst7
check_error "eq1"
"""
    
    # Add eq2 stages (all without restraints)
    for i in range(1, num_eq2 + 1):
        if i == 1:
            prev_rst = "eq1.rst7"
        else:
            prev_rst = f"eq2_{i-1:02d}.rst7"
        
        script += f"""
# Stage 5.{i}: Equilibration without restraints - sub-stage {i}/{num_eq2}
echo "Starting equilibration stage 2.{i}..."
$EXE -O -i eq2_{i:02d}.in -o eq2_{i:02d}.out -p $PRMTOP -c {prev_rst} -r eq2_{i:02d}.rst7 -x eq2_{i:02d}.nc
check_error "eq2_{i:02d}"
"""
    
    script += f"""
# Stage 6: Production MD (NPT)
echo "Starting production MD..."
$EXE -O -i prod.in -o prod.out -p $PRMTOP -c eq2_{num_eq2:02d}.rst7 -r prod.rst7 -x prod.nc
check_error "prod"

echo "=========================================="
echo "MD simulation completed successfully!"
echo "Date: $(date)"
echo "=========================================="
"""
    
    script_path = output_dir / "run_md.sh"
    script_path.write_text(script)
    script_path.chmod(0o755)
    
    return script_path


def setup_md_simulation(ligand_id: str, force: bool = False) -> bool:
    """
    Complete MD simulation setup.
    
    Args:
        ligand_id: Unique identifier for the ligand
        force: If True, overwrite existing files
        
    Returns:
        True if successful
    """
    system_dir = SYSTEMS_DIR / ligand_id
    
    # Check topology exists
    prmtop = system_dir / "complex_solv.prmtop"
    if not prmtop.exists():
        print(f"  ERROR: Topology not found for {ligand_id}")
        print(f"  Run step2_build_complex_topology.py first")
        return False
    
    # Check if already done
    run_script = system_dir / "run_md.sh"
    if not force and run_script.exists():
        print(f"  SKIP: MD files already exist for {ligand_id}")
        return True
    
    # Create input files (with Type II inhibitor distance restraint if applicable)
    print(f"  Creating minimization inputs...")
    create_minimization_inputs(system_dir, ligand_id)
    
    print(f"  Creating heating input...")
    create_heating_input(system_dir, ligand_id)
    
    print(f"  Creating equilibration inputs...")
    create_equilibration_inputs(system_dir, ligand_id)
    
    print(f"  Creating production input...")
    create_production_input(system_dir)
    
    print(f"  Creating SLURM script...")
    create_slurm_script(system_dir, ligand_id)
    
    print(f"  SUCCESS: MD simulation files created")
    return True


def main():
    """Standalone test."""
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python step3_setup_md_simulation.py <ligand_id>")
        sys.exit(1)
    
    ligand_id = sys.argv[1]
    force = "--force" in sys.argv or "-f" in sys.argv
    
    success = setup_md_simulation(ligand_id, force=force)
    
    if success:
        print(f"\nTo submit MD job:")
        print(f"  cd {SYSTEMS_DIR / ligand_id}")
        print(f"  sbatch run_md.sh")
    else:
        print("\nMD setup failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()
