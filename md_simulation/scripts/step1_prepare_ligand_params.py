"""
Step 1: Ligand Parameterization
Generate GAFF2 force field parameters and AM1-BCC charges for ligands.
Author: zhangshd
Date: 2024-12-17
"""

import os
import subprocess
from pathlib import Path
from typing import Optional, Tuple

from config import (
    AMBER_HOME, LIGAND_PARAMS_DIR, LIGAND_FF, CHARGE_MODEL,
    LIGAND_RESIDUE_NAME, ensure_directories
)


def validate_mol2(mol2_path: Path) -> bool:
    """
    Validate that a mol2 file is complete and has hydrogen atoms.
    
    Args:
        mol2_path: Path to mol2 file
        
    Returns:
        True if valid, False otherwise
    """
    if not mol2_path.exists():
        print(f"  ERROR: File not found: {mol2_path}")
        return False
    
    content = mol2_path.read_text()
    
    # Check for ATOM section
    if "@<TRIPOS>ATOM" not in content:
        print(f"  ERROR: No ATOM section in {mol2_path}")
        return False
    
    # Check for hydrogen atoms
    lines = content.split("\n")
    in_atom_section = False
    has_hydrogen = False
    atom_count = 0
    
    for line in lines:
        if "@<TRIPOS>ATOM" in line:
            in_atom_section = True
            continue
        if in_atom_section:
            if line.startswith("@"):
                break
            parts = line.split()
            if len(parts) >= 6:
                atom_count += 1
                atom_type = parts[5]
                if atom_type.startswith("H") or ".H" in atom_type:
                    has_hydrogen = True
    
    if atom_count == 0:
        print(f"  ERROR: No atoms found in {mol2_path}")
        return False
    
    if not has_hydrogen:
        print(f"  WARNING: No hydrogen atoms found in {mol2_path}")
        return False
    
    return True


def run_antechamber(
    input_mol2: Path,
    output_mol2: Path,
    net_charge: int = 0
) -> bool:
    """
    Run antechamber to generate GAFF2 atom types and AM1-BCC charges.
    
    Args:
        input_mol2: Input mol2 file
        output_mol2: Output mol2 file
        net_charge: Net charge of the ligand
        
    Returns:
        True if successful
    """
    cmd = (
        f"antechamber -i {input_mol2} -fi mol2 "
        f"-o {output_mol2} -fo mol2 "
        f"-c {CHARGE_MODEL} -at {LIGAND_FF} -nc {net_charge} "
        f"-rn {LIGAND_RESIDUE_NAME} -pf y"
    )
    
    result = subprocess.run(
        cmd,
        shell=True,
        cwd=output_mol2.parent,
        capture_output=True,
        text=True,
        executable="/bin/bash"
    )
    
    if result.returncode != 0:
        print(f"  ERROR: antechamber failed")
        print(f"  STDERR: {result.stderr}")
        return False
    
    return output_mol2.exists()


def run_parmchk2(mol2_file: Path, frcmod_file: Path) -> bool:
    """
    Run parmchk2 to generate missing force field parameters.
    
    Args:
        mol2_file: Input mol2 file (with GAFF2 atom types)
        frcmod_file: Output frcmod file
        
    Returns:
        True if successful
    """
    cmd = f"parmchk2 -i {mol2_file} -f mol2 -o {frcmod_file} -s {LIGAND_FF}"
    
    result = subprocess.run(
        cmd,
        shell=True,
        cwd=mol2_file.parent,
        capture_output=True,
        text=True,
        executable="/bin/bash"
    )
    
    if result.returncode != 0:
        print(f"  ERROR: parmchk2 failed")
        print(f"  STDERR: {result.stderr}")
        return False
    
    return frcmod_file.exists()


def parameterize_ligand(
    input_mol2: Path,
    ligand_id: str,
    net_charge: int = 0,
    force: bool = False
) -> Tuple[bool, Optional[Path], Optional[Path]]:
    """
    Complete ligand parameterization workflow.
    
    Args:
        input_mol2: Input mol2 file (must have hydrogens)
        ligand_id: Unique identifier for the ligand
        net_charge: Net charge of the ligand
        force: If True, overwrite existing parameters
        
    Returns:
        Tuple of (success, mol2_path, frcmod_path)
    """
    ensure_directories()
    
    output_dir = LIGAND_PARAMS_DIR / ligand_id
    output_mol2 = output_dir / f"{LIGAND_RESIDUE_NAME}.mol2"
    output_frcmod = output_dir / f"{LIGAND_RESIDUE_NAME}.frcmod"
    
    # Check if already done
    if not force and output_mol2.exists() and output_frcmod.exists():
        print(f"  SKIP: Parameters already exist for {ligand_id}")
        return True, output_mol2, output_frcmod
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Validate input
    if not validate_mol2(input_mol2):
        return False, None, None
    
    # Set up AMBER environment
    os.environ["AMBERHOME"] = AMBER_HOME
    os.environ["PATH"] = f"{AMBER_HOME}/bin:" + os.environ.get("PATH", "")
    
    # Step 1: Run antechamber
    print(f"  Running antechamber...")
    if not run_antechamber(input_mol2, output_mol2, net_charge):
        return False, None, None
    
    # Step 2: Run parmchk2
    print(f"  Running parmchk2...")
    if not run_parmchk2(output_mol2, output_frcmod):
        return False, None, None
    
    print(f"  SUCCESS: {output_mol2.name} and {output_frcmod.name} generated")
    return True, output_mol2, output_frcmod


def main():
    """Standalone test."""
    import sys
    
    if len(sys.argv) < 3:
        print("Usage: python step1_prepare_ligand_params.py <input.mol2> <ligand_id>")
        sys.exit(1)
    
    input_mol2 = Path(sys.argv[1])
    ligand_id = sys.argv[2]
    
    success, mol2, frcmod = parameterize_ligand(input_mol2, ligand_id)
    
    if success:
        print(f"\nOutput files:")
        print(f"  mol2: {mol2}")
        print(f"  frcmod: {frcmod}")
    else:
        print("\nParameterization failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()
