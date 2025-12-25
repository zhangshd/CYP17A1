"""
Step 2: Build Complex Topology
Generate AMBER topology files for protein-ligand complexes.
Uses MCPB.py-generated protein topology with Fe-Heme parameters.
Author: zhangshd
Date: 2024-12-17
"""

import os
import subprocess
from pathlib import Path
from typing import Optional, Tuple, List

from config import (
    AMBER_HOME, AMBER_SOURCE_CMD,
    MCPB_DIR, MCPB_FRCMOD, MCPB_PDB, MCPB_TEMPLATES_DIR,
    LIGAND_PARAMS_DIR, SYSTEMS_DIR,
    LIGAND_RESIDUE_NAME, LIGAND_FF, WATER_MODEL, BUFFER_SIZE,
    MCPB_BONDS, CUSTOM_ATOM_TYPES, LIGAND_FE_COORDINATION, ensure_directories
)


def read_mol2_coordinates(mol2_path: Path) -> List[Tuple[str, float, float, float, str]]:
    """
    Read atom names and coordinates from mol2 file.
    
    Returns:
        List of (atom_name, x, y, z, atom_type) tuples
    """
    atoms = []
    in_atom_section = False
    
    with open(mol2_path, 'r') as f:
        for line in f:
            if "@<TRIPOS>ATOM" in line:
                in_atom_section = True
                continue
            if in_atom_section:
                if line.startswith("@"):
                    break
                parts = line.split()
                if len(parts) >= 6:
                    atom_name = parts[1]
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    atom_type = parts[5]
                    atoms.append((atom_name, x, y, z, atom_type))
    
    return atoms


def read_pdb_coordinates(pdb_path: Path) -> List[Tuple[str, float, float, float]]:
    """
    Read atom names and coordinates from PDB file.
    
    Returns:
        List of (atom_name, x, y, z) tuples
    """
    atoms = []
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append((atom_name, x, y, z))
    
    return atoms


def create_complex_pdb(
    protein_pdb: Path,
    ligand_atoms: List[Tuple[str, float, float, float, str]],
    output_pdb: Path,
    lig_res_name: str = "LIG"
) -> bool:
    """
    Combine protein PDB and ligand atoms into a single PDB file.
    
    Args:
        protein_pdb: Path to protein PDB (from MCPB.py)
        ligand_atoms: List of (atom_name, x, y, z, atom_type) from parameterized mol2
        output_pdb: Output combined PDB path
        lig_res_name: Ligand residue name
        
    Returns:
        True if successful
    """
    # Read protein PDB
    with open(protein_pdb, 'r') as f:
        protein_lines = f.readlines()
    
    # Get atom lines only
    protein_atoms = [l for l in protein_lines if l.startswith(('ATOM', 'HETATM'))]
    
    # Get max atom and residue numbers
    max_atom_num = 0
    max_res_num = 0
    for line in protein_atoms:
        if len(line) >= 26:
            atom_num = int(line[6:11].strip())
            res_num = int(line[22:26].strip())
            max_atom_num = max(max_atom_num, atom_num)
            max_res_num = max(max_res_num, res_num)
    
    # Generate ligand PDB lines
    ligand_lines = []
    new_res_num = max_res_num + 1
    
    for i, (atom_name, x, y, z, atom_type) in enumerate(ligand_atoms):
        atom_num = max_atom_num + i + 1
        
        # Get element from atom type
        elem = ''.join(c for c in atom_name if c.isalpha())[:2]
        if len(elem) == 2:
            elem = elem[0].upper() + elem[1].lower()
        
        # Format atom name for PDB (4 characters, right justified for 1-char elements)
        if len(atom_name) < 4:
            if len(elem) == 1:
                atom_name_fmt = f" {atom_name:<3s}"
            else:
                atom_name_fmt = f"{atom_name:<4s}"
        else:
            atom_name_fmt = atom_name[:4]
        
        # Standard PDB format - CRITICAL: chain ID at column 22 (index 21)
        line = f"HETATM{atom_num:5d} {atom_name_fmt} {lig_res_name:<3s} A{new_res_num:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}\n"
        ligand_lines.append(line)
    
    # Write combined PDB
    with open(output_pdb, 'w') as f:
        for line in protein_atoms:
            f.write(line)
        f.write("TER\n")
        for line in ligand_lines:
            f.write(line)
        f.write("TER\n")
        f.write("END\n")
    
    return True


def generate_tleap_script(
    ligand_id: str,
    ligand_mol2: Path,
    ligand_frcmod: Path,
    complex_pdb: Path,
    output_dir: Path
) -> Path:
    """
    Generate tleap input script for topology building.
    
    Returns:
        Path to tleap script
    """
    tleap_script = output_dir / "tleap.in"
    
    # Template mol2 files from MCPB.py
    cm1_mol2 = MCPB_TEMPLATES_DIR / "CM1.mol2"
    hm1_mol2 = MCPB_TEMPLATES_DIR / "HM1.mol2"
    fe1_mol2 = MCPB_TEMPLATES_DIR / "FE1.mol2"
    
    content = f"""# tleap script for {ligand_id} complex
# Generated by step2_build_complex_topology.py

# Load force fields
source oldff/leaprc.ff14SB
source leaprc.{LIGAND_FF}
source leaprc.water.opc

# Add custom atom types for Fe-Heme
{CUSTOM_ATOM_TYPES}

# Load MCPB.py templates (residue types)
CM1 = loadmol2 {cm1_mol2}
HM1 = loadmol2 {hm1_mol2}
FE1 = loadmol2 {fe1_mol2}

# Load MCPB.py parameters for Fe-Heme
loadamberparams {MCPB_DIR}/HEM.frcmod
loadamberparams frcmod.ionslm_126_opc
loadamberparams {MCPB_FRCMOD}

# Load ligand parameters
{LIGAND_RESIDUE_NAME} = loadmol2 {ligand_mol2}
loadamberparams {ligand_frcmod}

# Load complex structure
mol = loadpdb {complex_pdb}

# Remove OXT atom from C-terminal (it's not in standard HIE template)
remove mol mol.473.OXT

# Create Fe-Heme-Cys bonds
"""
    
    # Add MCPB bond commands
    for res1_atom, res2_atom in MCPB_BONDS:
        content += f"bond mol.{res1_atom} mol.{res2_atom}\n"
    
    # Note: For Type II inhibitors with Fe coordination, we do NOT create a covalent bond
    # because MCPB.py parameters do not include Fe-ligand torsions.
    # Instead, use distance restraints in MD (handled by step3_setup_md_simulation.py)
    
    content += f"""
# Check structure
check mol

# Save dry topology
saveamberparm mol {output_dir}/complex_dry.prmtop {output_dir}/complex_dry.inpcrd
savepdb mol {output_dir}/complex_dry.pdb

# Solvate with OPC water (12 Å buffer)
solvatebox mol OPCBOX {BUFFER_SIZE}

# Add counter ions to neutralize
addions mol Na+ 0
addions mol Cl- 0

# Save solvated topology
saveamberparm mol {output_dir}/complex_solv.prmtop {output_dir}/complex_solv.inpcrd
savepdb mol {output_dir}/complex_solv.pdb

quit
"""
    
    tleap_script.write_text(content)
    return tleap_script


def run_tleap(script_path: Path) -> Tuple[bool, str]:
    """
    Execute tleap script.
    
    Returns:
        (success, log_output)
    """
    cmd = f"{AMBER_SOURCE_CMD} && tleap -f {script_path}"
    
    result = subprocess.run(
        cmd,
        shell=True,
        cwd=script_path.parent,
        capture_output=True,
        text=True,
        executable="/bin/bash"
    )
    
    log_output = result.stdout + result.stderr
    
    # Check for errors
    if "FATAL" in log_output or result.returncode != 0:
        return False, log_output
    
    return True, log_output


def build_complex_topology(
    ligand_id: str,
    ligand_mol2: Path,
    force: bool = False
) -> Tuple[bool, Optional[Path]]:
    """
    Complete topology building workflow.
    
    Args:
        ligand_id: Unique identifier for the ligand
        ligand_mol2: Path to ligand mol2 with docking coordinates
        force: If True, overwrite existing topology
        
    Returns:
        Tuple of (success, system_dir)
    """
    ensure_directories()
    
    # Check ligand parameters exist
    param_dir = LIGAND_PARAMS_DIR / ligand_id
    param_mol2 = param_dir / f"{LIGAND_RESIDUE_NAME}.mol2"
    param_frcmod = param_dir / f"{LIGAND_RESIDUE_NAME}.frcmod"
    
    if not param_mol2.exists() or not param_frcmod.exists():
        print(f"  ERROR: Ligand parameters not found for {ligand_id}")
        print(f"  Run step1_prepare_ligand_params.py first")
        return False, None
    
    # Create system directory
    system_dir = SYSTEMS_DIR / ligand_id
    system_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if already done
    prmtop_file = system_dir / "complex_solv.prmtop"
    if not force and prmtop_file.exists():
        print(f"  SKIP: Topology already exists for {ligand_id}")
        return True, system_dir
    
    # Read ligand coordinates from original file (docking pose)
    # but use atom names from parameterized mol2
    param_atoms = read_mol2_coordinates(param_mol2)
    
    # Support both mol2 and PDB input formats
    if ligand_mol2.suffix.lower() == '.pdb':
        print(f"  Reading coordinates from PDB file...")
        orig_atoms_pdb = read_pdb_coordinates(ligand_mol2)
        # Convert to mol2-like format (add dummy atom type)
        orig_atoms = [(name, x, y, z, "") for name, x, y, z in orig_atoms_pdb]
    else:
        orig_atoms = read_mol2_coordinates(ligand_mol2)
    
    if len(param_atoms) != len(orig_atoms):
        print(f"  ERROR: Atom count mismatch between param ({len(param_atoms)}) and original ({len(orig_atoms)})")
        return False, None
    
    # Use atom names from param_mol2, coordinates from original
    combined_atoms = []
    for (name, _, _, _, atype), (_, x, y, z, _) in zip(param_atoms, orig_atoms):
        combined_atoms.append((name, x, y, z, atype))
    
    # Create complex PDB
    complex_pdb = system_dir / "complex.pdb"
    print(f"  Creating complex PDB...")
    if not create_complex_pdb(MCPB_PDB, combined_atoms, complex_pdb):
        print(f"  ERROR: Failed to create complex PDB")
        return False, None
    
    # Generate tleap script
    print(f"  Generating tleap script...")
    tleap_script = generate_tleap_script(
        ligand_id, param_mol2, param_frcmod, complex_pdb, system_dir
    )
    
    # Run tleap
    print(f"  Running tleap...")
    success, log = run_tleap(tleap_script)
    
    # Save log
    (system_dir / "leap.log").write_text(log)
    
    if not success:
        print(f"  ERROR: tleap failed")
        print(f"  See {system_dir}/leap.log for details")
        return False, None
    
    # Verify outputs
    if not prmtop_file.exists():
        print(f"  ERROR: Topology file not created")
        return False, None
    
    # Report file sizes
    dry_size = (system_dir / "complex_dry.prmtop").stat().st_size / 1024
    solv_size = prmtop_file.stat().st_size / 1024
    print(f"  SUCCESS: dry={dry_size:.0f}KB, solv={solv_size:.0f}KB")
    
    return True, system_dir


def main():
    """Standalone test."""
    import sys
    
    if len(sys.argv) < 3:
        print("Usage: python step2_build_complex_topology.py <ligand.mol2> <ligand_id> [--force]")
        sys.exit(1)
    
    ligand_mol2 = Path(sys.argv[1])
    ligand_id = sys.argv[2]
    force = "--force" in sys.argv or "-f" in sys.argv
    
    success, system_dir = build_complex_topology(ligand_id, ligand_mol2, force=force)
    
    if success:
        print(f"\nSystem directory: {system_dir}")
    else:
        print("\nTopology building failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()
