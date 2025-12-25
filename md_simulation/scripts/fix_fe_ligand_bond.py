"""
Fix Fe-Ligand Coordination Bond
Manually add Fe-N22 bond to existing topology using parmed.
Author: zhangshd
Date: 2024-12-25
"""

import sys
from pathlib import Path

def add_fe_ligand_bond(
    prmtop_path: str,
    output_prmtop: str,
    fe_atom_index: int = 7668,
    ligand_n_index: int = 7691,
    bond_force_constant: float = 100.0,  # kcal/mol/Å²
    bond_length: float = 2.10,  # Å (typical Fe-N coordination)
) -> None:
    """
    Add Fe-N coordination bond to topology.
    
    Args:
        prmtop_path: Input topology file
        output_prmtop: Output topology file with added bond
        fe_atom_index: Fe atom index (1-based, from PDB)
        ligand_n_index: Ligand N atom index (1-based)
        bond_force_constant: Bond force constant (kcal/mol/Å²)
        bond_length: Equilibrium bond length (Å)
    """
    try:
        import parmed as pmd
    except ImportError:
        print("ERROR: parmed not available. Install with: conda install -c conda-forge parmed")
        sys.exit(1)
    
    # Load topology
    print(f"\n=== Loading topology: {prmtop_path} ===")
    parm = pmd.load_file(prmtop_path)
    
    # Convert to 0-based index for parmed
    fe_idx = fe_atom_index - 1
    n_idx = ligand_n_index - 1
    
    # Get atom objects
    fe_atom = parm.atoms[fe_idx]
    n_atom = parm.atoms[n_idx]
    
    print(f"\n=== Current atom information ===")
    print(f"Fe atom: index={fe_idx+1}, name={fe_atom.name}, type={fe_atom.type}, residue={fe_atom.residue}")
    print(f"N atom:  index={n_idx+1}, name={n_atom.name}, type={n_atom.type}, residue={n_atom.residue}")
    
    # Check if bond already exists
    existing_bond = None
    for bond in parm.bonds:
        if (bond.atom1 == fe_atom and bond.atom2 == n_atom) or \
           (bond.atom1 == n_atom and bond.atom2 == fe_atom):
            existing_bond = bond
            break
    
    if existing_bond:
        print(f"\n⚠️  Bond already exists: {existing_bond}")
        print(f"    Current parameters: k={existing_bond.type.k:.2f}, req={existing_bond.type.req:.3f}")
        
        # Update parameters
        existing_bond.type.k = bond_force_constant
        existing_bond.type.req = bond_length
        print(f"    Updated parameters: k={bond_force_constant:.2f}, req={bond_length:.3f}")
    else:
        print(f"\n✅ Adding new bond: {fe_atom.name}({fe_idx+1}) <-> {n_atom.name}({n_idx+1})")
        
        # Create bond type
        bond_type = pmd.BondType(bond_force_constant, bond_length)
        
        # Add bond
        new_bond = pmd.Bond(fe_atom, n_atom, type=bond_type)
        parm.bonds.append(new_bond)
        
        # Update bond partners
        if n_atom not in fe_atom.bond_partners:
            fe_atom.bond_partners.append(n_atom)
        if fe_atom not in n_atom.bond_partners:
            n_atom.bond_partners.append(fe_atom)
        
        print(f"    Bond parameters: k={bond_force_constant:.2f} kcal/mol/Å², req={bond_length:.3f} Å")
    
    # Save modified topology
    print(f"\n=== Saving modified topology: {output_prmtop} ===")
    parm.save(output_prmtop, overwrite=True)
    
    # Verify bond count
    fe_bonds = [b for b in parm.bonds if fe_atom in [b.atom1, b.atom2]]
    print(f"\n✅ Fe atom now has {len(fe_bonds)} bonds:")
    for bond in fe_bonds:
        other_atom = bond.atom2 if bond.atom1 == fe_atom else bond.atom1
        print(f"    {fe_atom.name} <-> {other_atom.name} (residue {other_atom.residue.name} {other_atom.residue.number})")
    
    print(f"\n=== Bond addition complete ===")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python fix_fe_ligand_bond.py <ligand_id> [fe_index] [n_index]")
        print("Example: python fix_fe_ligand_bond.py AER601 7668 7691")
        sys.exit(1)
    
    ligand_id = sys.argv[1]
    fe_idx = int(sys.argv[2]) if len(sys.argv) > 2 else 7668
    n_idx = int(sys.argv[3]) if len(sys.argv) > 3 else 7691
    
    # Paths
    md_base = Path("/home/zhangsd/repos/CYP17A1/md_simulation")
    system_dir = md_base / "complexes" / "systems_gpu" / ligand_id
    
    input_prmtop = system_dir / "complex_solv.prmtop"
    output_prmtop = system_dir / "complex_solv_bonded.prmtop"
    
    if not input_prmtop.exists():
        print(f"ERROR: Topology not found: {input_prmtop}")
        sys.exit(1)
    
    print(f"""
╔══════════════════════════════════════════════════════════════╗
║  Fix Fe-Ligand Coordination Bond                            ║
║  System: {ligand_id:49s} ║
╚══════════════════════════════════════════════════════════════╝
    """)
    
    add_fe_ligand_bond(
        prmtop_path=str(input_prmtop),
        output_prmtop=str(output_prmtop),
        fe_atom_index=fe_idx,
        ligand_n_index=n_idx,
        bond_force_constant=400.0,  # Increased from 100 (MCPB too weak at ~50)
        bond_length=2.10,  # Typical Fe-N coordination distance
    )
    
    print(f"\n⚠️  Next steps:")
    print(f"1. Backup original: mv {input_prmtop} {input_prmtop}.nobond")
    print(f"2. Use new topology: mv {output_prmtop} {input_prmtop}")
    print(f"3. Remove distance restraints from *.in files (nmropt=0)")
    print(f"4. Re-run MD simulation")
