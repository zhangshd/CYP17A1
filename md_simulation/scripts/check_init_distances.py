"""
Check Initial Ligand-Fe Distances
Verifies that docked conformations have ligand properly positioned near Fe.
Author: zhangshd
Date: 2024-12-19
"""

import math
from pathlib import Path
from typing import List, Tuple, Optional


def parse_pdb_coords(pdb_file: Path, residue_name: str, atom_name: str = None) -> List[Tuple[str, float, float, float]]:
    """Extract coordinates from PDB file."""
    coords = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                res = line[17:20].strip()
                aname = line[12:16].strip()
                if res == residue_name:
                    if atom_name is None or aname == atom_name:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append((aname, x, y, z))
    return coords


def calc_distance(c1: Tuple[float, float, float], c2: Tuple[float, float, float]) -> float:
    """Calculate Euclidean distance."""
    return math.sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2)


def get_min_distance(pdb_file: Path) -> Tuple[Optional[float], Optional[str], Optional[int]]:
    """
    Get minimum distance between any ligand heavy atom and Fe.
    
    Returns:
        Tuple of (min_distance, closest_atom_name, num_heavy_atoms)
    """
    # Get Fe coordinates - try multiple residue names
    # AMBER uses FE1 for the MCPB-derived heme iron
    fe_coords = None
    for fe_res in ["FE1", "HEC", "HEM"]:
        fe_coords = parse_pdb_coords(pdb_file, fe_res, "FE")
        if fe_coords:
            break
    
    if not fe_coords:
        return None, None, None
    
    fe_xyz = (fe_coords[0][1], fe_coords[0][2], fe_coords[0][3])
    
    # Get ligand coordinates - try multiple residue names
    # AMBER uses LIG for ligand, or MOL, or custom 3-letter codes
    lig_coords = None
    for lig_res in ["LIG", "MOL", "UNK", "L32", "L33", "L34", "L35"]:
        lig_coords = parse_pdb_coords(pdb_file, lig_res)
        if lig_coords:
            break
    
    # If still not found, look for residue 476 (standard ligand position)
    if not lig_coords:
        with open(pdb_file) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        resnum = int(line[22:26].strip())
                        if resnum == 476:
                            res = line[17:20].strip()
                            lig_coords = parse_pdb_coords(pdb_file, res)
                            if lig_coords:
                                break
                    except ValueError:
                        continue
    
    if not lig_coords:
        return None, None, None
    
    # Exclude hydrogens
    lig_heavy = [(a, x, y, z) for a, x, y, z in lig_coords if not a.startswith("H")]
    
    if not lig_heavy:
        return None, None, None
    
    # Find minimum distance
    min_dist = float('inf')
    closest_atom = None
    for aname, x, y, z in lig_heavy:
        d = calc_distance((x, y, z), fe_xyz)
        if d < min_dist:
            min_dist = d
            closest_atom = aname
    
    return min_dist, closest_atom, len(lig_heavy)


def main():
    """Check all systems."""
    import sys
    sys.path.insert(0, str(Path(__file__).parent))
    from config import SYSTEMS_DIR
    
    print("=" * 70)
    print("Initial Structure: Ligand-Fe Minimum Distance Check")
    print("=" * 70)
    print()
    print(f"{'System':<20} {'Min_Dist(Å)':<12} {'Closest_Atom':<15} {'Status':<15}")
    print("-" * 70)
    
    valid_count = 0
    invalid_count = 0
    
    for sys_dir in sorted(SYSTEMS_DIR.iterdir()):
        if not sys_dir.is_dir():
            continue
        
        pdb_file = sys_dir / "complex_dry.pdb"
        if not pdb_file.exists():
            continue
        
        min_dist, closest, n_atoms = get_min_distance(pdb_file)
        
        if min_dist is not None:
            # For Type II inhibitors, expected Fe-N distance is ~2.1-2.3 Å
            # For general binding, < 5 Å indicates proper active site placement
            if min_dist < 5.0:
                status = "✓ Good"
                valid_count += 1
            elif min_dist < 8.0:
                status = "⚠ Marginal"
                valid_count += 1
            else:
                status = "✗ Too far"
                invalid_count += 1
            
            print(f"{sys_dir.name:<20} {min_dist:<12.2f} {closest:<15} {status:<15}")
        else:
            print(f"{sys_dir.name:<20} {'N/A':<12} {'N/A':<15} {'ERROR':<15}")
            invalid_count += 1
    
    print("-" * 70)
    print(f"Valid (< 8 Å): {valid_count}")
    print(f"Invalid (>= 8 Å or error): {invalid_count}")
    print()
    print("Note: For Type II CYP inhibitors, Fe-N coordination distance should be ~2.1-2.3 Å")
    print("      Distances > 5 Å suggest the ligand is not in coordination distance")


if __name__ == "__main__":
    main()
