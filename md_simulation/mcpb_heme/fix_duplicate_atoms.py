"""
Fix duplicate atoms in MCPB.py generated Gaussian input files.
Removes atoms that are too close to each other (< 0.5 Å).

Author: zhangshd
Date: 2024-12-15
"""

import re
import math
from typing import List, Tuple

def parse_gaussian_com(filename: str) -> Tuple[str, List[Tuple[str, float, float, float]], str]:
    """Parse a Gaussian .com file and extract header, coordinates, and footer."""
    with open(filename, 'r') as f:
        content = f.read()
    
    lines = content.split('\n')
    
    # Find the charge/multiplicity line (format: "charge multiplicity")
    header_lines = []
    coord_start = 0
    for i, line in enumerate(lines):
        header_lines.append(line)
        # Charge/multiplicity line pattern: just two integers
        if re.match(r'^\s*-?\d+\s+\d+\s*$', line.strip()):
            coord_start = i + 1
            break
    
    # Parse coordinates until empty line
    atoms = []
    coord_end = coord_start
    for i in range(coord_start, len(lines)):
        line = lines[i].strip()
        if not line:
            coord_end = i
            break
        parts = line.split()
        if len(parts) >= 4:
            element = parts[0]
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            atoms.append((element, x, y, z))
        coord_end = i + 1
    
    header = '\n'.join(header_lines)
    footer = '\n'.join(lines[coord_end:])
    
    return header, atoms, footer


def distance(a1: Tuple[str, float, float, float], a2: Tuple[str, float, float, float]) -> float:
    """Calculate distance between two atoms."""
    return math.sqrt((a1[1]-a2[1])**2 + (a1[2]-a2[2])**2 + (a1[3]-a2[3])**2)


def remove_duplicate_atoms(atoms: List[Tuple[str, float, float, float]], 
                          min_dist: float = 0.5) -> List[Tuple[str, float, float, float]]:
    """Remove atoms that are too close to each other."""
    n = len(atoms)
    to_remove = set()
    
    # Find pairs that are too close
    for i in range(n):
        if i in to_remove:
            continue
        for j in range(i+1, n):
            if j in to_remove:
                continue
            d = distance(atoms[i], atoms[j])
            if d < min_dist:
                print(f"  Found close pair: atom {i+1} ({atoms[i][0]}) and atom {j+1} ({atoms[j][0]}), distance = {d:.4f} Å")
                # Remove the later one (usually the added capping H)
                to_remove.add(j)
    
    # Build cleaned list
    cleaned = [atoms[i] for i in range(n) if i not in to_remove]
    print(f"  Removed {len(to_remove)} duplicate atoms")
    return cleaned


def write_gaussian_com(filename: str, header: str, atoms: List[Tuple[str, float, float, float]], footer: str):
    """Write a Gaussian .com file."""
    with open(filename, 'w') as f:
        f.write(header)
        f.write('\n')
        for atom in atoms:
            f.write(f"{atom[0]:2s}  {atom[1]:12.6f} {atom[2]:12.6f} {atom[3]:12.6f}\n")
        f.write('\n')
        if footer.strip():
            f.write(footer)


def fix_file(filename: str, min_dist: float = 0.5):
    """Fix a single Gaussian .com file."""
    print(f"Processing {filename}...")
    header, atoms, footer = parse_gaussian_com(filename)
    print(f"  Original: {len(atoms)} atoms")
    
    cleaned = remove_duplicate_atoms(atoms, min_dist)
    print(f"  Cleaned: {len(cleaned)} atoms")
    
    # Backup original
    import shutil
    shutil.copy(filename, filename + '.bak')
    
    # Write cleaned file
    write_gaussian_com(filename, header, cleaned, footer)
    print(f"  Saved to {filename}")


if __name__ == '__main__':
    import sys
    import glob
    
    files = glob.glob('pig_CYP17A1_*.com')
    
    for f in files:
        fix_file(f, min_dist=0.5)
    
    print("\nDone! Original files backed up with .bak extension.")
