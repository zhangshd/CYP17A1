"""
Regenerate GPU version of MD input files for all systems.
Does not submit jobs - just updates the input files and SLURM scripts.
Author: zhangshd
Date: 2024-12-24
"""

from pathlib import Path
from typing import List
from step3_setup_md_simulation import setup_md_simulation
from config import SYSTEMS_DIR


def get_all_md_systems() -> List[str]:
    """Get list of all MD systems (excluding ligand_* which are already complete)."""
    systems = []
    
    for system_dir in sorted(SYSTEMS_DIR.iterdir()):
        if not system_dir.is_dir():
            continue
        
        ligand_id = system_dir.name
        
        # Skip ligand_* systems
        if ligand_id.startswith("ligand_"):
            continue
        
        # Check if topology exists
        prmtop = system_dir / "complex_solv.prmtop"
        if not prmtop.exists():
            continue
        
        systems.append(ligand_id)
    
    return systems


def regenerate_inputs_for_system(ligand_id: str) -> bool:
    """
    Regenerate GPU input files for a single system.
    Only updates input files (.in), SLURM script (run_md.sh), and restraint file.
    Does not remove MD output files (.out, .rst7, .nc).
    """
    system_dir = SYSTEMS_DIR / ligand_id
    
    print(f"  Regenerating GPU inputs for {ligand_id}...")
    
    try:
        # Only remove input files that will be regenerated
        for pattern in ["*.in", "run_md.sh", "dist_restraint.RST"]:
            for f in system_dir.glob(pattern):
                f.unlink()
        
        # Generate new GPU input files
        setup_md_simulation(system_dir, ligand_id)
        return True
    
    except Exception as e:
        print(f"    ERROR: {e}")
        return False


def main():
    """Main function."""
    print("=" * 80)
    print("Regenerating GPU version of MD input files")
    print("=" * 80)
    print("\nNote: This only updates input files (.in), SLURM scripts, and restraint files.")
    print("MD output files (.out, .rst7, .nc) are preserved.")
    print()
    
    systems = get_all_md_systems()
    print(f"Found {len(systems)} systems to process\n")
    
    success_count = 0
    failed_count = 0
    
    for i, ligand_id in enumerate(systems, 1):
        print(f"[{i}/{len(systems)}] {ligand_id}")
        if regenerate_inputs_for_system(ligand_id):
            success_count += 1
        else:
            failed_count += 1
    
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total systems: {len(systems)}")
    print(f"Successfully updated: {success_count}")
    print(f"Failed: {failed_count}")
    print("\nGPU input files are ready.")
    print("To submit GPU jobs, run:")
    print("  cd /home/zhangsd/repos/CYP17A1/md_simulation/scripts")
    print("  python3 resubmit_md_with_restraints.py")
    print()


if __name__ == "__main__":
    main()
