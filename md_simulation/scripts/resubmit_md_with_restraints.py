"""
Resubmit MD jobs with Fe-Ligand distance restraints.
Excludes ligand_* systems which already have valid results.
Author: zhangshd
Date: 2024-12-19
"""

import os
import subprocess
from pathlib import Path
from typing import List

from config import SYSTEMS_DIR
from step3_setup_md_simulation import setup_md_simulation


def get_systems_to_rerun() -> List[str]:
    """
    Get list of systems that need to be re-run.
    Excludes ligand_* systems (already have valid results from previous batch).
    """
    systems = []
    
    for system_dir in sorted(SYSTEMS_DIR.iterdir()):
        if not system_dir.is_dir():
            continue
        
        ligand_id = system_dir.name
        
        # Skip ligand_* systems - they already have valid results
        if ligand_id.startswith("ligand_"):
            print(f"  SKIP: {ligand_id} (previous batch with valid results)")
            continue
        
        # Check if topology exists
        prmtop = system_dir / "complex_solv.prmtop"
        if not prmtop.exists():
            print(f"  SKIP: {ligand_id} (no topology)")
            continue
        
        systems.append(ligand_id)
    
    return systems


def clean_previous_md_files(system_dir: Path) -> None:
    """Remove previous MD output files to allow fresh start."""
    patterns = [
        "*.out", "*.rst7", "*.nc", "md_*.log",
        "*.in", "run_md.sh", "dist_restraint.RST"
    ]
    
    for pattern in patterns:
        for f in system_dir.glob(pattern):
            # Keep topology and coordinate files
            if f.name in ["complex_solv.prmtop", "complex_solv.inpcrd"]:
                continue
            f.unlink()


def submit_md_job(system_dir: Path) -> bool:
    """Submit MD job using sbatch."""
    run_script = system_dir / "run_md.sh"
    
    if not run_script.exists():
        return False
    
    result = subprocess.run(
        ["sbatch", "run_md.sh"],
        cwd=system_dir,
        capture_output=True,
        text=True
    )
    
    if result.returncode == 0:
        # Extract job ID from output
        job_info = result.stdout.strip()
        print(f"    Submitted: {job_info}")
        return True
    else:
        print(f"    ERROR: {result.stderr}")
        return False


def main():
    """Main entry point."""
    print("=" * 70)
    print("Resubmitting MD jobs with Fe-Ligand distance restraints")
    print("=" * 70)
    
    print("\nIdentifying systems to re-run...")
    print("-" * 70)
    systems = get_systems_to_rerun()
    
    print(f"\nFound {len(systems)} systems to re-run")
    print("-" * 70)
    
    if not systems:
        print("No systems to process.")
        return
    
    # Confirm with user
    print("\nSystems to be re-run:")
    for s in systems:
        print(f"  - {s}")
    
    print(f"\nThis will clean previous MD files and resubmit {len(systems)} jobs.")
    response = input("Continue? [y/N]: ").strip().lower()
    
    if response != 'y':
        print("Aborted.")
        return
    
    # Process each system
    print("\n" + "=" * 70)
    print("Processing systems...")
    print("=" * 70)
    
    submitted = 0
    failed = 0
    
    for i, ligand_id in enumerate(systems, 1):
        print(f"\n[{i}/{len(systems)}] {ligand_id}")
        system_dir = SYSTEMS_DIR / ligand_id
        
        # Clean previous files
        print("  Cleaning previous MD files...")
        clean_previous_md_files(system_dir)
        
        # Setup new MD with restraints
        print("  Setting up MD with distance restraints...")
        success = setup_md_simulation(ligand_id, force=True)
        
        if not success:
            print("  FAILED: Could not setup MD")
            failed += 1
            continue
        
        # Submit job
        print("  Submitting job...")
        if submit_md_job(system_dir):
            submitted += 1
        else:
            failed += 1
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total systems: {len(systems)}")
    print(f"Successfully submitted: {submitted}")
    print(f"Failed: {failed}")
    
    if submitted > 0:
        print("\nMonitor jobs with: squeue -u $USER")
        print("Check individual logs in each system directory: md_*.log")


if __name__ == "__main__":
    main()
