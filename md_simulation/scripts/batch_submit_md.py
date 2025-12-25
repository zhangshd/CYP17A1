"""
Batch submit MD jobs for all prepared systems.
Simplified version - no CPU/GPU separation needed.
Author: zhangshd
Date: 2024-12-24
"""

import subprocess
import sys
from pathlib import Path
from config import SYSTEMS_DIR


def get_systems_to_submit():
    """Get list of systems ready for submission."""
    systems = []
    
    if not SYSTEMS_DIR.exists():
        print(f"ERROR: Systems directory not found: {SYSTEMS_DIR}")
        return systems
    
    for system_dir in sorted(SYSTEMS_DIR.iterdir()):
        if not system_dir.is_dir():
            continue
        
        ligand_id = system_dir.name
        run_script = system_dir / "run_md.sh"
        
        if not run_script.exists():
            print(f"  SKIP: {ligand_id} (no run_md.sh)")
            continue
        
        systems.append((ligand_id, system_dir))
    
    return systems


def submit_job(ligand_id: str, system_dir: Path) -> tuple[bool, str]:
    """
    Submit a single MD job.
    
    Returns:
        (success, message): Tuple of success status and message
    """
    try:
        result = subprocess.run(
            ["sbatch", "run_md.sh"],
            cwd=system_dir,
            capture_output=True,
            text=True,
            check=True
        )
        
        # Extract job ID from output
        job_id = result.stdout.strip().split()[-1]
        return True, job_id
    
    except subprocess.CalledProcessError as e:
        return False, e.stderr.strip()
    except Exception as e:
        return False, str(e)


def main():
    """Main function."""
    print("=" * 80)
    print("Batch MD Job Submission")
    print("=" * 80)
    print(f"Systems directory: {SYSTEMS_DIR}")
    print()
    
    # Get systems
    systems = get_systems_to_submit()
    
    if not systems:
        print("No systems ready for submission!")
        sys.exit(1)
    
    print(f"Found {len(systems)} systems ready for submission:")
    for ligand_id, _ in systems:
        print(f"  - {ligand_id}")
    print()
    
    # Confirm
    response = input(f"Submit all {len(systems)} jobs? [y/N]: ").strip().lower()
    if response != 'y':
        print("Cancelled.")
        sys.exit(0)
    
    print()
    print("Submitting jobs...")
    print("-" * 80)
    
    submitted = []
    failed = []
    
    for i, (ligand_id, system_dir) in enumerate(systems, 1):
        print(f"[{i}/{len(systems)}] {ligand_id:20s} ", end="", flush=True)
        
        success, message = submit_job(ligand_id, system_dir)
        
        if success:
            print(f"✓ Job {message}")
            submitted.append((ligand_id, message))
        else:
            print(f"✗ FAILED: {message}")
            failed.append(ligand_id)
    
    # Summary
    print()
    print("=" * 80)
    print("SUBMISSION SUMMARY")
    print("=" * 80)
    print(f"Total systems:        {len(systems)}")
    print(f"Successfully submitted: {len(submitted)}")
    print(f"Failed:                {len(failed)}")
    
    if failed:
        print(f"\nFailed systems:")
        for sys in failed:
            print(f"  - {sys}")
    
    if submitted:
        print(f"\nMonitor jobs with:")
        print(f"  squeue -u $USER")
        print(f"  squeue -u $USER --format='%.8i %.9P %.30j %.8T %.10M %.6D'")
    
    print()


if __name__ == "__main__":
    main()
