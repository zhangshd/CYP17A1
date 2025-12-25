#!/usr/bin/env python3
"""
Batch Setup and Submit MD Simulations
Clean old files, regenerate MD inputs, and submit jobs.

Author: zhangshd
Date: 2024-12-26
"""

import sys
from pathlib import Path
from typing import List
import subprocess
import time

# Add scripts directory to path
SCRIPTS_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPTS_DIR))

from config import SYSTEMS_DIR
from step3_setup_md_simulation import setup_md_simulation


def read_ligand_list(list_file: Path) -> List[str]:
    """Read ligand IDs from file."""
    ligand_ids = []
    with open(list_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                ligand_ids.append(line)
    return ligand_ids


def clean_old_files(system_dir: Path) -> None:
    """Clean old MD output files."""
    patterns = [
        "*.rst7", "*.out", "*.nc", "*.pdb",
        "min*.mdout", "heat.mdout", "eq*.mdout", "prod.mdout",
        "*_fe_lig_dist.dat", "*_fe_n22.dat"
    ]
    
    cleaned = []
    for pattern in patterns:
        for file in system_dir.glob(pattern):
            # Keep complex_dry.pdb and complex_solv.pdb
            if file.name in ("complex_dry.pdb", "complex_solv.pdb", "complex.pdb"):
                continue
            try:
                file.unlink()
                cleaned.append(file.name)
            except Exception as e:
                print(f"    Warning: Could not delete {file.name}: {e}")
    
    if cleaned:
        print(f"    Cleaned {len(cleaned)} old files")


def submit_job(system_dir: Path, ligand_id: str) -> bool:
    """Submit SLURM job."""
    run_script = system_dir / "run_md.sh"
    
    if not run_script.exists():
        print(f"    ERROR: run_md.sh not found")
        return False
    
    try:
        result = subprocess.run(
            ["sbatch", "run_md.sh"],
            cwd=system_dir,
            capture_output=True,
            text=True
        )
        
        if result.returncode == 0:
            # Extract job ID from output
            job_id = result.stdout.strip().split()[-1]
            print(f"    ✓ Submitted job {job_id}")
            return True
        else:
            print(f"    ERROR: sbatch failed: {result.stderr}")
            return False
            
    except Exception as e:
        print(f"    ERROR: Failed to submit job: {e}")
        return False


def main():
    """Main batch processing function."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Batch setup and submit MD simulations")
    parser.add_argument("ligand_list", type=Path, help="File containing ligand IDs (one per line)")
    parser.add_argument("--skip-setup", action="store_true", help="Skip MD setup, only submit jobs")
    parser.add_argument("--skip-submit", action="store_true", help="Skip job submission, only setup")
    parser.add_argument("--delay", type=int, default=0, help="Delay (seconds) between submissions")
    
    args = parser.parse_args()
    
    if not args.ligand_list.exists():
        print(f"ERROR: Ligand list file not found: {args.ligand_list}")
        sys.exit(1)
    
    ligand_ids = read_ligand_list(args.ligand_list)
    
    print("="*80)
    print(f"Batch MD Setup and Submission")
    print(f"Ligand list: {args.ligand_list}")
    print(f"Total ligands: {len(ligand_ids)}")
    print("="*80)
    print()
    
    stats = {
        "total": len(ligand_ids),
        "setup_success": 0,
        "setup_failed": 0,
        "submit_success": 0,
        "submit_failed": 0,
        "skipped": 0
    }
    
    for i, ligand_id in enumerate(ligand_ids, 1):
        print(f"[{i}/{len(ligand_ids)}] Processing {ligand_id}...")
        
        system_dir = SYSTEMS_DIR / ligand_id
        
        # Check if system exists
        if not system_dir.exists():
            print(f"    WARNING: System directory not found, skipping")
            stats["skipped"] += 1
            continue
        
        # Check if topology exists
        prmtop = system_dir / "complex_solv.prmtop"
        if not prmtop.exists():
            print(f"    WARNING: Topology not found, skipping")
            stats["skipped"] += 1
            continue
        
        # Clean old files
        if not args.skip_setup:
            print(f"    Cleaning old files...")
            clean_old_files(system_dir)
        
        # Setup MD
        if not args.skip_setup:
            print(f"    Setting up MD inputs...")
            success = setup_md_simulation(ligand_id, force=True)
            
            if success:
                stats["setup_success"] += 1
            else:
                print(f"    ERROR: MD setup failed")
                stats["setup_failed"] += 1
                continue
        
        # Submit job
        if not args.skip_submit:
            print(f"    Submitting job...")
            success = submit_job(system_dir, ligand_id)
            
            if success:
                stats["submit_success"] += 1
            else:
                stats["submit_failed"] += 1
            
            # Delay between submissions
            if args.delay > 0 and i < len(ligand_ids):
                time.sleep(args.delay)
        
        print()
    
    # Print summary
    print("="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Total ligands:     {stats['total']}")
    print(f"Skipped:           {stats['skipped']}")
    
    if not args.skip_setup:
        print(f"Setup successful:  {stats['setup_success']}")
        print(f"Setup failed:      {stats['setup_failed']}")
    
    if not args.skip_submit:
        print(f"Submit successful: {stats['submit_success']}")
        print(f"Submit failed:     {stats['submit_failed']}")
    
    print("="*80)
    
    if stats['setup_failed'] > 0 or stats['submit_failed'] > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()

