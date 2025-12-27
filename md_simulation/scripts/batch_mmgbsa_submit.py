#!/usr/bin/env python3
"""
Batch Setup and Submit MM/GBSA Calculations
Check MD completion, setup MM/GBSA inputs, and submit jobs.

Author: zhangshd
Date: 2024-12-26
"""

import sys
from pathlib import Path
from typing import List, Tuple
import subprocess
import time

# Add scripts directory to path
SCRIPTS_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPTS_DIR))

from config import SYSTEMS_DIR
from step4_run_mmgbsa import setup_mmgbsa_analysis


def read_ligand_list(list_file: Path) -> List[str]:
    """Read ligand IDs from file."""
    ligand_ids = []
    with open(list_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                ligand_ids.append(line)
    return ligand_ids


def check_md_completion(system_dir: Path) -> Tuple[bool, str]:
    """
    Check if MD simulation is complete.
    
    Returns:
        (is_complete, status_message)
    """
    prod_nc = system_dir / "prod.nc"
    
    if not prod_nc.exists():
        return False, "prod.nc not found"
    
    # Check file size (should be > 100MB for 50ns)
    size_mb = prod_nc.stat().st_size / (1024 * 1024)
    
    if size_mb < 10:
        return False, f"prod.nc too small ({size_mb:.1f} MB)"
    
    return True, f"Complete ({size_mb:.1f} MB)"


def submit_mmgbsa_job(system_dir: Path, ligand_id: str) -> bool:
    """Submit MM/GBSA SLURM job."""
    run_script = system_dir / "run_mmgbsa.sh"
    
    if not run_script.exists():
        print(f"    ERROR: run_mmgbsa.sh not found")
        return False
    
    try:
        result = subprocess.run(
            ["sbatch", "run_mmgbsa.sh"],
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
    
    parser = argparse.ArgumentParser(description="Batch setup and submit MM/GBSA calculations")
    parser.add_argument("ligand_list", type=Path, help="File containing ligand IDs (one per line)")
    parser.add_argument("--skip-setup", action="store_true", help="Skip setup, only submit jobs")
    parser.add_argument("--skip-submit", action="store_true", help="Skip submission, only setup")
    parser.add_argument("--check-only", action="store_true", help="Only check MD completion status")
    parser.add_argument("--delay", type=int, default=0, help="Delay (seconds) between submissions")
    parser.add_argument("--force", action="store_true", help="Overwrite existing MM/GBSA files")
    
    args = parser.parse_args()
    
    if not args.ligand_list.exists():
        print(f"ERROR: Ligand list file not found: {args.ligand_list}")
        sys.exit(1)
    
    ligand_ids = read_ligand_list(args.ligand_list)
    
    print("="*80)
    print(f"Batch MM/GBSA Setup and Submission")
    print(f"Ligand list: {args.ligand_list}")
    print(f"Total ligands: {len(ligand_ids)}")
    print("="*80)
    print()
    
    stats = {
        "total": len(ligand_ids),
        "md_complete": 0,
        "md_incomplete": 0,
        "setup_success": 0,
        "setup_failed": 0,
        "submit_success": 0,
        "submit_failed": 0,
        "skipped": 0
    }
    
    incomplete_list = []
    
    for i, ligand_id in enumerate(ligand_ids, 1):
        print(f"[{i}/{len(ligand_ids)}] Processing {ligand_id}...")
        
        system_dir = SYSTEMS_DIR / ligand_id
        
        # Check if system exists
        if not system_dir.exists():
            print(f"    WARNING: System directory not found, skipping")
            stats["skipped"] += 1
            continue
        
        # Check MD completion
        md_complete, md_status = check_md_completion(system_dir)
        print(f"    MD status: {md_status}")
        
        if not md_complete:
            stats["md_incomplete"] += 1
            incomplete_list.append((ligand_id, md_status))
            continue
        
        stats["md_complete"] += 1
        
        # If only checking, skip the rest
        if args.check_only:
            continue
        
        # Setup MM/GBSA
        if not args.skip_setup:
            print(f"    Setting up MM/GBSA inputs...")
            success = setup_mmgbsa_analysis(ligand_id, force=args.force)
            
            if success:
                stats["setup_success"] += 1
            else:
                stats["setup_failed"] += 1
                continue
        
        # Submit job
        if not args.skip_submit:
            print(f"    Submitting MM/GBSA job...")
            success = submit_mmgbsa_job(system_dir, ligand_id)
            
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
    print(f"Total ligands:      {stats['total']}")
    print(f"Skipped:            {stats['skipped']}")
    print(f"MD complete:        {stats['md_complete']}")
    print(f"MD incomplete:      {stats['md_incomplete']}")
    
    if not args.check_only and not args.skip_setup:
        print(f"Setup successful:   {stats['setup_success']}")
        print(f"Setup failed:       {stats['setup_failed']}")
    
    if not args.check_only and not args.skip_submit:
        print(f"Submit successful:  {stats['submit_success']}")
        print(f"Submit failed:      {stats['submit_failed']}")
    
    if incomplete_list:
        print()
        print("="*80)
        print("INCOMPLETE MD SIMULATIONS")
        print("="*80)
        for ligand_id, status in incomplete_list:
            print(f"  {ligand_id:20s} - {status}")
    
    print("="*80)
    
    if stats['setup_failed'] > 0 or stats['submit_failed'] > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()

