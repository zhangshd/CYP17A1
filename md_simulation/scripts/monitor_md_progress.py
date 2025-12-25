"""
Monitor MD Simulation Progress
Checks the current stage of all running MD simulations.
Author: zhangshd
Date: 2024-12-19
"""

from pathlib import Path
from typing import Dict, List, Tuple
import re

# Import configuration
from config import SYSTEMS_DIR


def extract_current_stage(log_file: Path) -> Tuple[str, bool]:
    """
    Extract the current stage from MD log file.
    
    Returns:
        Tuple of (stage_name, has_error)
    """
    if not log_file.exists():
        return ("Not started", False)
    
    try:
        with open(log_file, 'r') as f:
            lines = f.readlines()
    except Exception:
        return ("Error reading log", True)
    
    if not lines:
        return ("Empty log", True)
    
    # Check for errors (but ignore floating-point exceptions which are normal)
    has_error = False
    for line in lines:
        if "ERROR:" in line:
            # Only check for uppercase ERROR:
            has_error = True
            break
    
    # Find the last stage mentioned
    current_stage = "Starting"
    for line in reversed(lines):
        line = line.strip()
        
        # Match stage patterns
        if "Starting minimization stage 1" in line:
            current_stage = "min1"
            break
        elif "Starting minimization stage 2" in line:
            current_stage = "min2"
            break
        elif "Starting heating" in line:
            current_stage = "heat"
            break
        elif "Starting equilibration stage 1" in line:
            current_stage = "eq1"
            break
        elif re.search(r"Starting equilibration stage 2\.(\d+)", line):
            match = re.search(r"Starting equilibration stage 2\.(\d+)", line)
            current_stage = f"eq2_{match.group(1)}"
            break
        elif "Starting production MD" in line:
            current_stage = "prod"
            break
        elif "MD simulation completed successfully" in line:
            current_stage = "COMPLETED"
            break
        elif "Stage" in line and "completed successfully" in line:
            # Try to extract stage name from "Stage xxx completed successfully"
            match = re.search(r"Stage (\S+) completed successfully", line)
            if match:
                current_stage = f"{match.group(1)} (done)"
    
    return (current_stage, has_error)


def monitor_all_systems() -> Dict[str, Tuple[str, bool]]:
    """
    Monitor all MD systems.
    
    Returns:
        Dictionary mapping ligand_id to (stage, has_error)
    """
    results = {}
    
    for system_dir in sorted(SYSTEMS_DIR.glob("*")):
        if not system_dir.is_dir():
            continue
        
        ligand_id = system_dir.name
        
        # Find the most recent log file
        log_files = sorted(system_dir.glob("md_*.log"))
        if not log_files:
            results[ligand_id] = ("No log file", False)
            continue
        
        latest_log = log_files[-1]
        stage, has_error = extract_current_stage(latest_log)
        results[ligand_id] = (stage, has_error)
    
    return results


def print_summary(results: Dict[str, Tuple[str, bool]]):
    """Print a summary of all systems."""
    print("\n" + "=" * 80)
    print("MD Simulation Progress Summary")
    print("=" * 80)
    
    # Group by stage
    stage_counts = {}
    error_systems = []
    
    for ligand_id, (stage, has_error) in results.items():
        if has_error:
            error_systems.append(ligand_id)
        
        if stage not in stage_counts:
            stage_counts[stage] = []
        stage_counts[stage].append(ligand_id)
    
    # Print stage summary
    print(f"\nTotal systems: {len(results)}")
    print(f"Systems with errors: {len(error_systems)}")
    print("\nStage distribution:")
    for stage, systems in sorted(stage_counts.items()):
        print(f"  {stage:20s}: {len(systems):3d} systems")
    
    # Print detailed list
    print("\n" + "-" * 80)
    print(f"{'Ligand ID':<20s} {'Current Stage':<20s} {'Status'}")
    print("-" * 80)
    
    for ligand_id, (stage, has_error) in sorted(results.items()):
        status = "ERROR" if has_error else "OK"
        print(f"{ligand_id:<20s} {stage:<20s} {status}")
    
    print("=" * 80)
    
    # Print error details if any
    if error_systems:
        print("\n⚠ WARNING: The following systems have errors:")
        for ligand_id in error_systems:
            print(f"  - {ligand_id}")
        print("\nCheck individual log files for details.")
    
    print("\n")


def main():
    """Main function."""
    results = monitor_all_systems()
    print_summary(results)


if __name__ == "__main__":
    main()
