#!/usr/bin/env python3
"""
Extract and Compile MM/GBSA Results
Parse FINAL_RESULTS_MMGBSA.dat files and create summary table.

Author: zhangshd
Date: 2024-12-26
"""

import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import csv

# Add scripts directory to path
SCRIPTS_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPTS_DIR))

from config import SYSTEMS_DIR, LIGAND_FE_COORDINATION


def parse_mmgbsa_result(result_file: Path) -> Optional[Dict[str, float]]:
    """
    Parse FINAL_RESULTS_MMGBSA.dat to extract energy components.
    
    Returns:
        Dictionary with energy values or None if parsing fails
    """
    if not result_file.exists():
        return None
    
    content = result_file.read_text()
    
    results = {}
    in_delta_section = False
    
    for line in content.split('\n'):
        line = line.strip()
        
        if 'DELTA Energy Terms' in line or 'DELTA G' in line:
            in_delta_section = True
            continue
        
        if not in_delta_section:
            continue
        
        # Parse energy lines
        if 'VDWAALS' in line:
            parts = line.split()
            if len(parts) >= 4:
                results['vdw'] = float(parts[1])
                results['vdw_std'] = float(parts[2])
        elif 'EEL' in line and 'PEEL' not in line:
            parts = line.split()
            if len(parts) >= 4:
                results['eel'] = float(parts[1])
                results['eel_std'] = float(parts[2])
        elif 'EGB' in line:
            parts = line.split()
            if len(parts) >= 4:
                results['egb'] = float(parts[1])
                results['egb_std'] = float(parts[2])
        elif 'ESURF' in line:
            parts = line.split()
            if len(parts) >= 4:
                results['esurf'] = float(parts[1])
                results['esurf_std'] = float(parts[2])
        elif 'DELTA G gas' in line:
            parts = line.split()
            if len(parts) >= 5:
                results['g_gas'] = float(parts[3])
                results['g_gas_std'] = float(parts[4])
        elif 'DELTA G solv' in line:
            parts = line.split()
            if len(parts) >= 5:
                results['g_solv'] = float(parts[3])
                results['g_solv_std'] = float(parts[4])
        elif 'DELTA TOTAL' in line:
            parts = line.split()
            if len(parts) >= 4:
                results['total'] = float(parts[2])
                results['total_std'] = float(parts[3])
    
    if 'total' in results:
        return results
    else:
        return None


def check_mmgbsa_status(ligand_ids: List[str]) -> Dict[str, str]:
    """Check MM/GBSA completion status for each ligand."""
    status = {}
    
    for ligand_id in ligand_ids:
        system_dir = SYSTEMS_DIR / ligand_id
        result_file = system_dir / "FINAL_RESULTS_MMGBSA.dat"
        
        if result_file.exists():
            # Check if file is complete
            results = parse_mmgbsa_result(result_file)
            if results:
                status[ligand_id] = "Complete"
            else:
                status[ligand_id] = "Incomplete"
        else:
            # Check if job is running
            log_files = list(system_dir.glob("mmgbsa_*.log"))
            if log_files:
                status[ligand_id] = "Running"
            else:
                prod_nc = system_dir / "prod.nc"
                if prod_nc.exists():
                    status[ligand_id] = "Ready (not submitted)"
                else:
                    status[ligand_id] = "MD incomplete"
    
    return status


def is_type_ii_inhibitor(ligand_id: str) -> bool:
    """Check if ligand is Type II inhibitor."""
    return ligand_id in LIGAND_FE_COORDINATION


def main():
    """Main extraction function."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Extract MM/GBSA results")
    parser.add_argument("ligand_list", type=Path, help="File containing ligand IDs")
    parser.add_argument("--output", type=Path, default=None, 
                       help="Output CSV file (default: mmgbsa_results.csv)")
    parser.add_argument("--status-only", action="store_true",
                       help="Only check status, don't extract results")
    
    args = parser.parse_args()
    
    if not args.ligand_list.exists():
        print(f"ERROR: Ligand list file not found: {args.ligand_list}")
        sys.exit(1)
    
    # Read ligand list
    ligand_ids = []
    with open(args.ligand_list, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                ligand_ids.append(line)
    
    print("="*80)
    print("MM/GBSA Results Extraction")
    print(f"Ligand list: {args.ligand_list}")
    print(f"Total ligands: {len(ligand_ids)}")
    print("="*80)
    print()
    
    # Check status
    print("Checking MM/GBSA status...")
    status = check_mmgbsa_status(ligand_ids)
    
    # Count by status
    status_counts = {}
    for s in status.values():
        status_counts[s] = status_counts.get(s, 0) + 1
    
    print()
    print("Status Summary:")
    print("-" * 40)
    for s, count in sorted(status_counts.items()):
        print(f"  {s:20s}: {count:3d}")
    print("-" * 40)
    print()
    
    if args.status_only:
        # Print detailed status
        print("Detailed Status:")
        print("-" * 60)
        for ligand_id in ligand_ids:
            type_ii = "Type II" if is_type_ii_inhibitor(ligand_id) else "Type I"
            print(f"  {ligand_id:20s} {type_ii:8s} {status[ligand_id]}")
        print("-" * 60)
        return
    
    # Extract results
    print("Extracting binding energies...")
    results = []
    
    for ligand_id in ligand_ids:
        if status[ligand_id] != "Complete":
            continue
        
        system_dir = SYSTEMS_DIR / ligand_id
        result_file = system_dir / "FINAL_RESULTS_MMGBSA.dat"
        
        data = parse_mmgbsa_result(result_file)
        if data:
            data['ligand_id'] = ligand_id
            data['type'] = "Type II" if is_type_ii_inhibitor(ligand_id) else "Type I"
            results.append(data)
    
    if not results:
        print("No complete results found!")
        return
    
    print(f"Extracted {len(results)} results")
    print()
    
    # Sort by binding affinity (most negative first)
    results.sort(key=lambda x: x['total'])
    
    # Print summary table
    print("="*80)
    print("Binding Free Energy Results (sorted by affinity)")
    print("="*80)
    print(f"{'Rank':<6} {'Ligand ID':<20} {'Type':<10} {'ΔG (kcal/mol)':<15} {'Std Dev'}")
    print("-"*80)
    
    for i, data in enumerate(results, 1):
        print(f"{i:<6} {data['ligand_id']:<20} {data['type']:<10} "
              f"{data['total']:>8.2f} ± {data['total_std']:5.2f}")
    
    print("="*80)
    print()
    
    # Print top 5
    print("Top 5 Binders:")
    print("-"*60)
    for i, data in enumerate(results[:5], 1):
        print(f"{i}. {data['ligand_id']:<20} {data['type']:<10} "
              f"ΔG = {data['total']:>8.2f} ± {data['total_std']:5.2f} kcal/mol")
    print()
    
    # Statistics by type
    type_ii_results = [r for r in results if r['type'] == "Type II"]
    type_i_results = [r for r in results if r['type'] == "Type I"]
    
    if type_ii_results and type_i_results:
        avg_type_ii = sum(r['total'] for r in type_ii_results) / len(type_ii_results)
        avg_type_i = sum(r['total'] for r in type_i_results) / len(type_i_results)
        
        print("Statistics by Inhibitor Type:")
        print("-"*60)
        print(f"Type II (n={len(type_ii_results)}): Average ΔG = {avg_type_ii:.2f} kcal/mol")
        print(f"Type I  (n={len(type_i_results)}): Average ΔG = {avg_type_i:.2f} kcal/mol")
        print()
    
    # Save to CSV
    if args.output is None:
        args.output = Path("mmgbsa_results.csv")
    
    with open(args.output, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Header
        writer.writerow([
            "Rank", "Ligand_ID", "Type",
            "Delta_G_Total", "Total_Std",
            "Delta_G_Gas", "Gas_Std",
            "Delta_G_Solv", "Solv_Std",
            "VDW", "VDW_Std",
            "EEL", "EEL_Std",
            "EGB", "EGB_Std",
            "ESURF", "ESURF_Std"
        ])
        
        # Data
        for i, data in enumerate(results, 1):
            writer.writerow([
                i,
                data['ligand_id'],
                data['type'],
                f"{data['total']:.2f}",
                f"{data['total_std']:.2f}",
                f"{data.get('g_gas', 0):.2f}",
                f"{data.get('g_gas_std', 0):.2f}",
                f"{data.get('g_solv', 0):.2f}",
                f"{data.get('g_solv_std', 0):.2f}",
                f"{data.get('vdw', 0):.2f}",
                f"{data.get('vdw_std', 0):.2f}",
                f"{data.get('eel', 0):.2f}",
                f"{data.get('eel_std', 0):.2f}",
                f"{data.get('egb', 0):.2f}",
                f"{data.get('egb_std', 0):.2f}",
                f"{data.get('esurf', 0):.2f}",
                f"{data.get('esurf_std', 0):.2f}"
            ])
    
    print(f"Results saved to: {args.output}")
    print()


if __name__ == "__main__":
    main()

