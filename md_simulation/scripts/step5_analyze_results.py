"""
Step 5: Analyze MD Results and Collect MM/GBSA Data
Validates MD trajectories and summarizes binding free energy results.
Author: zhangshd
Date: 2024-12-19
"""

import os
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
import csv

from config import (
    SYSTEMS_DIR, RESULTS_DIR, AMBER_HOME,
    ensure_directories
)


@dataclass
class SystemAnalysis:
    """Analysis results for a single system."""
    ligand_id: str
    source: str
    prod_frames: int
    expected_frames: int
    fe_n_distance: Optional[float]  # Fe-N coordination distance (Type II inhibitors)
    centroid_distance: Optional[float]  # Ligand centroid to Fe distance
    has_nitrogen: bool  # Whether ligand contains nitrogen
    delta_g: Optional[float]
    delta_g_std: Optional[float]
    delta_g_sem: Optional[float]
    is_valid: bool
    validation_notes: str


# Distance threshold: Fe-N coordination distance should be ~2.1-2.3 Å for Type II inhibitors
# Using 3.5 Å as threshold to allow for thermal fluctuations
# If ligand has no nitrogen, use centroid distance with 8 Å threshold
FE_N_DISTANCE_THRESHOLD = 3.5  # Angstroms - for Fe-N coordination
LIGAND_CENTROID_DISTANCE_THRESHOLD = 8.0  # Angstroms - fallback for non-N ligands
EXPECTED_PROD_FRAMES = 5000


def get_prod_nc_frames(system_dir: Path) -> int:
    """Get number of frames in prod.nc using ncdump."""
    prod_nc = system_dir / "prod.nc"
    if not prod_nc.exists():
        return 0
    
    result = subprocess.run(
        ["ncdump", "-h", str(prod_nc)],
        capture_output=True, text=True
    )
    
    for line in result.stdout.split('\n'):
        if 'frame = ' in line:
            match = re.search(r'\d+', line)
            if match:
                return int(match.group())
    return 0


def calculate_ligand_fe_distance(system_dir: Path, sample_frames: int = 100) -> Tuple[Optional[float], Optional[float], bool]:
    """
    Calculate average ligand-Fe distance using cpptraj.
    Prioritizes Fe-N distance for Type II inhibitors.
    
    Args:
        system_dir: Path to system directory
        sample_frames: Number of frames to sample (evenly spaced)
        
    Returns:
        Tuple of (fe_n_distance, centroid_distance, has_nitrogen)
        - fe_n_distance: Closest N atom to Fe distance (None if no N)
        - centroid_distance: Ligand centroid to Fe distance
        - has_nitrogen: Whether ligand contains nitrogen atoms
    """
    prmtop = system_dir / "complex_solv.prmtop"
    prod_nc = system_dir / "prod.nc"
    
    if not prmtop.exists() or not prod_nc.exists():
        return None, None, False
    
    # Create cpptraj input - calculate both Fe-N and Fe-centroid distances
    cpptraj_in = system_dir / "check_dist.cpptraj"
    dist_n_out = system_dir / "fe_n_dist_check.dat"
    dist_centroid_out = system_dir / "fe_centroid_dist_check.dat"
    
    # Calculate stride to get ~sample_frames
    total_frames = get_prod_nc_frames(system_dir)
    if total_frames == 0:
        return None, None, False
    
    stride = max(1, total_frames // sample_frames)
    
    # :476@N* selects all nitrogen atoms in ligand (residue 476)
    # distance command with closest will find the minimum distance
    cpptraj_content = f"""parm {prmtop}
trajin {prod_nc} 1 last {stride}
# Fe-N distance (closest nitrogen to Fe)
distance fe_n :476@N* :475@FE out {dist_n_out}
# Fe-centroid distance (fallback)
distance fe_centroid :476 :475@FE out {dist_centroid_out}
go
"""
    cpptraj_in.write_text(cpptraj_content)
    
    # Run cpptraj
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = f"/opt/share/fftw/lib:{env.get('LD_LIBRARY_PATH', '')}"
    
    result = subprocess.run(
        [f"{AMBER_HOME}/bin/cpptraj", "-i", str(cpptraj_in)],
        capture_output=True, text=True, env=env
    )
    
    if result.returncode != 0:
        return None, None, False
    
    def parse_distance_file(filepath: Path) -> Optional[float]:
        """Parse cpptraj distance output and return average."""
        if not filepath.exists():
            return None
        distances = []
        with open(filepath) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        distances.append(float(parts[1]))
                    except ValueError:
                        continue
        if not distances:
            return None
        return sum(distances) / len(distances)
    
    # Parse Fe-N distance (may not exist if ligand has no nitrogen)
    fe_n_dist = parse_distance_file(dist_n_out)
    has_nitrogen = fe_n_dist is not None
    
    # Parse centroid distance (should always exist)
    centroid_dist = parse_distance_file(dist_centroid_out)
    
    return fe_n_dist, centroid_dist, has_nitrogen


def parse_mmgbsa_result(system_dir: Path) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """
    Parse FINAL_RESULTS_MMGBSA.dat to extract binding energy.
    
    Returns:
        Tuple of (delta_G, std_dev, std_err) or (None, None, None) if not found
    """
    result_file = system_dir / "FINAL_RESULTS_MMGBSA.dat"
    
    if not result_file.exists():
        return None, None, None
    
    with open(result_file) as f:
        for line in f:
            if "DELTA TOTAL" in line:
                parts = line.split()
                if len(parts) >= 5:
                    return float(parts[2]), float(parts[3]), float(parts[4])
    
    return None, None, None


def identify_source(ligand_id: str) -> str:
    """Identify the source database of a ligand."""
    if ligand_id == "AER601":
        return "Reference"
    elif ligand_id.startswith("ChMed_"):
        return "Ch_medicinal_DB"
    elif ligand_id.startswith("FeedAdd_"):
        return "feed_additive"
    elif ligand_id.startswith("GRAS_"):
        return "gras_new"
    elif ligand_id.startswith("ligand_"):
        return "gras_new"  # Legacy naming
    else:
        return "unknown"


def analyze_system(ligand_id: str) -> Optional[SystemAnalysis]:
    """Analyze a single system and return results."""
    system_dir = SYSTEMS_DIR / ligand_id
    
    if not system_dir.exists():
        return None
    
    # Check trajectory completeness
    prod_frames = get_prod_nc_frames(system_dir)
    
    # Calculate ligand-Fe distances (Fe-N for Type II, centroid as fallback)
    fe_n_dist, centroid_dist, has_nitrogen = calculate_ligand_fe_distance(system_dir)
    
    # Parse MMGBSA results
    delta_g, delta_g_std, delta_g_sem = parse_mmgbsa_result(system_dir)
    
    # Validation
    validation_notes = []
    is_valid = True
    
    if prod_frames < EXPECTED_PROD_FRAMES:
        validation_notes.append(f"Incomplete trajectory ({prod_frames}/{EXPECTED_PROD_FRAMES} frames)")
        is_valid = False
    
    # Distance validation based on ligand type
    # For Type II inhibitors (with nitrogen): check Fe-N coordination distance
    # For other ligands: check centroid distance
    if has_nitrogen and fe_n_dist is not None:
        if fe_n_dist > FE_N_DISTANCE_THRESHOLD:
            validation_notes.append(
                f"Fe-N coordination lost (dist={fe_n_dist:.2f} Å > {FE_N_DISTANCE_THRESHOLD} Å, expected 2.1-2.3 Å)"
            )
            is_valid = False
    elif centroid_dist is not None:
        if centroid_dist > LIGAND_CENTROID_DISTANCE_THRESHOLD:
            validation_notes.append(
                f"Ligand dissociated (centroid dist={centroid_dist:.1f} Å > {LIGAND_CENTROID_DISTANCE_THRESHOLD} Å)"
            )
            is_valid = False
    else:
        validation_notes.append("Could not calculate ligand-Fe distance")
        is_valid = False
    
    if delta_g is None:
        validation_notes.append("No MMGBSA result")
        is_valid = False
    
    return SystemAnalysis(
        ligand_id=ligand_id,
        source=identify_source(ligand_id),
        prod_frames=prod_frames,
        expected_frames=EXPECTED_PROD_FRAMES,
        fe_n_distance=fe_n_dist,
        centroid_distance=centroid_dist,
        has_nitrogen=has_nitrogen,
        delta_g=delta_g,
        delta_g_std=delta_g_std,
        delta_g_sem=delta_g_sem,
        is_valid=is_valid,
        validation_notes="; ".join(validation_notes) if validation_notes else "OK"
    )


def analyze_all_systems() -> List[SystemAnalysis]:
    """Analyze all systems in the SYSTEMS_DIR."""
    results = []
    
    if not SYSTEMS_DIR.exists():
        print(f"ERROR: Systems directory not found: {SYSTEMS_DIR}")
        return results
    
    systems = sorted([d.name for d in SYSTEMS_DIR.iterdir() if d.is_dir()])
    
    print(f"Analyzing {len(systems)} systems...")
    print("-" * 80)
    
    for i, ligand_id in enumerate(systems, 1):
        print(f"  [{i}/{len(systems)}] Analyzing {ligand_id}...", end=" ")
        analysis = analyze_system(ligand_id)
        if analysis:
            status = "✓ Valid" if analysis.is_valid else "✗ Invalid"
            print(status)
            results.append(analysis)
        else:
            print("✗ Failed")
    
    return results


def generate_summary_report(results: List[SystemAnalysis], output_dir: Path) -> Path:
    """Generate comprehensive summary report."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Generate detailed CSV
    csv_file = output_dir / "mmgbsa_analysis_full.csv"
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            "Ligand_ID", "Source", "Prod_Frames", "Expected_Frames",
            "Fe_N_Distance_A", "Centroid_Distance_A", "Has_Nitrogen",
            "Delta_G_kcal_mol", "Std_Dev", "Std_Err",
            "Is_Valid", "Validation_Notes"
        ])
        for r in results:
            # Use Fe-N distance for Type II inhibitors, centroid for others
            primary_dist = r.fe_n_distance if r.has_nitrogen else r.centroid_distance
            writer.writerow([
                r.ligand_id, r.source, r.prod_frames, r.expected_frames,
                f"{r.fe_n_distance:.2f}" if r.fe_n_distance is not None else "N/A",
                f"{r.centroid_distance:.2f}" if r.centroid_distance is not None else "N/A",
                r.has_nitrogen,
                f"{r.delta_g:.4f}" if r.delta_g else "N/A",
                f"{r.delta_g_std:.4f}" if r.delta_g_std else "N/A",
                f"{r.delta_g_sem:.4f}" if r.delta_g_sem else "N/A",
                r.is_valid, r.validation_notes
            ])
    
    # 2. Generate valid-only sorted CSV
    valid_results = [r for r in results if r.is_valid and r.delta_g is not None]
    valid_results.sort(key=lambda x: x.delta_g)
    
    valid_csv = output_dir / "mmgbsa_valid_ranked.csv"
    with open(valid_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            "Rank", "Ligand_ID", "Source", "Delta_G_kcal_mol", 
            "Std_Dev", "Fe_N_Distance_A", "Centroid_Distance_A", "Has_Nitrogen"
        ])
        for i, r in enumerate(valid_results, 1):
            writer.writerow([
                i, r.ligand_id, r.source,
                f"{r.delta_g:.2f}", f"{r.delta_g_std:.2f}",
                f"{r.fe_n_distance:.2f}" if r.fe_n_distance is not None else "N/A",
                f"{r.centroid_distance:.2f}" if r.centroid_distance is not None else "N/A",
                r.has_nitrogen
            ])
    
    # 3. Generate text report
    report_file = output_dir / "analysis_report.txt"
    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("CYP17A1 Virtual Screening - MM/GBSA Analysis Report\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Total systems analyzed: {len(results)}\n")
        f.write(f"Valid systems: {len(valid_results)}\n")
        f.write(f"Invalid systems: {len(results) - len(valid_results)}\n\n")
        
        f.write("-" * 80 + "\n")
        f.write("VALID RESULTS (Ranked by ΔG)\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'Rank':<6}{'Ligand ID':<20}{'Source':<18}{'ΔG (kcal/mol)':<15}{'Fe-N (Å)':<10}{'Centroid (Å)':<14}\n")
        f.write("-" * 80 + "\n")
        
        for i, r in enumerate(valid_results, 1):
            fe_n_str = f"{r.fe_n_distance:.2f}" if r.fe_n_distance is not None else "N/A"
            cent_str = f"{r.centroid_distance:.2f}" if r.centroid_distance is not None else "N/A"
            f.write(f"{i:<6}{r.ligand_id:<20}{r.source:<18}{r.delta_g:<15.2f}{fe_n_str:<10}{cent_str:<14}\n")
        
        f.write("\n")
        f.write("-" * 80 + "\n")
        f.write("INVALID RESULTS (Need Investigation)\n")
        f.write("-" * 80 + "\n")
        
        invalid_results = [r for r in results if not r.is_valid]
        for r in invalid_results:
            f.write(f"{r.ligand_id}: {r.validation_notes}\n")
        
        f.write("\n")
        f.write("=" * 80 + "\n")
        f.write("VALIDATION CRITERIA\n")
        f.write("=" * 80 + "\n")
        f.write(f"- Expected trajectory frames: {EXPECTED_PROD_FRAMES}\n")
        f.write(f"- Fe-N coordination distance threshold: {FE_N_DISTANCE_THRESHOLD} Å (for Type II inhibitors with N)\n")
        f.write(f"  (Expected Fe-N distance for stable coordination: 2.1-2.3 Å)\n")
        f.write(f"- Ligand centroid distance threshold: {LIGAND_CENTROID_DISTANCE_THRESHOLD} Å (fallback for non-N ligands)\n")
        f.write("- Distance exceeding threshold indicates ligand dissociation from active site\n")
        f.write("\n")
    
    return report_file


def print_summary(results: List[SystemAnalysis]):
    """Print summary to console."""
    valid_results = [r for r in results if r.is_valid and r.delta_g is not None]
    valid_results.sort(key=lambda x: x.delta_g)
    
    print("\n" + "=" * 80)
    print("ANALYSIS SUMMARY")
    print("=" * 80)
    print(f"\nTotal systems: {len(results)}")
    print(f"Valid systems: {len(valid_results)}")
    print(f"Invalid systems: {len(results) - len(valid_results)}")
    
    if valid_results:
        print("\n" + "-" * 80)
        print("TOP VALID CANDIDATES (by ΔG)")
        print("-" * 80)
        print(f"{'Rank':<6}{'Ligand ID':<20}{'Source':<18}{'ΔG (kcal/mol)':<15}{'Fe-N (Å)':<10}{'Centroid (Å)':<14}")
        print("-" * 80)
        
        for i, r in enumerate(valid_results[:10], 1):
            fe_n_str = f"{r.fe_n_distance:.2f}" if r.fe_n_distance is not None else "N/A"
            cent_str = f"{r.centroid_distance:.2f}" if r.centroid_distance is not None else "N/A"
            print(f"{i:<6}{r.ligand_id:<20}{r.source:<18}{r.delta_g:<15.2f}{fe_n_str:<10}{cent_str:<14}")
    
    invalid_results = [r for r in results if not r.is_valid]
    if invalid_results:
        print("\n" + "-" * 80)
        print("INVALID SYSTEMS (Summary)")
        print("-" * 80)
        for r in invalid_results[:5]:
            print(f"  {r.ligand_id}: {r.validation_notes}")
        if len(invalid_results) > 5:
            print(f"  ... and {len(invalid_results) - 5} more")


def main():
    """Main entry point."""
    import sys
    
    # Parse arguments
    if len(sys.argv) > 1 and sys.argv[1] in ['-h', '--help']:
        print("Usage: python step5_analyze_results.py [ligand_id ...]")
        print("\nAnalyze MD simulation results and validate MM/GBSA calculations.")
        print("\nOptions:")
        print("  No arguments    Analyze all systems in SYSTEMS_DIR")
        print("  ligand_id ...   Analyze specific systems")
        sys.exit(0)
    
    ensure_directories()
    
    # Analyze systems
    if len(sys.argv) > 1:
        # Analyze specific systems
        results = []
        for ligand_id in sys.argv[1:]:
            analysis = analyze_system(ligand_id)
            if analysis:
                results.append(analysis)
                status = "✓ Valid" if analysis.is_valid else "✗ Invalid"
                print(f"{ligand_id}: {status} - {analysis.validation_notes}")
    else:
        # Analyze all systems
        results = analyze_all_systems()
    
    if not results:
        print("No systems to analyze.")
        sys.exit(1)
    
    # Generate reports
    report_file = generate_summary_report(results, RESULTS_DIR)
    print(f"\nReports saved to: {RESULTS_DIR}")
    print(f"  - mmgbsa_analysis_full.csv (all results)")
    print(f"  - mmgbsa_valid_ranked.csv (valid results, ranked)")
    print(f"  - analysis_report.txt (text summary)")
    
    # Print summary
    print_summary(results)


if __name__ == "__main__":
    main()
