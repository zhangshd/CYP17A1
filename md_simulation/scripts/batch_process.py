"""
Batch Processing Pipeline for CYP17A1 Virtual Screening
Process multiple ligands through the complete MD simulation workflow.
Author: zhangshd
Date: 2024-12-17
"""

import os
import sys
import argparse
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
from enum import Enum

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent))

from config import (
    LIGAND_DATABASES, SYSTEMS_DIR, RESULTS_DIR, LIGAND_PARAMS_DIR,
    AMBER_HOME, ensure_directories
)
from step1_prepare_ligand_params import parameterize_ligand
from step2_build_complex_topology import build_complex_topology
from step3_setup_md_simulation import setup_md_simulation
from step4_run_mmgbsa import setup_mmgbsa_analysis, parse_mmgbsa_result


class Status(Enum):
    """Processing status for each step."""
    NOT_STARTED = "not_started"
    SUCCESS = "success"
    FAILED = "failed"
    SKIPPED = "skipped"


@dataclass
class LigandResult:
    """Result for a single ligand."""
    ligand_id: str
    source: str
    params_status: Status = Status.NOT_STARTED
    topology_status: Status = Status.NOT_STARTED
    md_setup_status: Status = Status.NOT_STARTED
    mmgbsa_status: Status = Status.NOT_STARTED
    delta_g: Optional[float] = None
    error_message: Optional[str] = None


def extract_ligand_from_database(
    db_name: str,
    ligand_index: int,
    output_path: Path
) -> bool:
    """
    Extract a single ligand from a multi-mol2 database file.
    
    Args:
        db_name: Database name (key in LIGAND_DATABASES)
        ligand_index: 1-based index of the ligand
        output_path: Where to save the extracted mol2
        
    Returns:
        True if successful
    """
    db_config = LIGAND_DATABASES.get(db_name)
    if not db_config:
        print(f"  ERROR: Unknown database: {db_name}")
        return False
    
    mol2_file = db_config["mol2_file"]
    if not mol2_file.exists():
        print(f"  ERROR: Database file not found: {mol2_file}")
        return False
    
    # Parse multi-mol2 file
    molecules = []
    current_mol = []
    in_molecule = False
    
    with open(mol2_file, 'r') as f:
        for line in f:
            if line.startswith("@<TRIPOS>MOLECULE"):
                if current_mol:
                    molecules.append(''.join(current_mol))
                current_mol = [line]
                in_molecule = True
            elif in_molecule:
                current_mol.append(line)
    
    # Don't forget the last molecule
    if current_mol:
        molecules.append(''.join(current_mol))
    
    if ligand_index < 1 or ligand_index > len(molecules):
        print(f"  ERROR: Invalid ligand index {ligand_index} (database has {len(molecules)} molecules)")
        return False
    
    # Write the selected molecule
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(molecules[ligand_index - 1])
    
    return True


def get_docked_mol2_path(db_name: str, mol_index: int) -> Optional[Path]:
    """
    Get path to the docked mol2 file for a specific molecule.
    
    Args:
        db_name: Database name (Ch_medicinal_DB, feed_additive, gras_new)
        mol_index: 1-based molecule index
        
    Returns:
        Path to docked mol2 file, or None if not found
    """
    db_config = LIGAND_DATABASES.get(db_name)
    if not db_config:
        return None
    
    docking_dir = db_config["docking_dir"]
    
    # Construct directory name based on database naming convention
    # Format: <db_prefix>_ligprep.sdf_<index>
    db_to_prefix = {
        "Ch_medicinal_DB": "Ch_medicinal_DB_ligprep.sdf",
        "feed_additive": "feed_additive_ligprep.sdf",
        "gras_new": "gras_new_ligprep.sdf",
    }
    
    prefix = db_to_prefix.get(db_name)
    if not prefix:
        return None
    
    mol_dir = docking_dir / f"{prefix}_{mol_index}"
    docked_mol2 = mol_dir / "GD2_HEME_fb.mol2"
    
    if docked_mol2.exists():
        return docked_mol2
    
    return None


def load_top_molecules_from_summary(db_name: str, top_n: int = 10) -> List[int]:
    """
    Load top N molecule indices from database summary CSV.
    
    The CSV is already sorted by docking score (best first).
    Molecule column format: "gras_new_ligprep.sdf_1015" -> extract 1015
    
    Args:
        db_name: Database name
        top_n: Number of top molecules to return
        
    Returns:
        List of molecule indices (1-based), deduplicated
    """
    import csv
    import re
    
    db_config = LIGAND_DATABASES.get(db_name)
    if not db_config:
        print(f"  ERROR: Unknown database: {db_name}")
        return []
    
    docking_dir = db_config["docking_dir"]
    
    # Find summary CSV file
    summary_pattern = list(docking_dir.glob("*_summary.csv"))
    if not summary_pattern:
        print(f"  ERROR: No summary CSV found in {docking_dir}")
        return []
    
    summary_file = summary_pattern[0]
    
    # Parse CSV and extract top N molecule indices (deduplicated)
    molecules = []
    seen = set()
    
    with open(summary_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Get the Molecule column (format: "xxx_ligprep.sdf_INDEX")
            mol_name = row.get("Molecule", "")
            
            # Extract index from molecule name using regex
            match = re.search(r'_(\d+)$', mol_name)
            if match:
                mol_idx = int(match.group(1))
                if mol_idx not in seen:
                    seen.add(mol_idx)
                    molecules.append(mol_idx)
                    if len(molecules) >= top_n:
                        break
    
    return molecules


def prepare_ligand_with_hydrogens(
    docked_mol2: Path,
    output_dir: Path
) -> Optional[Path]:
    """
    Prepare ligand mol2 file for parameterization using PyMOL.
    
    Uses PyMOL's h_add command to add hydrogens, which is more robust
    than OpenBabel for handling incomplete bond information from GalaxyDock2.
    
    Args:
        docked_mol2: Path to docked mol2
        output_dir: Directory to save prepared file
        
    Returns:
        Path to prepared mol2 with hydrogens, or None on failure
    """
    import subprocess
    import tempfile
    
    output_dir.mkdir(parents=True, exist_ok=True)
    output_mol2 = output_dir / "ligand_with_H.mol2"
    
    # Create PyMOL script for adding hydrogens
    pymol_script = f"""
# PyMOL script to add hydrogens to ligand
from pymol import cmd

# Load the docked mol2 (first state only)
cmd.load("{docked_mol2}", "ligand", state=1)

# Remove any existing hydrogens first to ensure clean addition
cmd.remove("hydro")

# Add hydrogens using PyMOL's h_add
cmd.h_add("ligand")

# Save as mol2
cmd.save("{output_mol2}", "ligand", state=1)

# Quit
cmd.quit()
"""
    
    # Write script to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as f:
        f.write(pymol_script)
        script_path = f.name
    
    # Run PyMOL in command-line mode
    cmd = f"pymol -cq {script_path}"
    
    result = subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
        text=True,
        executable="/bin/bash"
    )
    
    # Clean up script
    if Path(script_path).exists():
        Path(script_path).unlink()
    
    if not output_mol2.exists():
        print(f"  ERROR: PyMOL failed to create output file")
        if result.stderr:
            print(f"  STDERR: {result.stderr[:500]}")
        return None
    
    return output_mol2


def process_ligand(
    ligand_id: str,
    input_mol2: Path,
    docked_mol2: Optional[Path] = None,
    force: bool = False
) -> LigandResult:
    """
    Process a single ligand through the complete pipeline.
    
    Args:
        ligand_id: Unique identifier for the ligand
        input_mol2: Path to input mol2 (may or may not have hydrogens)
        docked_mol2: Path to docked mol2 (with binding pose coordinates)
        force: If True, overwrite existing files
        
    Returns:
        LigandResult with status for each step
    """
    result = LigandResult(ligand_id=ligand_id, source="batch")
    
    print(f"\n{'='*70}")
    print(f"Processing: {ligand_id}")
    print(f"{'='*70}")
    
    # Step 0: Prepare ligand (add hydrogens if needed)
    print("Step 0: Preparing ligand with hydrogens...")
    prep_dir = LIGAND_PARAMS_DIR / ligand_id
    prepared_mol2 = prepare_ligand_with_hydrogens(input_mol2, prep_dir)
    
    if prepared_mol2 is None:
        result.params_status = Status.FAILED
        result.error_message = "Failed to add hydrogens"
        print("  FAILED: Could not add hydrogens to ligand")
        return result
    print(f"  OK: {prepared_mol2.name}")
    
    # Use prepared mol2 for parameterization, keep original coords for docking pose
    coord_mol2 = docked_mol2 if docked_mol2 else prepared_mol2
    
    # Step 1: Parameterize ligand
    print("Step 1: Preparing GAFF2 parameters...")
    success, mol2_path, frcmod_path = parameterize_ligand(
        prepared_mol2, ligand_id, force=force
    )
    result.params_status = Status.SUCCESS if success else Status.FAILED
    
    if result.params_status == Status.FAILED:
        result.error_message = "Parameterization failed"
        return result
    print("  OK")
    
    # Step 2: Build topology
    print("Step 2: Building topology with tleap...")
    success, system_dir = build_complex_topology(ligand_id, prepared_mol2, force=force)
    result.topology_status = Status.SUCCESS if success else Status.FAILED
    
    if result.topology_status == Status.FAILED:
        result.error_message = "Topology build failed"
        return result
    print("  OK")
    
    # Step 3: Setup MD simulation
    print("Step 3: Setting up MD simulation...")
    success = setup_md_simulation(ligand_id, force=force)
    result.md_setup_status = Status.SUCCESS if success else Status.FAILED
    
    if result.md_setup_status == Status.FAILED:
        result.error_message = "MD setup failed"
        return result
    print("  OK")
    
    return result


def process_batch(
    ligands: List[Tuple[str, Path, Optional[Path]]],
    force: bool = False
) -> List[LigandResult]:
    """
    Process a batch of ligands.
    
    Args:
        ligands: List of (ligand_id, input_mol2, docked_mol2) tuples
        force: If True, overwrite existing files
        
    Returns:
        List of LigandResult
    """
    results = []
    
    for ligand_id, input_mol2, docked_mol2 in ligands:
        result = process_ligand(ligand_id, input_mol2, docked_mol2, force)
        results.append(result)
    
    return results


def print_summary(results: List[LigandResult]):
    """Print summary of batch processing."""
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    success_count = sum(1 for r in results if r.md_setup_status == Status.SUCCESS)
    print(f"Successful: {success_count}/{len(results)}")
    
    for r in results:
        if r.md_setup_status == Status.SUCCESS:
            print(f"  ✓ {r.ligand_id}")
        else:
            failed_step = "params" if r.params_status == Status.FAILED else \
                          "topology" if r.topology_status == Status.FAILED else \
                          "md_setup"
            print(f"  ✗ {r.ligand_id}: {failed_step}")


def submit_md_jobs(results: List[LigandResult], dry_run: bool = True):
    """
    Submit MD jobs for successfully processed ligands.
    
    Args:
        results: List of LigandResult
        dry_run: If True, only print commands without executing
    """
    print("\n" + "=" * 70)
    print("SUBMITTING MD JOBS")
    print("=" * 70)
    
    for r in results:
        if r.md_setup_status == Status.SUCCESS:
            system_dir = SYSTEMS_DIR / r.ligand_id
            if dry_run:
                print(f"  [DRY RUN] cd {system_dir} && sbatch run_md.sh")
            else:
                import subprocess
                result = subprocess.run(
                    "sbatch run_md.sh",
                    shell=True,
                    cwd=system_dir,
                    capture_output=True,
                    text=True
                )
                if result.returncode == 0:
                    print(f"  ✓ {r.ligand_id}: {result.stdout.strip()}")
                else:
                    print(f"  ✗ {r.ligand_id}: {result.stderr}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Batch process ligands for CYP17A1 MD simulation"
    )
    parser.add_argument(
        "--top-n", "-n",
        type=int,
        default=10,
        help="Number of top molecules to process from each database (default: 10)"
    )
    parser.add_argument(
        "--database", "-d",
        type=str,
        choices=list(LIGAND_DATABASES.keys()) + ["all"],
        default="all",
        help="Source database for ligands (default: all)"
    )
    parser.add_argument(
        "--indices", "-i",
        type=str,
        help="Comma-separated molecule indices (requires --database)"
    )
    parser.add_argument(
        "--force", "-f",
        action="store_true",
        help="Overwrite existing files"
    )
    parser.add_argument(
        "--submit", "-s",
        action="store_true",
        help="Submit MD jobs after setup (default: dry run)"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only print what would be done without executing"
    )
    
    args = parser.parse_args()
    
    # Set up AMBER environment
    os.environ["AMBERHOME"] = AMBER_HOME
    os.environ["PATH"] = f"{AMBER_HOME}/bin:" + os.environ.get("PATH", "")
    
    ensure_directories()
    
    # Determine which databases to process
    if args.database == "all":
        databases = list(LIGAND_DATABASES.keys())
    else:
        databases = [args.database]
    
    # Build ligand list: (ligand_id, db_name, mol_index, docked_mol2_path)
    ligands_to_process: List[Tuple[str, str, int, Path]] = []
    
    for db_name in databases:
        db_config = LIGAND_DATABASES[db_name]
        prefix = db_config["prefix"]
        
        if args.indices:
            # Use specified indices
            mol_indices = [int(x.strip()) for x in args.indices.split(",")]
        else:
            # Load from summary CSV
            print(f"\nLoading top-{args.top_n} molecules from {db_name}...")
            mol_indices = load_top_molecules_from_summary(db_name, args.top_n)
        
        for mol_idx in mol_indices:
            ligand_id = f"{prefix}_{mol_idx}"
            docked_mol2 = get_docked_mol2_path(db_name, mol_idx)
            
            if docked_mol2 is None:
                print(f"  WARNING: Docked mol2 not found for {ligand_id}, skipping")
                continue
            
            ligands_to_process.append((ligand_id, db_name, mol_idx, docked_mol2))
            print(f"  Found: {ligand_id} -> {docked_mol2.name}")
    
    print(f"\n{'='*70}")
    print(f"Total ligands to process: {len(ligands_to_process)}")
    print(f"{'='*70}")
    
    if args.dry_run:
        print("\n[DRY RUN] Would process the following ligands:")
        for ligand_id, db_name, mol_idx, docked_mol2 in ligands_to_process:
            print(f"  - {ligand_id} from {db_name}")
        return
    
    # Process each ligand
    results = []
    for ligand_id, db_name, mol_idx, docked_mol2 in ligands_to_process:
        result = process_ligand(
            ligand_id=ligand_id,
            input_mol2=docked_mol2,  # Use docked pose as input
            docked_mol2=docked_mol2,
            force=args.force
        )
        result.source = db_name
        results.append(result)
    
    # Print summary
    print_summary(results)
    
    # Optionally submit MD jobs
    if args.submit:
        submit_md_jobs(results, dry_run=False)
    else:
        submit_md_jobs(results, dry_run=True)


if __name__ == "__main__":
    main()
