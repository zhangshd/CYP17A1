"""
Configuration file for CYP17A1 MD Simulation Pipeline
Contains all paths, parameters, and settings.
Author: zhangshd
Date: 2024-12-17
"""

from pathlib import Path
from typing import Dict, Any
import os

# ============================================================================
# Directory Structure
# ============================================================================

# Base directories (relative paths)
PROJECT_ROOT = Path(os.path.abspath(__file__)).parent.parent.parent
MD_BASE = PROJECT_ROOT / "md_simulation"

# Input data directories
DATA_DIR = PROJECT_ROOT / "data"
DOCKING_RESULTS_DIR = PROJECT_ROOT / "docking_results"
LIGAND_DB_DIR = DATA_DIR / "mols"

# MCPB.py parameters (Fe-Heme)
MCPB_DIR = MD_BASE / "mcpb_heme"
MCPB_PRMTOP = MCPB_DIR / "pig_CYP17A1_dry.prmtop"
MCPB_INPCRD = MCPB_DIR / "pig_CYP17A1_dry.inpcrd"
MCPB_PDB = MCPB_DIR / "pig_CYP17A1_mcpbpy.pdb"

# MCPB.py templates for tleap (mol2 files in same directory)
MCPB_TEMPLATES_DIR = MCPB_DIR  # CM1.mol2, HM1.mol2, FE1.mol2 are in mcpb_heme/
MCPB_FRCMOD = MCPB_DIR / "pig_CYP17A1_mcpbpy.frcmod"
MCPB_LIB = None  # Not using lib file, loading mol2 directly

# Output directories
SCRIPTS_DIR = MD_BASE / "scripts"
COMPLEXES_DIR = MD_BASE / "complexes"  # Base directory for MD systems
LIGAND_PARAMS_DIR = COMPLEXES_DIR / "ligand_params"  # Under complexes/
SYSTEMS_DIR = COMPLEXES_DIR / "systems_round2"  # GPU-based systems (main working directory)
RESULTS_DIR = MD_BASE / "results"/"round2"

# ============================================================================
# AMBER Environment
# ============================================================================
AMBER_HOME = "/opt/share/Amber/amber22"
AMBER_SOURCE_CMD = f"source {AMBER_HOME}/amber.sh"

# SLURM configuration (GPU-based execution)
SLURM_CONFIG = {
    "partition": "GA30",
    "nodes": 1,
    "ntasks": 1,
    "cpus_per_task": 16,
    "gres": "gpu:1",
    "time_md": "72:00:00",
    "time_mmgbsa": "12:00:00",
    "modules": ["cuda/11.8", "fftw/fftw-3.3.8"],
}

# ============================================================================
# Force Field Parameters
# ============================================================================

# Ligand force field
LIGAND_FF = "gaff2"
CHARGE_MODEL = "bcc"  # AM1-BCC

# Water model
WATER_MODEL = "OPC"
BUFFER_SIZE = 12.0  # Angstrom

# ============================================================================
# System Definition
# ============================================================================

# Standard residue name for ligands
LIGAND_RESIDUE_NAME = "LIG"

# Residue numbering (for MM/GBSA masks)
# Protein: 1-473, Heme (CM1/HM1): 474, Fe (FE1): 475, Ligand: 476
PROTEIN_RESIDUES = "1-473"
HEME_RESIDUE = 474
FE_RESIDUE = 475
LIGAND_RESIDUE = 476

# Masks for MM/GBSA
RECEPTOR_MASK = ":1-475"  # Protein + Heme + Fe
LIGAND_MASK = ":476"
STRIP_MASK = ":WAT,Na+,Cl-,Na,K+,Mg2+"

# ============================================================================
# MCPB.py Bond Definitions
# ============================================================================

# Residue numbering from MCPB.py output
# CM1 (coordinating Cys442): 412
# HM1 (Heme): 474
# FE1 (Iron): 475
# Ligand: 476 (appended after protein)

# Custom atom types for Fe-Heme (required for tleap)
CUSTOM_ATOM_TYPES = """addAtomTypes {
        { "M1"  "Fe" "sp3" }
        { "Y1"  "S" "sp3" }
        { "Y2"  "N" "sp3" }
        { "Y3"  "N" "sp3" }
        { "Y4"  "N" "sp3" }
        { "Y5"  "N" "sp3" }
}"""

# Bond commands for tleap (Fe-Heme-Cys coordination)
MCPB_BONDS = [
    ("412.SG", "475.FE"),   # CM1(Cys)-FE
    ("474.NA", "475.FE"),   # HM1-FE
    ("474.NB", "475.FE"),
    ("474.NC", "475.FE"),
    ("474.ND", "475.FE"),
    ("411.C", "412.N"),     # Peptide bond before CM1
    ("412.C", "413.N"),     # Peptide bond after CM1
]

# ============================================================================
# Type II Inhibitor Configuration (Fe-Ligand Coordination)
# ============================================================================

# For Type II inhibitors that form coordination bonds with Fe
# We use distance restraints instead of covalent bonds (no Fe-ligand torsion parameters)
# Format: {ligand_id: {"atom": atom_name, "target_dist": distance_in_angstrom}}
# The restraint uses flat-bottom potential with strong force constant
LIGAND_FE_COORDINATION = {
    "AER601": {
        "atom": "N22",          # Coordinating nitrogen atom
        "target_dist": 2.1,     # Target Fe-N22 distance (Å)
        "force_const": 50.0,    # Force constant (kcal/mol/Å²) - strong restraint
        "tolerance": 0.3,       # Flat-bottom width (±Å)
    },
    # Add more Type II inhibitors here as needed:
    # "ligand_xxx": {"atom": "N1", "target_dist": 2.1, "force_const": 50.0, "tolerance": 0.3},
}

# Distance restraint parameters for Type II inhibitors
# Used during minimization and equilibration, gradually released during production
FE_COORDINATION_RESTRAINT = {
    "min1": 50.0,      # Strong restraint during minimization
    "min2": 50.0,
    "heat": 50.0,      # Maintain during heating
    "eq1": 30.0,       # Gradually reduce during equilibration
    "eq2": 10.0,
    "prod": 0.0,       # No restraint in production (let it equilibrate naturally)
}

# ============================================================================
# MD Simulation Protocol
# ============================================================================

MD_PROTOCOL = {
    "min1": {
        "description": "Minimization with heavy atom restraints",
        "maxcyc": 10000,
        "ncyc": 5000,
        "restraint_wt": 10.0,
        "restraintmask": "!@H=",
    },
    "min2": {
        "description": "Minimization without restraints",
        "maxcyc": 10000,
        "ncyc": 5000,
    },
    "heat": {
        "description": "Heating from 0K to 300K (NVT, 50ps)",
        "nstlim": 25000,  # 50 ps
        "dt": 0.002,
        "tempi": 0.0,
        "temp0": 300.0,
        "restraint_wt": 5.0,
        "restraintmask": "@CA,C,N,O",
    },
    "eq1": {
        "description": "NPT equilibration with restraints (100ps)",
        "nstlim": 50000,  # 100 ps
        "dt": 0.002,
        "temp0": 300.0,
        "restraint_wt": 2.0,
        "restraintmask": "@CA,C,N,O",
    },
    "eq2": {
        "description": "NPT equilibration without restraints (10 stages x 50ps)",
        "num_stages": 10,
        "nstlim_per_stage": 25000,  # 50 ps each
        "dt": 0.002,
        "temp0": 300.0,
    },
    "prod": {
        "description": "Production NPT (20 ns)",
        "nstlim": 10000000,  # 20 ns
        "dt": 0.002,
        "temp0": 300.0,
        "ntpr": 5000,  # 10 ps
        "ntwx": 5000,  # 10 ps
    },
}

# ============================================================================
# MM/GBSA Parameters
# ============================================================================

MMGBSA_PARAMS = {
    "igb": 5,           # GB-OBC2 model
    "saltcon": 0.15,    # 150 mM salt
    "startframe": 1,
    "endframe": 500,    # Last 500 frames (5 ns at 10 ps interval)
    "interval": 1,
}

# ============================================================================
# Ligand Source Databases
# ============================================================================

LIGAND_DATABASES: Dict[str, Dict[str, Any]] = {
    "Ch_medicinal_DB": {
        "mol2_file": LIGAND_DB_DIR / "Ch_medicinal_DB_ligprep-out.mol2",
        "docking_dir": DOCKING_RESULTS_DIR / "Ch_medicinal_DB_ligprep-out",
        "prefix": "ChMed",
    },
    "feed_additive": {
        "mol2_file": LIGAND_DB_DIR / "feed_additive_ligprep-out.mol2",
        "docking_dir": DOCKING_RESULTS_DIR / "feed_additive_ligprep-out",
        "prefix": "FeedAdd",
    },
    "gras_new": {
        "mol2_file": LIGAND_DB_DIR / "gras_new_ligprep-out.mol2",
        "docking_dir": DOCKING_RESULTS_DIR / "gras_new_ligprep-out",
        "prefix": "GRAS",
    },
}

# ============================================================================
# Utility Functions
# ============================================================================

def get_system_dir(ligand_id: str) -> Path:
    """Get the system directory path for a ligand."""
    return SYSTEMS_DIR / ligand_id


def get_ligand_param_dir(ligand_id: str) -> Path:
    """Get the ligand parameter directory path."""
    return LIGAND_PARAMS_DIR / ligand_id


def ensure_directories():
    """Create all necessary directories if they don't exist."""
    for d in [LIGAND_PARAMS_DIR, SYSTEMS_DIR, RESULTS_DIR, SCRIPTS_DIR]:
        d.mkdir(parents=True, exist_ok=True)
