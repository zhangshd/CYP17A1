# MM/GBSA Analysis Tools

Quick reference for MM/GBSA batch processing and results extraction.

---

## Setup and Submit MM/GBSA Calculations

### Batch Processing Script: `batch_mmgbsa_submit.py`

Process multiple systems at once - checks MD completion, sets up MM/GBSA inputs, and submits jobs.

#### Basic Usage

```bash
# Process all ligands in list
python3 batch_mmgbsa_submit.py ligand_list.txt

# Check MD completion status only
python3 batch_mmgbsa_submit.py ligand_list.txt --check-only

# Setup but don't submit jobs
python3 batch_mmgbsa_submit.py ligand_list.txt --skip-submit

# Submit existing setup files
python3 batch_mmgbsa_submit.py ligand_list.txt --skip-setup

# Force overwrite existing files
python3 batch_mmgbsa_submit.py ligand_list.txt --force

# Add delay between submissions
python3 batch_mmgbsa_submit.py ligand_list.txt --delay 2
```

#### Example

```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/scripts
python3 batch_mmgbsa_submit.py ../complexes/systems_gpu/ligand_list.txt --force
```

---

## Extract and Analyze Results

### Results Extraction Script: `extract_mmgbsa_results.py`

Parse `FINAL_RESULTS_MMGBSA.dat` files and create summary tables.

#### Basic Usage

```bash
# Check completion status
python3 extract_mmgbsa_results.py ligand_list.txt --status-only

# Extract results and save to CSV
python3 extract_mmgbsa_results.py ligand_list.txt

# Specify output file
python3 extract_mmgbsa_results.py ligand_list.txt --output my_results.csv
```

#### Example Output

```
================================================================================
Binding Free Energy Results (sorted by affinity)
================================================================================
Rank   Ligand ID            Type       ΔG (kcal/mol)   Std Dev
--------------------------------------------------------------------------------
1      GRAS_467             Type II      -45.23 ±  3.12
2      FeedAdd_570          Type II      -42.18 ±  2.98
3      ChMed_2470           Type II      -38.45 ±  3.45
...
```

#### CSV Output Columns

- `Rank` - Ranking by binding affinity
- `Ligand_ID` - Ligand identifier
- `Type` - Type I or Type II inhibitor
- `Delta_G_Total` - Total binding free energy (kcal/mol)
- `Total_Std` - Standard deviation
- `Delta_G_Gas` - Gas phase energy
- `Delta_G_Solv` - Solvation energy
- `VDW`, `EEL`, `EGB`, `ESURF` - Energy components

---

## Single System Processing

### Setup: `step4_run_mmgbsa.py`

Process a single ligand system.

```bash
# Setup MM/GBSA for one ligand
python3 step4_run_mmgbsa.py <LIGAND_ID>

# Force overwrite
python3 step4_run_mmgbsa.py <LIGAND_ID> --force

# Then submit manually
cd ../complexes/systems_gpu/<LIGAND_ID>
sbatch run_mmgbsa.sh
```

---

## Monitoring

### Check Job Queue

```bash
# All MM/GBSA jobs
squeue -u $(whoami) | grep mmgbsa

# Count running jobs
squeue -u $(whoami) | grep mmgbsa | wc -l
```

### Monitor Specific Job

```bash
# Real-time log
tail -f /path/to/system/mmgbsa_<JOB_ID>.log

# Check for completion
ls -lh /path/to/system/FINAL_RESULTS_MMGBSA.dat

# View binding energy
grep "DELTA TOTAL" /path/to/system/FINAL_RESULTS_MMGBSA.dat
```

### Check Status of All Systems

```bash
python3 extract_mmgbsa_results.py ligand_list.txt --status-only
```

---

## Typical Workflow

### 1. After MD Simulations Complete

```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/scripts

# Check which MD simulations are complete
python3 batch_mmgbsa_submit.py ../complexes/systems_gpu/ligand_list.txt --check-only
```

### 2. Submit MM/GBSA Calculations

```bash
# Setup and submit all completed MD systems
python3 batch_mmgbsa_submit.py ../complexes/systems_gpu/ligand_list.txt --force
```

### 3. Monitor Progress

```bash
# Check how many are running/complete
python3 extract_mmgbsa_results.py ../complexes/systems_gpu/ligand_list.txt --status-only

# Or check queue
squeue -u $(whoami) | grep mmgbsa
```

### 4. Extract Results (when complete)

```bash
# Extract and save results
python3 extract_mmgbsa_results.py ../complexes/systems_gpu/ligand_list.txt --output results.csv

# Results will be printed to terminal and saved to CSV
```

### 5. Analyze Results

```bash
# View CSV in Excel/LibreOffice
libreoffice results.csv

# Or view with column
column -t -s, results.csv | less -S
```

---

## File Locations

### Input Files (per system)
- `md_simulation/complexes/systems_gpu/<LIGAND_ID>/prod.nc` - Production trajectory

### MM/GBSA Setup Files (generated)
- `prep_topologies.cpptraj` - Topology preparation
- `strip_traj.cpptraj` - Trajectory stripping
- `mmgbsa.in` - MMPBSA.py input
- `run_mmgbsa.sh` - SLURM script

### Output Files (generated)
- `FINAL_RESULTS_MMGBSA.dat` - **Main results**
- `ENERGY_DECOMP.csv` - Per-frame energies
- `mmgbsa_<JOB_ID>.log` - Job log

---

## Troubleshooting

### MM/GBSA Job Failed

1. Check log file:
   ```bash
   tail -100 mmgbsa_<JOB_ID>.log
   ```

2. Common issues:
   - **"prod.nc not found"** → MD not complete
   - **"Box information"** → Should be fixed by parmed (check script)
   - **"Atom mismatch"** → Topology issue, check masks in `config.py`

3. Resubmit:
   ```bash
   python3 step4_run_mmgbsa.py <LIGAND_ID> --force
   cd ../complexes/systems_gpu/<LIGAND_ID>
   sbatch run_mmgbsa.sh
   ```

### Results Look Strange

1. Check trajectory quality:
   ```bash
   # Check Fe-ligand distance during MD
   cd /home/zhangsd/repos/CYP17A1/md_simulation/scripts
   bash check_md_progress.sh ../complexes/systems_gpu/<LIGAND_ID>
   ```

2. Check frame count:
   ```bash
   # Should show ~50000 frames
   cpptraj -p complex_solv.prmtop -y prod.nc -tl
   ```

3. Visualize trajectory:
   ```bash
   # Load in VMD/PyMOL to check for issues
   vmd complex_solv.prmtop prod.nc
   ```

---

## Configuration

Edit `config.py` to modify:

### MM/GBSA Parameters
```python
MMGBSA_PARAMS = {
    'startframe': 1,      # First frame to analyze
    'endframe': 50000,    # Last frame (50ns / 1ps)
    'interval': 5,        # Sample every N frames
    'igb': 2,             # GB model (2=OBC)
    'saltcon': 0.15       # Salt concentration (M)
}
```

### Masks
```python
RECEPTOR_MASK = ":1-475"      # Protein + heme + Fe
LIGAND_MASK = ":476"          # Ligand residue
STRIP_MASK = ":WAT,Na+,Cl-"   # Remove from trajectory
```

---

## Performance

- **Typical runtime:** 30-60 minutes per system
- **Frames analyzed:** ~10,000 per system (5 ps interval)
- **CPU usage:** 32 cores per job
- **Memory:** ~4-8 GB per job
- **Parallel capacity:** Limited by C9654 partition resources

---

## References

- **Amber MMPBSA.py manual:** https://ambermd.org/doc12/Amber22.pdf
- **MM/GBSA theory:** Kollman et al., Acc. Chem. Res. 2000
- **GB models:** Onufriev et al., Proteins 2004

---

**Created:** 2024-12-26  
**Author:** zhangshd

