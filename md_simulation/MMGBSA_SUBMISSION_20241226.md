# MM/GBSA Batch Submission Report

**Date:** 2024-12-26  
**Operator:** zhangshd  
**Script:** batch_mmgbsa_submit.py

---

## Summary

- **Total ligands:** 30
- **MD completed:** 26 (87%)
- **MD incomplete:** 4 (13%)
- **MM/GBSA setup successful:** 26 (100% of completed)
- **MM/GBSA jobs submitted:** 26 (100% of completed)
- **Job IDs:** 213503-213528

---

## MM/GBSA Calculations Submitted (26 systems)

| # | Ligand ID | Trajectory Size | Job ID | Status |
|---|-----------|-----------------|--------|--------|
| 1 | AER601 | 6105.7 MB | 213503 | Submitted |
| 2 | ChMed_1332 | 6107.1 MB | 213504 | Submitted |
| 3 | ChMed_1788 | 6107.1 MB | 213505 | Submitted |
| 4 | ChMed_1789 | 6107.2 MB | 213506 | Submitted |
| 5 | ChMed_2470 | 6107.7 MB | 213507 | Submitted |
| 6 | ChMed_2851 | 6107.1 MB | 213508 | Submitted |
| 7 | ChMed_3370 | 6107.1 MB | 213509 | Submitted |
| 8 | ChMed_5302 | 6106.4 MB | 213510 | Submitted |
| 9 | ChMed_5859 | 6107.6 MB | 213511 | Submitted |
| 10 | FeedAdd_146 | 6106.9 MB | 213512 | Submitted |
| 11 | FeedAdd_44 | 6106.6 MB | 213513 | Submitted |
| 12 | FeedAdd_45 | 6106.6 MB | 213514 | Submitted |
| 13 | FeedAdd_46 | 6106.6 MB | 213515 | Submitted |
| 14 | FeedAdd_47 | 6106.6 MB | 213516 | Submitted |
| 15 | FeedAdd_556 | 6106.8 MB | 213517 | Submitted |
| 16 | FeedAdd_570 | 6105.8 MB | 213518 | Submitted |
| 17 | FeedAdd_696 | 6105.4 MB | 213519 | Submitted |
| 18 | GRAS_1015 | 6106.2 MB | 213520 | Submitted |
| 19 | GRAS_2010 | 6106.4 MB | 213521 | Submitted |
| 20 | GRAS_2397 | 6106.1 MB | 213522 | Submitted |
| 21 | GRAS_325 | 6106.0 MB | 213523 | Submitted |
| 22 | GRAS_338 | 6105.9 MB | 213524 | Submitted |
| 23 | GRAS_467 | 6105.5 MB | 213525 | Submitted |
| 24 | GRAS_519 | 6106.4 MB | 213526 | Submitted |
| 25 | GRAS_573 | 6106.2 MB | 213527 | Submitted |
| 26 | GRAS_987 | 6105.8 MB | 213528 | Submitted |

**Average trajectory size:** ~6106 MB (~6.1 GB per system)

---

## MD Simulations Not Yet Complete (4 systems)

These systems are likely still running and can be processed later:

| # | Ligand ID | Status | Action Needed |
|---|-----------|--------|---------------|
| 1 | ChMed_5343 | prod.nc not found | Wait for MD completion |
| 2 | FeedAdd_237 | prod.nc not found | Wait for MD completion |
| 3 | FeedAdd_584 | prod.nc not found | Wait for MD completion |
| 4 | GRAS_2546 | prod.nc not found | Wait for MD completion |

### To Submit MM/GBSA for These Later:

```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/scripts
python3 step4_run_mmgbsa.py <LIGAND_ID> --force
cd ../complexes/systems_gpu/<LIGAND_ID>
sbatch run_mmgbsa.sh
```

Or batch submit when all complete:
```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/scripts
python3 batch_mmgbsa_submit.py ../complexes/systems_gpu/ligand_list.txt --force
```

---

## MM/GBSA Calculation Details

### SLURM Configuration
- **Partition:** C9654 (CPU)
- **Nodes:** 1
- **Tasks:** 32
- **CPUs per task:** 1
- **Expected runtime:** ~30-60 minutes per system

### Analysis Parameters

From `config.py::MMGBSA_PARAMS`:
```python
{
    'startframe': 1,      # Start from frame 1
    'endframe': 50000,    # Use all frames (50ns / 1ps = 50000)
    'interval': 5,        # Sample every 5 frames (every 5ps)
    'igb': 2,             # GB model (igb=2, OBC GBSA)
    'saltcon': 0.15       # Salt concentration (M)
}
```

**Total frames analyzed per system:** (50000 - 1) / 5 = ~10,000 frames  
**Time window:** 0-50 ns (production trajectory)

### Calculation Steps

For each system, the MM/GBSA job performs:

1. **Topology Preparation** (cpptraj + parmed)
   - Strip water/ions from `complex_solv.prmtop` → `complex_dry.prmtop`
   - Create receptor topology (protein + heme + Fe)
   - Create ligand topology
   - Remove periodic box information (required for GB)

2. **Trajectory Processing** (cpptraj)
   - Strip water/ions from `prod.nc` → `prod_dry.nc`

3. **MM/GBSA Calculation** (MMPBSA.py)
   - Calculate binding free energy using GB implicit solvent
   - ΔG_bind = G_complex - G_receptor - G_ligand
   - Output: `FINAL_RESULTS_MMGBSA.dat`, `ENERGY_DECOMP.csv`

---

## Output Files

For each successfully completed system, the following files will be generated:

### Input Files (created during setup)
- `prep_topologies.cpptraj` - Topology preparation script
- `strip_traj.cpptraj` - Trajectory stripping script
- `mmgbsa.in` - MMPBSA.py input parameters
- `run_mmgbsa.sh` - SLURM submission script

### Intermediate Files (created during calculation)
- `complex_dry.prmtop` - Complex topology without water/ions
- `receptor.prmtop` - Receptor-only topology
- `ligand.prmtop` - Ligand-only topology
- `prod_dry.nc` - Stripped production trajectory

### Final Results
- `FINAL_RESULTS_MMGBSA.dat` - **Main results file with ΔG values**
- `ENERGY_DECOMP.csv` - Energy decomposition per frame
- `mmgbsa_<JOB_ID>.log` - SLURM output log

---

## Monitoring MM/GBSA Progress

### Check job status
```bash
squeue -u $(whoami) | grep mmgbsa
```

### Monitor specific job
```bash
tail -f /home/zhangsd/repos/CYP17A1/md_simulation/complexes/systems_gpu/<LIGAND_ID>/mmgbsa_<JOB_ID>.log
```

### Check for completion
```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/complexes/systems_gpu/<LIGAND_ID>
ls -lh FINAL_RESULTS_MMGBSA.dat
```

### View binding energy
```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/complexes/systems_gpu/<LIGAND_ID>
grep "DELTA TOTAL" FINAL_RESULTS_MMGBSA.dat
```

---

## Expected Results Format

From `FINAL_RESULTS_MMGBSA.dat`:

```
DELTA Energy Terms (kcal/mol)

                         MEAN       STD       SEM
=========================================================
VDWAALS                 -XX.XX      X.XX      X.XX
EEL                     -XX.XX      X.XX      X.XX
EGB                      XX.XX      X.XX      X.XX
ESURF                    -X.XX      X.XX      X.XX

DELTA G gas            -XX.XX      X.XX      X.XX
DELTA G solv             XX.XX      X.XX      X.XX

DELTA TOTAL            -XX.XX      X.XX      X.XX  ← Binding Free Energy
```

**Key values:**
- **DELTA TOTAL (MEAN):** Binding free energy (kcal/mol)
  - More negative = stronger binding
  - Typical range: -20 to -50 kcal/mol for good inhibitors
- **STD:** Standard deviation across frames
- **SEM:** Standard error of the mean

---

## Next Steps

1. **Wait for completion** (~30-60 min per job, ~1-2 hours total with parallel execution)

2. **Check for failures**
   ```bash
   cd /home/zhangsd/repos/CYP17A1/md_simulation/scripts
   # Check which jobs completed
   for dir in ../complexes/systems_gpu/*/; do
       ligand=$(basename "$dir")
       if [ -f "$dir/FINAL_RESULTS_MMGBSA.dat" ]; then
           echo "✓ $ligand"
       else
           echo "✗ $ligand (not complete)"
       fi
   done
   ```

3. **Extract and compile results**
   - Create a summary script to extract ΔG values from all systems
   - Rank ligands by binding affinity
   - Compare Type II vs Type I inhibitors

4. **Process incomplete MD simulations**
   - Check status of ChMed_5343, FeedAdd_237, FeedAdd_584, GRAS_2546
   - Submit their MM/GBSA jobs when MD completes

---

## Performance Notes

- **Parallel execution:** 26 jobs running on C9654 partition
- **Resource utilization:** Each job uses 32 CPU cores
- **Expected completion:** ~1-2 hours (depending on queue)
- **Total compute time:** ~26-52 CPU-hours

---

## Files Generated by This Batch

Total files per system: ~10

**Total for 26 systems:** ~260 files

---

## Related Documentation

- **MM/GBSA script:** `md_simulation/scripts/step4_run_mmgbsa.py`
- **Batch script:** `md_simulation/scripts/batch_mmgbsa_submit.py`
- **Configuration:** `md_simulation/scripts/config.py`
- **MD submission:** `md_simulation/BATCH_SUBMISSION_20241226.md`
- **Type II support:** `md_simulation/TYPE_II_INHIBITOR_SUPPORT.md`

---

## Notes

1. MM/GBSA uses GB implicit solvent (igb=2) with salt concentration 0.15M
2. All trajectories are ~6.1 GB (50 ns at 2 fs timestep with 1 ps output interval)
3. Analysis uses ~10,000 frames per system (5 ps sampling interval)
4. Periodic box information is removed from topologies (required for GB calculations)
5. MMPBSA.py runs in serial mode (more stable than MPI for this system size)

---

**Status:** ✅ 26/26 MM/GBSA calculations successfully submitted and running

