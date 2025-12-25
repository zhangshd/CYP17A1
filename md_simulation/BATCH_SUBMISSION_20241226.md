# Batch MD Submission Report

**Date:** 2024-12-26  
**Operator:** zhangshd  
**Script:** batch_setup_and_submit.py

---

## Summary

- **Total ligands processed:** 29
- **Setup successful:** 29 (100%)
- **Jobs submitted:** 29 (100%)
- **Job IDs:** 213470-213498
- **Type II inhibitors auto-detected:** 11 (38%)

---

## Type II Inhibitor Detection Results

The pipeline automatically detected **11 Type II inhibitors** based on N/O atoms within 4.0 Å of Fe:

| # | Ligand ID | Coord Atom | Element | Fe Distance (Å) | Job ID |
|---|-----------|------------|---------|-----------------|--------|
| 1 | ChMed_2470 | O52 | O | 1.31 | 213473 |
| 2 | ChMed_5302 | O31 | O | 2.18 | 213476 |
| 3 | ChMed_5343 | O44 | O | 2.24 | 213477 |
| 4 | ChMed_5859 | O31 | O | 2.39 | 213478 |
| 5 | FeedAdd_44 | O6 | O | 2.82 | 213481 |
| 6 | FeedAdd_556 | O13 | O | 3.28 | 213485 |
| 7 | FeedAdd_570 | O4 | O | 1.85 | 213486 |
| 8 | FeedAdd_584 | O19 | O | 3.27 | 213487 |
| 9 | FeedAdd_696 | O28 | O | 3.73 | 213488 |
| 10 | GRAS_2010 | O19 | O | 2.29 | 213490 |
| 11 | GRAS_467 | N3 | N | 2.16 | 213495 |

**Note:** AER601 (N22, 2.19 Å) is configured manually in `LIGAND_FE_COORDINATION`.

---

## Non-Type II Ligands (Type I or Substrates)

These 18 ligands were processed without Fe-coordination restraints:

| # | Ligand ID | Job ID | Notes |
|---|-----------|--------|-------|
| 1 | ChMed_1332 | 213470 | No N/O within 4Å |
| 2 | ChMed_1788 | 213471 | No N/O within 4Å |
| 3 | ChMed_1789 | 213472 | No N/O within 4Å |
| 4 | ChMed_2851 | 213474 | No N/O within 4Å |
| 5 | ChMed_3370 | 213475 | No N/O within 4Å |
| 6 | FeedAdd_146 | 213479 | No N/O within 4Å |
| 7 | FeedAdd_237 | 213480 | No N/O within 4Å |
| 8 | FeedAdd_45 | 213482 | No N/O within 4Å |
| 9 | FeedAdd_46 | 213483 | No N/O within 4Å |
| 10 | FeedAdd_47 | 213484 | No N/O within 4Å |
| 11 | GRAS_1015 | 213489 | No N/O within 4Å |
| 12 | GRAS_2397 | 213491 | No N/O within 4Å |
| 13 | GRAS_2546 | 213492 | No N/O within 4Å |
| 14 | GRAS_325 | 213493 | No N/O within 4Å |
| 15 | GRAS_338 | 213494 | No N/O within 4Å |
| 16 | GRAS_519 | 213496 | No N/O within 4Å |
| 17 | GRAS_573 | 213497 | No N/O within 4Å |
| 18 | GRAS_987 | 213498 | No N/O within 4Å |

---

## Processing Steps

For each ligand:

1. **Clean old files**
   - Removed `*.rst7`, `*.out`, `*.nc`, MD output files
   - Kept topology files (`complex_*.prmtop`, `complex_*.inpcrd`)

2. **Setup MD inputs**
   - Generated minimization, heating, equilibration, production inputs
   - **Auto-detected Type II inhibitors** (new feature!)
   - Applied Fe-coordination distance restraints for Type II inhibitors:
     - Minimization/Heating: 50.0 kcal/mol/Å²
     - Equilibration 1: 30.0 kcal/mol/Å²
     - Equilibration 2: 10.0 kcal/mol/Å²
     - Production: No restraint (natural equilibration)

3. **Submit SLURM jobs**
   - Partition: G4090
   - Nodes: 1, GPUs: 1
   - Expected runtime: ~3 days per ligand

---

## Distance Restraint Details (Type II only)

### Restraint Parameters
```
Target distance: Detected Fe-coordination distance
Tolerance: ±0.3 Å (flat-bottom)
Force constant schedule:
  - min1, min2, heat: 50.0 kcal/mol/Å²
  - eq1: 30.0 kcal/mol/Å²
  - eq2: 10.0 kcal/mol/Å²
  - prod: 0.0 (no restraint)
```

### Restraint File Format
```
 &rst
  iat=<Fe_atom_num>,<ligand_atom_num>,
  r1=0.0, r2=<target-0.3>, r3=<target+0.3>, r4=999.0,
  rk2=<force_const>, rk3=<force_const>,
 &end
```

---

## Monitoring

### Check job status
```bash
squeue -u $(whoami)
```

### Monitor specific ligand
```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/complexes/systems_gpu/<LIGAND_ID>
tail -f md_<JOB_ID>.log
```

### Check Fe-coordination distance (Type II only)
```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/scripts
bash check_md_progress.sh /path/to/system_dir
```

---

## Expected Outcomes

### Type II Inhibitors
- Fe-coordination distance should remain within 1.5-3.0 Å range
- Ligand should NOT drift away during simulation
- Smooth progression through all MD stages

### Type I/Substrates
- Ligand can move freely in active site
- May approach or move away from Fe naturally
- No artificial restraints applied

---

## Notes

1. **Auto-detection threshold:** N/O atoms within 4.0 Å of Fe
2. **Priority:** Manual configuration in `LIGAND_FE_COORDINATION` overrides auto-detection
3. **Validation:** Check first few frames of Type II simulations to ensure restraints are working
4. **Performance:** Batch submission with 1-second delay between jobs to avoid overwhelming scheduler

---

## Files Generated

For each ligand system:
- `min1.in`, `min2.in` - Minimization inputs
- `heat.in` - Heating input
- `eq1.in` - Equilibration stage 1 input
- `eq2_01.in` to `eq2_10.in` - Equilibration stage 2 inputs (10 sub-stages)
- `prod.in` - Production MD input
- `fe_coord_restraint_*.RST` - Distance restraint files (Type II only)
- `run_md.sh` - SLURM submission script

---

## References

- Type II inhibitor detection: `step3_setup_md_simulation.py::detect_type_ii_inhibitor()`
- Distance restraint generation: `step3_setup_md_simulation.py::create_distance_restraint_file()`
- Batch processing script: `batch_setup_and_submit.py`
- Documentation: `TYPE_II_INHIBITOR_SUPPORT.md`

