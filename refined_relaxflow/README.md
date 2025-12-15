# Refined Electride Workflow - High-Precision VASP Calculations

Refined electride workflow for high-precision structure optimization of promising electride candidates identified from the initial screening workflow.

---

## Quick Start

```bash
# 1. Ensure prerequisites are met (see Prerequisites section)
#    - workflow_manager.py completed
#    - analyze.py completed
#    - filter_comb_db.py completed

# 2. Submit refined workflow
bash run_refineflow.sh \
    --input electride_candidates.db \
    --vasp-jobs-dir VASP_JOBS/ \
    --max-concurrent 10

# 3. Monitor progress
tail -f refine_workflow_<JOBID>.out

# 4. Check workflow status
python3 workflow_status.py REFINE_VASP_JOBS/workflow.json
```

---

## Overview

The **Refined Electride Workflow** is a high-precision DFT relaxation pipeline designed to perform careful 3-step structure optimization on electride candidates that passed initial screening and filtering. It starts from **original unrelaxed structures** (from mattergen CIF files) and applies progressive relaxation with strict electronic convergence checks at each step.

### Relationship to Original Workflow

```
┌─────────────────────────────────────────────────────────────────────┐
│                     ORIGINAL WORKFLOW                               │
│  (workflow_manager.py - Initial Screening)                          │
├─────────────────────────────────────────────────────────────────────┤
│  1. Relax (NSW=30, ISIF=3, EDIFFG=-0.01)                            │
│  2. SC (self-consistent static calculation)                          │
│  3. PARCHG (5 energy windows: e0025, e05, e10, band0, band1)       │
│  4. ELF (electron localization function)                             │
├─────────────────────────────────────────────────────────────────────┤
│  → analyze.py: Identify electride candidates                        │
│  → filter_comb_db.py: Apply strict filtering + symmetrization       │
│                                                                      │
│  Filtering criteria:                                                 │
│    - max(e0025, e05, e10) >= 20 Å³                                  │
│    - max(band0, band1) >= 20 Å³                                     │
│    - PyXtal symmetrization                                           │
│    - Duplicate removal                                               │
│                                                                      │
│  Output: electride_candidates.db / electride_candidates.csv         │
└─────────────────────────────────────────────────────────────────────┘
                                ↓
┌─────────────────────────────────────────────────────────────────────┐
│                     REFINED WORKFLOW                                │
│  (refine_electrideflow.py - High-Precision Refinement)             │
├─────────────────────────────────────────────────────────────────────┤
│  Load: Relaxed CONTCARs from original VASP_JOBS                    │
│  Symmetrize: PyXtal with progressive tolerance (same as original)   │
│  Filter: Only high-symmetry candidates (space group > 15)           │
│                                                                      │
│  3-Step High-Precision Relaxation:                                  │
│    Step 1: NSW=30, ISIF=3, EDIFFG=-0.005, POTIM=0.3               │
│      → Check electronic convergence (OUTCAR)                        │
│    Step 2: NSW=30, ISIF=3, EDIFFG=-0.005, POTIM=0.3               │
│      → Check electronic convergence (OUTCAR)                        │
│    Step 3: NSW=40, ISIF=3, EDIFFG=-0.002, POTIM=0.2 (final)       │
│      → Check electronic convergence (OUTCAR)                        │
│                                                                      │
│  Timeout handling: Continue if CONTCAR exists + electronic converged│
│  Output: High-quality refined structures in REFINE_VASP_JOBS/       │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Key Differences from Original Workflow

| Aspect | Original Workflow | Refined Workflow |
|--------|------------------|------------------|
| **Input structures** | Original CIFs from mattergen | Relaxed CONTCARs from VASP_JOBS |
| **Symmetrization** | PyXtal (progressive tolerance) | PyXtal (progressive tolerance, same) |
| **Relaxation steps** | 1 step (NSW=30) | 3 steps (30+30+40 ionic steps) |
| **Force convergence** | EDIFFG=-0.01 eV/Å | -0.005 → -0.005 → -0.002 eV/Å |
| **POTIM** | 0.3 (default) | 0.3 → 0.3 → 0.2 (careful) |
| **Electronic check** | End only (vasprun.xml) | After each step (OUTCAR grep) |
| **Structure count** | All generated (~thousands) | Filtered candidates (~hundreds) |
| **Space group filter** | None | Only space group > 15 |
| **Timeout handling** | Mark as failed or timeout | Continue if electronic converged |
| **Purpose** | Initial screening | High-precision refinement |
| **Computational cost** | ~30 min/structure | ~1.5 hours/structure |

---

## Features

### 1. Progressive Relaxation Strategy
- **Step 1**: Initial relaxation with moderate convergence (EDIFFG=-0.005)
- **Step 2**: Continued relaxation with same settings
- **Step 3**: Final high-precision relaxation (EDIFFG=-0.002, POTIM=0.2)

### 2. Electronic Convergence Validation
After each step, checks OUTCAR for electronic SCF convergence:
```bash
if grep -q "aborting loop because EDIFF is reached" OUTCAR; then
    echo "  Electronic SCF converged in step X"
else
    echo "ERROR: Electronic SCF did not converge in step X"
    exit 1
fi
```
**Benefit**: Fails fast if electronic SCF doesn't converge, saving 1-2 hours per structure

### 3. Graceful Timeout Handling
For SLURM timeout (exit codes 140/143):
- Checks if CONTCAR exists
- Checks if electronic SCF converged (OUTCAR)
- If both: marks as `RELAX_TMOUT` and `VASP_DONE` → usable partial result
- If not: marks as `VASP_FAILED`

### 4. High-Symmetry Filtering
Only processes orthorhombic and higher symmetry structures:
- **Excluded**: Triclinic (1-2), Monoclinic (3-15)
- **Included**: Orthorhombic (16-74), Tetragonal (75-142), Trigonal (143-167), Hexagonal (168-194), Cubic (195-230)

**Rationale**: High-symmetry structures are more likely to be synthesizable and stable

### 5. Debug-Friendly Structure Tracking
Saves `POSCAR-1`, `POSCAR-2`, `POSCAR-3` before each step for failure diagnosis

### 6. Comprehensive Cleanup
Removes large intermediate files at each step:
- Between steps: Removes WAVECAR, CHGCAR, CHG, WFULL, TMPCAR, AECCAR*
- Keeps for debugging: OUTCAR, OSZICAR, vasprun.xml
- Final cleanup: Removes all large files, keeps essential outputs

---

## Prerequisites

### 1. Complete Original Workflow
```bash
# Step 1: Run initial screening workflow
python3 workflow_manager.py \
    --results-dir mattergen_results/ternary_csp_electrides/ \
    --output-dir /scratch/$USER/VASP_JOBS \
    --max-concurrent 10

# Wait for completion, then:
```

### 2. Run Analysis and Filtering
```bash
# Step 2: Analyze ELF calculations to identify electride candidates
sbatch analyze.sh

# Step 3: Apply strict filtering criteria
python3 filter_comb_db.py \
    --input electride_data.db \
    --csv electride_analysis.csv \
    --output electride_candidates.db \
    --output-csv electride_candidates.csv \
    --min-energy 20 \
    --min-band 20

# This creates:
#   - electride_candidates.db (PyXtal database with symmetrized structures)
#   - electride_candidates.csv (CSV for quick viewing, sorted by space group)
```

### 3. Required Files
- `electride_candidates.db` or `electride_candidates.csv` (output from filter_comb_db.py)
- Original VASP_JOBS directory with relaxed CONTCARs (from workflow_manager.py)
- Refined workflow scripts:
  - `refine_electrideflow.py` (Python workflow manager)
  - `run_refineflow.sh` (submission wrapper)
  - `submit_refineflow.sh` (SLURM job script)

---

## Submission Scripts

### `run_refineflow.sh` (Wrapper Script)
- User-facing submission script
- Validates input files and paths
- Counts electride candidates
- Exports configuration as environment variables
- Submits `submit_refineflow.sh` as SLURM job
- Provides job monitoring commands

### `submit_refineflow.sh` (SLURM Job Script)
- Runs as SLURM job with 30-day time limit
- Activates conda environment
- Validates dependencies (pymatgen, PyXtal)
- Executes `refine_electrideflow.py` with configuration
- Logs all output to `refine_workflow_<JOBID>.out`

### Workflow Flow
```
User → run_refineflow.sh → sbatch submit_refineflow.sh → refine_electrideflow.py
                                     (SLURM job)              (Python manager)
                                          ↓
                                 Monitors and submits
                                   VASP jobs (SLURM)
```

---

## Installation

No additional installation needed beyond the original workflow requirements:
- Python 3.8+
- pymatgen
- PyXtal (optional, for database loading)
- pandas

---

## Usage

### Recommended: Using Submission Scripts

The refined workflow manager should be submitted as a SLURM job for long-running execution:

```bash
# Basic usage (uses defaults)
bash run_refineflow.sh

# Custom settings
bash run_refineflow.sh \
    --input electride_candidates.db \
    --vasp-jobs-dir VASP_JOBS/ \
    --output-dir ./REFINE_VASP_JOBS \
    --max-concurrent 10

# Using CSV file
bash run_refineflow.sh \
    --input electride_candidates.csv \
    --vasp-jobs-dir VASP_JOBS/
```

**What it does**:
- Submits `submit_refineflow.sh` as a SLURM job (30-day time limit)
- Workflow manager runs in background, monitoring and submitting VASP jobs
- Automatically resumes if interrupted

**Monitoring the submitted job**:
```bash
# Check SLURM job status
squeue -u $USER

# View workflow manager log
tail -f refine_workflow_<JOBID>.out

# Check workflow database
python3 workflow_status.py REFINE_VASP_JOBS/workflow.json

# Cancel workflow manager
scancel <JOBID>
```

### Alternative: Direct Python Execution

For testing or debugging, you can run the workflow manager directly:

```bash
# Using database file
python3 refine_electrideflow.py \
    --input electride_candidates.db \
    --vasp-jobs-dir VASP_JOBS/ \
    --output-dir ./REFINE_VASP_JOBS \
    --db workflow.json \
    --max-concurrent 5

# Using CSV file
python3 refine_electrideflow.py \
    --input electride_candidates.csv \
    --vasp-jobs-dir VASP_JOBS/ \
    --output-dir ./REFINE_VASP_JOBS \
    --db workflow.json \
    --max-concurrent 5
```

**Note**: Direct execution is not recommended for production runs as it requires keeping your terminal session active for days/weeks.

### Common Options (run_refineflow.sh)

| Option | Description | Default |
|--------|-------------|---------|
| `--input` | Path to electride_candidates.db or .csv | `./electride_candidates.db` |
| `--vasp-jobs-dir` | Path to original VASP_JOBS directory | `./VASP_JOBS` |
| `--output-dir` | Output directory for VASP jobs | `./REFINE_VASP_JOBS` |
| `--max-concurrent` | Max concurrent structures running | 10 |
| `--max-structures` | Limit number of structures to process | 0 (all) |
| `--check-interval` | Status check interval (seconds) | 60 |
| `--conda-env` | Conda environment name | `vaspflow` |
| `--help` | Show help message | - |

### Additional Options (refine_electrideflow.py only)

| Option | Description | Default |
|--------|-------------|---------|
| `--db` | JSON database filename | `workflow.json` |
| `--init-only` | Only initialize database, don't start monitoring | False |

### Examples

**1. Basic submission with defaults**
```bash
bash run_refineflow.sh
```

**2. Process first 50 high-symmetry candidates**
```bash
bash run_refineflow.sh \
    --max-structures 50 \
    --max-concurrent 10
```

**3. Custom paths and settings**
```bash
bash run_refineflow.sh \
    --input /path/to/electride_candidates.db \
    --vasp-jobs-dir /scratch/$USER/VASP_JOBS/ \
    --output-dir /scratch/$USER/REFINE_VASP_JOBS \
    --max-concurrent 20
```

**4. Using CSV file instead of database**
```bash
bash run_refineflow.sh \
    --input electride_candidates.csv \
    --vasp-jobs-dir VASP_JOBS/ \
    --max-concurrent 15
```

**5. Resume interrupted workflow**
```bash
# Simply re-run the same command - it will automatically resume
bash run_refineflow.sh

# The database (workflow.json) tracks progress
# No need to specify resume - it detects existing database
```

**6. Test run (limit structures and concurrency)**
```bash
bash run_refineflow.sh \
    --max-structures 5 \
    --max-concurrent 2
```

**7. View help and all options**
```bash
bash run_refineflow.sh --help
```

---

## Workflow States

The refined workflow uses a JSON database (`workflow.json`) to track job states:

```
PENDING → RELAX_RUNNING → RELAX_DONE (final state)
                       ↓
                   RELAX_TMOUT (usable partial result, final state)
                       ↓
                   RELAX_FAILED
```

### State Descriptions

| State | Description | Next Action |
|-------|-------------|-------------|
| `PENDING` | Structure loaded, waiting for submission | Submit 3-step relaxation job |
| `RELAX_RUNNING` | VASP relaxation job running | Monitor SLURM status |
| `RELAX_DONE` | All 3 steps completed successfully (final state) | Analysis complete |
| `RELAX_TMOUT` | Timeout but usable result (final state) | Usable for analysis |
| `RELAX_FAILED` | VASP failed or electronic not converged | Manual inspection needed |

---

## Directory Structure

```
REFINE_VASP_JOBS/
├── Ba2N_s001/
│   └── Relax/
│       ├── INCAR          # VASP input settings
│       ├── POSCAR         # Original structure
│       ├── POTCAR         # Pseudopotentials
│       ├── KPOINTS        # K-point mesh
│       ├── CONTCAR        # Final relaxed structure
│       ├── OUTCAR         # Detailed output (for convergence check)
│       ├── OSZICAR        # Convergence trajectory
│       ├── vasprun.xml    # Full VASP output
│       ├── POSCAR-1       # Structure before step 1 (debug)
│       ├── POSCAR-2       # Structure before step 2 (debug)
│       ├── POSCAR-3       # Structure before step 3 (debug)
│       ├── job.sh         # SLURM submission script
│       ├── vasp_*.out     # SLURM stdout
│       ├── vasp_*.err     # SLURM stderr
│       ├── VASP_DONE      # Success marker
│       ├── VASP_FAILED    # Failure marker
│       └── RELAX_TMOUT    # Timeout marker (if applicable)
├── Ca5P3_s002/
│   └── Relax/
│       └── ...
└── ...
```

---

## Monitoring Progress

### 1. Real-time Console Output

The workflow manager prints status updates:
```
[2024-12-12 14:30:15] Checking job status...
Currently running: 5/5
  Ba2N_s001: Relax completed
  Ca5P3_s002: Relax completed

Statistics:
  PENDING: 45
  RELAX_RUNNING: 5
  RELAX_DONE: 32
  RELAX_TMOUT: 3
  RELAX_FAILED: 5

Sleeping for 60s...
```

### 2. Query Database Status

```bash
# View all structures and their states
python3 workflow_status.py --db workflow.json

# View only failed structures
python3 workflow_status.py --db workflow.json --state RELAX_FAILED

# View structures by composition
python3 workflow_status.py --db workflow.json --composition Ba2N
```

### 3. Check SLURM Queue

```bash
# Check running jobs
squeue -u $USER | grep refine

# Check job details
squeue -j <job_id> -l

# Check job output
tail -f REFINE_VASP_JOBS/Ba2N_s001/Relax/vasp_*.out
```

### 4. Inspect Individual Jobs

```bash
# Check which step is running
ls -lh REFINE_VASP_JOBS/Ba2N_s001/Relax/POSCAR-*

# POSCAR-1, POSCAR-2, POSCAR-3 indicate steps 1, 2, 3 were attempted

# Check electronic convergence
grep "aborting loop because EDIFF is reached" \
    REFINE_VASP_JOBS/Ba2N_s001/Relax/OUTCAR

# Check ionic convergence
tail -20 REFINE_VASP_JOBS/Ba2N_s001/Relax/OSZICAR
```

---

## Understanding Results

### Successful Completion (`RELAX_DONE`)

Files to analyze:
- **`CONTCAR`**: Final relaxed structure (use for further analysis)
- **`OUTCAR`**: Electronic and ionic convergence information
- **`OSZICAR`**: Energy and force convergence trajectory
- **`vasprun.xml`**: Full calculation details

Check convergence:
```bash
# Check final forces (should be < 0.002 eV/Å for step 3)
grep "TOTAL-FORCE" OUTCAR | tail -1

# Check energy convergence
tail -10 OSZICAR

# View structure
ase gui CONTCAR
```

### Timeout with Usable Result (`RELAX_TMOUT`)

The structure timed out but:
- Electronic SCF converged (checked via OUTCAR)
- CONTCAR exists (partial ionic relaxation)

**Action**: Acceptable for most purposes. If critical, can extend time limit and re-run.

### Failed Jobs (`RELAX_FAILED`)

Common causes:
1. **Electronic SCF not converged**: Difficult electronic structure (narrow band gap, magnetic, etc.)
2. **CONTCAR missing**: VASP crashed before writing output
3. **Structure unstable**: Large forces or energy not decreasing

**Diagnosis**:
```bash
# Check error messages
cat REFINE_VASP_JOBS/Ba2N_s001/Relax/vasp_*.err

# Check which step failed (count POSCAR-* files)
ls REFINE_VASP_JOBS/Ba2N_s001/Relax/POSCAR-*

# Check OSZICAR for convergence behavior
cat REFINE_VASP_JOBS/Ba2N_s001/Relax/OSZICAR
```

**Solutions**:
- Adjust POTIM (reduce if oscillating)
- Increase EDIFF (if electronic convergence is tight)
- Change ALGO (try All or VeryFast)
- Check structure validity (no overlapping atoms)

---

## Output Files

### Primary Outputs

After successful completion, each structure directory contains:

| File | Description | Size |
|------|-------------|------|
| `CONTCAR` | Final relaxed structure | ~10 KB |
| `POSCAR` | Original input structure | ~10 KB |
| `OUTCAR` | Detailed VASP output | ~10-50 MB |
| `OSZICAR` | Convergence trajectory | ~10-100 KB |
| `vasprun.xml` | Full calculation data | ~10-50 MB |
| `POSCAR-1/2/3` | Debug snapshots | ~10 KB each |

### Database Files

- **`workflow.json`**: Job tracking database
  - Structure metadata
  - Job states and IDs
  - Timestamps
  - Error messages (if any)

---

## Performance and Resource Usage

### Computational Cost

| Resource | Original Workflow | Refined Workflow |
|----------|------------------|------------------|
| **Time/structure** | ~30 minutes | ~1.5 hours |
| **CPU cores** | 16 per job | 16 per job |
| **Memory** | 32 GB per job | 32 GB per job |
| **Disk/structure** | ~200 MB | ~100 MB (with cleanup) |
| **Total structures** | ~2000-5000 | ~100-300 (filtered) |

### Estimated Total Time

For 200 filtered candidates with 10 concurrent jobs:
- Wall time: ~30 hours (20 batches × 1.5 hours/batch)
- Total core-hours: 200 × 1.5 × 16 = 4,800 core-hours
- Disk usage: 200 × 100 MB = 20 GB

---

## Troubleshooting

### Common Issues

**1. "ERROR: VASP_JOBS directory not found"**
```
Solution: Ensure VASP_JOBS/ path is correct
Check: ls VASP_JOBS/
```

**2. "CONTCAR not found for structure_id"**
```
Cause: Structure ID in candidates file doesn't exist in VASP_JOBS
Solution: Verify the structure was successfully relaxed in the original workflow
Check: ls VASP_JOBS/Ba2N/Ba2N_s001/Relax/CONTCAR
```

**3. "Electronic SCF did not converge in step X"**
```
Cause: Electronic structure convergence issues
Solutions:
  - Check OUTCAR for convergence behavior
  - May need to adjust ALGO, EDIFF, or electronic mixing parameters
  - Structure may have metallic character or narrow band gap
```

**4. "Job timed out without producing CONTCAR"**
```
Cause: VASP crashed or didn't finish before time limit
Solutions:
  - Increase --time in SLURM script (currently 00:30:00 per step)
  - Check vasp_*.err for error messages
  - May need to reduce POTIM if oscillating
```

**5. "Too many structures, disk space running out"**
```
Solution: Use cleanup script to remove large intermediate files
python3 cleanup_vaspfail.py --db workflow.json --clean
```

### Reset Failed Jobs

To resubmit failed structures:

```bash
# Reset failed jobs to PENDING
python3 reset_failed_jobs.py \
    --db workflow.json \
    --state RELAX_FAILED \
    --reset-to PENDING \
    --clean  # Optional: clean directories
```

---

## Advanced Usage

### Custom VASP Settings

Edit `refine_electrideflow.py` to modify INCAR settings for each step:

```python
# Step 1 settings (line ~400)
user_incar_settings={
    'NSW': 30,
    'EDIFFG': -0.005,
    'POTIM': 0.3,
    # Add custom settings here
}

# Step 3 settings (line ~430)
user_incar_settings={
    'NSW': 40,
    'EDIFFG': -0.002,
    'POTIM': 0.2,
    # Add custom settings here
}
```

### Custom Space Group Filter

By default, only processes space groups > 15. To modify:

```python
# In load_electride_candidates() function
if space_group_num is not None and space_group_num <= 15:  # Change 15 to your threshold
    filtered_count += 1
    continue
```

### Batch Processing by Composition

Process specific compositions:

```bash
# Extract composition-specific structures from CSV
grep "^Ba2N" electride_candidates.csv > Ba2N_candidates.csv

# Run refined workflow for Ba2N only
bash run_refineflow.sh \
    --input Ba2N_candidates.csv \
    --output-dir ./REFINE_VASP_JOBS_Ba2N \
    --max-concurrent 5
```

---

## Best Practices

### 1. Start Small
```bash
# Test with first 10 structures
bash run_refineflow.sh --max-structures 10 --max-concurrent 2

# Monitor the test run
tail -f refine_workflow_<JOBID>.out
```

### 2. Monitor Disk Space
```bash
# Check disk usage regularly
du -sh REFINE_VASP_JOBS/

# Clean up failed jobs periodically
python3 cleanup_vaspfail.py --db workflow.json --clean
```

### 3. Use Appropriate Concurrency
- HPC cluster: `--max-concurrent 10-20`
- Shared resources: `--max-concurrent 3-5`
- Test run: `--max-concurrent 1-2`

### 4. Check First Completed Structure
Before running all structures, verify first completion:
```bash
# Wait for first structure to complete
# Then inspect convergence
grep "TOTAL-FORCE" REFINE_VASP_JOBS/*/Relax/OUTCAR | head -20
tail REFINE_VASP_JOBS/*/Relax/OSZICAR
```

### 5. Regular Database Backups
```bash
# Backup database regularly
cp workflow.json workflow_backup_$(date +%Y%m%d).json
```

---

## Integration with Original Workflow

### Workflow Pipeline

```
1. Initial Screening (workflow_manager.py)
   → VASP_JOBS/
   
2. Analysis (analyze.py)
   → electride_data.db
   → electride_analysis.csv
   
3. Filtering (filter_comb_db.py)
   → electride_candidates.db
   → electride_candidates.csv
   
4. Refined Relaxation (refine_electrideflow.py)  ← THIS WORKFLOW
   → REFINE_VASP_JOBS/
   
5. Final Analysis
   - Extract refined CONTCARs
   - Re-run ELF/PARCHG on refined structures (if needed)
   - Compare with initial screening results
```

### Data Flow

```
Original CIFs → Initial Relax → SC/PARCHG/ELF → Analysis → Filtering
                     ↓                                         ↓
                 VASP_JOBS/                    electride_candidates.db/csv
                  (POSCARs)                           (structure IDs)
                     └────────────────┬─────────────────────┘
                                      ↓
                        Refined 3-Step Relaxation
                         (using VASP_JOBS POSCARs)
                                      ↓
                        High-quality refined structures
```

---

## Support and Contact

For issues or questions:
1. Check this README for common solutions
2. Inspect VASP error files (`vasp_*.err`)
3. Review OUTCAR and OSZICAR for convergence issues
4. Consult VASP manual for parameter adjustments

---

## Citation

If you use this workflow, please cite:
- VASP: [Kresse & Furthmüller, Phys. Rev. B 54, 11169 (1996)]
- Pymatgen: [Ong et al., Comput. Mater. Sci. 68, 314 (2013)]
- PyXtal: [Fredericks et al., Comput. Phys. Commun. 261, 107810 (2021)]

---

## License

Same as the main VASP workflow package.

---

**Last Updated**: December 2024

