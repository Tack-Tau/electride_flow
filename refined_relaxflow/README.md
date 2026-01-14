# Refined Electride Workflow - High-Precision VASP Calculations

Refined electride workflow for high-precision structure optimization of promising electride candidates identified from the initial screening workflow.

---

## Quick Start

```bash
# 1. Ensure prerequisites are met (see Prerequisites section)
#    - workflow_manager.py completed
#    - analyze.py completed
#    - filter_comb_db.py completed

# 2. Submit refined relaxation workflow
bash run_refineflow.sh \
    --input electride_candidates.db \
    --vasp-jobs-dir VASP_JOBS/ \
    --max-concurrent 10

# 3. Monitor progress
tail -f refine_workflow_<JOBID>.out

# 4. Check workflow status
python3 workflow_status.py REFINE_VASP_JOBS/workflow.json

# 5. (Optional) Run electronic structure analysis
bash run_electronic_flow.sh \
    --refine-jobs REFINE_VASP_JOBS/ \
    --output-dir ELECTRONIC_JOBS/ \
    --max-concurrent 10
```

---

## Overview

The **Refined Electride Workflow** is a high-precision DFT pipeline designed to perform careful structure optimization and electronic structure analysis on electride candidates. It consists of two main components:

1. **Refined Relaxation** (`refine_electrideflow.py`): Progressive 3-step structure optimization with tight convergence
2. **Electronic Structure Analysis** (`electronic_flow_manager.py`): SCF/PARCHG/ELF/BAND/PDOS calculations for detailed electronic properties

The workflow starts from relaxed structures (from the initial screening) and applies progressive relaxation with strict electronic convergence checks, followed by optional detailed electronic structure calculations.

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
│  3-Step Progressive Relaxation:                                     │
│    Step 1: NSW=100, ISIF=2, EDIFFG=-0.02, IBRION=2, POTIM=0.3     │
│      → Quick ionic relaxation (CG optimizer), fixed cell            │
│    Step 2: NSW=100, ISIF=3, EDIFFG=-0.01, IBRION=1, POTIM=0.2     │
│      → Full relaxation with RMM-DIIS (robust optimizer)             │
│    Step 3: NSW=100, ISIF=3, EDIFFG=-0.005, IBRION=1, POTIM=0.1    │
│      → Final high-precision with RMM-DIIS                           │
│                                                                      │
│  Timeout handling: Continue if CONTCAR exists + electronic converged│
│  Output: High-quality refined structures in REFINE_VASP_JOBS/       │
└─────────────────────────────────────────────────────────────────────┘
                                ↓ (optional)
┌─────────────────────────────────────────────────────────────────────┐
│              ELECTRONIC STRUCTURE WORKFLOW                          │
│  (electronic_flow_manager.py - Detailed Electronic Analysis)       │
├─────────────────────────────────────────────────────────────────────┤
│  Load: Refined CONTCARs from REFINE_VASP_JOBS/*/Relax/CONTCAR     │
│  Symmetrize: PyXtal → Conventional + Primitive cells               │
│                                                                      │
│  Sequential Electronic Structure Calculations:                      │
│    SPE Job (Conventional Cell):                                     │
│      → SC: Self-consistent field (250 k-pts/Å³)                    │
│      → PARCHG: 5 energy windows (visualization)                     │
│      → ELF: Electron localization function                          │
│    BAND Job (Primitive Cell):                                       │
│      → Band structure along high-symmetry k-path                    │
│    DOS Job (Primitive Cell):                                        │
│      → Projected density of states (NEDOS=3000)                     │
│                                                                      │
│  Output: Detailed electronic structure in ELECTRONIC_JOBS/         │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Key Differences from Original Workflow

| Aspect | Original Workflow | Refined Workflow |
|--------|------------------|------------------|
| **Input structures** | Original CIFs from mattergen | Relaxed CONTCARs from VASP_JOBS |
| **Symmetrization** | PyXtal (progressive tolerance) | PyXtal (progressive tolerance, same) |
| **Relaxation steps** | 1 step (NSW=30) | 3 progressive steps (NSW=100 each) |
| **ISIF** | 3 (constant) | 2 → 3 → 3 (progressive) |
| **IBRION** | 2 (CG) | 2 → 1 → 1 (CG then RMM-DIIS) |
| **Force convergence** | EDIFFG=-0.01 eV/Å | -0.02 → -0.01 → -0.005 eV/Å |
| **POTIM** | 0.3 (default) | 0.3 → 0.2 → 0.1 (progressive) |
| **Electronic check** | End only (vasprun.xml) | After each step (OUTCAR grep) |
| **Structure count** | All generated (~thousands) | Filtered candidates (~hundreds) |
| **Space group filter** | None | Only space group > 15 |
| **Timeout handling** | Mark as failed or timeout | Continue if electronic converged |
| **Purpose** | Initial screening | High-precision refinement |
| **Computational cost** | ~30 min/structure | ~2 hours/structure |
| **SLURM wall time** | 20 minutes | 2 hours |

---

## Features

### 1. Progressive Relaxation Strategy
Three-step relaxation with progressively tighter convergence and robust optimizers:

- **Step 1**: Quick ionic relaxation (CG optimizer)
  - `ISIF=2` (relax ions only, fixed cell volume and shape)
  - `IBRION=2` (Conjugate Gradient optimizer)
  - `EDIFFG=-0.02` eV/Å (moderate force convergence)
  - `POTIM=0.3` (standard step size for CG)
  - Purpose: Quick initial relaxation of ionic positions

- **Step 2**: Full relaxation (RMM-DIIS optimizer)
  - `ISIF=3` (relax ions + cell volume and shape)
  - `IBRION=1` (RMM-DIIS optimizer - more robust than CG)
  - `EDIFFG=-0.01` eV/Å (tighter forces)
  - `POTIM=0.2` (smaller step size for stability)
  - Purpose: Full cell+ion relaxation with robust optimizer

- **Step 3**: Final high-precision (RMM-DIIS optimizer)
  - `ISIF=3` (relax ions + cell volume and shape)
  - `IBRION=1` (RMM-DIIS optimizer)
  - `EDIFFG=-0.005` eV/Å (tight forces)
  - `POTIM=0.1` (small steps for precise minimum)
  - Purpose: Achieve high-precision structure for accurate ELF analysis

**Benefits**:
- ISIF=2 in step 1 avoids cell shape changes during initial relaxation
- **IBRION=1 (RMM-DIIS) in steps 2-3 avoids "ZBRENT: can't locate minimum" errors**
- RMM-DIIS is more robust than CG for cell relaxation (ISIF=3)
- Progressive tightening balances speed and accuracy
- Conservative POTIM values prevent oscillations
- Final EDIFFG=-0.005 ensures forces < 5 meV/Å

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
- Electronic structure workflow scripts:
  - `electronic_flow_manager.py` (Python workflow manager for SCF/PARCHG/ELF/BAND/PDOS)
  - `run_electronic_flow.sh` (submission wrapper)
  - `submit_electronic_flow.sh` (SLURM job script)
- MatterSim e_hull scripts:
  - `compute_mattersim_e_hull.py` (MatterSim relaxation and hull calculation)
  - `run_mattersim_e_hull.sh` (submission wrapper)
  - `submit_mattersim_e_hull.sh` (SLURM job script)
- Stable electride extraction scripts:
  - `get_stable_ele_db.py` (extract stable candidates to CSV/database)
  - `run_get_stable_db.sh` (submission wrapper)
  - `submit_get_stable_db.sh` (SLURM job script)

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
# Check final forces (should be < 0.005 eV/Å for step 3)
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
- Use IBRION=1 (RMM-DIIS) instead of IBRION=2 (CG) to avoid ZBRENT errors

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
| **Time/structure** | ~30 minutes | ~2 hours |
| **SLURM wall time** | 20 minutes | 2 hours |
| **CPU cores** | 16 per job | 16 per job |
| **Memory** | 32 GB per job | 32 GB per job |
| **Disk/structure** | ~200 MB | ~100 MB (with cleanup) |
| **Total structures** | ~2000-5000 | ~100-300 (filtered) |

### Estimated Total Time

For 200 filtered candidates with 10 concurrent jobs:
- Wall time: ~40 hours (20 batches × 2 hours/batch)
- Total core-hours: 200 × 2 × 16 = 6,400 core-hours
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
  - Increase --time in SLURM script (currently 02:00:00 total)
  - Check vasp_*.err for error messages
  - May need to reduce POTIM if oscillating
```

**5. "ZBRENT: can't locate minimum, use default step"**
```
Cause: Conjugate Gradient (IBRION=2) line search failed during cell relaxation
Solution: Use IBRION=1 (RMM-DIIS) instead of IBRION=2 for ISIF=3 steps
Note: This is already implemented in the refined workflow (IBRION=1 for steps 2-3)
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
# Step 1 settings (line ~367-370)
if step == 1:
    ediffg = -0.02   # Moderate force convergence
    isif = 2         # Fixed cell, relax ions only
    ibrion = 2       # Conjugate Gradient
    potim = 0.3      # Standard CG step size

# Step 2 settings (line ~371-374)
elif step == 2:
    ediffg = -0.01   # Tighter forces
    isif = 3         # Full relaxation
    ibrion = 1       # RMM-DIIS (robust, avoids ZBRENT)
    potim = 0.2      # Smaller step size

# Step 3 settings (line ~375-378)
elif step == 3:
    ediffg = -0.005  # Tight forces
    isif = 3         # Full relaxation
    ibrion = 1       # RMM-DIIS (robust)
    potim = 0.1      # Small step for precision
```

**Why IBRION=1 for cell relaxation (ISIF=3)?**
- Conjugate Gradient (IBRION=2) can fail with "ZBRENT: can't locate minimum" when relaxing cell shape
- RMM-DIIS (IBRION=1) is more robust for variable-cell optimizations
- This is particularly important for structures with soft modes or anisotropic cells

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

## Electronic Structure Workflow

After completing the refined relaxation workflow, you can perform detailed electronic structure calculations (SCF/PARCHG/ELF/BAND/PDOS) on the refined structures for in-depth analysis and visualization.

### Purpose

High-precision electronic structure analysis workflow:
- **SCF (Self-Consistent Field)**: Static calculation with tight convergence
- **PARCHG (Partial Charge Density)**: 5 energy windows for interstitial electron visualization
- **ELF (Electron Localization Function)**: Electronic localization analysis
- **BAND (Band Structure)**: Electronic band structure along high-symmetry k-paths
- **PDOS (Projected Density of States)**: Orbital-projected density of states

**Key Features**:
- **Dual cell strategy**: Conventional cell for visualization (SPE), primitive cell for efficiency (BAND/PDOS)
- **Sequential workflow**: SPE → BAND+DOS (BAND and DOS can run in parallel after SPE completes)
- **Separate job tracking**: Individual status markers for SC, PARCHG, ELF, BAND, and PDOS stages
- **High k-point density**: 250 k-points/Å³ for accurate electronic properties
- **Automatic compression**: PARCHG files are compressed to save disk space

### Workflow Stages

```
┌─────────────────────────────────────────────────────────────────────┐
│                  ELECTRONIC STRUCTURE WORKFLOW                      │
├─────────────────────────────────────────────────────────────────────┤
│  Input: Refined CONTCARs from REFINE_VASP_JOBS/*/Relax/CONTCAR    │
│                                                                      │
│  Stage 1: Symmetrization (PyXtal)                                  │
│    → Conventional cell (for SPE: visualization)                     │
│    → Primitive cell (for BAND/PDOS: efficiency)                     │
│                                                                      │
│  Stage 2: SPE Job (Conventional Cell)                              │
│    → SC: Self-consistent calculation (LWAVE=True, LCHARG=True)     │
│    → PARCHG: 5 energy windows (e0025, e05, e10, band0, band1)     │
│    → ELF: Electron localization function                            │
│    Status: SC_DONE → PARCHG_DONE → ELF_DONE                        │
│                                                                      │
│  Stage 3: BAND Job (Primitive Cell, parallel with DOS)            │
│    → Band structure along high-symmetry k-path                      │
│    → ICHARG=2 (atomic charge density, independent of SPE)          │
│    → 40 k-points between each high-symmetry point                   │
│    Status: BAND_DONE                                                │
│                                                                      │
│  Stage 4: PDOS Job (Primitive Cell, parallel with BAND)           │
│    → Projected density of states (LORBIT=11)                        │
│    → ICHARG=2 (atomic charge density, independent of SPE)          │
│    → NEDOS=3000 for smooth DOS                                      │
│    Status: PDOS_DONE                                                │
│                                                                      │
│  Final State: COMPLETE (when both BAND_DONE and PDOS_DONE)        │
└─────────────────────────────────────────────────────────────────────┘
```

### Scripts

1. **`electronic_flow_manager.py`** - Python workflow manager
   - Loads refined CONTCARs and symmetrizes with PyXtal
   - Generates VASP inputs for each calculation type
   - Submits three separate SLURM jobs per structure (SPE, BAND, DOS)
   - Monitors job status and advances through stages
   - Tracks individual stage completion (SC/PARCHG/ELF/BAND/PDOS)

2. **`run_electronic_flow.sh`** - User-facing submission wrapper
   - Validates input paths and structure IDs
   - Exports configuration as environment variables
   - Submits `submit_electronic_flow.sh` as SLURM job

3. **`submit_electronic_flow.sh`** - SLURM job script for workflow manager
   - Runs workflow manager with 15-day time limit
   - Activates conda environment
   - Monitors and submits VASP jobs

### Usage

#### Basic Usage

```bash
cd refined_relaxflow/

# Process all RELAX_DONE structures
bash run_electronic_flow.sh \
    --refine-jobs ./REFINE_VASP_JOBS \
    --output-dir ./ELECTRONIC_JOBS \
    --max-concurrent 10

# Process specific structure IDs
bash run_electronic_flow.sh \
    --refine-jobs ./REFINE_VASP_JOBS \
    --structure-ids Ba2N_s001 Ca5P3_s002 Y3N2_s003 \
    --max-concurrent 5

# Test run with 5 structures
bash run_electronic_flow.sh \
    --refine-jobs ./REFINE_VASP_JOBS \
    --max-concurrent 2 \
    --check-interval 120
```

#### Options (run_electronic_flow.sh)

| Option | Description | Default |
|--------|-------------|---------|
| `--refine-jobs` | Path to REFINE_VASP_JOBS directory | `./REFINE_VASP_JOBS` |
| `--output-dir` | Output directory for electronic jobs | `./ELECTRONIC_JOBS` |
| `--structure-ids` | Specific structure IDs to process | All RELAX_DONE |
| `--max-concurrent` | Max concurrent structures | 10 |
| `--check-interval` | Status check interval (seconds) | 60 |
| `--conda-env` | Conda environment name | `vaspflow` |

### Directory Structure

```
ELECTRONIC_JOBS/
├── Ba2N_s001/
│   ├── SPE/                  # Conventional cell calculations
│   │   ├── POSCAR            # Conventional cell structure
│   │   ├── INCAR-SC          # SC calculation settings
│   │   ├── INCAR-ELF         # ELF calculation settings
│   │   ├── INCAR-PARCHG-*    # PARCHG INCARs (5 energy windows)
│   │   ├── KPOINTS           # Gamma-centered Monkhorst-Pack
│   │   ├── POTCAR            # Pseudopotentials
│   │   ├── vasprun.xml-SC    # SC results
│   │   ├── OSZICAR-SC        # SC convergence
│   │   ├── ELFCAR            # ELF output
│   │   ├── PARCHG.tar.gz     # Compressed PARCHG files
│   │   ├── job.sh            # SLURM script (SPE workflow)
│   │   ├── vasp_*.out/err    # SLURM output
│   │   ├── SC_DONE           # SC completion marker
│   │   ├── PARCHG_DONE       # PARCHG completion marker
│   │   ├── ELF_DONE          # ELF completion marker
│   │   └── VASP_DONE         # SPE workflow complete
│   ├── BAND/                 # Primitive cell band structure
│   │   ├── POSCAR            # Primitive cell structure
│   │   ├── INCAR             # Band calculation settings
│   │   ├── KPOINTS           # Line-mode k-path
│   │   ├── POTCAR            # Pseudopotentials
│   │   ├── vasprun.xml       # Band structure data
│   │   ├── EIGENVAL          # Eigenvalues
│   │   ├── job.sh            # SLURM script (BAND)
│   │   ├── vasp_*.out/err    # SLURM output
│   │   └── BAND_DONE         # Completion marker
│   └── DOS/                  # Primitive cell PDOS
│       ├── POSCAR            # Primitive cell structure
│       ├── INCAR             # DOS calculation settings
│       ├── KPOINTS           # Dense Monkhorst-Pack
│       ├── POTCAR            # Pseudopotentials
│       ├── vasprun.xml       # DOS data
│       ├── DOSCAR            # Density of states
│       ├── job.sh            # SLURM script (DOS)
│       ├── vasp_*.out/err    # SLURM output
│       └── PDOS_DONE         # Completion marker
├── Ca5P3_s002/
│   └── ...
└── workflow.json             # Job tracking database
```

### Workflow States

The electronic workflow tracks structures through multiple stages:

```
PENDING → SC_RUNNING → SC_DONE
                     ↓
               PARCHG_RUNNING → PARCHG_DONE
                              ↓
                        ELF_RUNNING → ELF_DONE
                                    ↓
                              ┌─────┴─────┐
                              ↓           ↓
                      BAND_RUNNING    PDOS_RUNNING
                              ↓           ↓
                         BAND_DONE    PDOS_DONE
                              └─────┬─────┘
                                    ↓
                                COMPLETE
```

**State Descriptions**:

| State | Description | Next Stage |
|-------|-------------|------------|
| `PENDING` | Structure loaded, waiting for submission | Submit SPE job |
| `SC_RUNNING` | SC calculation running | PARCHG stage |
| `PARCHG_RUNNING` | PARCHG calculations running | ELF stage |
| `ELF_RUNNING` | ELF calculation running | BAND+DOS stages |
| `ELF_DONE` | SPE workflow complete | Submit BAND job |
| `BAND_RUNNING` | Band structure calculation running | - |
| `BAND_DONE` | Band structure complete | Submit DOS job |
| `PDOS_RUNNING` | PDOS calculation running | - |
| `PDOS_DONE` | PDOS complete | Check for COMPLETE |
| `COMPLETE` | All stages complete | Analysis ready |
| `*_FAILED` | Stage failed | Manual inspection |

### VASP Input Settings

**SC Calculation (Conventional Cell)**:
- `PREC = Accurate`
- `ALGO = Normal`
- `EDIFF = 1e-6`
- `ISMEAR = 0, SIGMA = 0.05`
- `LWAVE = True, LCHARG = True` (for PARCHG/ELF)
- `reciprocal_density = 250` (high k-point density)

**PARCHG Calculations (Conventional Cell)**:
- 5 energy windows: e0025 (E_F - 0.025), e05 (E_F - 0.5), e10 (E_F - 1.0), band0, band1
- `ICHARG = 11` (read CHGCAR/WAVECAR from SC)
- `LPARD = True` (write PARCHG)
- Generated automatically from SC vasprun.xml

**ELF Calculation (Conventional Cell)**:
- `LELF = True` (compute ELF)
- `ICHARG = 11` (read CHGCAR/WAVECAR from SC)
- `ISMEAR = -5` (tetrahedron method)
- `NEDOS = 1000`

**Band Structure (Primitive Cell)**:
- `ICHARG = 2` (atomic charge density, independent calculation)
- `LORBIT = 11` (orbital-projected)
- `ISMEAR = 0, SIGMA = 0.05`
- Line-mode k-path: 40 k-points between each high-symmetry point
- High-symmetry path determined automatically by pymatgen

**PDOS (Primitive Cell)**:
- `ICHARG = 2` (atomic charge density, independent calculation)
- `LORBIT = 11` (orbital-projected DOS)
- `ISMEAR = -5` (tetrahedron method)
- `NEDOS = 3000` (smooth DOS)
- `reciprocal_density = 250` (dense k-point sampling)

### Key Design Choices

**Why Conventional Cell for SPE?**
- **Visualization**: ELF and PARCHG are typically visualized in the conventional cell for clarity
- **Symmetry**: Easier to identify interstitial sites in conventional cell

**Why Primitive Cell for BAND/PDOS?**
- **Efficiency**: Fewer atoms → faster calculations
- **Standard practice**: Band structures are typically computed with primitive cells
- **Compatibility**: Easier comparison with literature and databases

**Why ICHARG=2 for BAND/PDOS?**
- **Cell mismatch**: BAND/PDOS use primitive cell, but SC uses conventional cell
- **Independence**: Cannot reuse CHGCAR from different cell type
- **Consistency**: Ensures VASP starts from atomic charge density for the correct cell

**Why Separate Jobs for SPE/BAND/DOS?**
- **Modularity**: Independent failure handling for each stage
- **Flexibility**: BAND and DOS can run in parallel after SPE
- **Resource efficiency**: Different time requirements for each stage

### Monitoring Progress

#### Real-time Output

```bash
# View workflow manager log
tail -f electronic_flow_*.out

# Example output:
# [2024-12-12 14:30:15] Checking status...
# Currently running: 8/10
#   Ba2N_s001: SC completed, PARCHG running
#   Ca5P3_s002: ELF completed (SPE workflow done)
#   Y3N2_s003: Band structure running
#
# Statistics:
#   PENDING: 25
#   SC_RUNNING: 3
#   PARCHG_RUNNING: 2
#   ELF_RUNNING: 1
#   ELF_DONE: 5
#   BAND_RUNNING: 2
#   PDOS_RUNNING: 1
#   COMPLETE: 12
```

#### Query Database

```bash
# Check workflow status
python3 ../workflow_status.py ELECTRONIC_JOBS/workflow.json

# View specific state
python3 ../workflow_status.py ELECTRONIC_JOBS/workflow.json --state COMPLETE

# Check individual structure
grep "Ba2N_s001" ELECTRONIC_JOBS/workflow.json
```

#### Check Individual Jobs

```bash
# Check SPE job status
tail -f ELECTRONIC_JOBS/Ba2N_s001/SPE/vasp_*.out

# Check which stage SPE is at
ls -lh ELECTRONIC_JOBS/Ba2N_s001/SPE/*_DONE

# Check band structure job
tail -f ELECTRONIC_JOBS/Ba2N_s001/BAND/vasp_*.out

# Check PDOS job
tail -f ELECTRONIC_JOBS/Ba2N_s001/DOS/vasp_*.out
```

### Output Analysis

#### Band Structure

```python
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.plotter import BSPlotter

# Load band structure
vr = Vasprun("ELECTRONIC_JOBS/Ba2N_s001/BAND/vasprun.xml")
bs = vr.get_band_structure(line_mode=True)

# Plot
plotter = BSPlotter(bs)
plotter.get_plot(ylim=[-5, 5]).savefig("band_structure.png")
```

#### Density of States

```python
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.plotter import DosPlotter

# Load DOS
vr = Vasprun("ELECTRONIC_JOBS/Ba2N_s001/DOS/vasprun.xml")
dos = vr.complete_dos

# Plot
plotter = DosPlotter()
plotter.add_dos("Total", dos)
plotter.get_plot(xlim=[-5, 5]).savefig("dos.png")
```

#### ELF Visualization

```bash
# Extract ELFCAR and visualize with VESTA or pymatgen
# ELFCAR is in ELECTRONIC_JOBS/Ba2N_s001/SPE/ELFCAR

# Using VESTA (GUI)
vesta ELECTRONIC_JOBS/Ba2N_s001/SPE/ELFCAR

# Or extract isosurface data
python3 -c "
from pymatgen.io.vasp import VolumetricData
elf = VolumetricData.from_file('ELECTRONIC_JOBS/Ba2N_s001/SPE/ELFCAR')
# Process ELF data...
"
```

#### PARCHG Visualization

```bash
# Extract PARCHG files from archive
cd ELECTRONIC_JOBS/Ba2N_s001/SPE/
tar -xzf PARCHG.tar.gz

# Visualize with VESTA
vesta PARCHG-e0025
vesta PARCHG-band0
```

### Performance and Resources

**Per Structure**:
- **SPE job**: ~30-60 minutes (depends on system size and PARCHG calculations)
- **BAND job**: ~15-30 minutes
- **PDOS job**: ~15-30 minutes
- **Total time**: ~1-2 hours per structure
- **Disk space**: ~200-500 MB per structure (with PARCHG compression)

**For 50 structures with 10 concurrent**:
- **Wall time**: ~5-10 hours (5 batches × 1-2 hours/batch)
- **Total core-hours**: 50 × 1.5 × 16 × 3 jobs = 3,600 core-hours
- **Disk usage**: 50 × 300 MB = 15 GB

### Troubleshooting

**"CONTCAR not found for structure"**
```bash
# Ensure structure completed refined relaxation
grep "structure_id" REFINE_VASP_JOBS/workflow.json | grep RELAX_DONE

# Check CONTCAR exists
ls -lh REFINE_VASP_JOBS/*/structure_id/Relax/CONTCAR
```

**"SC calculation failed"**
```bash
# Check VASP error
cat ELECTRONIC_JOBS/Ba2N_s001/SPE/vasp_*.err

# Check electronic convergence
grep "reaching required accuracy" ELECTRONIC_JOBS/Ba2N_s001/SPE/OUTCAR
```

**"PARCHG generation failed"**
```bash
# Check if vasprun.xml-SC is valid
grep "finalpos" ELECTRONIC_JOBS/Ba2N_s001/SPE/vasprun.xml-SC

# Check generate_parchg_incars.py output
cat ELECTRONIC_JOBS/Ba2N_s001/SPE/vasp_*.out | grep "PARCHG"
```

**"Band structure calculation stalled"**
```bash
# Check if POSCAR exists (primitive cell)
ls -lh ELECTRONIC_JOBS/Ba2N_s001/BAND/POSCAR

# Verify k-path was generated
cat ELECTRONIC_JOBS/Ba2N_s001/BAND/KPOINTS
```

**"ICHARG error in BAND/DOS"**
```bash
# This should not happen if ICHARG=2 is set correctly
# Check INCAR
grep "ICHARG" ELECTRONIC_JOBS/Ba2N_s001/BAND/INCAR
# Should show: ICHARG = 2
```

### Resume Workflow

If interrupted, simply re-run the same command:

```bash
# The workflow manager will:
# - Load existing workflow.json database
# - Resume from current states
# - Skip completed structures
# - Continue monitoring running jobs

bash run_electronic_flow.sh \
    --refine-jobs ./REFINE_VASP_JOBS \
    --output-dir ./ELECTRONIC_JOBS
```

### Prerequisites

1. **Completed refined relaxation** with structures in `RELAX_DONE` state
2. **VASP modules** and conda environment (`vaspflow`)
3. **Disk space**: ~300-500 MB per structure
4. **Computational resources**: CPU nodes (Orion, Apus)

### Best Practices

1. **Test first**: Start with 5-10 structures to verify settings
2. **Monitor disk**: PARCHG files are large before compression
3. **Check convergence**: Verify first completed structure before processing all
4. **Backup database**: Copy `workflow.json` regularly

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
   
4. Refined Relaxation (refine_electrideflow.py)
   → REFINE_VASP_JOBS/
   
5. Electronic Structure Analysis (electronic_flow_manager.py)
   → ELECTRONIC_JOBS/ (SCF/PARCHG/ELF/BAND/PDOS)
   
6. MatterSim Validation (compute_mattersim_e_hull.py)
   → REFINE_VASP_JOBS/mattersim_stability_results.json
   
7. DFT Hull Comparison (compute_dft_e_hull.py)
   → REFINE_VASP_JOBS/dft_stability_results.json
   → REFINE_VASP_JOBS/hull_comparison.json + plots
   
8. Extract Stable Electrides (get_stable_ele_db.py)
   → REFINE_VASP_JOBS/stable_electrides.csv
   → REFINE_VASP_JOBS/stable_electrides.db
   
9. Final Analysis
   - Analyze stable_electrides.csv/db for publication
   - Electronic structure analysis from ELECTRONIC_JOBS/
   - Band structure and DOS plots
   - ELF and PARCHG visualization
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
                              (REFINE_VASP_JOBS/)
                                      ↓
                    ┌─────────────────┴─────────────────┐
                    ↓                                   ↓
          Electronic Structure              Stability Validation
          (ELECTRONIC_JOBS/)                (MatterSim + DFT)
          - SCF/PARCHG/ELF                           ↓
          - Band Structure                  Stable Electrides
          - PDOS                            Database + Analysis
```

---

## MatterSim Energy Above Hull Comparison

After completing the refined relaxation workflow, you can relax the refined structures with MatterSim and compare energy above hull calculations between MatterSim and VASP-DFT.

### Purpose

Two-step validation and comparison workflow:

**Step 1: MatterSim Relaxation** (`compute_mattersim_e_hull.py`)
- Query MP API for GGA/GGA+U reference phases
- Relax MP phases with MatterSim using tight convergence (`fmax=0.001`, `max_steps=800`)
- Relax refined VASP structures with MatterSim (same convergence)
- Compute MatterSim energy_above_hull using MatterSim-relaxed MP phases
- **Consistent convergence**: Both MP phases and structures use `fmax=0.001`, matching refined VASP

**Step 2: DFT Hull Comparison** (`../compute_dft_e_hull.py`)
- Extract VASP-DFT energies from refined relaxations
- Compute DFT energy_above_hull using MP DFT reference phases
- Compare MatterSim vs DFT hulls with statistics and plots
- Quantify agreement between MatterSim and high-precision DFT

### Scripts

**Refined workflow scripts (in `refined_relaxflow/`)**:

1. **`compute_mattersim_e_hull.py`** - Query MP, relax MP phases & structures with MatterSim
   - Input: 
     - `CONTCAR` files from `REFINE_VASP_JOBS/*/Relax/`
     - MP API key (for GGA/GGA+U reference phases)
   - Process:
     - Query MP API for GGA/GGA+U phases
     - Relax MP phases with MatterSim (`fmax=0.001`, `max_steps=800`)
     - Relax refined structures with MatterSim (same convergence)
     - Compute energy_above_hull
   - Output: 
     - `REFINE_VASP_JOBS/mattersim_stability_results.json` (MatterSim energies & hulls)
     - `REFINE_VASP_JOBS/mp_mattersim.json` (MP phases relaxed with MatterSim)

2. **`submit_mattersim_e_hull.sh`** - SLURM job script for MatterSim relaxation

3. **`run_mattersim_e_hull.sh`** - Wrapper to validate inputs and submit MatterSim job

**Existing script from parent directory**:

4. **`../compute_dft_e_hull.py`** - Compute DFT hulls and compare with MatterSim
   - Input: 
     - VASP energies from `REFINE_VASP_JOBS/*/Relax/vasprun.xml` (or OUTCAR)
     - MatterSim results from `REFINE_VASP_JOBS/mattersim_stability_results.json`
   - Process: 
     - Query MP API for DFT reference phases (GGA/GGA+U)
     - Compute DFT energy_above_hull
     - Compare DFT vs MatterSim hulls
   - Output:
     - `REFINE_VASP_JOBS/dft_stability_results.json` (DFT energies & hulls)
     - `REFINE_VASP_JOBS/mp_vaspdft.json` (MP DFT cache)
     - `REFINE_VASP_JOBS/hull_comparison.json` (statistics)
     - `REFINE_VASP_JOBS/hull_comparison_*.png` (plots)

### Usage

**Step 1: Run MatterSim relaxation**

```bash
cd refined_relaxflow/

# Export MP API key
export MP_API_KEY="your_32_character_api_key"

# Submit MatterSim job (uses GPU by default)
bash run_mattersim_e_hull.sh \
    --refine-jobs ./REFINE_VASP_JOBS

# For CPU mode (not recommended, 3x slower)
# bash run_mattersim_e_hull.sh --device cpu

# Monitor job
tail -f mattersim_e_hull_*.out
```

**Step 2: Compute DFT hulls and compare with MatterSim**

After MatterSim job completes, use the existing `run_dft_e_hull.sh` wrapper script from the parent directory:

```bash
cd ..  # Back to main vaspflow directory

# Submit DFT energy_above_hull calculation as SLURM job
bash run_dft_e_hull.sh \
    --vasp-jobs REFINE_VASP_JOBS/ \
    --prescreen-results REFINE_VASP_JOBS/mattersim_stability_results.json \
    --pure-pbe
```

**What this does**:
1. Reads VASP-DFT energies from `REFINE_VASP_JOBS/*/Relax/vasprun.xml` (or OUTCAR)
2. Queries MP API for DFT reference phases (GGA/GGA+U, cached to `mp_vaspdft.json`)
3. Computes DFT energy_above_hull for each refined structure
4. Compares DFT vs MatterSim hulls (since `--prescreen-results` points to MatterSim results)
5. Generates comparison statistics and plots

**Output files**:
- `REFINE_VASP_JOBS/dft_stability_results.json` - VASP-DFT energies and hulls
- `REFINE_VASP_JOBS/mp_vaspdft.json` - MP DFT reference phases (cached)
- `REFINE_VASP_JOBS/hull_comparison.json` - Comparison statistics (correlation, MAE, RMSE, etc.)
- `REFINE_VASP_JOBS/hull_comparison_scatter.png` - Scatter plot (DFT vs MatterSim)
- `REFINE_VASP_JOBS/hull_comparison_residuals.png` - Residual plot

### Key Features

- **Consistent convergence**: Both MP phases and structures relaxed with `fmax=0.001 eV/Å`, `max_steps=800`
- **Tighter than prescreen**: Prescreen used `fmax=0.01`, `max_steps=500`
- **Matching refined VASP**: MatterSim convergence matches refined VASP force tolerance (EDIFFG=-0.001)
- **Independent validation**: MP phases freshly queried and relaxed (not reused from prescreen)
- **Publication-quality comparison**: Correlation, MAE, RMSE, precision/recall metrics, density plots

### Options (run_mattersim_e_hull.sh)

| Option | Description | Default |
|--------|-------------|---------|
| `--refine-jobs` | Refined VASP jobs directory | `./REFINE_VASP_JOBS` |
| `--device` | MatterSim device: `cpu` or `cuda` | `cuda` |
| `--conda-env` | Conda environment name | `mattersim` |

**Environment Variables**:
- `MP_API_KEY` - Materials Project API key (required)

**SLURM Resources**:
- Partition: GPU
- GPUs: 1x (A40/A100/H200/L40S)
- Memory: 128 GB
- Time limit: 2 days

### Prerequisites

1. **Completed refined relaxation workflow** (structures in `RELAX_DONE` state)
   ```bash
   grep -c '"state": "RELAX_DONE"' REFINE_VASP_JOBS/workflow.json
   ```

2. **MatterSim environment**
   ```bash
   conda activate mattersim
   python3 -c "from mattersim.forcefield import MatterSimCalculator"
   ```

3. **Materials Project API key** (for both MatterSim and DFT steps)
   ```bash
   export MP_API_KEY="your_32_character_api_key"
   # Get key from: https://next-gen.materialsproject.org/api
   ```

### Performance

- **MatterSim relaxation (GPU)**: ~10-20 seconds per structure
- **MatterSim relaxation (CPU)**: ~30-60 seconds per structure (not recommended)
- **MP phase relaxation**: Variable depending on number of phases (~10-50 per chemsys)
- **DFT comparison**: 1-5 seconds per structure
- **Typical workflow**: 
  - 50 refined structures across 10 chemsys
  - MP phases: ~200-300 phases total (~30-60 minutes on GPU)
  - Structure relaxation: ~15-20 minutes on GPU
  - DFT comparison: ~2 minutes
  - **Total: ~1-1.5 hours on GPU**

### Example Workflow

```bash
# 1. Complete refined relaxation workflow
bash run_refineflow.sh --vasp-jobs-dir VASP_JOBS/

# Wait for all structures to reach RELAX_DONE state...

# 2. Set MP API key
export MP_API_KEY="your_32_character_api_key"

# 3. Run MatterSim relaxation on refined structures (GPU)
cd refined_relaxflow/
bash run_mattersim_e_hull.sh

# This:
#   - Queries MP API for GGA/GGA+U reference phases
#   - Relaxes MP phases with MatterSim (fmax=0.001, max_steps=800)
#   - Relaxes refined structures with MatterSim (same convergence)
#   - Computes MatterSim energy_above_hull
#
# Output:
#   REFINE_VASP_JOBS/mattersim_stability_results.json  (MatterSim energies & hulls)
#   REFINE_VASP_JOBS/mp_mattersim.json                 (MP phases relaxed with MatterSim)

# 4. After MatterSim job completes, compute DFT hulls and compare
cd ..
bash run_dft_e_hull.sh \
    --vasp-jobs REFINE_VASP_JOBS/ \
    --prescreen-results REFINE_VASP_JOBS/mattersim_stability_results.json \
    --pure-pbe

# This:
#   - Extracts VASP-DFT energies from vasprun.xml
#   - Queries MP API for DFT reference phases (cached to mp_vaspdft.json)
#   - Computes DFT energy_above_hull
#   - Compares DFT vs MatterSim hulls
#   - Generates statistics and plots
#
# Output:
#   REFINE_VASP_JOBS/dft_stability_results.json        (VASP-DFT energies & hulls)
#   REFINE_VASP_JOBS/mp_vaspdft.json                   (MP DFT reference cache)
#   REFINE_VASP_JOBS/hull_comparison.json              (comparison statistics)
#   REFINE_VASP_JOBS/hull_comparison_scatter.png       (scatter plot)
#   REFINE_VASP_JOBS/hull_comparison_residuals.png     (residual plot)

# 5. View comparison results
cat REFINE_VASP_JOBS/hull_comparison.json

# 6. Extract stable electride candidates (both MatterSim and DFT stable)
cd refined_relaxflow/
bash run_get_stable_db.sh \
    --refine-jobs ./REFINE_VASP_JOBS \
    --electride-csv ../electride_candidates.csv

# This:
#   - Filters structures with mattersim_e_hull < 0.005 AND dft_e_hull < 0.005
#   - Loads refined CONTCARs and extracts spacegroups via PyXtal
#   - Merges electride volumes (e0025, e05, e10, band0, band1) from electride_candidates.csv
#   - Merges energy data from mattersim_stability_results.json and dft_stability_results.json
#
# Output:
#   REFINE_VASP_JOBS/stable_electrides.csv  (sorted by spacegroup, descending)
#   REFINE_VASP_JOBS/stable_electrides.db   (PyXtal database with symmetrized structures)
```

---

## Extract Stable Electride Database

After completing both MatterSim and DFT hull calculations, extract the final list of thermodynamically stable electride candidates that are confirmed by both methods.

### Purpose

Filter and consolidate stable electride candidates:
- **Dual validation**: Only includes structures stable by both MatterSim AND VASP-DFT
- **High-precision structures**: Uses refined 3-step relaxed CONTCARs
- **Complete metadata**: Combines electride volumes, energies, and e_hull values
- **PyXtal database**: Symmetrized structures ready for further analysis

### Scripts

1. **`get_stable_ele_db.py`** - Main Python script
   - Reads `hull_comparison.json` for e_hull values
   - Reads `mattersim_stability_results.json` for MatterSim energies
   - Reads `dft_stability_results.json` for VASP-DFT energies
   - Reads `electride_candidates.csv` for volume data (e0025, e05, e10, band0, band1)
   - Loads CONTCARs from `REFINE_VASP_JOBS/*/Relax/`
   - Extracts spacegroups via PyXtal symmetrization
   - Outputs CSV and PyXtal database

2. **`submit_get_stable_db.sh`** - SLURM job script (CPU-only)

3. **`run_get_stable_db.sh`** - Wrapper to validate inputs and submit job

### Usage

```bash
cd refined_relaxflow/

# Basic usage (uses defaults)
bash run_get_stable_db.sh

# Custom settings
bash run_get_stable_db.sh \
    --refine-jobs ./REFINE_VASP_JOBS \
    --electride-csv ../electride_candidates.csv \
    --threshold 0.005
```

### Options (run_get_stable_db.sh)

| Option | Description | Default |
|--------|-------------|---------|
| `--refine-jobs` | Refined VASP jobs directory | `./REFINE_VASP_JOBS` |
| `--electride-csv` | Electride candidates CSV with volume data | `./electride_candidates.csv` |
| `--threshold` | Energy above hull threshold (eV/atom) | `0.005` |
| `--conda-env` | Conda environment name | `vaspflow` |

### Input Files

| File | Source | Contents |
|------|--------|----------|
| `hull_comparison.json` | `run_dft_e_hull.sh` | `mattersim_e_hull`, `dft_e_hull` |
| `mattersim_stability_results.json` | `run_mattersim_e_hull.sh` | `mattersim_energy_per_atom` |
| `dft_stability_results.json` | `run_dft_e_hull.sh` | `vasp_energy_per_atom` |
| `electride_candidates.csv` | `filter_comb_db.py` | `e0025`, `e05`, `e10`, `band0`, `band1` |
| `workflow.json` | `run_refineflow.sh` | Structure paths and states |
| `*/Relax/CONTCAR` | Refined workflow | Relaxed structures |

### Output Files

| File | Description |
|------|-------------|
| `stable_electrides.csv` | CSV sorted by spacegroup (highest first) |
| `stable_electrides.db` | PyXtal database with symmetrized structures |

### CSV Columns

| Column | Description |
|--------|-------------|
| `formula` | Structure ID (e.g., `Ba2N_s001`) |
| `composition` | Chemical composition |
| `e0025` | Interstitial volume at E_F-0.025 eV |
| `e05` | Interstitial volume at E_F-0.5 eV |
| `e10` | Interstitial volume at E_F-1.0 eV |
| `band0` | Interstitial volume for band 0 |
| `band1` | Interstitial volume for band 1 |
| `spacegroup` | Space group number from PyXtal |
| `mattersim_energy_per_atom` | MatterSim energy (eV/atom) |
| `mattersim_e_hull` | MatterSim energy above hull (eV/atom) |
| `vasp_energy_per_atom` | VASP-DFT energy (eV/atom) |
| `dft_e_hull` | VASP-DFT energy above hull (eV/atom) |

### Database Fields

Same as CSV columns, with:
- `space_group_number` (instead of `spacegroup`)
- `symmetrized` (bool: True if PyXtal symmetrization succeeded)
- Full atomic structure stored

### Filtering Criteria

Structures pass if **both**:
- `mattersim_e_hull < threshold` (default: 0.005 eV/atom)
- `dft_e_hull < threshold` (default: 0.005 eV/atom)

### SLURM Resources

- **Partition**: Orion, Apus (CPU-only, no GPU needed)
- **CPUs**: 8 cores
- **Memory**: 16 GB
- **Time limit**: 1 hour

### Example Output

```
======================================================================
Summary
======================================================================
Total matched structures: 150
Stable (both e_hull < 0.005): 42
Successfully processed: 42
Failed to process: 0
Saved to database: 42
Database save failed: 0

Output CSV: REFINE_VASP_JOBS/stable_electrides.csv
Output database: REFINE_VASP_JOBS/stable_electrides.db
======================================================================

Stable Electride Candidates (sorted by spacegroup, descending):
+-------------------+-------------+-------+------+------+-------+-------+------------+-------------------+------------+
| formula           | composition | e0025 | e05  | e10  | band0 | band1 | spacegroup | mattersim_e_hull  | dft_e_hull |
+-------------------+-------------+-------+------+------+-------+-------+------------+-------------------+------------+
| Ca12Al14_s001     | Ca12Al14    | 45.2  | 38.1 | 28.3 | 52.1  | 48.7  |        225 |          0.000000 |   0.001234 |
| Ba2N_s003         | Ba2N        | 32.5  | 28.9 | 22.1 | 35.6  | 33.2  |        166 |          0.000000 |   0.002156 |
...
```

### Troubleshooting

**"mattersim_energy_per_atom: None" errors**
```bash
# Verify MatterSim results file exists and has data
cat REFINE_VASP_JOBS/mattersim_stability_results.json | head -50

# Check if structure is in MatterSim results
grep "structure_id_here" REFINE_VASP_JOBS/mattersim_stability_results.json
```

**"Could not save to database" warnings**
```bash
# Usually means PyXtal symmetrization failed for that structure
# Structure is still saved to CSV, just not to database
# Check if structure has unusual geometry or very low symmetry
```

**Empty output (no stable structures)**
```bash
# Try increasing threshold
bash run_get_stable_db.sh --threshold 0.01

# Or check hull comparison for available structures
jq '.matched_structures | length' REFINE_VASP_JOBS/hull_comparison.json
jq '[.matched_structures[] | select(.mattersim_e_hull < 0.01 and .dft_e_hull < 0.01)] | length' REFINE_VASP_JOBS/hull_comparison.json
```

---

### Troubleshooting

**MP API key not set**
```bash
# Check if set
echo $MP_API_KEY

# Set it
export MP_API_KEY="your_32_character_api_key"

# Verify it's valid (should be 32 characters)
echo ${#MP_API_KEY}
```

**No RELAX_DONE structures**
```bash
# Check workflow status
grep -c '"state": "RELAX_DONE"' REFINE_VASP_JOBS/workflow.json

# If zero, complete refined workflow first
bash run_refineflow.sh
```

**MP query/relaxation failures**
```bash
# Check log for specific errors
grep "ERROR" mattersim_e_hull_*.out
grep "Warning: Failed to relax" mattersim_e_hull_*.out

# Common causes:
# - Network issues (retry)
# - MP API rate limiting (wait and retry)
# - MatterSim relaxation failures for specific phases (acceptable, will skip)

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


