# VASPflow - High-Throughput VASP Workflow Manager

**No MongoDB. No FireWorks. Just VASP + SLURM + JSON + AI Pre-screening.**

Intelligent workflow manager for high-throughput VASP calculations with batch submission, MatterSim pre-screening, and dynamic monitoring. Perfect for electride screening on HPC clusters.

## What's New (v3.0)

**Workflow Improvements**:
- **Unified SPE workflow**: 2 SLURM jobs per structure (vs 8)
- **Simplified directories**: Relax/ and SPE/ (vs 7+ subdirectories)
- **Two-step filtering**: Initial analysis → stricter filtering for high-quality candidates
- **Fast symmetry analysis**: PyXtal-based (faster than AFLOW)
- **Query-able database**: `ase db` commands for candidate exploration

**Performance Optimizations**:
- **Batch processing**: MatterSimCalculator reused across structures (3-5× faster prescreening)
- **Parallel prescreening**: Split across multiple GPU nodes (linear speedup with GPUs)
- **Legacy MP API**: Uses pymatgen.ext.matproj.MPRester for complete GGA coverage (modern API misses phases)
- **Strict GGA filtering**: Only accepts entries with `-GGA` or `-GGA+U` suffix (no r2SCAN/SCAN)
- **Uncorrected energies**: Uses raw DFT energies (no MP corrections) for accurate hull calculations
- **Single cache files**: Global `mp_mattersim.json` and `mp_vaspdft.json` with automatic subsystem handling
- **PyXtal symmetrization**: Automatic structure symmetrization for better VASP convergence
- **PDEntry phase diagrams**: More accurate convex hull analysis with decomposition products

---

## Table of Contents

- [Features](#features)
- [Quick Start](#quick-start)
- [Complete Setup Guide](#complete-setup-guide)
- [Workflow Overview](#workflow-overview)
- [Monitoring](#monitoring)
- [Troubleshooting](#troubleshooting)
- [Configuration](#configuration)
- [Performance](#performance)
- [Comparison with FireWorks](#comparison-with-fireworksatomate)
- [MatterSim Pre-screening](#mattersim-pre-screening)
- [SLURM Job Management](#slurm-job-management-reference)

---

##  Features

 **MatterSim pre-screening** - filter structures by energy_above_hull before VASP (saves 40-60% computation)  
 **Batch processing** - reuses MatterSimCalculator across structures for GPU memory efficiency  
 **Parallel prescreening** - split prescreening across multiple GPU nodes for faster throughput  
 **Stable phases only** - queries only on-hull MP phases for accurate, fast convex hull analysis  
 **Smart caching** - single global cache files (mp_mattersim.json, mp_vaspdft.json) handle overlapping systems  
 **PyXtal symmetrization** - automatic structure symmetrization for better VASP convergence  
 **Unified SPE workflow** - SC-PARCHG-ELF in single job (75% fewer SLURM jobs)  
 **Batch submission** - submit n structures at a time, auto-submit more as jobs complete  
 **Dynamic monitoring** - check status and submit new jobs automatically  
 **Local JSON database** - no network dependencies, works on HPC  
 **Smart workflow** - Pre-screen → Relax → SPE (SC + PARCHG + ELF sequential)  
 **PARCHG analysis** - band-edge charge density for accurate electride identification  
 **Two-step filtering** - initial analysis saves all, stricter filtering for candidates  
 **Fast symmetry** - PyXtal-based (14× faster than AFLOW)  
 **Query-able database** - `ase db` commands for candidate exploration  
 **Pymatgen input generation** - INCAR, POSCAR, POTCAR, KPOINTS  
 **Bader analysis** - identify electride candidates with topological analysis  
 **Interrupt-safe** - resume from any interruption  
 **HPC-friendly** - workflow manager runs as SLURM job, no MongoDB needed  

---

##  Quick Start

### 1. Setup Environment

```bash
# On HPC login node
cd /scratch/$USER
mkdir -p vaspflow
cd vaspflow

# Copy files from your local machine:
# - prescreen.py                     # Pre-screening script
# - run_prescreen.sh                 # Submit pre-screening
# - submit_prescreen.sh              # SLURM script for pre-screening
# - plot_e_above_hull.py             # Visualize pre-screening results
# - extract_stable_structs.py        # Extract stable structure CIFs
# - workflow_manager.py              # VASP workflow orchestration
# - run_workflow.sh                  # Submit VASP workflow
# - submit_workflow_manager.sh       # SLURM script for workflow
# - workflow_status.py               # Monitor workflow progress
# - reset_failed_jobs.py             # Reset failed jobs for retry
# - Electride.py                     # Bader analysis
# - bader (executable, optional)

# Create conda environment (with MatterSim for pre-screening)
conda create -n mattersim python=3.10
conda activate mattersim
conda install -c conda-forge pymatgen mp-api
pip install mattersim ase

# Set environment variables in ~/.bashrc
echo 'export PMG_VASP_PSP_DIR=$HOME/apps/PBE64' >> ~/.bashrc
echo 'export MP_API_KEY=your_materials_project_api_key' >> ~/.bashrc
source ~/.bashrc
```

### 2. Pre-screening with MatterSim (NEW Step)

```bash
# Option A: Single GPU job (for small datasets)
cd /scratch/$USER/vaspflow

bash run_prescreen.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --batch-size 32 \
  --max-structures 0

# Option B: Parallel multi-GPU (recommended for large datasets)
bash submit_parallel_prescreen.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --batch-size 32 \
  --hull-threshold 0.05 \
  --max-structures 0 \
  --device cuda \
  --compositions-per-job 400

# Optional: Add --pure-pbe flag to filter MP entries to pure GGA-PBE only
# (default is mixed PBE/PBE+U which is recommended for best accuracy)

# After all parallel jobs complete, merge results (JSON + database)
python3 merge_prescreen_batches.py --output-dir ./VASP_JOBS
# Automatically splits compositions across multiple GPU nodes
# Each batch processes a different subset (batch 0: comps 0-399, batch 1: 400-799, etc.)

# Monitor pre-screening
squeue -u $USER | grep prescreen
tail -f prescreen_*.out

# When done, check results
cat ./VASP_JOBS/prescreening_stability.json | jq '.summary'

# View cache files (single global files with chemsys field)
ls -lh ./VASP_JOBS/mp_mattersim.json
cat ./VASP_JOBS/mp_mattersim.json | jq '.[0]'

# (Optional) Visualize results and extract CIF files
python3 plot_e_above_hull.py --input ./VASP_JOBS/prescreening_stability.json
python3 extract_stable_structs.py \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./stable_structures
```

### 3. VASP Workflow (Run After Pre-screening)

```bash
# Submit VASP workflow manager (processes structures that passed pre-screening)
cd /scratch/$USER/vaspflow

bash run_workflow.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --max-concurrent 10 \
  --max-structures 0

# Check job status
squeue -u $USER
tail -f workflow_manager_<JOBID>.out
```

### 4. Monitor Progress

```bash
cd /scratch/$USER/vaspflow
conda activate mattersim

# Check workflow manager job
squeue -u $USER

# View workflow manager log
tail -f workflow_manager_<JOBID>.out

# Overall workflow status
python workflow_status.py ./VASP_JOBS/workflow.json

# Composition breakdown
python workflow_status.py ./VASP_JOBS/workflow.json --compositions

# Running VASP jobs
python workflow_status.py ./VASP_JOBS/workflow.json --running

# Failed structures
python workflow_status.py ./VASP_JOBS/workflow.json --failed
```

### 5. Analyze Results

```bash
# Submit electride analysis as SLURM job
bash submit_analysis.sh \
  --vasp-jobs ./VASP_JOBS \
  --bader-exe ~/apps/Bader/bader \
  --threshold 0.6 \
  --output electride_analysis.csv

# Monitor analysis job
squeue -u $USER | grep analyze
tail -f electride_analysis_*.out

# When complete, view all analyzed structures
cat electride_analysis.csv

# Filter with stricter criteria to get high-quality candidates
python3 filter_comb_db.py \
  --input electride_data.db \
  --csv electride_analysis.csv \
  --min-energy 20 \
  --min-band 20

# View filtered candidates
cat electride_candidates.csv

# Extract CIF files for candidates
python3 extract_electride_struct.py \
  --db electride_candidates.db \
  --output-dir electride_CIF

# Query database for specific criteria
ase db electride_candidates.db 'e_above_hull<0.05'
```

---

## Complete Setup Guide

### Initial Setup (One-time)

```bash
# 1. SSH to HPC
ssh your_cluster

# 2. Create conda environment with all dependencies
conda create -n mattersim python=3.10
conda activate mattersim
conda install -c conda-forge pymatgen mp-api
pip install mattersim ase

# 3. Setup workspace
cd /scratch/$USER
mkdir -p vaspflow
cd vaspflow

# 4. Upload workflow scripts from your local machine
# - prescreen.py                     # Pre-screening script
# - run_prescreen.sh                 # Submit pre-screening
# - submit_prescreen.sh              # SLURM script for pre-screening
# - plot_e_above_hull.py             # Visualize pre-screening results
# - extract_stable_structs.py        # Extract stable structure CIFs
# - workflow_manager.py              # VASP workflow
# - run_workflow.sh                  # Submit VASP workflow
# - submit_workflow_manager.sh       # SLURM script for workflow
# - workflow_status.py               # Status monitoring
# - reset_failed_jobs.py             # Reset failed jobs for retry
# - Electride.py                     # Bader analysis
# - bader (executable, optional)

# 5. Make scripts executable
chmod +x *.py *.sh

# 6. Set environment variables
echo 'export PMG_VASP_PSP_DIR=$HOME/apps/PBE64' >> ~/.bashrc
echo 'export MP_API_KEY=your_materials_project_api_key' >> ~/.bashrc
source ~/.bashrc

# 7. Verify setup
ls $PMG_VASP_PSP_DIR/POT_GGA_PAW_PBE/Li/POTCAR
echo "MP API Key: $MP_API_KEY"
```

### Starting a Workflow

```bash
# 1. SSH to HPC
ssh your_cluster
cd /scratch/$USER/vaspflow

# 2. Run pre-screening first (MatterSim + MP)
bash run_prescreen.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --max-structures 0

# Monitor pre-screening
squeue -u $USER | grep prescreen
tail -f prescreen_*.out

# Wait for pre-screening to complete (~1-2 hours for 100 structures)
# Check results:
cat ./VASP_JOBS/prescreening_stability.json | jq '.summary'

# 2b. (Optional) Visualize and extract stable structures
python3 plot_e_above_hull.py --input ./VASP_JOBS/prescreening_stability.json
python3 extract_stable_structs.py \
  --json ./VASP_JOBS/prescreening_stability.json \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./stable_structures

# 3. Submit VASP workflow manager (processes passing structures)
bash run_workflow.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --max-concurrent 20 \
  --max-structures 0

# 4. Monitor (optional)
squeue -u $USER                    # Check job status
tail -f workflow_manager_*.out     # View log

# You can now safely log out - workflow runs as SLURM job!
```

### Resuming After Interruption

```bash
# SSH back to HPC
ssh your_cluster
cd /scratch/$USER/vaspflow

# Check if workflow manager is still running
squeue -u $USER | grep vaspflow_manager

# If not running, resubmit (database persists!)
bash run_workflow.sh --max-concurrent 20
# It automatically resumes from saved state in workflow.json
```

---

## Workflow Overview

### Why Run Workflow Manager as SLURM Job?

The workflow manager itself runs as a SLURM job on a compute node (not login node) because:
- **HPC policy**: Login nodes prohibit multi-threaded processes
- **Python threading**: Pymatgen and underlying libraries use multiple threads
- **Long-running**: Workflow can run for days/weeks without SSH connection
- **Lightweight**: Manager uses minimal resources (1 CPU, 2GB RAM)
- **Can submit jobs**: SLURM allows jobs to submit other jobs

### Architecture

```
MatterGen Structures
    ↓
┌──────────────────────────────────────────────────────────────┐
│ Step 1: Pre-screening (prescreen.py + MatterSim)            │
│ - SLURM job on compute node                                 │
│ - Relax structures with MatterSim ML potential              │
│ - Query MP for competing phases (cached per chemical system)│
│ - Compute energy_above_hull                                 │
│ - Output: prescreening_stability.json                       │
│ - Filter: only E_hull < 0.1 eV/atom proceed to VASP        │
└──────────────────────────────────────────────────────────────┘
    ↓
prescreening_stability.json
    ↓
┌──────────────────────────────────────────────────────────────┐
│ Step 2: VASP Workflow (workflow_manager.py)                 │
│ - SLURM job on compute node                                 │
│ - Initialize workflow database (JSON)                       │
│ - Batch submission loop:                                    │
│   * Submit up to max_concurrent VASP jobs                   │
│   * Each structure: Relax → SC → PARCHG* → ELF (sequential)│
│     (* PARCHG only for semiconductors, band gap > 0.01 eV)  │
│   * Monitor job completion every 60s                        │
│   * Submit new structures as slots free up                  │
└──────────────────────────────────────────────────────────────┘
    ↓
Analyze ELF + PARCHG Results (Bader analysis)
    ↓
Identify Electride Candidates
```

### Job Chain Per Structure

1. **Pre-screening** (MatterSim): Quick stability check → `energy_above_hull`
2. **Relax** (VASP): Optimize geometry → `CONTCAR`
3. **SC** (Self-Consistent): Accurate electronic structure → `CHGCAR`, `WAVECAR`, `vasprun.xml`
4. **PARCHG** (all materials): Band-decomposed charge density → `PARCHG-band0`, `PARCHG-band1`, `PARCHG-e0025`, `PARCHG-e05`, `PARCHG-e10`
5. **ELF**: Calculate electron localization → `ELFCAR` (PARCHG files copied here for analysis)

### Directory Structure

```
/scratch/$USER/VASP_JOBS/
├── Li10B1N3/
│   ├── Li10B1N3_s001/
│   │   ├── Relax/
│   │   │   ├── INCAR, POSCAR, POTCAR, KPOINTS
│   │   │   ├── job.sh
│   │   │   ├── CONTCAR (after completion)
│   │   │   ├── VASP_DONE (or VASP_FAILED)
│   │   │   └── vasp_*.out/err
│   │   └── SPE/                           ← Unified SC-PARCHG-ELF directory
│   │       ├── job.sh                     ← Single SLURM job for all stages
│   │       ├── INCAR, POSCAR, POTCAR, KPOINTS
│   │       ├── generate_parchg_incars.py ← Helper script (copied by workflow)
│   │       ├── SC_DONE, PARCHG_DONE       ← Stage markers
│   │       ├── VASP_DONE (or VASP_FAILED) ← Final completion marker
│   │       ├── INCAR-SC, OSZICAR-SC, vasprun.xml-SC  ← SC outputs
│   │       ├── INCAR-PARCHG-e0025, INCAR-PARCHG-e05, INCAR-PARCHG-e10
│   │       ├── INCAR-PARCHG-band0, INCAR-PARCHG-band1
│   │       ├── PARCHG.tar.gz              ← Compressed PARCHG files
│   │       ├── INCAR-ELF
│   │       ├── ELFCAR                     ← ELF output
│   │       ├── ELF_DONE                   ← ELF stage marker
│   │       └── vasp_*.out/err
│   └── Li10B1N3_s002/
├── workflow.json                          # Job tracking database
├── prescreening_stability.json            ← Pre-screening results
├── mp_mattersim.json                      ← Cached MP stable phases (MatterSim energies)
├── mp_vaspdft.json                        ← Cached MP stable phases (DFT energies)
└── electride_candidates.json              # Analysis results
```

### Workflow States

Each structure progresses through these states in workflow_manager.py:

```
PENDING (structure passed pre-screening, ready for VASP)
   ↓
   RELAX_RUNNING → RELAX_DONE → SC_RUNNING → PARCHG_RUNNING → ELF_RUNNING → ELF_DONE
      ↓               ↓             ↓              ↓               ↓             ↓
   RELAX_FAILED   (skip)      SC_FAILED    PARCHG_FAILED   ELF_FAILED    (success)
                                                                                
   Note: SC_RUNNING, PARCHG_RUNNING, and ELF_RUNNING are all part of one SPE job
         State transitions happen automatically as stage markers are detected
```

**Note**: Pre-screening (MatterSim + MP hull analysis) is done separately by `prescreen.py`. Only structures that pass pre-screening (E_hull < threshold) are loaded into workflow_manager.py.

**Key behaviors**:
- Only structures that passed pre-screening are processed
- **Reduced SLURM jobs**: 2 jobs per structure (1 Relax + 1 SPE) instead of 8 jobs
- **SPE job runs sequentially**: SC -> PARCHG (5 energy windows) -> ELF in one script
- **Stage markers**: SC_DONE, PARCHG_DONE, VASP_DONE track progress within SPE job
- **State transitions**: Workflow manager detects markers and transitions states automatically
- **PARCHG for ALL materials**: Following HT-electride methodology
- **ISYM=0**: VASP internal symmetrization disabled (PyXtal already applied)
- Workflow manager automatically submits next structure when slots free up
- Failed jobs don't block other structures
- Max concurrent limit applies to total running structures

### SPE Workflow Architecture

The unified SPE (SC-PARCHG-ELF) workflow uses a clean separation of concerns:

**Python (`workflow_manager.py`) handles**:
- INCAR/POSCAR/POTCAR/KPOINTS for Relax, SC, ELF jobs
- Structure loading from CIF files and symmetrization
- Database state tracking and job orchestration
- Directory structure creation

**Helper script (`generate_parchg_incars.py`) handles**:
- PARCHG INCAR generation (requires parsing vasprun.xml-SC)
- Band index calculation from NELECT and SOC parameters
- FFT grid extraction (NGXF, NGYF, NGZF)
- Called from bash after SC completes

**Bash (`job.sh` in SPE/) handles**:
- Sequential execution: SC → PARCHG (5 types) → ELF
- File copying between stages (CHGCAR-SC, WAVECAR-SC reused)
- PARCHG output renaming (PARCHG → PARCHG-e0025, etc.)
- PARCHG compression (PARCHG-* → PARCHG.tar.gz) to save disk space
- Cleanup between stages to prevent VASP issues
- Stage marker creation (SC_DONE, PARCHG_DONE, ELF_DONE, VASP_DONE)

**Why this design**:
- PARCHG is the ONLY calculation requiring runtime VASP output parsing
- Keeps Python code IDE-friendly (no embedded heredocs)
- Helper script is testable standalone
- Other INCARs are straightforward pymatgen calls (no external scripts needed)

---

## Monitoring

### Option 1: Workflow Status Tool (Recommended)

```bash
# Overall summary
python workflow_status.py ./VASP_JOBS/workflow.json
```

Shows:
- Total structures
- Status breakdown (PENDING, RELAX_RUNNING, SC_RUNNING, etc.)
- Currently running count
- Completion percentage

**Additional views**:

```bash
# Running jobs with SLURM job IDs
python workflow_status.py workflow.json --running

# Failed structures (for debugging)
python workflow_status.py workflow.json --failed

# Composition-level summary
python workflow_status.py workflow.json --compositions

# Recently completed structures
python workflow_status.py workflow.json --completed --limit 20

# Filter by specific state
python workflow_status.py workflow.json --details --state SC_RUNNING
```

### Option 2: Reattach to Tmux

```bash
# See live workflow output
tmux attach -t vaspflow

# View, then detach: Ctrl+B, then D
```

### Option 3: Check SLURM Queue

```bash
# Check running jobs
squeue -u $USER

# Check job history
sacct -S $(date --date="7 days ago" -I) -u $USER

# Detailed job info
scontrol show job <JOBID>
```

### Option 4: Manual Inspection

```bash
# View database directly
cat ./VASP_JOBS/workflow.json | jq '.structures | to_entries | map({id: .key, state: .value.state})'

# Check specific job directory
cd ./VASP_JOBS/Li10B1N3/Li10B1N3_s001/Relax
ls -lh
tail vasp_*.out
tail OSZICAR
```

---

## Analyzing Results for Electride Candidates

After ELF calculations complete, analyze structures for electride characteristics using Bader topological analysis on both ELFCAR and PARCHG files.

**New in v2.0**: Incremental analysis with MatterSim e_above_hull integration, PyXtal database storage, and stricter candidate filtering.

### Step 1: Initial Analysis (All Structures)

**Recommended**: Run analysis as a SLURM job (not on login node):

```bash
# Submit analysis job with defaults
bash submit_analysis.sh

# Custom settings
bash submit_analysis.sh \
  --vasp-jobs ./VASP_JOBS \
  --bader-exe ~/apps/Bader/bader \
  --threshold 0.6 \
  --output electride_analysis.csv

# Monitor analysis job
squeue -u $USER | grep analyze
tail -f electride_analysis_*.out

# View results when complete
cat electride_analysis.csv
```

**Output**:
- `electride_analysis.csv` - All analyzed structures with interstitial volumes
- `electride_data.db` - PyXtal database with all structures

### Step 2: Filter High-Quality Candidates

Apply stricter volume thresholds to identify high-quality electride candidates:

```bash
# Filter with default thresholds (20 Å³)
python3 filter_comb_db.py \
  --input electride_data.db \
  --csv electride_analysis.csv

# Custom thresholds (stricter)
python3 filter_comb_db.py \
  --input electride_data.db \
  --csv electride_analysis.csv \
  --min-energy 25 \
  --min-band 25
```

**Filtering Criteria**:
- `max(e0025, e05, e10) >= min-energy` (default: 20 Å³)
- `max(band0, band1) >= min-band` (default: 20 Å³)
- Removes duplicate structures
- Adds PyXtal symmetry information (space groups)

**Output**:
- `electride_candidates.db` - Filtered high-quality candidates with symmetry
- `electride_candidates.csv` - CSV table for quick viewing

### Step 3: Extract CIF Files

Extract CIF files for the filtered candidates:

```bash
# Extract from filtered database
python3 extract_electride_struct.py \
  --db electride_candidates.db \
  --output-dir electride_CIF
```

**Output**: `electride_CIF/*.cif` - One CIF file per candidate

### Step 4: Query Database

Use `ase db` to query and filter candidates:

```bash
# View all candidates
ase db electride_candidates.db

# Sort by space group (high symmetry first)
ase db electride_candidates.db \
  formula,space_group_number,e_above_hull \
  -s space_group_number-

# Filter by stability (< 50 meV/atom)
ase db electride_candidates.db 'e_above_hull<0.05'

# Show interstitial volumes
ase db electride_candidates.db \
  formula,e0025,e05,e10,band0,band1,e_above_hull

# Cubic structures (space group >= 195)
ase db electride_candidates.db 'space_group_number>=195'

# High interstitial volume (e0025 > 50 Å³)
ase db electride_candidates.db 'e0025>50'
```

### What the Analysis Does

**Initial Analysis** (`analyze.py`):
1. **Reads ELFCAR** - electron localization function from VASP SPE directory
2. **Extracts PARCHG files** from `PARCHG.tar.gz` (5 types):
   - `PARCHG-band0`: Valence band maximum (VBM)
   - `PARCHG-band1`: Second highest occupied band (VBM-1)
   - `PARCHG-e0025`: Energy window [-0.025, 0.025] eV around Fermi level
   - `PARCHG-e05`: Energy window [-0.5, 0.025] eV
   - `PARCHG-e10`: Energy window [-1.0, 0.025] eV
3. **Runs Bader analysis** on ELFCAR and each PARCHG file
4. **Filters interstitial electrons** - sites far from atoms with high localization
5. **Reports volumes** of interstitial electron density (Å³) for each type
6. **Adds stability data** - reads `e_above_hull` from `prescreening_stability.json`
7. **Extracts spacegroup** - from CONTCAR using PyXtal (adaptive tolerance)
8. **Saves to database** - stores all structures in `electride_data.db`

**Filtering** (`filter_comb_db.py`):
1. **Applies volume thresholds** - only keeps structures with substantial interstitial volumes
2. **Symmetrizes structures** - uses PyXtal with adaptive tolerance (coarse→fine)
3. **Removes duplicates** - identifies and filters redundant structures
4. **Adds space group info** - from PyXtal symmetry analysis
5. **Creates unified database** - `electride_candidates.db` with structures + volumes
6. **Exports CSV** - `electride_candidates.csv` for quick viewing

### Key Features (v2.0)

**Incremental Analysis**:
- Automatically skips already-analyzed structures from CSV and database
- Run repeatedly as new ELF calculations complete without re-processing
- Perfect for monitoring ongoing workflows

**Enhanced Electride Criteria**:
- **Old criteria**: ANY of (e0025 > 0 OR e05 > 0 OR e10 > 0 OR band0 > 0)
- **New criteria**: ALL of (e0025 > 0 AND e05 > 0 AND e10 > 0 AND band0 > 0)
- More stringent to reduce false positives

**DFT Integration**:
- Reads DFT energy_above_hull from `compute_dft_e_hull.py` output
- Results sorted by stability (lowest e_above_hull first)
- Most stable electride candidates rise to the top

**PyXtal Database**:
- Saves all analyzed structures to `electride_data.db`
- Adaptive symmetrization (tolerances: 1e-5 to 0.5)
- P1 fallback ensures 100% save success
- Query-able for advanced analysis

### Analysis Parameters

**analyze.py**:
- `--threshold`: Minimum ELF value for electride detection (default: 0.6)
- `--bader-exe`: Path to Bader executable (default: `bader` in PATH)
- `--output`: Output CSV file name (default: `electride_analysis.csv`)
- `--pyxtal-db`: PyXtal database file (default: `electride_data.db`)
- `--prescreening`: MatterSim stability JSON (default: `./VASP_JOBS/prescreening_stability.json`)
- `--workers`: Number of parallel workers (default: CPU count)

**filter_comb_db.py**:
- `--min-energy`: Minimum for max(e0025, e05, e10) in Å³ (default: 20)
- `--min-band`: Minimum for max(band0, band1) in Å³ (default: 20)
- `--output`: Output database (default: `electride_candidates.db`)
- `--output-csv`: Output CSV (default: `electride_candidates.csv`)

### Understanding Results

**Initial Analysis** (`electride_analysis.csv`):
```
formula,composition,e0025,e05,e10,band0,band1,spacegroup,e_above_hull
Li10B1N4_s001,Li10B1N4,32.62,31.82,31.31,35.11,0.00,166,0.023
Li10B1N4_s002,Li10B1N4,15.23,12.45,18.67,22.34,0.00,139,0.045
Li10B1N4_s003,Li10B1N4,0.00,0.00,0.00,5.28,0.00,12,0.068
```

**Columns**:
- `formula`: Structure ID
- `composition`: Chemical formula
- `e0025, e05, e10`: Interstitial volumes (Å³) from energy windows
- `band0, band1`: Interstitial volumes (Å³) from band-specific charges
- `spacegroup`: PyXtal space group (1 = P1 fallback)
- `e_above_hull`: MatterSim energy above hull (eV/atom)

**Electride Criteria**:
- `(e0025 > 0 OR e05 > 0) AND (e10 > 0 OR band0 > 0)`

**Filtered Candidates** (`electride_candidates.csv`):
```
formula,composition,e0025,e05,e10,band0,band1,spacegroup,e_above_hull
Li10B1N4_s001,Li10B1N4,32.62,31.82,31.31,35.11,0.00,166,0.023
Li10B1N4_s002,Li10B1N4,15.23,12.45,18.67,22.34,0.00,139,0.045
```

**Additional Filtering**:
- `max(e0025, e05, e10) >= 20 Å³` (substantial energy window volume)
- `max(band0, band1) >= 20 Å³` (substantial band-specific volume)
- Duplicates removed
- Space group from PyXtal symmetry analysis

**Interpreting volumes** (in Å³):
- **> 20**: High-quality electride (passes filtering)
- **0-20**: Weak electride (filtered out by default)
- **High across multiple columns**: Strong, robust electride character

**Example interpretation**:
- `Li10B1N4_s001`: **High-quality electride** (all volumes > 20 Å³, stable)
- `Li10B1N4_s002`: **Moderate electride** (volumes > 12 Å³, borderline)
- `Li10B1N4_s003`: **Not electride** (very low volumes)

**Sorting**:
- Both CSVs sorted by `e_above_hull` (low → high)
- Most thermodynamically stable candidates appear first

Console output during analysis:
```
Total ELF_DONE: 100
  Semiconductors (with PARCHG): 15
  Metals (ELFCAR only): 85
  (All structures analyzed with PARCHG in SPE workflow)

Analyzing: Li10B1N4_s001 (VASP_JOBS/Li10B1N4/Li10B1N4_s001/SPE)
  Spacegroup: 166
  *** ELECTRIDE CANDIDATE *** (e0025=32.62, e05=31.82, e10=31.31, band0=35.11, band1=0.00)

Analyzing: Li10B1N4_s003 (VASP_JOBS/Li10B1N4/Li10B1N4_s003/SPE)
  Spacegroup: 139
  Not electride (e0025=0.00, e05=0.00, e10=0.00, band0=5.28, band1=0.00)
```

Filter output:
```
======================================================================
Filtering Electride Database with Stricter Criteria
======================================================================
Criteria:
  max(e0025, e05, e10) >= 20.0 Å³
  max(band0, band1) >= 20.0 Å³

Processing structures...
  Added 100 structures...
  Using P1 fallback for Ba8S7_s024
  Added 200 structures...

Filtering Summary:
  Input structures: 6134
  Filtered out (below thresholds): 4512
  Duplicates removed: 89
  Failed to process: 0
  Added to output database: 1533
```

### Incremental Analysis Workflow

```bash
# Day 1: 50 ELF calculations done
bash submit_analysis.sh
# Output: 50 structures analyzed
#   → electride_analysis.csv (50 structures)
#   → electride_data.db (50 structures)

# Filter to get high-quality candidates
python3 filter_comb_db.py --input electride_data.db --csv electride_analysis.csv
# Output: electride_candidates.db (e.g., 25 high-quality)
#         electride_candidates.csv

# Day 2: 30 more ELF calculations done (total 80)
bash submit_analysis.sh
# Output: Only 30 new analyzed (skips existing 50)
#   → electride_analysis.csv updated (80 total, re-sorted)
#   → electride_data.db updated (80 total)

# Re-filter with updated data
python3 filter_comb_db.py --input electride_data.db --csv electride_analysis.csv
# Output: electride_candidates.db (e.g., 40 high-quality)

# Day 3: 20 more ELF calculations done (total 100)
bash submit_analysis.sh
# Output: Only 20 new analyzed
#   → Final electride_analysis.csv (100 structures)

# Final filtering
python3 filter_comb_db.py --input electride_data.db --csv electride_analysis.csv
# Output: Final electride_candidates.db (e.g., 50 high-quality)

# Extract CIF files for candidates
python3 extract_electride_struct.py \
  --db electride_candidates.db \
  --output-dir electride_CIF
```

**Benefits**:
- Incremental analysis saves time (only analyzes new structures)
- Filtering is fast (~2 min) and can be adjusted without re-running Bader
- Database is query-able at any stage

### Direct Analysis (Advanced)

The analysis uses `Electride.py` internally via `analyze.py`. For manual testing:

```bash
# Test on single structure (SPE workflow)
cd VASP_JOBS/Li10B1N4/Li10B1N4_s001/SPE
python3 ../../../../Electride.py . --bader-exe ~/apps/Bader/bader

# Run full analysis with all options
python3 analyze.py \
  --db VASP_JOBS/workflow.json \
  --bader-exe ~/apps/Bader/bader \
  --threshold 0.6 \
  --output electride_analysis.csv \
  --pyxtal-db electride_data.db \
  --prescreening VASP_JOBS/prescreening_stability.json \
  --workers 32

# Filter high-quality candidates
python3 filter_comb_db.py \
  --input electride_data.db \
  --csv electride_analysis.csv \
  --min-energy 20 \
  --min-band 20
```

**Requirements**:
- ELFCAR file (generated by VASP with `LELF=True` in SPE/)
- PARCHG.tar.gz (5 files compressed, automatically extracted during analysis)
- Bader executable (download from: https://theory.cm.utexas.edu/henkelman/code/bader/)
- (Optional) `prescreening_stability.json` for MatterSim e_above_hull values

### Output Files

**Initial Analysis**:
1. **`electride_analysis.csv`**: All analyzed structures with volumes
2. **`electride_data.db`**: PyXtal database with all structures
3. **`electride_analysis_detailed.log`**: Full console output

**Filtered Candidates**:
1. **`electride_candidates.csv`**: High-quality candidates (CSV table)
2. **`electride_candidates.db`**: PyXtal database with filtered structures + volumes
3. **`electride_CIF/`**: Directory with CIF files for candidates

### Performance

**Analysis** (`analyze.py`):
- **Incremental mode**: Only analyzes new structures (minutes for 30 new structures)
- **Full analysis**: ~5-10 minutes per structure (PARCHG extraction + Bader analysis)
- **Database saves**: ~1 second per structure (with adaptive tolerance)
- **Parallel processing**: Scales with CPU cores (32 cores recommended)

**Filtering** (`filter_comb_db.py`):
- **Time**: ~2 minutes for 6000 structures
- **PyXtal symmetry**: ~0.02 sec per structure (vs AFLOW: ~0.45 sec)
- **Speedup**: 14× faster than AFLOW-based workflow
- **Duplicate removal**: Automatic during processing

---

## DFT Energy Above Hull Calculation

After VASP relaxations complete, compute **DFT-level thermodynamic stability** using your VASP energies combined with Materials Project DFT energies for competing phases.

### Why Compute DFT Hull?

**Pre-screening uses MatterSim** (fast, approximate):
- Target structures: MatterSim energies
- Competing phases: MatterSim energies (re-relaxed from MP)
- Purpose: Filter ~60% of structures before expensive VASP

**DFT hull uses VASP + MP** (accurate, publication-quality):
- Target structures: **VASP-PBE energies** (from your relaxations)
- Competing phases: **MP DFT-PBE energies** (original MP data)
- Purpose: Accurate thermodynamic stability for validation and publication

**Key advantage**: VASP and MP both use GGA-PBE functional with PAW pseudopotentials → energies are directly comparable!

### Usage

```bash
# Submit DFT hull calculation (after VASP relaxations complete)
bash run_dft_e_hull.sh \
  --vasp-jobs ./VASP_JOBS \
  --output dft_stability_results.json \
  --prescreen-results ./VASP_JOBS/prescreening_stability.json

# Monitor job
squeue -u $USER | grep dft_e_hull
tail -f dft_e_hull_*.out

# Hull comparison (MatterSim vs DFT) is automatically generated
# if --prescreen-results is provided to compute_dft_e_hull.py
# Results will be in: VASP_JOBS/hull_comparison.json
# Plots: hull_comparison_scatter.png, hull_comparison_residuals.png
```

### What It Does

1. **Scans for completed VASP relaxations** (RELAX_DONE or later states)
2. **Extracts VASP energies** from `vasprun.xml`
3. **Queries MP for stable phases only** (`is_stable=True`) - much faster than querying all phases
4. **Uses single cache file** (`mp_vaspdft.json`) with automatic subsystem handling
5. **Computes DFT energy_above_hull** for each structure using PDEntry
6. **Outputs** `dft_stability_results.json` with decomposition products

### Output Format

```json
{
  "summary": {
    "total_structures": 100,
    "processed_successfully": 95,
    "failed": 5,
    "energy_reference": "DFT uncorrected (VASP-PBE + MP stable phases only)",
    "mp_phase_selection": "stable_only (is_stable=True)",
    "mp_functional_filter": "mixed_pbe_pbeU"
  },
  "results": [
    {
      "structure_id": "Li10B1N4_s001",
      "composition": "Li10B1N4",
      "chemsys": "B-Li-N",
      "vasp_energy_per_atom": -3.45678,
      "dft_energy_above_hull": 0.023,
      "decomposition_products": "Li3N (0.714) + BN (0.286)",
      "num_mp_stable_phases": 18,
      "error": null
    }
  ]
}
```

### Validation: MatterSim vs DFT Hull Comparison

The `compute_dft_e_hull.py` script automatically validates pre-screening accuracy by comparing MatterSim vs DFT hull values when `--prescreen-results` is provided:

**Statistical metrics**:
- Pearson correlation (how well MatterSim predicts DFT)
- Mean Absolute Error (MAE)
- Root Mean Square Error (RMSE)

**Pre-screening performance**:
- True Positives: Correctly passed (stable in both)
- False Positives: Incorrectly passed (unstable in DFT, wasted VASP)
- True Negatives: Correctly filtered (unstable in both)
- False Negatives: Incorrectly filtered (stable in DFT, missed!)
- Precision, Recall, F1 Score
- Computational savings (% VASP avoided)

**Output**: `hull_comparison.json` with detailed analysis

### When to Use

- **After some VASP relaxations complete**: Validate pre-screening accuracy early
- **Before publication**: Get accurate DFT stability data
- **To identify outliers**: Find structures that were incorrectly filtered/passed

### Performance

- **Time**: ~1-2 hours for 100 structures (mostly MP API queries)
- **Resources**: 8 CPUs, 32GB RAM (for large systems)
- **Caching**: Uses MP IDs from pre-screening cache to minimize API calls

---

## Troubleshooting

### Workflow Manager Not Running

**Symptom**: Jobs complete but database still shows `RUNNING` state

**Check if manager job is running**:
```bash
squeue -u $USER | grep vaspflow_manager
```

**Solution**: Restart workflow manager if stopped
```bash
bash run_workflow.sh --max-concurrent 10
# It will detect completed jobs and update states
```

### Jobs Not Progressing

**Check**:
```bash
# Check workflow manager job
squeue -u $USER | grep vaspflow_manager
tail -f workflow_manager_*.out

# See what's completed
find ./VASP_JOBS -name "VASP_DONE" | wc -l
find ./VASP_JOBS -name "VASP_FAILED" | wc -l

# Check VASP jobs in queue
squeue -u $USER

# Check database state
python workflow_status.py ./VASP_JOBS/workflow.json
```

**Fix**: Restart workflow manager if it stopped
```bash
bash run_workflow.sh --max-concurrent 10
```

### Failed Jobs

```bash
# List all failures
python workflow_status.py workflow.json --failed

# Investigate specific failure
cd ./VASP_JOBS/Li10B1N3/Li10B1N3_s003/Relax
tail -50 vasp_*.err
tail -50 OSZICAR
grep -i "error\|killed\|fail" vasp_*.out
```

**Common failure causes**:
- Walltime limit hit (24 hours default)
- Convergence issues (check OSZICAR)
- Memory issues (check vasp_*.err)
- POTCAR problems (check first few lines of vasp_*.out)

### Resetting Failed Jobs for Retry

Use `reset_failed_jobs.py` to reset failed structures and retry them:

```bash
# Step 1: List all failed structures to understand what went wrong
python3 reset_failed_jobs.py --list

# Output shows failures grouped by stage:
#   RELAX_FAILED: 72 structures
#   SC_FAILED: 4 structures
#   PARCHG_FAILED: 4 structures

# Step 2: Dry run to preview changes (recommended)
python3 reset_failed_jobs.py --dry-run --clean

# Step 3: Reset all failed jobs
python3 reset_failed_jobs.py --clean

# OR reset specific stage only
python3 reset_failed_jobs.py --stage RELAX --clean
python3 reset_failed_jobs.py --stage SC --clean
python3 reset_failed_jobs.py --stage PARCHG --clean

# Step 4: Resume workflow to retry
bash run_workflow.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --max-concurrent 10
```

**What it does**:
- `RELAX_FAILED` → `PENDING` (retry from scratch)
- `SC_FAILED` → `RELAX_DONE` (retry SC using existing relaxed structure)
- `PARCHG_FAILED` → `SC_DONE` (retry PARCHG using existing SC results)
- `--clean` flag removes `VASP_FAILED` markers for clean restart
- Creates automatic backup: `workflow.json.bak`

**Important**: Investigate failure causes before blindly retrying. Some structures may have fundamental issues (bad geometry, unconverged, etc.) and will fail repeatedly.

### POTCAR Not Found

```bash
# Check environment
echo $PMG_VASP_PSP_DIR

# Verify structure
ls $PMG_VASP_PSP_DIR/POT_GGA_PAW_PBE/Li/POTCAR

# Set if needed
export PMG_VASP_PSP_DIR=$HOME/apps/PBE64
# Add to ~/.bashrc to persist
```

### Workflow Manager Job Issues

```bash
# Check if manager job is running
squeue -u $USER | grep vaspflow_manager

# View manager log
tail -100 workflow_manager_*.out

# Cancel stuck manager job
scancel <JOBID>

# Restart (database persists!)
bash run_workflow.sh --max-concurrent 10
```

### Jobs Stuck in SLURM Queue

```bash
# Check job status
squeue -u $USER -o "%.18i %.9P %.30j %.8T %.10M %.6D %R"

# Check why pending
squeue -u $USER -j <JOBID> --start

# Cancel stuck job if needed
scancel <JOBID>
# Workflow will mark it as failed on next check
```

### Database Corruption

```bash
# Backup database
cp workflow.json workflow.json.backup

# Reinitialize (keeps existing directories)
python workflow_manager.py --init-only \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS

# Manually merge completed jobs from backup if needed
```

### Parallel Prescreening Results Missing

**Symptom**: After running `submit_parallel_prescreen.sh`, the `prescreening_stability.json` file is empty or only contains results from one batch.

**Cause**: Multiple parallel jobs writing to the same file cause race conditions and data loss.

**Solution**: 
1. Each batch processes a **different subset of compositions** using `--start-composition` (batch 0: comps 0-399, batch 1: 400-799, etc.)
2. Each batch writes to its own file (`prescreening_stability_batch0.json`, etc.)
3. After all jobs complete, merge them with `merge_prescreen_batches.py`:

```bash
# Check that all batch jobs completed
squeue -u $USER | grep prescreen  # Should be empty

# List batch result files
ls ./VASP_JOBS/prescreening_stability_batch*.json
ls ./VASP_JOBS/prescreening_structures_batch*.db
ls ./VASP_JOBS/mp_mattersim.json  # Shared cache (not batch-specific)

# Merge all batch results (JSON + database)
python3 merge_prescreen_batches.py --output-dir ./VASP_JOBS

# OR skip database merge (faster, JSON only)
python3 merge_prescreen_batches.py --output-dir ./VASP_JOBS --skip-database

# Verify merged results
cat ./VASP_JOBS/prescreening_stability.json | jq '.metadata.total_structures'
cat ./VASP_JOBS/mp_mattersim.json | jq 'length'  # Check shared MP cache
ls -lh ./VASP_JOBS/prescreening_structures.db  # Check database exists
```

**Prevention**: Always run `merge_prescreen_batches.py` after parallel prescreening jobs complete.

---

## Configuration

### Key Parameters

**Pre-screening (prescreen.py)**:
```bash
# Single GPU job
python prescreen.py \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --max-structures 0 \         # Max structures per composition
  --batch-size 32 \            # Structures per batch (calculator reuse)
  --mp-api-key YOUR_KEY \      # Materials Project API key
  --hull-threshold 0.05 \       # E_hull threshold (eV/atom)
  --device cuda \              # MatterSim device: cpu or cuda
  --pure-pbe                   # OPTIONAL: Filter MP to pure GGA-PBE only (exclude PBE+U)

# Parallel multi-GPU (recommended for large datasets)
bash submit_parallel_prescreen.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --batch-size 32 \
  --hull-threshold 0.05 \
  --compositions-per-job 400 \
  --pure-pbe                   # OPTIONAL: Filter MP to pure GGA-PBE only (exclude PBE+U)
# Automatically splits compositions across GPU nodes
```

**VASP Workflow (workflow_manager.py)**:
```bash
python workflow_manager.py \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --max-concurrent 10 \        # Max structures running simultaneously
  --max-compositions 20 \      # Limit number of compositions (optional)
  --max-structures 0 \         # Max structures per composition
  --check-interval 60          # Status check interval (seconds)
```

**Recommendations**:
- **Pre-screening**:
  - `--batch-size`: 32 for A40/A100 GPUs (48GB VRAM), reduce to 16 for smaller GPUs
  - `--hull-threshold`: 0.1 eV/atom filters most unstable structures, use 0.05 for stricter filtering
  - `--device`: Use `cuda` for GPU acceleration (5-10× faster than CPU)
  - `--mp-api-key`: Get your free API key at https://next-gen.materialsproject.org/api
  - `--pure-pbe`: OPTIONAL flag to exclude PBE+U (recommended to OMIT for best accuracy)
    - Default (no flag): Mixed PBE/PBE+U (matches MP phase diagram methodology)
    - With flag: Pure GGA-PBE only (may reduce accuracy for transition metal systems)
  - **MP API**: Uses legacy pymatgen.ext.matproj.MPRester for complete GGA coverage
  - **Large datasets**: Use `submit_parallel_prescreen.sh` to split across multiple GPUs
  - **GPU memory**: Set `export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True` in SLURM script
- **VASP Workflow**:
- `--max-concurrent`: 10-20 for typical HPC cluster
- `--check-interval`: 60s is good balance
- `--max-compositions`: Omit to process all
  - **PyXtal**: Structures automatically symmetrized for better convergence

### SLURM Settings

Edit `create_slurm_script()` in `workflow_manager.py` (lines 163-176):

```python
script = f"""#!/bin/bash
#SBATCH --partition=Apus      # Change partition
#SBATCH --ntasks=32            # Change cores
#SBATCH --time=1-00:00:00      # Change walltime (1 day)
#SBATCH --mem=64G              # Change memoery
#SBATCH --output={job_dir}/vasp_%j.out
#SBATCH --error={job_dir}/vasp_%j.err
```

### VASP Input Settings

Edit `create_vasp_inputs()` in `workflow_manager.py` (lines 129-161):

**Relax** (line 135-142):
```python
vis = MPRelaxSet(structure, user_incar_settings={
    'EDIFF': 1e-6,      # Energy convergence
    'EDIFFG': -0.01,    # Force convergence (eV/Å)
    'NSW': 200,         # Max ionic steps
    'PREC': 'Accurate',
    'ISMEAR': 0,        # Gaussian smearing
    'SIGMA': 0.05,
})
```

**SC/Static** (line 143-148):
```python
vis = MPStaticSet(structure, user_incar_settings={
    'PREC': 'Accurate',
    'LAECHG': True,     # Write core charge
    'LCHARG': True,     # Write CHGCAR
})
```

**ELF** (line 149-156):
```python
vis = MPStaticSet(structure, user_incar_settings={
    'PREC': 'Accurate',
    'LELF': True,       # Calculate ELF
    'ADDGRID': True,    # Finer FFT grid
    'NSW': 0,           # No ionic steps
    'IBRION': -1,
})
```

### Electride Detection Parameters

```bash
python analyze_completed.py \
  --threshold 0.6 \           # ELF threshold (0.5-0.7 typical)
  --min-distance 1.5 \        # Min distance from atoms (Å)
  --volume-threshold 0.5      # Min interstitial volume (Å³)
```

---

## Performance

### Per Structure (typical, 32 cores)

**Pre-screening (MatterSim)**:
- **Relaxation**: 10-30s per structure (CPU) or 5-10s (GPU)
- **MP query**: One-time per chemical system (~2-5 min, cached)
- **Hull analysis**: <1s per structure

**VASP calculations**:
- **Relax**: 2-6 hours
- **SPE job** (SC + PARCHG + ELF sequential): 5-10 hours
  - SC: 1-3 hours
  - PARCHG (5 sequential): 2.5-5 hours (5 × 30-60 min)
  - ELF: 1-2 hours
- **Total VASP time**: 7-16 hours/structure
- **SLURM overhead**: Reduced by 75% (2 jobs instead of 8 jobs per structure)

**With pre-screening** (filters ~60% structures):
- **Effective time**: 3-6 hours/structure (accounting for filtered structures)
- **Savings**: 40-60% of total VASP computation time

### Example Workflow (100 structures, max_concurrent=10, with pre-screening)

| Time | Action | Running | Completed | Notes |
|------|--------|---------|-----------|-------|
| T=0 | Pre-screen 100 structures | 0 | 0 | MatterSim + MP |
| T=1h | Pre-screening done | 0 | 0 | 60 pass, 40 filtered out |
| T=1h | Submit first 10 Relax jobs | 10 | 0 | Only stable structures |
| T=5h | 10 Relax done, submit 10 SPE | 10 | 0 | SPE = SC+PARCHG+ELF |
| T=13h | 10 SPE done, submit new | 10 | 10 | ~8h per SPE job |
| ... | Repeat | 10 | ... | |
| T=53h | All 60 complete | 0 | 60 | 40 filtered by pre-screening |

**Wall time**: ~53 hours for 100 structures (60 pass pre-screening)  
**Core-hours**: ~24,000 (60 structures × 12.5h avg × 32 cores / 2 efficiency)  
**SLURM jobs**: 120 total (vs 480 with old workflow)  
**Savings**: 40-60% computational time from pre-screening, 75% fewer SLURM jobs

---

## Comparison with FireWorks/Atomate

| Feature | VASPflow (this) | FireWorks/Atomate |
|---------|-----------------|-------------------|
| Setup | 4 Python files | 10+ files + MongoDB |
| Dependencies | SLURM + JSON | MongoDB + network |
| Reliability |  High | Network issues |
| Batch control |  Built-in | Manual/complex |
| Monitoring | Python script | Web GUI / database |
| Debugging | Simple logs | Complex traces |
| HPC-friendly |  Login node only | Needs network |
| Resume workflow |  Automatic | Manual restart |
| Scalability | Good (10-200) | Excellent (1000s+) |

**Recommendation**:
- **10-200 structures on HPC**: Use VASPflow (simpler, more reliable)
- **1000+ structures**: Consider FireWorks if you have proper infrastructure

---

## MatterSim Pre-screening

### What is Pre-screening?

Pre-screening uses the MatterSim ML potential to quickly evaluate thermodynamic stability before running expensive VASP calculations. This is done in a **separate step** using `prescreen.py`, which generates `prescreening_stability.json` with structures that pass the stability threshold.

### How It Works

1. **MatterSim Relaxation**: Structures relaxed in batches with calculator reuse for GPU memory efficiency
2. **MP Reference Phases Query**: Uses **legacy pymatgen.ext.matproj.MPRester** for complete GGA coverage
   - The modern mp_api.client misses many stable GGA phases needed for accurate hull calculations
   - Strict filtering: only accepts entries with `-GGA` or `-GGA+U` suffix in entry_id
   - Uses uncorrected_energy (raw DFT, no MP corrections) to match MatterSim energies
   - Optional `--pure-pbe` flag to exclude PBE+U (use pure GGA-PBE only)
3. **Re-relax MP Structures**: MP GGA structures are re-relaxed with MatterSim for consistent energy reference
4. **Convex Hull**: Build phase diagram using MatterSim energies with PDEntry
5. **Output**: Generate `prescreening_stability.json` with structures that pass threshold
6. **VASP Workflow Filter**: `workflow_manager.py` reads `prescreening_stability.json` and only processes structures with `E_hull < threshold` (default: 0.1 eV/atom)

### Key Optimizations (NEW)

**Batch Processing**:
- MatterSimCalculator reused across 32 structures per batch (configurable with `--batch-size`)
- Dramatically reduces GPU memory overhead and initialization time
- Automatic GPU memory cleanup between batches
- Supports CUDA memory expansion: `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True`

**Complete GGA Coverage with Legacy API**:
- Uses **legacy pymatgen.ext.matproj.MPRester** instead of modern mp_api.client
- Modern API misses many stable GGA phases critical for accurate hull calculations
- Queries ALL GGA/GGA+U entries via `get_entries_in_chemsys()`
- **Strict suffix filtering**: Only accepts entries with `-GGA` or `-GGA+U` suffix in entry_id
- Excludes r2SCAN, SCAN, and entries without explicit functional suffix
- Uses **uncorrected_energy** (raw DFT, no MP anion/composition corrections)
- Optional `--pure-pbe` flag to exclude PBE+U (for pure GGA-PBE calculations only)
- More accurate phase diagrams with comprehensive GGA reference phases

**Single Global Cache Files**:
- `mp_mattersim.json`: All MatterSim-relaxed MP phases with `chemsys` field
- `mp_vaspdft.json`: All MP DFT energies with `chemsys` field
- Automatically handles overlapping chemical systems (e.g., B-Li-N includes B, Li, N subsystems)
- No more separate cache directories with duplicate data

**Parallel Multi-GPU Processing**:
- `submit_parallel_prescreen.sh`: Splits compositions across multiple GPU nodes
- Each batch job writes to unique files:
  - JSON: `prescreening_stability_batch0.json`, `prescreening_stability_batch1.json`, ...
  - Database: `prescreening_structures_batch0.db`, `prescreening_structures_batch1.db`, ...
- **MP cache is SHARED** across all batches using dual-level file locking:
  - Single file: `mp_mattersim.json` (all batches read/write safely)
  - **Cache I/O locking**: `mp_mattersim.lock` prevents file corruption
  - **Per-chemsys locking**: `mp_cache_{chemsys}.lock` prevents duplicate MP queries
  - When Batch 0 queries Ba-O, other batches wait for it to finish
  - After Batch 0 completes, other batches read from cache (no duplicate API calls)
  - Lock files auto-cleaned after merge
- After all jobs complete, run `merge_prescreen_batches.py` to combine results
- Command-line configurable (use `--help` for full options)
- Configurable compositions per job (default: 400)
- Linear speedup with number of GPUs
- Example: 2000 compositions → 5 GPU nodes = 5x faster
- **Critical**: File locking ensures safe concurrent access to shared MP cache

**PyXtal Symmetrization**:
- Structures automatically symmetrized before MatterSim relaxation
- Same symmetrization applied to VASP initial structures in workflow_manager.py
- Improves convergence and consistency between pre-screening and VASP

### Configuration

```bash
# Option A: Single GPU job (for small-medium datasets)
bash run_prescreen.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --batch-size 32 \
  --hull-threshold 0.05 \
  --device cuda \
  --max-structures 0

# Option B: Parallel multi-GPU (recommended for large datasets)
bash submit_parallel_prescreen.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --batch-size 32 \
  --hull-threshold 0.05 \
  --max-structures 0 \
  --compositions-per-job 400 \
  --pure-pbe                   # OPTIONAL: Pure GGA-PBE only (omit for best accuracy)

# View all options
bash submit_parallel_prescreen.sh --help

# Monitor parallel jobs
squeue -u $USER | grep prescreen

# After all parallel jobs complete, merge batch results
python3 merge_prescreen_batches.py --output-dir ./VASP_JOBS
# This combines:
#   - prescreening_stability_batch*.json → prescreening_stability.json
#   - prescreening_structures_batch*.db → prescreening_structures.db
#   - mp_mattersim.json is already shared (no merging needed)
# Options:
#   --keep-batches: preserve individual batch files (default: delete after merge)
#   --skip-database: skip database merge, only merge JSON (faster)

# Custom options for single GPU
bash run_prescreen.sh \
  --batch-size 64 \           # Larger batch (requires more GPU memory)
  --hull-threshold 0.05 \     # Stricter threshold (fewer pass)
  --device cuda               # Use GPU (much faster than CPU)

# Check results and cache
cat ./VASP_JOBS/prescreening_stability.json | jq '.summary'
cat ./VASP_JOBS/mp_mattersim.json | jq 'length'  # Number of cached MP phases
cat ./VASP_JOBS/mp_mattersim.json | jq 'group_by(.chemsys) | map({chemsys: .[0].chemsys, count: length})'  # By chemical system

# Step 2: Run VASP workflow (automatically filters based on pre-screening)
bash run_workflow.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --max-concurrent 10 \
  --max-structures 0

# Note: VASP structures are automatically symmetrized with PyXtal before relaxation
# This ensures consistency with pre-screening and improves convergence

# By default, workflow_manager.py looks for:
#   ./VASP_JOBS/prescreening_stability.json
# 
# To use a different file:
bash run_workflow.sh --prescreen-results /path/to/custom_results.json

# To skip filtering (process all structures):
# (Not recommended - remove prescreening_stability.json first)
bash run_workflow.sh  # Will warn if file not found
```

### Performance Impact

**Overhead** (with optimizations):
- **Single GPU**: ~30-60 minutes for 100 structures across 10 chemical systems
- **Parallel (5 GPUs)**: ~10-15 minutes for same dataset (5x speedup)
- **Batch processing**: 3-5x faster than previous per-structure calculator creation
- **Stable phases only**: 10-20x faster MP queries (20-50 phases vs 500+ phases)
- First structure in each chemical system takes longer (MP query + cache)
- Subsequent structures in same system benefit from shared cache

**Cache Efficiency**:
- Single `mp_mattersim.json` file (~1-5 MB) replaces multiple per-chemsys files
- Automatic deduplication of overlapping subsystems
- B-Li-N system includes: B, Li, N, B-Li, B-N, Li-N, B-Li-N (7 subsystems in 1 cache)

**Savings**:
- Filters 50-80% of structures (typical)
- Saves 4-10 VASP hours per filtered structure
- Net savings: 40-60% of total computation time
- **Additional**: Pre-screening itself is 3-5x faster with batch processing

### Checking Pre-screening Results

```bash
# View pre-screening summary
cat VASP_JOBS/prescreening_stability.json | jq '.summary'

# View individual results
cat VASP_JOBS/prescreening_stability.json | jq '.results[] | {
  id: .structure_id,
  chemsys: .chemsys,
  E_hull: .energy_above_hull,
  passed: .passed_prescreening
}' | less

# Count passed vs filtered
cat VASP_JOBS/prescreening_stability.json | jq '.results[] | select(.passed_prescreening==true)' | jq -s 'length'
cat VASP_JOBS/prescreening_stability.json | jq '.results[] | select(.passed_prescreening==false)' | jq -s 'length'
```

### Pre-screening Utilities

#### Visualizing Energy Distribution

Use `plot_e_above_hull.py` to visualize the energy distribution and verify your threshold:

```bash
# Plot histogram of energy_above_hull values
python3 plot_e_above_hull.py \
  --input ./VASP_JOBS/prescreening_stability.json \
  --output histogram.png \
  --bins 100 \
  --max-energy 2.0
```

#### Extracting Stable Structures

Use `extract_stable_structs.py` to extract CIF files for structures that passed pre-screening:

```bash
# Extract CIF files for structures that passed pre-screening
python3 extract_stable_structs.py \
  --json ./VASP_JOBS/prescreening_stability.json \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./stable_structures
```

### Troubleshooting

**All structures filtered out**:
- Threshold may be too strict (try `--hull-threshold 0.2`)
- Check if MP API key is working
- Verify MatterSim is installed correctly

**Pre-screening taking too long**:
- Use GPU: `--device cuda`
- MP cache should build once per chemical system
- Check MP API rate limits (built-in backoff)

---

## SLURM Job Management Reference

### Why Submit Workflow Manager as SLURM Job?

- **Login node restrictions**: HPC login nodes prohibit multi-threaded processes
- **Persistent execution**: Runs independently of SSH session
- **Resource control**: Dedicated CPU and memory allocation
- **Job management**: Standard SLURM tools for monitoring and control

### Basic SLURM Commands

| Action | Command |
|--------|---------|
| Submit workflow manager | `bash run_workflow.sh --max-concurrent 20` |
| Check job queue | `squeue -u $USER` |
| Check specific job | `squeue -j <JOBID>` |
| View job log | `tail -f workflow_manager_<JOBID>.out` |
| Cancel job | `scancel <JOBID>` |
| Job history | `sacct -u $USER` |
| Job details | `scontrol show job <JOBID>` |

### Typical Usage Pattern

```bash
# Day 1: Start workflow
ssh your_cluster
cd /scratch/$USER/vaspflow
bash run_workflow.sh --max-concurrent 20
# Note the job ID, then log out

# Day 2: Check progress (from anywhere)
ssh your_cluster
cd /scratch/$USER/vaspflow

# Check if manager is still running
squeue -u $USER | grep vaspflow_manager

# Check workflow progress
python workflow_status.py ./VASP_JOBS/workflow.json --compositions

# View manager log
tail -f workflow_manager_*.out
```

---

## Files

| File | Purpose |
|------|---------|
| `prescreen.py` | MatterSim pre-screening with batch processing |
| `run_prescreen.sh` | Submit single-GPU pre-screening wrapper |
| `submit_prescreen.sh` | SLURM script for single-GPU pre-screening |
| `submit_parallel_prescreen.sh` | **NEW**: Parallel multi-GPU prescreening across nodes |
| `merge_prescreen_batches.py` | **NEW**: Merge parallel batch results into single file |
| `plot_e_above_hull.py` | Visualize pre-screening energy distribution |
| `extract_stable_structs.py` | Extract CIF files for stable structures |
| `workflow_manager.py` | VASP workflow with PyXtal symmetrization |
| `run_workflow.sh` | Submit workflow wrapper (user-facing) |
| `submit_workflow_manager.sh` | SLURM script for workflow manager |
| `workflow_status.py` | Status checking and reporting |
| `reset_failed_jobs.py` | Reset failed jobs to retry (RELAX/SC/PARCHG) |
| `compute_dft_e_hull.py` | Compute DFT energy_above_hull + MatterSim comparison |
| `run_dft_e_hull.sh` | Submit DFT hull wrapper (user-facing) |
| `submit_dft_e_hull.sh` | SLURM script for DFT hull calculation |
| `analyze.py` | Orchestrates electride analysis (calls Electride.py) |
| `analyze.sh` | SLURM script for analysis job |
| `Electride.py` | Bader analysis on ELFCAR + PARCHG files |
| `submit_analysis.sh` | Submit analysis wrapper (user-facing) |
| `filter_comb_db.py` | **NEW**: Filter candidates with stricter volume criteria + PyXtal symmetry |
| `extract_electride_struct.py` | Extract CIF files from filtered database |
| `generate_parchg_incars.py` | Helper script to generate PARCHG INCARs from vasprun.xml-SC |

**Key Features**:
- **Unified SPE workflow**: SC-PARCHG-ELF run sequentially in one SLURM job (2 jobs per structure vs 8)
- **Simplified directory structure**: Relax/ and SPE/ subdirectories (vs 7+ subdirectories)
- **Stage markers**: SC_DONE, PARCHG_DONE, ELF_DONE, VASP_DONE track progress within SPE job
- **Reduced SLURM overhead**: 75% fewer jobs submitted to queue
- **Clean code separation**: Python for VASP inputs, bash for workflow orchestration
- **Helper script for PARCHG**: Standalone script generates PARCHG INCARs (band indices from vasprun.xml-SC)
- **ISYM=0**: VASP internal symmetrization disabled (PyXtal already applied during relaxation)
- **Batch processing**: `prescreen.py` reuses MatterSimCalculator across structures (3-5× faster)
- **Parallel prescreening**: `submit_parallel_prescreen.sh` splits across GPU nodes (linear speedup)
- **Complete GGA coverage**: Uses legacy pymatgen.ext.matproj.MPRester for comprehensive GGA/GGA+U phases
- **Strict functional filtering**: Only accepts entries with `-GGA` or `-GGA+U` suffix (excludes r2SCAN/SCAN)
- **Uncorrected energies**: Uses raw DFT energies (no MP corrections) for accurate hull calculations
- **Optional pure-pbe**: `--pure-pbe` flag to filter MP entries to pure GGA-PBE only (exclude PBE+U)
- **Single cache files**: `mp_mattersim.json` and `mp_vaspdft.json` in `VASP_JOBS/` directory
- **PyXtal symmetrization**: Automatic structure symmetrization in both prescreening and VASP workflow
- **Two-step electride filtering**: Initial analysis saves all data, stricter filtering for high-quality candidates
- **Fast symmetry analysis**: PyXtal-based filtering (14× faster than AFLOW)
- **Query-able database**: `ase db` commands for filtering candidates by stability, symmetry, or volumes
- **PDEntry for hulls**: Uses `PDEntry` instead of `ComputedEntry` for accurate phase diagram analysis
- `reset_failed_jobs.py` resets failed VASP jobs to retry them without data loss
- **Integrated hull validation**: `compute_dft_e_hull.py` automatically validates pre-screening accuracy with plots and statistics
- User-friendly wrapper scripts provide consistent interface across all stages

---

##  Best Practices

1. **Run pre-screening first** with `prescreen.py` to save 40-60% of VASP computation
2. **Use parallel prescreening** (`submit_parallel_prescreen.sh`) for large datasets (>500 compositions)
3. **Merge batch results** with `merge_prescreen_batches.py` after all parallel jobs complete
4. **Optimize batch size** based on GPU memory: 32 for A40/A100 (48GB), 16 for smaller GPUs
5. **Wait for pre-screening to complete** before starting VASP workflow
6. **Check pre-screening results** (`prescreening_stability.json`) to verify threshold is appropriate
7. **Inspect cache files** (`mp_mattersim.json`, `mp_vaspdft.json`) to verify MP data quality
8. **Submit both jobs as SLURM jobs** (required on most HPC systems)
9. **Monitor from login node** using `workflow_status.py` instead of viewing logs constantly
10. **Set reasonable max_concurrent** (10-20) to avoid queue congestion
11. **Monitor your quota** on /scratch periodically (cache files are now much smaller: ~1-5 MB total)
12. **Check failed jobs** early to catch systematic issues
13. **Start small** (test with 2-3 structures first using `--max-compositions 2 --max-structures 2`)
14. **Trust PyXtal symmetrization** - structures are automatically symmetrized for consistency

---

## Common Workflows

### Production Run (Full Scale)

```bash
# Step 1: Pre-screening

# Option A: Single GPU (for small-medium datasets)
bash run_prescreen.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --device cuda \
  --batch-size 32 \
  --hull-threshold 0.05 \
  --max-structures 0

# Option B: Parallel multi-GPU (recommended for large datasets >500 compositions)
bash submit_parallel_prescreen.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --batch-size 32 \
  --hull-threshold 0.05 \
  --max-structures 0 \
  --compositions-per-job 400

# After all parallel jobs complete, merge results (JSON + database)
python3 merge_prescreen_batches.py --output-dir ./VASP_JOBS
# Auto-splits compositions across multiple GPU nodes
# Each batch processes a different subset (batch 0: comps 0-399, batch 1: 400-799, etc.)
# Monitor with: squeue -u $USER | grep prescreen

# Wait for completion, check and visualize results
cat VASP_JOBS/prescreening_stability.json | jq '.summary'
cat VASP_JOBS/mp_mattersim.json | jq 'group_by(.chemsys) | map({chemsys: .[0].chemsys, count: length})'

jq '.results[] | select(.passed_prescreening == true) | {structure_id, energy_above_hull}' VASP_JOBS/prescreening_stability.json

python3 plot_e_above_hull.py --bins 100

# (Optional) Extract stable structures for inspection
python3 extract_stable_structs.py \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./stable_structures

# Step 2: VASP workflow (structures automatically symmetrized with PyXtal)
bash run_workflow.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --max-concurrent 10 \
  --max-structures 0 \
  --check-interval 120

# Step 3: Analyze electride candidates
bash submit_analysis.sh --vasp-jobs ./VASP_JOBS

# Step 4: Filter high-quality candidates
python3 filter_comb_db.py \
  --input electride_data.db \
  --csv electride_analysis.csv \
  --min-energy 20 \
  --min-band 20

# Step 5: Extract CIF files
python3 extract_electride_struct.py \
  --db electride_candidates.db \
  --output-dir electride_CIF

# Step 6: Query and explore
ase db electride_candidates.db
ase db electride_candidates.db 'e_above_hull<0.05' -s space_group_number-
```

### Resume After Interruption

```bash
# Check if still running
squeue -u $USER | grep vaspflow_manager

# If stopped, restart - database persists
bash run_workflow.sh --max-concurrent 20
```

---

## Integration with MatterGen

```bash
# 1. Generate structures with MatterGen
cd ~/mattergen_test/ternary_electride
bash generate_ternary_csp.sh

# 2. Transfer to HPC scratch
cd /scratch/$USER
cp -r ~/mattergen_test/results/ternary_csp_electrides mattergen_results/

# 3. Run pre-screening (as SLURM job)
cd /scratch/$USER/vaspflow
bash run_prescreen.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --max-structures 0

# Wait for pre-screening to complete
squeue -u $USER | grep prescreen
cat ./VASP_JOBS/prescreening_stability.json | jq '.summary'

# 3b. Visualize and extract stable structures (optional)
python3 plot_e_above_hull.py --input ./VASP_JOBS/prescreening_stability.json
python3 extract_stable_structs.py \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./stable_structures

# 4. Start VASP workflow (as SLURM job)
bash run_workflow.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --max-concurrent 20

# 5. Monitor and analyze
python workflow_status.py ./VASP_JOBS/workflow.json --compositions
bash submit_analysis.sh --vasp-jobs ./VASP_JOBS

# 6. Filter high-quality candidates
python3 filter_comb_db.py \
  --input electride_data.db \
  --csv electride_analysis.csv

# 7. Extract CIF files
python3 extract_electride_struct.py \
  --db electride_candidates.db \
  --output-dir electride_CIF

# 8. Query candidates
ase db electride_candidates.db 'e_above_hull<0.05'
```

---

## Credits & References

**Workflow Design**:
- DopeFlow methodology (simple monitoring)
- HT-electride PARCHG analysis approach
- This implementation by Tony/SOFT team

**Tools & Libraries**:
- [Pymatgen](https://pymatgen.org/) - VASP input generation & MP API
- [MatterSim](https://github.com/microsoft/mattersim) - ML potential for pre-screening
- [Bader](https://theory.cm.utexas.edu/henkelman/code/bader/) - Topological analysis
- [ASE](https://wiki.fysik.dtu.dk/ase/) - Atomic simulation environment
- [Materials Project](https://next-gen.materialsproject.org/) - Competing phase data

**Key Papers**:
- MatterSim: Wang et al., "MatterSim: A Deep Learning Atomistic Model Across Elements, Temperatures and Pressures" (2024)
- Electride identification: Topological analysis of ELF and PARCHG at band edges

---

**Smart, fast, reliable - perfect for HPC electride screening with AI-powered pre-screening!**
