# VASPflow - High-Throughput VASP Workflow Manager

**No MongoDB. No FireWorks. Just VASP + SLURM + JSON + AI Pre-screening.**

Intelligent workflow manager for high-throughput VASP calculations with batch submission, MatterSim pre-screening, and dynamic monitoring. Perfect for electride screening on HPC clusters.

**New**: MatterSim pre-screening saves 40-60% of VASP computation by filtering unstable structures. PARCHG analysis for semiconductors improves electride identification accuracy.

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
 **Batch submission** - submit n structures at a time, auto-submit more as jobs complete  
 **Dynamic monitoring** - check status and submit new jobs automatically  
 **Local JSON database** - no network dependencies, works on HPC  
 **Smart workflow** - Pre-screen → Relax → SC → PARCHG (semiconductors) → ELF  
 **PARCHG analysis** - band-edge charge density for accurate electride identification  
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
# - prescreen.py                    # Pre-screening script
# - run_prescreen.sh                 # Submit pre-screening
# - submit_prescreen.sh              # SLURM script for pre-screening
# - workflow_manager.py              # VASP workflow orchestration
# - run_workflow.sh                  # Submit VASP workflow
# - submit_workflow_manager.sh       # SLURM script for workflow
# - workflow_status.py               # Monitor workflow progress
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
# Submit pre-screening job (fast, runs first)
cd /scratch/$USER/vaspflow

bash run_prescreen.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --max-structures 5

# Monitor pre-screening
squeue -u $USER | grep prescreen
tail -f prescreen_*.out

# When done, check results
cat ./VASP_JOBS/prescreening_stability.json | jq '.summary'
```

### 3. VASP Workflow (Run After Pre-screening)

```bash
# Submit VASP workflow manager (processes structures that passed pre-screening)
cd /scratch/$USER/vaspflow

bash run_workflow.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --max-concurrent 10 \
  --max-structures 5

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
  --output electride_candidates.csv

# Monitor analysis job
squeue -u $USER | grep analyze
tail -f electride_analysis_*.out

# View results when complete
cat electride_candidates.csv
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
# - prescreen.py                    # Pre-screening script
# - run_prescreen.sh                 # Submit pre-screening
# - submit_prescreen.sh              # SLURM script for pre-screening
# - workflow_manager.py              # VASP workflow
# - run_workflow.sh                  # Submit VASP workflow
# - submit_workflow_manager.sh       # SLURM script for workflow
# - workflow_status.py               # Status monitoring
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
  --max-structures 5

# Monitor pre-screening
squeue -u $USER | grep prescreen
tail -f prescreen_*.out

# Wait for pre-screening to complete (~1-2 hours for 100 structures)
# Check results:
cat ./VASP_JOBS/prescreening_stability.json | jq '.summary'

# 3. Submit VASP workflow manager (processes passing structures)
bash run_workflow.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --max-concurrent 20 \
  --max-structures 5

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
│   │   ├── SC/
│   │   │   ├── POSCAR (copied from Relax/CONTCAR)
│   │   │   ├── CHGCAR, WAVECAR (after completion)
│   │   │   ├── vasprun.xml (band gap info)
│   │   │   └── ...
│   │   ├── PARCHG/                        ← Band-decomposed charge density
│   │   │   ├── band0/                     ← Valence band maximum (VBM)
│   │   │   │   ├── PARCHG-band0
│   │   │   │   └── VASP_DONE
│   │   │   ├── band1/                     ← VBM-1 (one band below VBM)
│   │   │   ├── e0025/                     ← E-window: [-0.025, 0] eV
│   │   │   ├── e05/                       ← E-window: [-0.5, -0.475] eV
│   │   │   └── e10/                       ← E-window: [-1.0, -0.975] eV
│   │   └── ELF/
│   │       ├── CHGCAR (copied from SC)
│   │       ├── WAVECAR (copied from SC)
│   │       ├── ELFCAR (after completion)
│   │       ├── PARCHG-band0 (copied from PARCHG/band0/)
│   │       ├── PARCHG-band1 (copied from PARCHG/band1/)
│   │       ├── PARCHG-e0025 (copied from PARCHG/e0025/)
│   │       ├── PARCHG-e05 (copied from PARCHG/e05/)
│   │       ├── PARCHG-e10 (copied from PARCHG/e10/)
│   │       ├── BCF.dat, ACF.dat (Bader analysis outputs)
│   │       └── ...
│   └── Li10B1N3_s002/
├── workflow.json                          # Job tracking database
├── prescreening_stability.json            ← NEW: Pre-screening results
├── mp_mattersim_cache/                    ← NEW: Cached MP structures
│   ├── mp_B-Li-N_mattersim.json
│   └── ...
└── electride_candidates.json              # Analysis results
```

### Workflow States

Each structure progresses through these states in workflow_manager.py:

```
PENDING (structure passed pre-screening, ready for VASP)
   ↓
   RELAX_RUNNING → RELAX_DONE → SC_RUNNING → SC_DONE → PARCHG_RUNNING (5 parallel jobs)
      ↓               ↓             ↓            ↓              ↓
   RELAX_FAILED   (skip)      SC_FAILED    (skip)      PARCHG_DONE/SKIPPED
                                                              ↓
                                                        ELF_RUNNING
                                                              ↓
                                                        ELF_DONE (success)
                                                              ↓
                                                        ELF_FAILED
```

**Note**: Pre-screening (MatterSim + MP hull analysis) is done separately by `prescreen.py`. Only structures that pass pre-screening (E_hull < threshold) are loaded into workflow_manager.py.

**Key behaviors**:
- Only structures that passed pre-screening are processed
- Only one stage runs at a time per structure (sequential, except PARCHG)
- **PARCHG runs for ALL materials** (metals and semiconductors) following HT-electride methodology
- PARCHG jobs run in parallel (5 jobs: band0, band1, e0025, e05, e10) but tracked as one stage
- PARCHG files are automatically copied to ELF directory before ELF job submission
- Workflow manager automatically submits next stage when previous completes
- Failed jobs don't block other structures
- Max concurrent limit applies to total running structures (any stage)

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

### Submit Analysis Job

**Recommended**: Run analysis as a SLURM job (not on login node):

```bash
# Submit analysis job with defaults
bash submit_analysis.sh

# Custom settings
bash submit_analysis.sh \
  --vasp-jobs ./VASP_JOBS \
  --bader-exe ~/apps/Bader/bader \
  --threshold 0.6 \
  --output electride_candidates.csv

# Monitor analysis job
squeue -u $USER | grep analyze
tail -f electride_analysis_*.out

# View results when complete
cat electride_candidates.csv
```

### What the Analysis Does

For each completed ELF calculation, the workflow:
1. **Reads ELFCAR** - electron localization function from VASP
2. **Reads PARCHG files** (5 types) - band-decomposed partial charge density:
   - `PARCHG-band0`: Valence band maximum (VBM)
   - `PARCHG-band1`: Second highest occupied band (VBM-1)
   - `PARCHG-e0025`: Energy window [-0.025, 0] eV below Fermi level
   - `PARCHG-e05`: Energy window [-0.5, -0.475] eV below Fermi level
   - `PARCHG-e10`: Energy window [-1.0, -0.975] eV below Fermi level
3. **Runs Bader analysis** on each file to identify electron localization maxima
4. **Filters interstitial electrons** - sites far from atoms with high localization
5. **Reports volumes** of interstitial electron density for each type

### Detection Parameters

- `--threshold`: Minimum ELF value for electride detection (default: 0.6)
- `--bader-exe`: Path to Bader executable (default: `bader` in PATH)
- `--output`: Output CSV file name (default: `electride_results.csv`)

### Understanding Results

The analysis generates a **CSV file** with columns:
```
formula,e0025,e05,e10,band0,band1
Li10B1N4_s001,32.62,31.82,31.31,35.11,0.00
Li10B1N4_s004,0.00,0.00,0.00,33.28,0.00
```

**Interpreting volumes** (in Å³):
- **> 0**: Interstitial electrons detected (potential electride)
- **High values across multiple columns**: Strong electride character
- **Only band0 or band1 non-zero**: Band-selective electride
- **All zeros**: Not an electride

**Example cases**:
- `Li10B1N4_s001`: Strong electride (all energy windows + VBM show interstitial character)
- `Li10B1N4_s004`: Band-selective electride (only VBM shows interstitial electrons)

Console output during analysis:
```
Analyzing: Li10B1N4_s001 (VASP_JOBS/Li10B1N4/Li10B1N4_s001/ELF)
  *** ELECTRIDE CANDIDATE *** (e0025=32.62, e05=31.82, e10=31.31, band0=35.11, band1=0.00)

Analyzing: Li10B1N4_s003 (VASP_JOBS/Li10B1N4/Li10B1N4_s003/ELF)
  Not electride (no interstitial electrons)
```

### Direct Analysis (Advanced)

The analysis uses `Electride.py` internally via `analyze.py`. For manual testing:

```bash
# Test on single structure
cd VASP_JOBS/Li10B1N4/Li10B1N4_s001/ELF
python3 ../../../../Electride.py . --bader-exe ~/apps/Bader/bader
```

**Note**: The analysis requires:
- ELFCAR file (generated by VASP with `LELF=True`)
- PARCHG files (5 types, automatically copied to ELF directory by workflow)
- Bader executable (download from: https://theory.cm.utexas.edu/henkelman/code/bader/)

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

---

## Configuration

### Key Parameters

**Pre-screening (prescreen.py)**:
```bash
python prescreen.py \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --max-structures 5 \         # Max structures per composition
  --mp-api-key YOUR_KEY \      # Materials Project API key
  --hull-threshold 0.1 \       # E_hull threshold (eV/atom)
  --device cpu                 # MatterSim device: cpu or cuda
```

**VASP Workflow (workflow_manager.py)**:
```bash
python workflow_manager.py \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --max-concurrent 10 \        # Max structures running simultaneously
  --max-compositions 20 \      # Limit number of compositions (optional)
  --max-structures 5 \         # Max structures per composition
  --check-interval 60          # Status check interval (seconds)
```

**Recommendations**:
- **Pre-screening**:
  - `--hull-threshold`: 0.1 eV/atom filters most unstable structures, use 0.05 for stricter filtering
  - `--device`: Use `cuda` if GPU available for 2-3× faster MatterSim relaxation
  - `--mp-api-key`: Get your free API key at https://next-gen.materialsproject.org/api
- **VASP Workflow**:
  - `--max-concurrent`: 10-20 for typical HPC cluster
  - `--check-interval`: 60s is good balance
  - `--max-compositions`: Omit to process all

### SLURM Settings

Edit `create_slurm_script()` in `workflow_manager.py` (lines 163-176):

```python
script = f"""#!/bin/bash
#SBATCH --partition=Apus      # Change partition
#SBATCH --ntasks=32            # Change cores
#SBATCH --time=1-00:00:00      # Change walltime (1 day)
#SBATCH --exclusive            # Exclusive node access
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
- **SC**: 1-3 hours  
- **PARCHG** (all materials): 5 × 30-60 min (parallel) = ~1 hour total
- **ELF**: 1-2 hours
- **Total VASP time**: 6-13 hours/structure

**With pre-screening** (filters ~60% structures):
- **Effective time**: 2-5 hours/structure (accounting for filtered structures)
- **Savings**: 40-60% of total VASP computation time

### Example Workflow (100 structures, max_concurrent=10, with pre-screening)

| Time | Action | Running | Completed | Notes |
|------|--------|---------|-----------|-------|
| T=0 | Pre-screen 100 structures | 0 | 0 | MatterSim + MP |
| T=1h | Pre-screening done | 0 | 0 | 60 pass, 40 filtered out |
| T=1h | Submit first 10 Relax jobs | 10 | 0 | Only stable structures |
| T=5h | 10 Relax done, submit 10 SC | 10 | 0 | |
| T=7h | 10 SC done, 5 semiconductors | 10 | 0 | 5 → PARCHG, 5 → ELF |
| T=8h | PARCHG+ELF done, submit new | 10 | 10 | |
| ... | Repeat | 10 | ... | |
| T=35h | All 60 complete | 0 | 60 | 40 filtered by pre-screening |

**Wall time**: ~35 hours for 100 structures (60 pass pre-screening)  
**Core-hours**: ~18,000 (60 structures × 9h avg × 32 cores / 2 efficiency)  
**Savings**: ~50% compared to running all 100 structures without pre-screening

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

1. **MatterSim Relaxation**: Each structure is relaxed with MatterSim (~10-30s per structure)
2. **MP Phase Diagram**: Query Materials Project for competing phases in the chemical system
3. **Re-relax MP Structures**: All MP structures are re-relaxed with MatterSim for consistent energy reference
4. **Convex Hull**: Build phase diagram using MatterSim energies
5. **Output**: Generate `prescreening_stability.json` with structures that pass threshold
6. **VASP Workflow Filter**: `workflow_manager.py` reads `prescreening_stability.json` and only processes structures with `E_hull < threshold` (default: 0.1 eV/atom)

### Key Features

**Separate Script**:
- Pre-screening is handled by `prescreen.py` (separate from `workflow_manager.py`)
- Run pre-screening first, then start VASP workflow with passing structures
- Cleaner separation of concerns and easier to debug

**Smart Caching**:
- MP data is cached per **chemical system** (e.g., `B-Li-N`), not per formula
- Li10B1N3, Li3B2N1, Li5B3N2 all share the same B-Li-N hull
- Dramatically reduces MP API calls and re-computation
- Cache is saved in `OUTPUT_DIR/mp_mattersim_cache/`

**Consistent Energy Reference**:
- All MP competing phases are re-relaxed with MatterSim
- Ensures apples-to-apples energy comparison
- Avoids issues with mixing DFT and ML potential energies

### Configuration

```bash
# Step 1: Run pre-screening (required first step)
bash run_prescreen.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --max-structures 5

# Custom threshold
bash run_prescreen.sh --hull-threshold 0.05  # Stricter (fewer pass)

# Use GPU for faster pre-screening
bash run_prescreen.sh --device cuda

# Check results
cat ./VASP_JOBS/prescreening_stability.json | jq '.summary'

# Step 2: Run VASP workflow (automatically filters based on pre-screening)
bash run_workflow.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --max-concurrent 10 \
  --max-structures 5

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

**Overhead**:
- ~1-2 hours for 100 structures across 10 chemical systems
- First structure in each chemical system takes longer (MP query + cache)
- Subsequent structures in same system are very fast

**Savings**:
- Filters 50-80% of structures (typical)
- Saves 4-10 VASP hours per filtered structure
- Net savings: 40-60% of total computation time

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
| `prescreen.py` | MatterSim pre-screening (separate step) |
| `run_prescreen.sh` | Submit pre-screening wrapper (user-facing) |
| `submit_prescreen.sh` | SLURM script for pre-screening |
| `workflow_manager.py` | VASP workflow orchestration + PARCHG |
| `run_workflow.sh` | Submit workflow wrapper (user-facing) |
| `submit_workflow_manager.sh` | SLURM script for workflow manager |
| `workflow_status.py` | Status checking and reporting |
| `analyze.py` | Orchestrates electride analysis (calls Electride.py) |
| `analyze.sh` | SLURM script for analysis job |
| `Electride.py` | Bader analysis on ELFCAR + PARCHG files |
| `submit_analysis.sh` | Submit analysis wrapper (user-facing) |

**Note**: 
- Pre-screening runs separately (`prescreen.py`) before VASP workflow (`workflow_manager.py`)
- Analysis uses `analyze.py` to orchestrate `Electride.py` for batch processing
- User-friendly wrapper scripts provide consistent interface across all stages

---

##  Best Practices

1. **Run pre-screening first** with `prescreen.py` to save 40-60% of VASP computation
2. **Wait for pre-screening to complete** before starting VASP workflow
3. **Check pre-screening results** (`prescreening_stability.json`) to verify threshold is appropriate
4. **Submit both jobs as SLURM jobs** (required on most HPC systems)
5. **Monitor from login node** using `workflow_status.py` instead of viewing logs constantly
6. **Set reasonable max_concurrent** (10-20) to avoid queue congestion
7. **Monitor your quota** on /scratch periodically (pre-screening cache adds ~50MB per chemical system)
8. **Check failed jobs** early to catch systematic issues
9. **Start small** (test with 2-3 structures first using `--max-compositions 2 --max-structures 2`)

---

## Common Workflows

### Test Run (Small Scale)

```bash
# Step 1: Pre-screening
bash run_prescreen.sh \
  --max-compositions 2 \
  --max-structures 2
# Total: 4 structures for testing

# Wait for pre-screening to complete
squeue -u $USER | grep prescreen
cat ./VASP_JOBS/prescreening_stability.json | jq '.summary'

# Step 2: VASP workflow
bash run_workflow.sh \
  --max-compositions 2 \
  --max-structures 2 \
  --max-concurrent 4
```

### Production Run (Full Scale)

```bash
# Step 1: Pre-screening (all structures)
bash run_prescreen.sh --max-structures 5

# Wait for completion, then start VASP
bash run_workflow.sh \
  --max-structures 5 \
  --max-concurrent 20
# Processes all compositions that passed pre-screening
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
  --max-structures 5

# Wait for pre-screening to complete
squeue -u $USER | grep prescreen
cat ./VASP_JOBS/prescreening_stability.json | jq '.summary'

# 4. Start VASP workflow (as SLURM job)
bash run_workflow.sh \
  --results-dir ./mattergen_results/ternary_csp_electrides \
  --output-dir ./VASP_JOBS \
  --max-concurrent 20

# 5. Monitor and analyze
python workflow_status.py ./VASP_JOBS/workflow.json --compositions
bash submit_analysis.sh --vasp-jobs ./VASP_JOBS
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
