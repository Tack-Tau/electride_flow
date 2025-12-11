# MP Phase DFT Relaxation Workflow

Validates DFT energies by running VASP relaxations on MP phases used during pre-screening.

## Purpose

During pre-screening, MP competing phases are relaxed with MatterSim. This workflow:
1. Downloads the original MP crystal structures
2. Relaxes them with VASP using the **exact same settings** as the main workflow
3. Compares VASP energies with MP database energies

This validates energy consistency between your VASP setup and MP's DFT methodology.

## Directory Structure

```
mp_phase_dft_relax/
├── get_mp_struct.py            # Download MP structures
├── run_get_mp.sh               # Wrapper to submit get_mp_struct.py
├── submit_get_mp.sh            # SLURM job script for get_mp_struct.py
├── mp_dft_relaxflow.py         # Workflow manager
├── run_mp_dft_relax.sh         # Wrapper to submit workflow manager
├── submit_mp_dft_relax.sh      # SLURM job script for workflow manager
├── compare_mp_vasp_energies.py # Energy comparison analysis
├── run_compare_energies.sh     # Wrapper to submit energy comparison
├── submit_compare_energies.sh  # SLURM job script for energy comparison
├── README.md                   # This file
├── mp_cache_structs/           # Downloaded structures
│   ├── mp-*.cif                # Flat structure (CIF files)
│   └── mp_metadata.json        # Single metadata file for all structures
├── VASP_JOBS/                  # VASP calculations
│   └── mp_B-Li-N/
│       ├── mp-12345/Relax/
│       └── mp-23456/Relax/
├── mp_relax_workflow.json      # Workflow database
├── mp_vasp_comparison.json     # Energy comparison results
├── mp_vasp_comparison_scatter.png    # Scatter plot
└── mp_vasp_comparison_residuals.png  # Residual plots
```

## Workflow Steps

### Step 1: Download MP Structures

Extract MP IDs from cached prescreening files and download structures:

```bash
# Option 1: Submit as SLURM job (recommended for large downloads)
bash run_get_mp.sh

# Option 2: Run interactively (for testing)
python3 get_mp_struct.py \
    --cache-dir ./mp_mattersim_cache \
    --output-dir ./mp_cache_structs
```

**For specific chemical system:**
```bash
# SLURM job
bash run_get_mp.sh --chemsys B-Li-N

# Interactive
python3 get_mp_struct.py \
    --cache-dir ./mp_mattersim_cache \
    --output-dir ./mp_cache_structs \
    --chemsys B-Li-N
```

**For pure PBE filtering:**
```bash
bash run_get_mp.sh --pure-pbe
```

**Output:**
- Structures saved to `mp_cache_structs/mp_*/mp-*/`
- Metadata saved to `mp_cache_structs/mp_metadata.json` (single file for all structures)
- Overlapping MP IDs between chemical systems are automatically detected and reused (symlinks)

### Step 2: Submit VASP Workflow

```bash
# For testing with B-Li-N system (5 concurrent jobs)
bash run_mp_dft_relax.sh --chemsys B-Li-N --max-concurrent 5

# For all chemical systems (10 concurrent jobs)
bash run_mp_dft_relax.sh --max-concurrent 10
```

**What it does:**
1. Loads structures from `mp_cache_structs/`
2. Creates VASP input files (INCAR, POSCAR, POTCAR, KPOINTS)
3. Submits VASP relaxation jobs to SLURM
4. Monitors job completion
5. Extracts final energies

**VASP Settings:**
- **Same as main workflow** (`workflow_manager.py` lines 237-255)
- Spin-polarized (ISPIN=2)
- PBE functional
- High accuracy (PREC=Accurate, EDIFF=1e-6, EDIFFG=-0.01)

### Step 3: Monitor Progress

```bash
# Check SLURM jobs
squeue -u $USER

# Monitor workflow manager log
tail -f mp_dft_relax_*.out

# Check workflow database
cat mp_relax_workflow.json | jq '.structures | to_entries | map({mp_id: .key, state: .value.state})'

# Count by state
cat mp_relax_workflow.json | jq '.structures | to_entries | map(.value.state) | group_by(.) | map({state: .[0], count: length})'
```

### Step 4: Compare MP and VASP Energies

Use the automated comparison script to analyze energy differences:

```bash
# Make sure MP_API_KEY is set
export MP_API_KEY=your_api_key

# Submit as SLURM job (recommended)
bash run_compare_energies.sh

# With debug output
bash run_compare_energies.sh --debug

# Or run interactively (for testing)
python3 compare_mp_vasp_energies.py \
    --db mp_relax_workflow.json \
    --output mp_vasp_comparison.json
```

**What it does:**
1. Fetches MP GGA-PBE DFT energies via API (filters by functional!)
2. Verifies VASP convergence (reads vasprun.xml to check `vr.converged`)
3. Filters out non-converged VASP calculations
4. Detects and removes statistical outliers (using IQR method + absolute threshold)
5. Compares filtered VASP energies with MP GGA-PBE energies
6. Calculates statistics (MAE, RMSE) on clean data
7. Generates scatter plots and residual analysis (outliers excluded)

**Output:**
- `mp_vasp_comparison.json`: Detailed statistics and all structure data
- `mp_vasp_comparison_scatter.png`: Scatter plot (MP vs VASP energies with minor ticks)
- `mp_vasp_comparison_residuals.png`: Residual plot (energy differences vs VASP energy)
- Minimal summary statistics printed to stdout

**Example output:**
```
Analyzing energy differences...
  (Verifying VASP convergence - use --no-convergence-check to disable)
  (Outlier threshold: |Δ| > 0.50 eV/atom)
  WARNING: Skipping mp-12345 - VASP calculation not converged
  WARNING: Excluding outlier mp-67890 (Li2O3): Δ = +11.6250 eV/atom (|Δ| = 11.6250)

  Filtered out 1 non-converged structures
  Filtered out 1 statistical outliers (|Δ| > 0.50 eV/atom or > Q3+3*IQR)
  Analyzing 148 structures after filtering

======================================================================
MP vs VASP Energy Comparison (Pure GGA-PBE Functional)
======================================================================
Structures analyzed: 148
Skipped (not converged): 1
Skipped (statistical outliers): 1
MAE:  0.0234 eV/atom
RMSE: 0.0312 eV/atom
Mean: +0.0023 eV/atom
======================================================================

Generating comparison plots...
Scatter plot saved to: mp_vasp_comparison_scatter.png
Residual plot saved to: mp_vasp_comparison_residuals.png

======================================================================
Analysis complete!
  JSON results: mp_vasp_comparison.json
  Scatter plot: mp_vasp_comparison_scatter.png
  Residual plot: mp_vasp_comparison_residuals.png
======================================================================
```

**Plots:**
- **Scatter plot**: Shows MP (Y-axis) vs VASP (X-axis) energies for **converged, non-outlier calculations only** with perfect agreement line (y=x), ±0.05 eV/atom error bands, and minor ticks for precision. Title shows number of structures analyzed and number excluded (non-converged + outliers).
- **Residual plot**: Shows energy residuals (VASP - MP) vs VASP energy for **converged, non-outlier calculations only** with focused Y-axis range, reference lines at zero and ±0.05 eV/atom, mean residual line, and minor ticks. Automatically scales based on filtered data for clear visualization.

**Manual inspection:**

The workflow database (`mp_relax_workflow.json`) contains:

```json
{
  "structures": {
    "mp-12345": {
      "state": "COMPLETED",
      "mp_id": "mp-12345",
      "chemsys": "B-Li-N",
      "formula": "Li3BN2",
      "relax_dir": "VASP_JOBS/mp_B-Li-N/mp-12345/Relax",
      "slurm_id": "7891011",
      "mp_energy_per_atom": null,
      "vasp_energy_per_atom": -5.257
    }
  }
}
```

Note: `mp_energy_per_atom` is fetched by the comparison script, not stored in the workflow database.

## Command Line Options

### get_mp_struct.py

```bash
python3 get_mp_struct.py [OPTIONS]

OPTIONS:
  --cache-dir DIR       Directory with mp_*_mattersim.json files
                        (default: ./mp_mattersim_cache)
  --output-dir DIR      Output directory for structures
                        (default: ./mp_cache_structs)
  --mp-api-key KEY      Materials Project API key
                        (default: from MP_API_KEY environment)
  --chemsys CHEMSYS     Process specific chemical system (e.g., 'B-Li-N')
                        (default: process all)
  --pure-pbe            Filter to pure GGA-PBE only (exclude PBE+U, R2SCAN, SCAN)
                        Use only if your VASP uses pure PBE without +U
                        (default: accept all functionals - recommended)

Note:
  - Uses uncorrected_energy_per_atom from MP (raw DFT, no corrections)
  - Functional filtering matches compute_dft_e_hull.py behavior
  - Default (no --pure-pbe) is recommended for validation
```

### run_mp_dft_relax.sh

```bash
bash run_mp_dft_relax.sh [OPTIONS]

OPTIONS:
  --cache-dir DIR         Directory with cached MP structures
                          (default: ./mp_cache_structs)
  --vasp-jobs-dir DIR     Directory for VASP jobs
                          (default: ./VASP_JOBS)
  --db FILE               Workflow database file
                          (default: ./mp_relax_workflow.json)
  --max-concurrent N      Maximum concurrent VASP jobs (default: 10)
  --check-interval N      Job status check interval in seconds (default: 300)
  --chemsys CHEMSYS       Process specific chemical system only
  --help                  Show help message

Note: Arguments are passed to the SLURM job via environment variables.
```

### compare_mp_vasp_energies.py

```bash
python3 compare_mp_vasp_energies.py [OPTIONS]

OPTIONS:
  --db FILE                  Workflow database file
                             (default: ./mp_relax_workflow.json)
  --mp-api-key KEY           Materials Project API key
                             (default: from MP_API_KEY environment)
  --output FILE              Output file prefix for results
                             (default: mp_vasp_comparison.json)
                             Also creates: *_scatter.png, *_residuals.png
  --debug                    Print detailed debugging information
                             (default: False)
  --no-convergence-check     Skip VASP convergence verification
                             (default: False - convergence is checked)
  --outlier-threshold FLOAT  Absolute energy difference threshold (eV/atom) 
                             for outlier detection (default: 0.5)

Note:
  - Requires MP_API_KEY environment variable to be set
  - Filters for pure GGA-PBE entries (excludes +U, r2SCAN, SCAN)
  - Uses uncorrected_energy_per_atom field from materials.thermo
  - Automatically excludes non-converged VASP calculations (reads vasprun.xml)
  - Automatically removes statistical outliers (|Δ| > threshold or Q3+3*IQR)
  - Typical MAE should be < 0.05 eV/atom for good agreement
  - Generates scatter plots and residual analysis (outliers excluded)
```

### run_compare_energies.sh (Wrapper Script)

```bash
bash run_compare_energies.sh [OPTIONS]

OPTIONS:
  --db FILE          Workflow database file (default: ./mp_relax_workflow.json)
  --output FILE      Output file for results (default: ./mp_vasp_comparison.json)
  --debug            Enable detailed debugging output
  --help             Show help message

ENVIRONMENT:
  MP_API_KEY         Materials Project API key (required)

EXAMPLES:
  # Basic usage
  bash run_compare_energies.sh

  # With debug output
  bash run_compare_energies.sh --debug

  # Custom database
  bash run_compare_energies.sh --db /path/to/workflow.json

Note:
  - Wrapper script that submits submit_compare_energies.sh to SLURM
  - Automatically checks for MP_API_KEY before submitting
  - Passes configuration via environment variables
```

### submit_compare_energies.sh (SLURM Script)

```bash
# Do not run directly - use run_compare_energies.sh instead

SLURM OPTIONS:
  --job-name=mp_compare
  --partition=Apus,Orion
  --cpus-per-task=8
  --mem=8G
  --time=4:00:00

Note:
  - Activated conda environment: mattersim
  - Requires pymatgen and matplotlib
  - Automatically validates MP_API_KEY
  - Output saved to compare_energies_<jobid>.out
```

## Resuming

The workflow automatically resumes from the database:

```bash
# If workflow manager stops, just rerun
bash run_mp_dft_relax.sh --chemsys B-Li-N

# It will:
# - Skip COMPLETED structures
# - Resubmit FAILED structures as PENDING
# - Continue monitoring RUNNING structures
```

## Troubleshooting

### No structures found

```bash
# Make sure you ran get_mp_struct.py first
ls mp_cache_structs/

# Should see mp_A-B-C directories with mp-* subdirectories
```

### VASP job failed

```bash
# Check VASP error log
cat VASP_JOBS/mp_B-Li-N/mp-12345/Relax/vasp.err

# Common issues:
# - POTCAR missing: Check VASP_PSP_DIR environment variable
# - Out of memory: Increase mem-per-cpu in submit script
# - Not converged: Increase NSW or relax convergence criteria
```

### Reset failed jobs

Edit `mp_relax_workflow.json` and change `"state": "FAILED"` to `"state": "PENDING"`, then resubmit.

## Expected Results

### Energy Differences

With **correct functional filtering** (MP GGA vs VASP PBE), typical energy differences:

| System Type | Expected Δ (eV/atom) | Reason |
|-------------|----------------------|--------|
| Non-magnetic | 0.001-0.02 | Numerical precision, k-points |
| Magnetic | 0.01-0.05 | Magnetic ordering, ISPIN settings |
| +U systems | 0.05-0.1 | If MP used +U and your VASP doesn't |

**Large differences (>0.1 eV/atom) after functional filtering:**
- Check ISPIN settings (should match MP: typically ISPIN=2)
- Check if MP used +U corrections (your VASP may need matching +U)
- Check POTCAR versions (MP uses specific PAW versions)
- Verify structure is identical (check POSCAR vs MP CIF)

**CRITICAL:** If you see huge differences (~5 eV/atom), you're likely comparing different functionals (e.g., VASP PBE vs MP r2SCAN). The script now filters for GGA automatically.

### Validation Criteria

  **Excellent:** Δ < 0.02 eV/atom - Setup perfectly matches MP
  **Good:** Δ = 0.02-0.05 eV/atom - Normal variation, setup is correct
  **Check:** Δ = 0.05-0.1 eV/atom - May be +U, magnetic, or k-point differences
  **Problem:** Δ > 0.1 eV/atom - Configuration mismatch or missing +U corrections

## Files Generated

| File | Description |
|------|-------------|
| `mp_cache_structs/` | Downloaded MP structures |
| `mp_metadata.json` | Metadata for all downloaded structures |
| `VASP_JOBS/mp_*/` | VASP calculations |
| `mp_relax_workflow.json` | Workflow database |
| `mp_vasp_comparison.json` | Energy comparison statistics (JSON) |
| `mp_vasp_comparison_scatter.png` | Scatter plot: MP vs VASP energies |
| `mp_vasp_comparison_residuals.png` | Residual plot: energy differences vs VASP energy |
| `mp_dft_relax_*.out` | Workflow manager log |
| `mp_dft_relax_*.err` | Workflow manager errors |
| `compare_energies_*.out` | Energy comparison log |
| `compare_energies_*.err` | Energy comparison errors |

## Notes

**Critical: MP Functional Filtering**

MP stores calculations with **multiple functionals** (GGA, r2SCAN, SCAN). For comparison with VASP PBE calculations, you **MUST** filter for GGA entries only!

Example: mp-23309 (ScCl3)
- MP GGA (PBE):    `uncorrected_energy_per_atom = -5.088 eV/atom`   **Use this!**
- MP r2SCAN:       `uncorrected_energy_per_atom = -9.807 eV/atom`   Different functional!
- Your VASP PBE:   `-5.066 eV/atom` → Δ = 0.022 eV/atom (excellent!)

**Implementation (2025-11-20):**
The `compare_mp_vasp_energies.py` script:
1. Uses `materials.thermo` endpoint (stores multiple functionals)
2. Filters for `energy_type='GGA'` entries
3. Uses `uncorrected_energy_per_atom` (raw DFT, no composition corrections)
4. This matches the "Final Energy/Atom" on MP website calculation summary

Without functional filtering, you would compare VASP PBE (-5.066 eV) against MP r2SCAN (-9.807 eV) and incorrectly conclude there's a huge discrepancy!

**Conda Environments:**
- `get_mp_struct.py`: Uses `mattersim` (requires mp-api)
- `mp_dft_relaxflow.py`: Uses `vaspflow` (requires pymatgen)
- `compare_mp_vasp_energies.py`: Uses `mattersim` (requires mp-api, matplotlib)

**Cluster Configuration:**
- CPU-only jobs (no GPU needed)
- Partition: Apus,Orion (same as main workflow)
- Each VASP job: 1 node, 32 cores, exclusive, 24 hours
- Workflow manager: 1 core, 8GB, 72 hours
- MP structure download: 4 cores, 16GB, 6 hours
- Energy comparison: 1 core, 4GB, 2 hours
- Same VASP settings as main workflow for consistency

## Quick Start (B-Li-N Test)

```bash
# Step 1: Download structures
cd mp_phase_dft_relax
python3 get_mp_struct.py --cache-dir ./mp_mattersim_cache --chemsys B-Li-N

# Step 2: Submit workflow
bash run_mp_dft_relax.sh --chemsys B-Li-N --max-concurrent 5

# Step 3: Monitor
tail -f mp_dft_relax_*.out

# Step 4: Compare energies (after completion)
export MP_API_KEY=your_api_key_here
bash run_compare_energies.sh

# Or with debug output
bash run_compare_energies.sh --debug

# Step 5: View results
# View plots
open mp_vasp_comparison_scatter.png
open mp_vasp_comparison_residuals.png

# View JSON statistics
cat mp_vasp_comparison.json | jq '.stats'
cat mp_vasp_comparison.json | jq '.structures[] | select(.abs_diff > 0.05)'

# Expected: MAE should be < 0.05 eV/atom with GGA functional filtering
# Example: mp-23309 should show Δ ≈ 0.02 eV/atom (excellent!)
```

## Integration with Main Workflow

After validation, you can use these results to:
1. Verify energy consistency with `compute_dft_e_hull.py`
2. Identify problematic phases (large energy differences)
3. Validate phase diagram construction
4. Generate correction factors if needed

The validated VASP energies can be compared with `dft_stability_results.json` from the main workflow to ensure consistent DFT energy references.

