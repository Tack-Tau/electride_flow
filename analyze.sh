#!/bin/bash
#SBATCH --job-name=analyze_electrides
#SBATCH --partition=Apus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --output=electride_analysis_%j.out
#SBATCH --error=electride_analysis_%j.err

# SLURM script to analyze completed ELF calculations for electride candidates

# Configuration (can be overridden by environment variables)
VASP_JOBS_DIR=${VASP_JOBS_DIR:-"VASP_JOBS"}
BADER_EXE=${BADER_EXE:-"$HOME/apps/Bader/bader"}
OUTPUT_CSV=${OUTPUT_CSV:-"electride_analysis.csv"}
PYXTAL_DB=${PYXTAL_DB:-"electride_data.db"}
DFT_STABILITY=${DFT_STABILITY:-"${VASP_JOBS_DIR}/dft_stability_results.json"}
ELF_THRESHOLD=${ELF_THRESHOLD:-0.6}
CONDA_ENV=${CONDA_ENV:-"vaspflow"}

echo "========================================================================"
echo "Electride Analysis Job"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""
echo "Configuration:"
echo "  VASP jobs directory: $VASP_JOBS_DIR"
echo "  Bader executable: $BADER_EXE"
echo "  Output CSV: $OUTPUT_CSV"
echo "  PyXtal database: $PYXTAL_DB"
echo "  DFT stability: $DFT_STABILITY"
echo "  ELF threshold: $ELF_THRESHOLD"
echo "  Conda environment: $CONDA_ENV"
echo ""
echo "Features:"
echo "  - Incremental analysis (skips already-analyzed structures)"
echo "  - Electride criteria: e0025 > 0 AND e05 > 0 AND e10 > 0 AND band0 > 0"
echo "  - Adds e_above_hull from DFT calculations"
echo "  - Adds spacegroup from CONTCAR"
echo "  - Saves to PyXtal database"
echo "========================================================================"
echo ""

# Activate conda environment
source ~/.bashrc
conda activate $CONDA_ENV

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to activate conda environment: $CONDA_ENV"
    exit 1
fi

echo "Conda environment activated: $(which python3)"
echo ""

# Set number of threads for CPU parallelization
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "Parallelization settings:"
echo "  OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "  MKL_NUM_THREADS: $MKL_NUM_THREADS"
echo ""

# Check if VASP_JOBS directory exists
if [ ! -d "$VASP_JOBS_DIR" ]; then
    echo "ERROR: VASP jobs directory not found: $VASP_JOBS_DIR"
    exit 1
fi

# Check if bader executable exists
if [ ! -x "$BADER_EXE" ]; then
    echo "WARNING: Bader executable not found or not executable: $BADER_EXE"
    echo "         Will try 'bader' from PATH"
    BADER_CMD="bader"
else
    BADER_CMD="$BADER_EXE"
fi

# Check if workflow database exists
WORKFLOW_DB="${VASP_JOBS_DIR}/workflow.json"
if [ ! -f "$WORKFLOW_DB" ]; then
    echo "ERROR: Workflow database not found: $WORKFLOW_DB"
    exit 1
fi

# Run analysis script
echo "Running electride analysis..."
echo "Note: Metals use ELFCAR only, Semiconductors use ELFCAR + PARCHG"
echo ""

python3 analyze.py \
    --db "$WORKFLOW_DB" \
    --bader-exe "$BADER_CMD" \
    --threshold "$ELF_THRESHOLD" \
    --output "electride_analysis.csv" \
    --pyxtal-db "$PYXTAL_DB" \
    --dft-stability "$DFT_STABILITY" \
    2>&1 | tee electride_analysis_detailed.log

EXIT_CODE=$?

# Move output files to the specified locations
if [ -f "electride_analysis.csv" ]; then
    mv electride_analysis.csv "${OUTPUT_CSV}"
    echo ""
    echo "Final CSV output: ${OUTPUT_CSV}"
else
    echo "WARNING: No output CSV generated"
fi

if [ -f "electride_data.db" ]; then
    if [ "$PYXTAL_DB" != "electride_data.db" ]; then
        mv electride_data.db "${PYXTAL_DB}"
    fi
    echo "PyXtal database: ${PYXTAL_DB}"
fi

echo ""
echo "========================================================================"
echo "Job completed"
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"
echo "========================================================================"

exit $EXIT_CODE

