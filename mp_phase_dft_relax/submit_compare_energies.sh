#!/bin/bash
#SBATCH --job-name=mp_compare
#SBATCH --partition=Apus,Orion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=compare_energies_%j.out
#SBATCH --error=compare_energies_%j.err

# MP vs VASP Energy Comparison (SLURM Job)
# Compares MP GGA-PBE energies (from JSON) with VASP PBE energies

set -e

# Default parameters (can be overridden by environment variables from run_compare_energies.sh)
SCRIPT_DIR=${SCRIPT_DIR:-"$(pwd)"}
DB_FILE=${DB_FILE:-"./mp_relax_workflow.json"}
OUTPUT_FILE=${OUTPUT_FILE:-"./mp_vasp_comparison.json"}
MP_ENERGIES_FILE=${MP_ENERGIES_FILE:-""}

echo "========================================================================"
echo "MP vs VASP Energy Comparison (SLURM Job)"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo ""

# Load environment
source ~/.bashrc
conda activate vaspflow

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to activate conda environment 'vaspflow'"
    exit 1
fi

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Python: $(which python3)"
echo ""

# Check required packages
if ! python3 -c "import pymatgen" 2>/dev/null; then
    echo "ERROR: pymatgen not found"
    exit 1
fi

if ! python3 -c "import matplotlib" 2>/dev/null; then
    echo "ERROR: matplotlib not found"
    exit 1
fi

# Set parallelization
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "Parallelization settings:"
echo "  OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "  MKL_NUM_THREADS: $MKL_NUM_THREADS"
echo ""

# Expand paths
DB_FILE=$(eval echo "$DB_FILE")
OUTPUT_FILE=$(eval echo "$OUTPUT_FILE")

# Check database file
if [ ! -f "$DB_FILE" ]; then
    echo "ERROR: Database file not found: $DB_FILE"
    exit 1
fi

# Change to script directory
cd "$SCRIPT_DIR"
echo "Working directory: $(pwd)"
echo ""

# Build command
CMD="python3 ${SCRIPT_DIR}/compare_mp_vasp_energies.py"
CMD="$CMD --db $DB_FILE"
CMD="$CMD --output $OUTPUT_FILE"

if [ -n "$MP_ENERGIES_FILE" ]; then
    CMD="$CMD --mp-energies-file $MP_ENERGIES_FILE"
    echo "MP energies file: $MP_ENERGIES_FILE"
else
    echo "MP energies file: Auto-detect"
fi

# Print configuration
echo "Configuration:"
echo "  Script directory: $SCRIPT_DIR"
echo "  Database: $DB_FILE"
echo "  Output: $OUTPUT_FILE"
echo ""

echo "========================================================================"
echo "Running comparison script..."
echo "========================================================================"
echo ""
echo "Command: $CMD"
echo ""

# Run comparison
$CMD

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "========================================================================"
    echo "Comparison completed successfully!"
    echo "========================================================================"
    echo "Results saved to:"
    echo "  JSON: $OUTPUT_FILE"
    echo "  Scatter plot: ${OUTPUT_FILE%.json}_scatter.png"
    echo "  Residual plots: ${OUTPUT_FILE%.json}_residuals.png"
    echo "========================================================================"
else
    echo ""
    echo "========================================================================"
    echo "ERROR: Comparison failed with exit code $EXIT_CODE"
    echo "========================================================================"
fi

echo ""
echo "End time: $(date)"
echo ""

exit $EXIT_CODE

