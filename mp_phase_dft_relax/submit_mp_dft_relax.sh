#!/bin/bash
#SBATCH --job-name=mp_dft_relax
#SBATCH --partition=Apus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=14-00:00:00
#SBATCH --output=mp_dft_relax_%j.out
#SBATCH --error=mp_dft_relax_%j.err

# MP Phase DFT Relaxation Workflow Manager
# Submits and monitors VASP relaxation jobs for MP phases

set -e

# Default parameters (can be overridden by environment variables)
CACHE_DIR=${CACHE_DIR:-"./mp_cache_structs"}
VASP_JOBS_DIR=${VASP_JOBS_DIR:-"./VASP_JOBS"}
DB_FILE=${DB_FILE:-"./mp_relax_workflow.json"}
MAX_CONCURRENT=${MAX_CONCURRENT:-10}
CHECK_INTERVAL=${CHECK_INTERVAL:-120}
CHEMSYS=${CHEMSYS:-""}

echo "========================================================================"
echo "MP Phase DFT Relaxation Workflow (SLURM Job)"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo ""

# Load environment
source ~/.bashrc
conda activate vaspflow

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate conda environment 'vaspflow'"
    exit 1
fi

echo "Conda environment: vaspflow"
echo "Python: $(which python3)"
echo ""

# Check pymatgen
if ! python3 -c "import pymatgen" 2>/dev/null; then
    echo "Error: pymatgen not found"
    exit 1
fi

# Set environment variables
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Expand paths
CACHE_DIR=$(eval echo "$CACHE_DIR")
VASP_JOBS_DIR=$(eval echo "$VASP_JOBS_DIR")
DB_FILE=$(eval echo "$DB_FILE")

# Check cache directory
if [ ! -d "$CACHE_DIR" ]; then
    echo "Error: Cache directory not found: $CACHE_DIR"
    exit 1
fi

# Create VASP jobs directory
mkdir -p "$VASP_JOBS_DIR"

# Build command
CMD="python3 mp_dft_relaxflow.py"
CMD="$CMD --cache-dir $CACHE_DIR"
CMD="$CMD --vasp-jobs-dir $VASP_JOBS_DIR"
CMD="$CMD --db $DB_FILE"
CMD="$CMD --max-concurrent $MAX_CONCURRENT"
CMD="$CMD --check-interval $CHECK_INTERVAL"

if [ -n "$CHEMSYS" ]; then
    CMD="$CMD --chemsys $CHEMSYS"
fi

# Print configuration
echo "Configuration:"
echo "  Cache dir: $CACHE_DIR"
echo "  VASP jobs dir: $VASP_JOBS_DIR"
echo "  Database: $DB_FILE"
echo "  Max concurrent: $MAX_CONCURRENT"
echo "  Check interval: ${CHECK_INTERVAL}s"
if [ -n "$CHEMSYS" ]; then
    echo "  Chemical system: $CHEMSYS"
fi
echo ""

# Check if resuming
if [ -f "$DB_FILE" ]; then
    echo "Database exists - resuming from previous state"
    echo ""
fi

echo "========================================================================"
echo "Starting workflow manager..."
echo "========================================================================"
echo ""
echo "Command: $CMD"
echo ""

# Run workflow manager
$CMD

EXIT_CODE=$?

echo ""
echo "========================================================================"
echo "Workflow manager finished"
echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"
echo "========================================================================"

exit $EXIT_CODE

