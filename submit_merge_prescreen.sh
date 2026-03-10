#!/bin/bash
#SBATCH --job-name=merge_prescreen
#SBATCH --partition=Apus,Orion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=merge_prescreen_%j.out
#SBATCH --error=merge_prescreen_%j.err

# Merge parallel prescreening batch results
# Runs on CPU nodes: merges JSON + PyXtal DB + post-relaxation deduplication

set -e

# Default parameters
OUTPUT_DIR=${OUTPUT_DIR:-"./VASP_JOBS"}
CONDA_ENV=${CONDA_ENV:-"mattersim"}
SKIP_DATABASE=${SKIP_DATABASE:-""}
SKIP_DEDUP=${SKIP_DEDUP:-""}
CLEAN_BATCHES=${CLEAN_BATCHES:-""}
OUTPUT_FILE=${OUTPUT_FILE:-""}

echo "========================================================================"
echo "Merge Prescreening Batch Results"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo ""

# Load environment
source ~/.bashrc
conda activate $CONDA_ENV

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate conda environment '$CONDA_ENV'"
    exit 1
fi

echo "Conda environment: $CONDA_ENV"
echo "Python: $(which python3)"
echo ""

# Set number of threads for CPU parallelization
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
echo "Parallelization settings:"
echo "  OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "  MKL_NUM_THREADS: $MKL_NUM_THREADS"
echo ""

# Configuration
echo "Configuration:"
echo "  Output directory: $OUTPUT_DIR"
echo "  Skip database merge: $([ -n "$SKIP_DATABASE" ] && echo "yes" || echo "no")"
echo "  Skip post-relax dedup:$([ -n "$SKIP_DEDUP" ] && echo "yes" || echo "no")"
echo "  Clean batch files:    $([ -n "$CLEAN_BATCHES" ] && echo "yes" || echo "no")"
echo "  Output file: ${OUTPUT_FILE:-prescreening_stability.json}"
echo "========================================================================"
echo ""

# Build command
CMD="python3 merge_prescreen_batches.py --output-dir $OUTPUT_DIR"

if [ -n "$SKIP_DATABASE" ]; then
    CMD="$CMD --skip-database"
fi

if [ -n "$SKIP_DEDUP" ]; then
    CMD="$CMD --skip-dedup"
fi

if [ -n "$CLEAN_BATCHES" ]; then
    CMD="$CMD --clean-batches"
fi

if [ -n "$OUTPUT_FILE" ]; then
    CMD="$CMD --output-file $OUTPUT_FILE"
fi

echo "Running: $CMD"
echo ""

$CMD

EXIT_CODE=$?

echo ""
echo "========================================================================"
echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"
echo "========================================================================"

exit $EXIT_CODE
