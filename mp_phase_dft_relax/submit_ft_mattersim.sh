#!/bin/bash
#SBATCH --job-name=ft_msim
#SBATCH --partition=GPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=ft_mattersim_%j.out
#SBATCH --error=ft_mattersim_%j.err

# MatterSim Fine-tuning on VASP-relaxed MP phases (SLURM Job)

set -e

# Default parameters (overridden by environment variables from run_ft_mattersim.sh)
SCRIPT_DIR=${SCRIPT_DIR:-"$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"}
WORKFLOW_DB=${WORKFLOW_DB:-"${SCRIPT_DIR}/mp_relax_workflow.json"}
COMPARISON_JSON=${COMPARISON_JSON:-""}
OUTPUT_DIR=${OUTPUT_DIR:-"${SCRIPT_DIR}/ft_mattersim"}
MODEL=${MODEL:-"mattersim-v1.0.0-5m"}
EPOCHS=${EPOCHS:-300}
BATCH_SIZE=${BATCH_SIZE:-4}
LR=${LR:-2e-5}
VAL_FRACTION=${VAL_FRACTION:-0.1}
SEED=${SEED:-42}
DEVICE=${DEVICE:-"cuda"}
INCLUDE_STRESSES=${INCLUDE_STRESSES:-"--include-stresses"}
DATASET_ONLY=${DATASET_ONLY:-""}
SKIP_FIRST=${SKIP_FIRST:-0}
MAX_FORCE=${MAX_FORCE:-50.0}
TEST_FRACTION=${TEST_FRACTION:-0.1}
PATIENCE=${PATIENCE:-50}
RE_NORMALIZE=${RE_NORMALIZE:-""}
CONDA_ENV=${CONDA_ENV:-"mattersim"}

echo "========================================================================"
echo "MatterSim Fine-tuning (SLURM Job)"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "GPUs: ${SLURM_GPUS:-0} (${CUDA_VISIBLE_DEVICES:-none})"
echo "Start time: $(date)"
echo ""

# Load environment
source ~/.bashrc
conda activate $CONDA_ENV

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to activate conda environment '$CONDA_ENV'"
    exit 1
fi

# Load CUDA module for shared libraries
module load cuda/11.8

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Python: $(which python3)"
echo ""

# Check GPU
if command -v nvidia-smi &> /dev/null; then
    echo "GPU Information:"
    nvidia-smi --query-gpu=index,name,memory.total,memory.free --format=csv,noheader
    echo ""

    if [ "$DEVICE" = "cuda" ]; then
        if ! python3 -c "import torch; assert torch.cuda.is_available()" 2>/dev/null; then
            echo "Warning: CUDA requested but not available, falling back to CPU"
            DEVICE="cpu"
        fi
    fi
else
    echo "No GPU detected"
    if [ "$DEVICE" = "cuda" ]; then
        echo "Warning: Falling back to CPU"
        DEVICE="cpu"
    fi
fi

# Set parallelization
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

# Suppress harmless Python warnings (pkg_resources deprecation, torch.load FutureWarning)
export PYTHONWARNINGS="ignore:pkg_resources is deprecated:UserWarning,ignore::FutureWarning"

echo "Device: $DEVICE"
echo "OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo ""

# Check dependencies
if ! python3 -c "import mattersim" 2>/dev/null; then
    echo "ERROR: MatterSim not found in '$CONDA_ENV' environment"
    exit 1
fi

if ! python3 -c "import pymatgen" 2>/dev/null; then
    echo "ERROR: pymatgen not found"
    exit 1
fi

# Expand paths
WORKFLOW_DB=$(eval echo "$WORKFLOW_DB")
OUTPUT_DIR=$(eval echo "$OUTPUT_DIR")

# Check workflow database
if [ ! -f "$WORKFLOW_DB" ]; then
    echo "ERROR: Workflow database not found: $WORKFLOW_DB"
    echo "Run mp_dft_relaxflow.py first, then run_compare_energies.sh"
    exit 1
fi

cd "$SCRIPT_DIR"
echo "Working directory: $(pwd)"
echo ""

# Build command
CMD="python3 ${SCRIPT_DIR}/ft_mattersim_mp.py"
CMD="$CMD --workflow-db $WORKFLOW_DB"
CMD="$CMD --output-dir $OUTPUT_DIR"
CMD="$CMD --model $MODEL"
CMD="$CMD --epochs $EPOCHS"
CMD="$CMD --batch-size $BATCH_SIZE"
CMD="$CMD --lr $LR"
CMD="$CMD --val-fraction $VAL_FRACTION"
CMD="$CMD --seed $SEED"
CMD="$CMD --device $DEVICE"

if [ -n "$COMPARISON_JSON" ] && [ -f "$COMPARISON_JSON" ]; then
    CMD="$CMD --comparison-json $COMPARISON_JSON"
    echo "Using comparison JSON: $COMPARISON_JSON"
fi

if [ -n "$INCLUDE_STRESSES" ]; then
    CMD="$CMD $INCLUDE_STRESSES"
fi

if [ -n "$DATASET_ONLY" ]; then
    CMD="$CMD --dataset-only"
fi

CMD="$CMD --skip-first $SKIP_FIRST --max-force $MAX_FORCE --test-fraction $TEST_FRACTION"
CMD="$CMD --patience $PATIENCE"

if [ -n "$RE_NORMALIZE" ]; then
    CMD="$CMD --re-normalize"
fi

echo "Configuration:"
echo "  Workflow DB: $WORKFLOW_DB"
echo "  Output dir: $OUTPUT_DIR"
echo "  Model: $MODEL"
echo "  Epochs: $EPOCHS"
echo "  Batch size: $BATCH_SIZE"
echo "  Learning rate: $LR"
echo "  Val fraction: $VAL_FRACTION"
echo "  Seed: $SEED"
echo ""

echo "========================================================================"
echo "Starting MatterSim fine-tuning..."
echo "========================================================================"
echo ""
echo "Command: $CMD"
echo ""

$CMD

EXIT_CODE=$?

echo ""
echo "========================================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "Fine-tuning completed successfully!"
    echo "Model saved to: ${OUTPUT_DIR}/model"
else
    echo "ERROR: Fine-tuning failed with exit code $EXIT_CODE"
fi
echo "End time: $(date)"
echo "========================================================================"

exit $EXIT_CODE
