#!/bin/bash
#SBATCH --job-name=prescreen
#SBATCH --partition=Apus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=prescreen_%j.out
#SBATCH --error=prescreen_%j.err

# VASPflow Pre-screening with MatterSim
# Fast thermodynamic stability screening before expensive VASP calculations

set -e

# Default parameters
RESULTS_DIR=${RESULTS_DIR:-"./mattergen_results/ternary_csp_electrides"}
OUTPUT_DIR=${OUTPUT_DIR:-"./VASP_JOBS"}
MAX_COMPOSITIONS=${MAX_COMPOSITIONS:-""}
MAX_STRUCTURES=${MAX_STRUCTURES:-0}
CONDA_ENV=${CONDA_ENV:-"mattersim"}
MP_API_KEY=${MP_API_KEY:-""}
HULL_THRESHOLD=${HULL_THRESHOLD:-0.1}
DEVICE=${DEVICE:-"cuda"}

echo "========================================================================"
echo "VASPflow Pre-screening (MatterSim)"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "GPUs: ${SLURM_GPUS:-0} (${CUDA_VISIBLE_DEVICES:-none})"
echo "Memory: 32 GB"
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

# Check GPU availability and auto-detect device
if command -v nvidia-smi &> /dev/null; then
    echo "GPU Information:"
    nvidia-smi --query-gpu=index,name,memory.total,memory.free --format=csv,noheader
    echo ""
    
    # Auto-detect: use cuda if GPU available, otherwise fallback to cpu
    if [ "$DEVICE" = "cuda" ]; then
        if ! python3 -c "import torch; assert torch.cuda.is_available()" 2>/dev/null; then
            echo "Warning: CUDA requested but not available in PyTorch, falling back to CPU"
            DEVICE="cpu"
        else
            echo "Using CUDA device for MatterSim"
        fi
    fi
else
    echo "No GPU detected (nvidia-smi not found)"
    if [ "$DEVICE" = "cuda" ]; then
        echo "Warning: CUDA requested but no GPU available, falling back to CPU"
        DEVICE="cpu"
    fi
fi
echo ""

# Set number of threads for CPU parallelization
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "Computation device: $DEVICE"
echo "Parallelization settings:"
echo "  OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "  MKL_NUM_THREADS: $MKL_NUM_THREADS"
echo ""

# Check dependencies
if ! python3 -c "import mattersim" 2>/dev/null; then
    echo "Error: MatterSim not found"
    exit 1
fi

if ! python3 -c "import pymatgen" 2>/dev/null; then
    echo "Error: pymatgen not found"
    exit 1
fi

# Expand paths
RESULTS_DIR=$(eval echo "$RESULTS_DIR")
OUTPUT_DIR=$(eval echo "$OUTPUT_DIR")

# Check results directory
if [ ! -d "$RESULTS_DIR" ]; then
    echo "Error: Results directory not found: $RESULTS_DIR"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Build command
CMD="python3 prescreen.py"
CMD="$CMD --results-dir $RESULTS_DIR"
CMD="$CMD --output-dir $OUTPUT_DIR"
CMD="$CMD --max-structures $MAX_STRUCTURES"
CMD="$CMD --hull-threshold $HULL_THRESHOLD"
CMD="$CMD --device $DEVICE"

if [ -n "$MAX_COMPOSITIONS" ]; then
    CMD="$CMD --max-compositions $MAX_COMPOSITIONS"
fi

if [ -n "$MP_API_KEY" ]; then
    CMD="$CMD --mp-api-key $MP_API_KEY"
fi

# Print configuration
echo "Configuration:"
echo "  Results dir: $RESULTS_DIR"
echo "  Output dir: $OUTPUT_DIR"
echo "  Max structures: $MAX_STRUCTURES"
echo "  Max compositions: ${MAX_COMPOSITIONS:-all}"
echo "  Hull threshold: ${HULL_THRESHOLD} eV/atom"
echo "  Device: $DEVICE"
echo "  MP API key: ${MP_API_KEY:+[SET]}${MP_API_KEY:-[NOT SET]}"
echo ""

echo "========================================================================"
echo "Starting pre-screening..."
echo "========================================================================"
echo ""
echo "Command: $CMD"
echo ""

# Run pre-screening
$CMD

EXIT_CODE=$?

echo ""
echo "========================================================================"
echo "Pre-screening finished"
echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"
echo "========================================================================"

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "Next step: Run workflow manager"
    echo "  sbatch submit_workflow_manager.sh"
fi

exit $EXIT_CODE

