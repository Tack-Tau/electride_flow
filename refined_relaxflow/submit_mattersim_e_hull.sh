#!/bin/bash
#SBATCH --job-name=mattersim_e_hull
#SBATCH --partition=GPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=mattersim_e_hull_%j.out
#SBATCH --error=mattersim_e_hull_%j.err

# MatterSim Energy Above Hull Calculation for Refined Structures
# Queries MP API for GGA/GGA+U phases and relaxes with MatterSim
# Relaxes VASP-refined structures with MatterSim (fmax=0.001, max_steps=800)

echo "========================================================================"
echo "MatterSim Energy Above Hull - Refined Structures"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "========================================================================"
echo ""

# Load conda environment
echo "Activating conda environment..."
source ~/.bashrc
conda activate ${CONDA_ENV:-mattersim}

# Verify MatterSim is available
echo "Verifying MatterSim..."
python3 -c "from mattersim.forcefield import MatterSimCalculator; print('MatterSim OK')" || {
    echo "ERROR: MatterSim not available in conda environment"
    exit 1
}

echo ""

# Check GPU availability
if command -v nvidia-smi &> /dev/null; then
    echo "GPU Information:"
    echo "========================================================================"
    nvidia-smi --query-gpu=index,name,memory.total,memory.free --format=csv,noheader
    
    GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
    GPU_MEM_TOTAL=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -1)
    GPU_MEM_FREE=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
    
    echo ""
    echo "Allocated GPU: $GPU_NAME"
    echo "Total VRAM: $GPU_MEM_TOTAL MiB (~$((GPU_MEM_TOTAL/1024)) GB)"
    echo "Free VRAM: $GPU_MEM_FREE MiB (~$((GPU_MEM_FREE/1024)) GB)"
    echo "========================================================================"
    echo ""
    
    # Auto-detect CUDA availability
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

# Set threading for CPU parallelization
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

# PyTorch CUDA memory allocator settings
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

# Configuration from environment variables
REFINE_JOBS_DIR=${REFINE_JOBS_DIR:-"./REFINE_VASP_JOBS"}
DEVICE=${DEVICE:-"cuda"}
MP_API_KEY=${MP_API_KEY:-""}

# Verify MP API key
if [ -z "$MP_API_KEY" ]; then
    echo "ERROR: MP_API_KEY environment variable not set"
    echo "Set it with: export MP_API_KEY=your_32_character_key"
    echo "Get key from: https://next-gen.materialsproject.org/api"
    exit 1
fi

echo "Configuration:"
echo "  Refined VASP jobs: $REFINE_JOBS_DIR"
echo "  Device: $DEVICE"
if [ -n "$PURE_PBE" ]; then
    echo "  Functional filter: Pure GGA-PBE only (PBE+U excluded)"
else
    echo "  Functional filter: Mixed PBE/PBE+U (recommended)"
fi
echo "  MP API key: ${MP_API_KEY:0:8}... (${#MP_API_KEY} chars)"
echo ""
echo "Parallelization settings:"
echo "  OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "  MKL_NUM_THREADS: $MKL_NUM_THREADS"
echo "  PYTORCH_CUDA_ALLOC_CONF: $PYTORCH_CUDA_ALLOC_CONF"
echo ""

# Verify required files exist
if [ ! -d "$REFINE_JOBS_DIR" ]; then
    echo "ERROR: Refined VASP jobs directory not found: $REFINE_JOBS_DIR"
    exit 1
fi

if [ ! -f "$REFINE_JOBS_DIR/workflow.json" ]; then
    echo "ERROR: Workflow database not found: $REFINE_JOBS_DIR/workflow.json"
    exit 1
fi

echo "All required files found"
echo ""

# Run MatterSim relaxation
echo "========================================================================"
echo "Starting MatterSim Relaxation (fmax=0.001, max_steps=800)"
echo "========================================================================"
echo ""

python3 compute_mattersim_e_hull.py \
    --refine-jobs "$REFINE_JOBS_DIR" \
    --mp-api-key "$MP_API_KEY" \
    --device "$DEVICE" \
    $PURE_PBE

EXIT_CODE=$?

echo ""
echo "========================================================================"
echo "Job Completed"
echo "========================================================================"
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"

if [ $EXIT_CODE -eq 0 ]; then
    echo "Status: SUCCESS"
    echo ""
    echo "Output files:"
    echo "  $REFINE_JOBS_DIR/mattersim_stability_results.json (MatterSim energies & hulls)"
    echo "  $REFINE_JOBS_DIR/mp_mattersim.json (MP phases relaxed with MatterSim)"
else
    echo "Status: FAILED"
fi

echo "========================================================================"

exit $EXIT_CODE

