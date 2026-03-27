#!/bin/bash

# Wrapper script to submit MatterSim fine-tuning on VASP-relaxed MP phases

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default values
WORKFLOW_DB="${SCRIPT_DIR}/mp_relax_workflow.json"
COMPARISON_JSON="${SCRIPT_DIR}/mp_vasp_comparison.json"
OUTPUT_DIR="${SCRIPT_DIR}/ft_mattersim"
MODEL="mattersim-v1.0.0-5m"
EPOCHS=300
BATCH_SIZE=4
LR=2e-5
VAL_FRACTION=0.1
SEED=42
DEVICE="cuda"
INCLUDE_STRESSES="--include-stresses"
DATASET_ONLY=""
SKIP_FIRST=0
MAX_FORCE=10.0
MAX_FRAMES_PER_ID=""
TEST_FRACTION=0.1
PATIENCE=50
RE_NORMALIZE=""
EVAL_ONLY=""
MODEL_PATH=""
OUTLIER_WEIGHT=1.0
CONDA_ENV="mattersim"

show_help() {
    cat << EOF
Usage: bash run_ft_mattersim.sh [OPTIONS]

Fine-tune MatterSim on VASP-relaxed MP phase structures.

Workflow:
  1. Run mp_dft_relaxflow.py to relax MP phases with VASP
  2. Run run_compare_energies.sh to validate energies
  3. Run this script to fine-tune MatterSim on the validated structures

OPTIONS:
    --workflow-db FILE       Workflow database (default: ./mp_relax_workflow.json)
    --comparison-json FILE   Comparison results (default: ./mp_vasp_comparison.json)
    --output-dir DIR         Output directory (default: ./ft_mattersim)
    --model NAME             Pre-trained model (default: mattersim-v1.0.0-5m)
                             Options: mattersim-v1.0.0-1m, mattersim-v1.0.0-5m
    --epochs N               Training epochs (default: 300)
    --batch-size N           Batch size (default: 4)
    --lr RATE                Learning rate (default: 1e-4)
    --val-fraction F         Validation fraction (default: 0.2)
    --seed N                 Random seed (default: 42)
    --no-stresses            Disable stress in training
    --dataset-only           Only prepare train/val/test XYZ files, skip training
    --test-fraction F        Fraction for test set (default: 0.1)
    --skip-first N           Skip the first N ionic steps per trajectory (default: 0)
    --max-force F            Skip steps with max |force| > F eV/A (default: 10.0)
    --max-frames-per-id N    Cap frames per MP-ID by evenly sampling (default: all)
    --patience N             Early stopping patience in epochs (default: 50)
    --re-normalize           Re-normalize energy/force scales to finetuning data
    --eval-only              Only evaluate an existing finetuned model (skip training)
    --model-path FILE        Path to finetuned model checkpoint for evaluation
    --outlier-weight F       Outlier oversampling weight (default: 1.0 = balanced)
    --help                   Show this help message

EXAMPLES:
    # Standard fine-tuning (uses comparison JSON if available)
    bash run_ft_mattersim.sh

    # Fine-tune with 1M model (faster)
    bash run_ft_mattersim.sh --model mattersim-v1.0.0-1m --epochs 500

    # Skip first 2 ionic steps, filter large forces
    bash run_ft_mattersim.sh --skip-first 2 --max-force 30.0

    # Custom test fraction
    bash run_ft_mattersim.sh --test-fraction 0.15 --val-fraction 0.15

    # Only prepare dataset (no GPU needed)
    bash run_ft_mattersim.sh --dataset-only

    # Custom output and learning rate
    bash run_ft_mattersim.sh --output-dir ./ft_model_v2 --lr 5e-5 --epochs 500

    # Evaluate an existing finetuned model
    bash run_ft_mattersim.sh --eval-only --model-path ft_mattersim/model/mattersim-v1.0.0-5m_20260317_ft.pth

EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --workflow-db)
            WORKFLOW_DB="$2"
            shift 2
            ;;
        --comparison-json)
            COMPARISON_JSON="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --model)
            MODEL="$2"
            shift 2
            ;;
        --epochs)
            EPOCHS="$2"
            shift 2
            ;;
        --batch-size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --lr)
            LR="$2"
            shift 2
            ;;
        --val-fraction)
            VAL_FRACTION="$2"
            shift 2
            ;;
        --seed)
            SEED="$2"
            shift 2
            ;;
        --no-stresses)
            INCLUDE_STRESSES=""
            shift
            ;;
        --dataset-only)
            DATASET_ONLY="1"
            shift
            ;;
        --test-fraction)
            TEST_FRACTION="$2"
            shift 2
            ;;
        --skip-first)
            SKIP_FIRST="$2"
            shift 2
            ;;
        --max-force)
            MAX_FORCE="$2"
            shift 2
            ;;
        --max-frames-per-id)
            MAX_FRAMES_PER_ID="$2"
            shift 2
            ;;
        --patience)
            PATIENCE="$2"
            shift 2
            ;;
        --re-normalize)
            RE_NORMALIZE="1"
            shift
            ;;
        --eval-only)
            EVAL_ONLY="1"
            shift
            ;;
        --model-path)
            MODEL_PATH="$2"
            shift 2
            ;;
        --outlier-weight)
            OUTLIER_WEIGHT="$2"
            shift 2
            ;;
        --help)
            show_help
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Check workflow database
if [ ! -f "${WORKFLOW_DB}" ]; then
    echo "ERROR: Workflow database not found: ${WORKFLOW_DB}"
    echo ""
    echo "Run the MP relaxation workflow first:"
    echo "  bash run_mp_dft_relax.sh"
    echo "Then run the comparison:"
    echo "  bash run_compare_energies.sh"
    exit 1
fi

# Check comparison JSON (optional but preferred)
if [ -f "${COMPARISON_JSON}" ]; then
    echo "Found comparison JSON: ${COMPARISON_JSON}"
    echo "  (will use validated structures only)"
else
    echo "No comparison JSON found at: ${COMPARISON_JSON}"
    echo "  (will read all RELAX_DONE/RELAX_TMOUT from workflow DB)"
    COMPARISON_JSON=""
fi

echo ""
echo "========================================================================"
echo "Submitting MatterSim Fine-tuning"
echo "========================================================================"
echo "Workflow DB: ${WORKFLOW_DB}"
echo "Comparison JSON: ${COMPARISON_JSON:-none}"
echo "Output dir: ${OUTPUT_DIR}"
echo "Model: ${MODEL}"
echo "Epochs: ${EPOCHS}"
echo "Batch size: ${BATCH_SIZE}"
echo "Learning rate: ${LR}"
echo "Val fraction: ${VAL_FRACTION}"
echo "Include stresses: $([ -n "$INCLUDE_STRESSES" ] && echo yes || echo no)"
echo "Dataset only: $([ -n "$DATASET_ONLY" ] && echo yes || echo no)"
echo "Test fraction: ${TEST_FRACTION}"
echo "Skip first: ${SKIP_FIRST}"
echo "Max force: ${MAX_FORCE} eV/A"
echo "Max frames/ID: ${MAX_FRAMES_PER_ID:-all}"
echo "Patience: ${PATIENCE}"
echo "Re-normalize: $([ -n "$RE_NORMALIZE" ] && echo yes || echo no)"
echo "Eval only: $([ -n "$EVAL_ONLY" ] && echo yes || echo no)"
echo "Model path: ${MODEL_PATH:-auto}"
echo "Outlier weight: ${OUTLIER_WEIGHT}"
echo "========================================================================"
echo ""

# Export variables for SLURM script
export SCRIPT_DIR
export WORKFLOW_DB
export COMPARISON_JSON
export OUTPUT_DIR
export MODEL
export EPOCHS
export BATCH_SIZE
export LR
export VAL_FRACTION
export SEED
export DEVICE
export INCLUDE_STRESSES
export DATASET_ONLY
export SKIP_FIRST
export MAX_FORCE
export MAX_FRAMES_PER_ID
export TEST_FRACTION
export PATIENCE
export RE_NORMALIZE
export EVAL_ONLY
export MODEL_PATH
export OUTLIER_WEIGHT
export CONDA_ENV

# Submit to SLURM
cd "${SCRIPT_DIR}"
sbatch submit_ft_mattersim.sh

if [ $? -eq 0 ]; then
    echo ""
    echo "Job submitted successfully!"
    echo ""
    echo "Monitor with:"
    echo "  squeue -u \$USER"
    echo "  tail -f ft_mattersim_*.out"
    echo ""
    echo "After completion:"
    echo "  ls ${OUTPUT_DIR}/model/"
    echo "  cat ${OUTPUT_DIR}/dataset_meta.json | jq '.n_train, .n_val'"
    echo ""
else
    echo "ERROR: Job submission failed"
    exit 1
fi
