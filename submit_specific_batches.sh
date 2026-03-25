#!/bin/bash
# Submit specific batch ranges for prescreening
# Usage: ./submit_specific_batches.sh START_BATCH END_BATCH [OPTIONS]

if [ $# -lt 2 ]; then
    echo "Usage: $0 START_BATCH END_BATCH [OPTIONS]"
    echo ""
    echo "Example:"
    echo "  $0 10 16 --results-dir ./mattergen_results/ternary_csp_electrides --batch-size 64"
    echo ""
    echo "This will submit batches 10 through 16 (inclusive)"
    exit 1
fi

START_BATCH=$1
END_BATCH=$2
shift 2

# Default configuration (must match your original submission)
COMPOSITIONS_PER_JOB=400
RESULTS_DIR="./mattergen_results/ternary_csp_electrides"
OUTPUT_DIR="./VASP_JOBS"
BATCH_SIZE=32
HULL_THRESHOLD=0.1
DEVICE="cuda"
MAX_STRUCTURES=0
MAX_ATOMS_GPU=2048
MODEL_PATH=""
PURE_PBE=""

# Parse remaining arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --results-dir)
            RESULTS_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --batch-size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --hull-threshold)
            HULL_THRESHOLD="$2"
            shift 2
            ;;
        --device)
            DEVICE="$2"
            shift 2
            ;;
        --compositions-per-job)
            COMPOSITIONS_PER_JOB="$2"
            shift 2
            ;;
        --max-structures)
            MAX_STRUCTURES="$2"
            shift 2
            ;;
        --max-atoms-gpu)
            MAX_ATOMS_GPU="$2"
            shift 2
            ;;
        --model-path)
            MODEL_PATH="$2"
            shift 2
            ;;
        --pure-pbe)
            PURE_PBE="--pure-pbe"
            shift
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            exit 1
            ;;
    esac
done

echo "========================================"
echo "Submitting Specific Batch Range"
echo "========================================"
echo "Batch range: $START_BATCH to $END_BATCH"
echo "Results directory: $RESULTS_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Compositions per job: $COMPOSITIONS_PER_JOB"
echo "Batch size: $BATCH_SIZE"
echo "Hull threshold: $HULL_THRESHOLD eV/atom"
echo "Max structures: $MAX_STRUCTURES"
echo "Max atoms on GPU: $MAX_ATOMS_GPU"
echo "Device: $DEVICE"
echo "Model: ${MODEL_PATH:-MatterSim-v1.0.0-5M.pth (default)}"
echo "Functional filter: ${PURE_PBE:+Pure PBE only}${PURE_PBE:-Mixed PBE/PBE+U}"
echo "========================================"
echo ""

# Submit each batch
for i in $(seq $START_BATCH $END_BATCH); do
    START_IDX=$((i * COMPOSITIONS_PER_JOB))
    
    TEMP_SCRIPT="submit_prescreen_batch${i}.sh"
    
    # Create SLURM script
    cat > $TEMP_SCRIPT <<EOF_HEADER
#!/bin/bash
#SBATCH --job-name=prescreen_b${i}
#SBATCH --partition=GPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=prescreen_batch${i}_%j.out
#SBATCH --error=prescreen_batch${i}_%j.err

# Load conda
eval "\$(conda shell.bash hook)"
conda activate mattersim

if [ \$? -ne 0 ]; then
    echo "ERROR: Failed to activate mattersim conda environment"
    exit 1
fi

# Print GPU info
nvidia-smi

# Set CUDA memory configuration
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

# Batch parameters
EOF_HEADER
    
    cat >> $TEMP_SCRIPT <<EOF_PARAMS
START_INDEX=$START_IDX
MAX_COMPS=$COMPOSITIONS_PER_JOB
BATCH_ID=$i
RESULTS_DIR="$RESULTS_DIR"
OUTPUT_DIR="$OUTPUT_DIR"
BATCH_SIZE=$BATCH_SIZE
MAX_ATOMS_GPU=$MAX_ATOMS_GPU
HULL_THRESHOLD=$HULL_THRESHOLD
MAX_STRUCTURES=$MAX_STRUCTURES
DEVICE="$DEVICE"
MODEL_PATH="$MODEL_PATH"
PURE_PBE="$PURE_PBE"

EOF_PARAMS
    
    cat >> $TEMP_SCRIPT << 'EOF_FOOTER'
echo "========================================"
echo "Prescreening Batch Job"
echo "========================================"
echo "Batch ID: $BATCH_ID"
echo "Start composition index: $START_INDEX"
echo "Max compositions: $MAX_COMPS"
echo "Composition range: $START_INDEX to $((START_INDEX + MAX_COMPS - 1))"
echo "Results directory: $RESULTS_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Batch size: $BATCH_SIZE"
echo "Hull threshold: $HULL_THRESHOLD eV/atom"
echo "Max structures per composition: $MAX_STRUCTURES"
echo "Max atoms on GPU: $MAX_ATOMS_GPU"
echo "Device: $DEVICE"
echo "Model: ${MODEL_PATH:-MatterSim-v1.0.0-5M.pth (default)}"
echo "========================================"

# Run prescreening
CMD="python3 prescreen.py"
CMD="$CMD --results-dir $RESULTS_DIR"
CMD="$CMD --output-dir $OUTPUT_DIR"
CMD="$CMD --hull-threshold $HULL_THRESHOLD"
CMD="$CMD --device $DEVICE"
CMD="$CMD --batch-size $BATCH_SIZE"
CMD="$CMD --max-atoms-gpu $MAX_ATOMS_GPU"
CMD="$CMD --start-composition $START_INDEX"
CMD="$CMD --max-compositions $MAX_COMPS"
CMD="$CMD --max-structures $MAX_STRUCTURES"
CMD="$CMD --batch-id $BATCH_ID"

if [ -n "$MODEL_PATH" ]; then
    CMD="$CMD --model-path $MODEL_PATH"
fi

if [ -n "$PURE_PBE" ]; then
    CMD="$CMD $PURE_PBE"
fi

$CMD

EXIT_CODE=$?

echo ""
echo "========================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "Batch completed successfully"
else
    echo "Batch failed with exit code $EXIT_CODE"
fi
echo "========================================"

exit $EXIT_CODE
EOF_FOOTER
    
    # Submit the job
    JOB_ID=$(sbatch $TEMP_SCRIPT | awk '{print $NF}')
    
    echo "Submitted batch $i (compositions $START_IDX-$((START_IDX + COMPOSITIONS_PER_JOB - 1))): Job ID $JOB_ID"
    
    # Clean up temporary script
    rm $TEMP_SCRIPT
    
    # Small delay
    sleep 0.5
done

echo ""
echo "========================================"
echo "Batch range $START_BATCH-$END_BATCH submitted!"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "========================================"

