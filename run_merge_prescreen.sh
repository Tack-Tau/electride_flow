#!/bin/bash
# VASPflow - Merge prescreening batch results via SLURM
# Usage: bash run_merge_prescreen.sh [options]

set -e

# Default values
OUTPUT_DIR="./VASP_JOBS"
CONDA_ENV="mattersim"
SKIP_DATABASE=""
SKIP_DEDUP=""
CLEAN_BATCHES=""
OUTPUT_FILE=""

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo "========================================================================"
echo "VASPflow - Merge Prescreening Batches"
echo "========================================================================"
echo ""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        --skip-database)
            SKIP_DATABASE="1"
            shift
            ;;
        --skip-dedup)
            SKIP_DEDUP="1"
            shift
            ;;
        --clean-batches)
            CLEAN_BATCHES="1"
            shift
            ;;
        --output-file)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        --help)
            echo "Usage: bash run_merge_prescreen.sh [options]"
            echo ""
            echo "Submits a SLURM job to merge parallel prescreening batch results."
            echo "Merges JSON results, PyXtal databases, and runs post-relaxation deduplication."
            echo ""
            echo "Options:"
            echo "  --output-dir DIR       Directory containing batch files (default: ./VASP_JOBS)"
            echo "  --conda-env NAME       Conda environment name (default: mattersim)"
            echo "  --skip-database        Skip merging PyXtal database files (JSON only)"
            echo "  --skip-dedup           Skip post-relaxation StructureMatcher deduplication"
            echo "  --clean-batches        Delete batch files after successful merge"
            echo "                         (default: keep them for safety)"
            echo "  --output-file NAME     Output JSON filename (default: prescreening_stability.json)"
            echo ""
            echo "Examples:"
            echo "  bash run_merge_prescreen.sh"
            echo "  bash run_merge_prescreen.sh --output-dir ./VASP_JOBS"
            echo "  bash run_merge_prescreen.sh --skip-dedup"
            echo "  bash run_merge_prescreen.sh --clean-batches"
            echo ""
            echo "Monitoring:"
            echo "  Check status:    squeue -u \$USER | grep merge"
            echo "  View log:        tail -f merge_prescreen_<JOBID>.out"
            echo ""
            exit 0
            ;;
        *)
            echo -e "${RED}Error: Unknown option $1${NC}"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Validate inputs
if [ ! -d "$OUTPUT_DIR" ]; then
    echo -e "${RED}Error: Output directory not found: $OUTPUT_DIR${NC}"
    exit 1
fi

# Check batch files exist
BATCH_COUNT=$(ls "$OUTPUT_DIR"/prescreening_stability_batch*.json 2>/dev/null | wc -l)
if [ "$BATCH_COUNT" -eq 0 ]; then
    echo -e "${RED}Error: No batch files found in $OUTPUT_DIR${NC}"
    echo "  Expected: prescreening_stability_batch*.json"
    exit 1
fi

# Check if merged file already exists
if [ -f "$OUTPUT_DIR/prescreening_stability.json" ] && [ -z "$OUTPUT_FILE" ]; then
    echo -e "${YELLOW}Warning: Merged results already exist:${NC}"
    echo "  $OUTPUT_DIR/prescreening_stability.json"
    echo ""
    read -p "Overwrite existing results? (y/N) " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborting."
        exit 0
    fi
fi

# Print configuration
echo "Configuration:"
echo "  Output directory:     $OUTPUT_DIR"
echo "  Batch files found:    $BATCH_COUNT"
echo "  Conda environment:    $CONDA_ENV"
echo "  Skip database merge:  $([ -n "$SKIP_DATABASE" ] && echo "yes" || echo "no")"
echo "  Skip post-relax dedup:$([ -n "$SKIP_DEDUP" ] && echo "yes" || echo "no")"
echo "  Clean batch files:    $([ -n "$CLEAN_BATCHES" ] && echo "yes" || echo "no")"
echo "  Output file:          ${OUTPUT_FILE:-prescreening_stability.json}"
echo ""

# Submit job
echo "Submitting merge job to SLURM..."
echo ""

export OUTPUT_DIR
export CONDA_ENV
export SKIP_DATABASE
export SKIP_DEDUP
export CLEAN_BATCHES
export OUTPUT_FILE

JOB_ID=$(sbatch submit_merge_prescreen.sh | awk '{print $NF}')

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Merge job submitted successfully!${NC}"
    echo ""
    echo "Job ID: $JOB_ID"
    echo ""
    echo "Monitor progress:"
    echo "  squeue -u $USER | grep merge"
    echo "  tail -f merge_prescreen_${JOB_ID}.out"
    echo ""
    echo "When complete, check results:"
    echo "  cat $OUTPUT_DIR/prescreening_stability.json | jq '.summary'"
    echo ""
else
    echo -e "${RED}Error: Failed to submit job${NC}"
    exit 1
fi
