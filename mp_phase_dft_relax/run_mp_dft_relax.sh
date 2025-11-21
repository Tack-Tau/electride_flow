#!/bin/bash

# Wrapper script to submit MP Phase DFT Relaxation Workflow

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default values
CACHE_DIR="${SCRIPT_DIR}/mp_cache_structs"
VASP_JOBS_DIR="${SCRIPT_DIR}/VASP_JOBS"
DB_FILE="${SCRIPT_DIR}/mp_relax_workflow.json"
MAX_CONCURRENT=10
CHECK_INTERVAL=120
CHEMSYS=""

# Help message
show_help() {
    cat << EOF
Usage: bash run_mp_dft_relax.sh [OPTIONS]

Submit MP Phase DFT Relaxation Workflow Manager

OPTIONS:
    --cache-dir DIR         Directory with cached MP structures (default: ./mp_cache_structs)
    --vasp-jobs-dir DIR     Directory for VASP jobs (default: ./VASP_JOBS)
    --db FILE               Workflow database file (default: ./mp_relax_workflow.json)
    --max-concurrent N      Maximum concurrent VASP jobs (default: 10)
    --check-interval N      Job status check interval in seconds (default: 300)
    --chemsys CHEMSYS       Process specific chemical system only (e.g., 'B-Li-N')
    --help                  Show this help message

EXAMPLES:
    # Process all chemical systems
    bash run_mp_dft_relax.sh --max-concurrent 10

    # Process only B-Li-N system (for testing)
    bash run_mp_dft_relax.sh --chemsys B-Li-N --max-concurrent 5

    # Resume with existing database
    bash run_mp_dft_relax.sh --db ./mp_relax_workflow.json

EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --cache-dir)
            CACHE_DIR="$2"
            shift 2
            ;;
        --vasp-jobs-dir)
            VASP_JOBS_DIR="$2"
            shift 2
            ;;
        --db)
            DB_FILE="$2"
            shift 2
            ;;
        --max-concurrent)
            MAX_CONCURRENT="$2"
            shift 2
            ;;
        --check-interval)
            CHECK_INTERVAL="$2"
            shift 2
            ;;
        --chemsys)
            CHEMSYS="$2"
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

# Check if cache directory exists
if [ ! -d "${CACHE_DIR}" ]; then
    echo "ERROR: Cache directory not found: ${CACHE_DIR}"
    echo ""
    echo "Run get_mp_struct.py first to download MP structures:"
    echo "  python3 get_mp_struct.py --cache-dir ./mp_mattersim_cache --chemsys B-Li-N"
    exit 1
fi

# Print configuration
echo "========================================================================"
echo "Submitting MP Phase DFT Relaxation Workflow"
echo "========================================================================"
echo "Cache directory: ${CACHE_DIR}"
echo "VASP jobs directory: ${VASP_JOBS_DIR}"
echo "Database file: ${DB_FILE}"
echo "Max concurrent jobs: ${MAX_CONCURRENT}"
echo "Check interval: ${CHECK_INTERVAL}s"
if [ -n "${CHEMSYS}" ]; then
    echo "Chemical system: ${CHEMSYS}"
fi
echo "========================================================================"
echo ""

# Export variables for SLURM script (same pattern as run_workflow.sh)
export CACHE_DIR
export VASP_JOBS_DIR
export DB_FILE
export MAX_CONCURRENT
export CHECK_INTERVAL
export CHEMSYS

# Submit to SLURM
cd "${SCRIPT_DIR}"
sbatch submit_mp_dft_relax.sh

if [ $? -eq 0 ]; then
    echo ""
    echo "Job submitted successfully!"
    echo ""
    echo "Monitor with:"
    echo "  squeue -u \$USER"
    echo "  tail -f mp_dft_relax_*.out"
    echo ""
    echo "Check database:"
    echo "  cat ${DB_FILE} | jq '.structures | to_entries | map({mp_id: .key, state: .value.state})'"
    echo ""
else
    echo "ERROR: Job submission failed"
    exit 1
fi

