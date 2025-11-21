#!/bin/bash

# Wrapper script to submit MP vs VASP Energy Comparison

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default values
DB_FILE="${SCRIPT_DIR}/mp_relax_workflow.json"
OUTPUT_FILE="${SCRIPT_DIR}/mp_vasp_comparison.json"
DEBUG=""

# Help message
show_help() {
    cat << EOF
Usage: bash run_compare_energies.sh [OPTIONS]

Submit MP vs VASP Energy Comparison

OPTIONS:
    --db FILE           Workflow database file (default: ./mp_relax_workflow.json)
    --output FILE       Output file for results (default: ./mp_vasp_comparison.json)
    --debug             Enable detailed debugging output
    --help              Show this help message

EXAMPLES:
    # Basic usage
    bash run_compare_energies.sh

    # With debug output
    bash run_compare_energies.sh --debug

    # Custom database location
    bash run_compare_energies.sh --db /path/to/workflow.json

ENVIRONMENT:
    MP_API_KEY          Materials Project API key (required)

EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --db)
            DB_FILE="$2"
            shift 2
            ;;
        --output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        --debug)
            DEBUG="--debug"
            shift
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

# Check if database exists
if [ ! -f "${DB_FILE}" ]; then
    echo "ERROR: Database file not found: ${DB_FILE}"
    echo ""
    echo "Run the MP relaxation workflow first:"
    echo "  bash run_mp_dft_relax.sh --chemsys B-Li-N"
    exit 1
fi

# Check MP_API_KEY
if [ -z "$MP_API_KEY" ]; then
    echo "ERROR: MP_API_KEY environment variable is not set"
    echo ""
    echo "Set it with:"
    echo "  export MP_API_KEY=your_api_key"
    exit 1
fi

# Print configuration
echo "========================================================================"
echo "Submitting MP vs VASP Energy Comparison"
echo "========================================================================"
echo "Database file: ${DB_FILE}"
echo "Output file: ${OUTPUT_FILE}"
if [ -n "${DEBUG}" ]; then
    echo "Debug mode: enabled"
fi
echo "========================================================================"
echo ""

# Export variables for SLURM script
export SCRIPT_DIR
export DB_FILE
export OUTPUT_FILE
export DEBUG

# Submit to SLURM
cd "${SCRIPT_DIR}"
sbatch submit_compare_energies.sh

if [ $? -eq 0 ]; then
    echo ""
    echo "Job submitted successfully!"
    echo ""
    echo "Monitor with:"
    echo "  squeue -u \$USER"
    echo "  tail -f compare_energies_*.out"
    echo ""
    echo "After completion, view results:"
    echo "  cat ${OUTPUT_FILE} | jq '.stats'"
    echo "  open mp_vasp_comparison_scatter.png"
    echo ""
else
    echo "ERROR: Job submission failed"
    exit 1
fi
