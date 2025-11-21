#!/bin/bash
# Wrapper script to submit MP structure download as SLURM job

set -e

# Default values
CACHE_FILE="./VASP_JOBS/mp_mattersim.json"
OUTPUT_DIR="./mp_cache_structs"
CHEMSYS=""
FORCE=""
SUBMIT_SCRIPT="submit_get_mp.sh"

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --cache-file)
            CACHE_FILE="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --chemsys)
            CHEMSYS="$2"
            shift 2
            ;;
        --force)
            FORCE="--force"
            shift 1
            ;;
        -h|--help)
            cat << EOF
Usage: $0 [OPTIONS]

Download MP structures from cached MP IDs in mp_mattersim.json.

Options:
    --cache-file PATH     Path to mp_mattersim.json file
                         (default: ./VASP_JOBS/mp_mattersim.json)
    --output-dir PATH     Output directory for structures
                         (default: ./mp_cache_structs)
    --chemsys SYSTEM     Process specific chemical system (e.g., B-Li-N)
                         (default: process all chemical systems in cache)
    --force              Force re-download even if CIF files already exist
    -h, --help           Show this help message

Example:
    # Download all cached chemical systems (recommended)
    bash run_get_mp.sh

    # Download specific chemical system
    bash run_get_mp.sh --chemsys B-Li-N

    # Force re-download (overwrite existing CIF files)
    bash run_get_mp.sh --force

    # Custom cache file path
    bash run_get_mp.sh --cache-file /path/to/mp_mattersim.json

Notes:
    - Requires mattersim conda environment (for mp-api)
    - Requires MP_API_KEY environment variable (32 characters)
    - Reads single mp_mattersim.json file (new cache strategy)
    - Downloads structure CIF files only (no metadata)
    - Saves flat: mp_cache_structs/mp-XXXXX.cif (no subdirectories)
    - MP energies are fetched separately by compare_mp_vasp_energies.py

Output Structure:
    mp_cache_structs/
    ├── mp-12345.cif
    ├── mp-67890.cif
    └── ...

EOF
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Check if submit script exists
if [ ! -f "$SUBMIT_SCRIPT" ]; then
    echo "ERROR: Submit script not found: $SUBMIT_SCRIPT"
    exit 1
fi

# Export variables for SLURM script
export CACHE_FILE
export OUTPUT_DIR
export CHEMSYS
export FORCE

echo "======================================================================"
echo "Submitting MP Structure Download"
echo "======================================================================"
echo "Cache file: $CACHE_FILE"
echo "Output directory: $OUTPUT_DIR"
if [ -n "$CHEMSYS" ]; then
    echo "Chemical system: $CHEMSYS"
else
    echo "Chemical system: All (from cache file)"
fi
if [ -n "$FORCE" ]; then
    echo "Force re-download: Yes"
else
    echo "Force re-download: No (skip existing)"
fi
echo "======================================================================"
echo ""

# Submit job
JOB_ID=$(sbatch "$SUBMIT_SCRIPT" | awk '{print $4}')

if [ -n "$JOB_ID" ]; then
    echo "Job submitted successfully!"
    echo "Job ID: $JOB_ID"
    echo ""
    echo "Monitor with:"
    echo "  squeue -j $JOB_ID"
    echo "  tail -f get_mp_struct_${JOB_ID}.out"
    echo ""
    echo "When complete:"
    echo "  - Structures: $OUTPUT_DIR/mp-*.cif"
else
    echo "ERROR: Job submission failed"
    exit 1
fi
