#!/bin/bash
#SBATCH --job-name=get_mp_struct
#SBATCH --partition=Apus,Orion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --output=get_mp_struct_%j.out
#SBATCH --error=get_mp_struct_%j.err

# MP Structure Download Script (Legacy MPRester)
# Downloads crystal structures and energies from Materials Project

echo "======================================================================"
echo "MP Structure Download"
echo "======================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Running on: $(hostname)"
echo "======================================================================"

# Activate conda environment
source ~/.bashrc
conda activate mattersim

if [ $? -ne 0 ]; then
    echo "ERROR: Could not activate mattersim environment"
    exit 1
fi

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Python: $(which python3)"
echo ""

# Check MP API key
if [ -z "$MP_API_KEY" ]; then
    echo "ERROR: MP_API_KEY environment variable not set"
    echo "Set it in ~/.bashrc: export MP_API_KEY=your_32_character_key"
    exit 1
fi

echo "MP_API_KEY: ${MP_API_KEY:0:8}..." # Show first 8 chars only
echo "======================================================================"
echo ""

# Default parameters (can be overridden by environment variables)
CACHE_FILE="${CACHE_FILE:-./VASP_JOBS/mp_mattersim.json}"
OUTPUT_DIR="${OUTPUT_DIR:-./mp_cache_structs}"
CHEMSYS="${CHEMSYS:-}"
FORCE="${FORCE:-}"
PURE_PBE="${PURE_PBE:-}"

# Build command
CMD="python3 get_mp_struct.py --cache-file \"$CACHE_FILE\" --output-dir \"$OUTPUT_DIR\""

if [ -n "$CHEMSYS" ]; then
    CMD="$CMD --chemsys $CHEMSYS"
fi

if [ -n "$FORCE" ]; then
    CMD="$CMD --force"
fi

if [ -n "$PURE_PBE" ]; then
    CMD="$CMD --pure-pbe"
fi

# Display configuration
echo "Configuration:"
echo "  Cache file: $CACHE_FILE"
echo "  Output directory: $OUTPUT_DIR"
if [ -n "$CHEMSYS" ]; then
    echo "  Chemical system: $CHEMSYS"
else
    echo "  Chemical system: All (from cache file)"
fi
if [ -n "$FORCE" ]; then
    echo "  Force re-download: Yes"
else
    echo "  Force re-download: No (skip existing)"
fi
if [ -n "$PURE_PBE" ]; then
    echo "  Functional filtering: Pure GGA-PBE only (PBE+U excluded)"
else
    echo "  Functional filtering: Mixed PBE/PBE+U (recommended)"
fi
echo "======================================================================"
echo ""

# Run the download
echo "Downloading MP structures..."
eval $CMD

EXIT_CODE=$?

echo ""
echo "======================================================================"
echo "Job completed with exit code: $EXIT_CODE"
echo "End time: $(date)"
echo "======================================================================"

exit $EXIT_CODE

