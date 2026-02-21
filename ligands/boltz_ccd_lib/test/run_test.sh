#!/bin/bash
#SBATCH --job-name=cif_pkl_test
#SBATCH --account=nn1003k
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=logs/cif_pkl_test_%j.out
#SBATCH --error=logs/cif_pkl_test_%j.err

# Test runner for CIF to PKL conversion using boltz_ccd_lib

set -euo pipefail

# Define absolute paths BEFORE PATH modification
TEST_DIR="/cluster/work/projects/nn1003k/eirik/Masteroppgave/ligands/boltz_ccd_lib/test"
LIB_DIR="/cluster/work/projects/nn1003k/eirik/Masteroppgave/ligands/boltz_ccd_lib"
CIF_DIR="$TEST_DIR/cif"
PKL_DIR="$TEST_DIR/pkl"

# Environment setup
export PATH="/cluster/work/projects/nn1003k/eirik/conda/boltz_ccd_env/bin:$PATH"

# Create output directory
mkdir -p "$PKL_DIR"

echo "======================================================================"
echo "CIF to PKL Conversion Test"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Start time: $(date)"
echo "Environment: boltz_ccd_env"
echo "Python: $(which python)"
echo "CIF directory: $CIF_DIR"
echo "PKL directory: $PKL_DIR"
echo "======================================================================"
echo ""

# Verify input directory
if [ ! -d "$CIF_DIR" ]; then
    echo "ERROR: CIF directory not found: $CIF_DIR"
    exit 1
fi

# Clean output
rm -f "$PKL_DIR"/*.pkl 2>/dev/null || true

echo "=== Conversion Phase ==="
echo ""

# Convert each CIF file
CONVERSION_SUCCESS=0
CONVERSION_TOTAL=0
for cif_file in "$CIF_DIR"/*.cif; do
    if [ -f "$cif_file" ]; then
        basename=$(basename "$cif_file")
        CONVERSION_TOTAL=$((CONVERSION_TOTAL + 1))
        
        echo "Processing: $basename"
        echo "---"
        
        if python "$LIB_DIR/cif_to_pkl.py" "$cif_file" "$PKL_DIR"; then
            CONVERSION_SUCCESS=$((CONVERSION_SUCCESS + 1))
        else
            echo "✗ Conversion failed for $basename"
        fi
        
        echo ""
    fi
done

echo "Conversion results: $CONVERSION_SUCCESS/$CONVERSION_TOTAL succeeded"
echo ""

echo "=== Validation Phase ==="
echo ""

# Validate each PKL file
VALIDATION_SUCCESS=0
VALIDATION_TOTAL=0

for pkl_file in "$PKL_DIR"/*.pkl; do
    if [ -f "$pkl_file" ]; then
        basename=$(basename "$pkl_file")
        VALIDATION_TOTAL=$((VALIDATION_TOTAL + 1))
        
        echo "Validating: $basename"
        echo "---"
        
        if python "$LIB_DIR/validate_pkl.py" "$pkl_file"; then
            VALIDATION_SUCCESS=$((VALIDATION_SUCCESS + 1))
        else
            echo "✗ Validation failed for $basename"
        fi
        
        echo ""
    fi
done

echo "Validation results: $VALIDATION_SUCCESS/$VALIDATION_TOTAL succeeded"
echo ""

echo "======================================================================"
if [ $CONVERSION_SUCCESS -eq $CONVERSION_TOTAL ] && [ $VALIDATION_SUCCESS -eq $VALIDATION_TOTAL ]; then
    echo "✓ All tests completed successfully!"
    EXIT_CODE=0
else
    echo "✗ Some tests failed"
    EXIT_CODE=1
fi

echo "End time: $(date)"
echo "======================================================================"

exit $EXIT_CODE
