#!/bin/bash
#SBATCH --job-name=mod2_new_test
#SBATCH --account=nn1003k
#SBATCH --time=03:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=logs/mod2_new_cazy_%j.out
#SBATCH --error=logs/mod2_new_cazy_%j.err

# Test runner for the reimplemented Module 2 (scripts/module_2_new)
# using the 'cazy' mode and test data found in pipe_test/data/test_data

set -euo pipefail

# Environment
export PATH=/cluster/work/projects/nn1003k/eirik/conda/lpmo_pipe_env/bin:$PATH

PROJECT_DIR="/cluster/work/projects/nn1003k/eirik/Masteroppgave/pipe_test"
SCRIPT_PATH="scripts/module_2_new/main_driver.py"

# Input Configuration
MODE="cazy"
INPUT="data/test_data/AA15.txt"  # Using CAZy family name for test (downsized data)
# INPUT="AA15" # Alternativ if whole family name is to be used directly
OUTPUT_DIR="data"

cd "$PROJECT_DIR"
mkdir -p logs
mkdir -p "$OUTPUT_DIR"

echo "======================================================================"
echo "Module 2 (New) Test Run"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Start time: $(date)"
echo "Mode: $MODE"
echo "Input: $INPUT"
echo "Output: $OUTPUT_DIR"
echo "Python: $(which python)"
echo "Script: $SCRIPT_PATH"
echo "======================================================================"

# Note: In cazy mode with a family name (e.g. "AA15"),
# the script handles the download itself, so we don't check for file existence here unless we provide a path.


python "$SCRIPT_PATH" \
  --mode "$MODE" \
  --input "$INPUT" \
  --output_dir "$OUTPUT_DIR"


EXIT_CODE=$?

echo "======================================================================"
if [ $EXIT_CODE -eq 0 ]; then
  echo "✓ MODULE 2 (New) completed successfully"
  echo "Results in: $OUTPUT_DIR"
else
  echo "✗ MODULE 2 (New) failed (exit $EXIT_CODE)"
fi

echo "End time: $(date)"
echo "======================================================================"

exit $EXIT_CODE
