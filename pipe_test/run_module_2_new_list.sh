#!/bin/bash
#SBATCH --job-name=mod2_list_test
#SBATCH --account=nn1003k
#SBATCH --time=01:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=logs/mod2_new_list_%j.out
#SBATCH --error=logs/mod2_new_list_%j.err

# Test runner for the reimplemented Module 2 (scripts/module_2_new)
# using the 'list' mode with CAZy characterized CSV format

set -euo pipefail

# Environment
export PATH=/cluster/work/projects/nn1003k/eirik/conda/lpmo_pipe_env/bin:$PATH

PROJECT_DIR="/cluster/work/projects/nn1003k/eirik/Masteroppgave/pipe_test"
SCRIPT_PATH="scripts/module_2_new/main_driver.py"

# Input Configuration
MODE="list"
INPUT_FILE="data/run/failed_id_fix/aa10_failed_new_code.txt"
OUTPUT_DIR="data"

cd "$PROJECT_DIR"
mkdir -p logs
mkdir -p "$OUTPUT_DIR"

echo "======================================================================"
echo "Module 2 (New) Test Run - List Mode with CAZy Characterized failed Input"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Start time: $(date)"
echo "Mode: $MODE"
echo "Input: $INPUT_FILE"
echo "Output: $OUTPUT_DIR"
echo "Python: $(which python)"
echo "Script: $SCRIPT_PATH"
echo "======================================================================"

# Ensure input file exists before running
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found!"
    ls -l data/test_data/ || echo "Directory data/test_data/ not found"
    exit 1
fi


python "$SCRIPT_PATH" \
  --mode "$MODE" \
  --input "$INPUT_FILE" \
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
