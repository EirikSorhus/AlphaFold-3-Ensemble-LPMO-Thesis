#!/bin/bash
#SBATCH --job-name=module_2_fasta_test
#SBATCH --account=nn1003k
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=logs/module_2_fasta_%j.out
#SBATCH --error=logs/module_2_fasta_%j.err

# Batch runner for scripts/module_2/main_driver.py using test FASTA input.
# Uses the conda-containerized env at lpmo_pipe_env/bin per user request.

set -euo pipefail

# Environment
export PATH=/cluster/work/projects/nn1003k/eirik/conda/lpmo_pipe_env/bin:$PATH

CONTACT_EMAIL="eirik.sorhus@nmbu.no"
PROJECT_DIR="/cluster/work/projects/nn1003k/eirik/Masteroppgave/pipe_test"
FASTA_INPUT="${PROJECT_DIR}/data/test_data/test_fasta.fasta"

cd "$PROJECT_DIR"
mkdir -p logs

echo "======================================================================"
echo "Module 2 FASTA run"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Node: ${SLURM_NODELIST:-local}"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo "Python: $(which python)"
echo "======================================================================"

timeout 900 \
python scripts/module_2/main_driver.py \
  --input-fasta "$FASTA_INPUT" \
  --project-dir "$PROJECT_DIR" \
  --contact-email "$CONTACT_EMAIL" \
  --allow-ncbi-fallback

EXIT_CODE=$?

echo "======================================================================"
if [ $EXIT_CODE -eq 0 ]; then
  echo "✓ MODULE 2 completed"
else
  echo "✗ MODULE 2 failed (exit $EXIT_CODE)"
fi

echo "End time: $(date)"
echo "======================================================================"

exit $EXIT_CODE
