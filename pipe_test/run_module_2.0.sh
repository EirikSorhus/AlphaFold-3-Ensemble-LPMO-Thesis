#!/bin/bash
#SBATCH --job-name=module_2_0_ref     # Job name in SLURM queue
#SBATCH --account=nn1003k              # CHANGE: Your project ID
#SBATCH --time=02:00:00                # Max runtime
#SBATCH --mem=4G                       # Memory per node
#SBATCH --cpus-per-task=2              # CPU cores
#SBATCH --ntasks=1
#SBATCH --output=logs/module_2_0_%j.out  # Output file (%j = job ID)
#SBATCH --error=logs/module_2_0_%j.err   # Error file

################################################################################
# SLURM BATCH SCRIPT FOR MODULE 2.0 (Refactored)
################################################################################
#
# Kjører den refaktorerte Module 2.0 pipelines:
#  - Per-family CAZy fetch
#  - Multi-source enrichment (UniProt/NCBI/PDB/Taxonomy)
#  - Per-family CSV/FASTA
#  - Merge + run metadata
#
################################################################################

# ==========================================================================
# 1. JOBBINFO
# ==========================================================================
set -euo pipefail

echo "======================================================================"
echo "Module 2.0 – Refactored LPMO Enrichment"
echo "======================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo ""

# ==========================================================================
# 2. MILJØ
# ==========================================================================
module purge
module load NRIS/CPU
module load hpc-container-wrapper
module load Python/3.11.5-GCCcore-13.2.0

# Prefer conda-containerized env if available
export PATH=/cluster/projects/nn1003k/eirik/conda/lpmo_pipe_env/bin:$PATH

# Activate local venv
source /cluster/projects/nn1003k/eirik/Masteroppgave/.venv/bin/activate

echo "Python version: $(python --version)"
echo "Python location: $(which python)"
echo ""

# ==========================================================================
# 3. PARAMS
# ==========================================================================
# NCBI krever gyldig e-post
export NCBI_EMAIL="eirik.sorhus@nmbu.no"   # CHANGE if needed

# Families (default: defined in code). Override by CLI below if desired.
FAMILIES=(AA13)

# Generate RUN_ID for unique outputs
RUN_ID=$(date +"%Y%m%d_%H%M%S")

# ==========================================================================
# 4. ARBEIDSMAPPE
# ==========================================================================
cd /cluster/projects/nn1003k/eirik/Masteroppgave/pipe_test
mkdir -p data/metadata data/sequences logs

echo "Project directory: $(pwd)"
echo "Run ID: ${RUN_ID}"
echo "Families: ${FAMILIES[*]}"
echo ""

# ==========================================================================
# 5. KJØRING
# ==========================================================================
echo "======================================================================" 
echo "Starting Module 2.0 pipeline..."
echo "======================================================================"

# You can pass families and email explicitly; otherwise defaults are used.
# Example: --families AA13 AA9 --ncbi-email $NCBI_EMAIL

timeout 600 \
python -m scripts.module_2_0.run_module2_refactored \
    --families ${FAMILIES[@]} \
    --output-dir "$(pwd)/data" \
    --ncbi-email "$NCBI_EMAIL" \
    --run-id "$RUN_ID" \
    --log-level INFO
PIPELINE_EXIT=$?

if [ $PIPELINE_EXIT -eq 124 ]; then
    echo ""
    echo "======================================================================" 
    echo "ERROR: Pipeline timed out (10 minutes)"
    echo "Check network/API responsiveness or increase timeout."
    echo "======================================================================"
    exit 1
fi

# ==========================================================================
# 6. EXIT SUMMARY
# ==========================================================================
EXIT_CODE=$PIPELINE_EXIT

echo ""
echo "======================================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ MODULE 2.0 COMPLETED SUCCESSFULLY"
    echo "======================================================================"
    echo ""
    echo "Output files:"
    echo "  - data/metadata/m2_sequence_3d_metadata_${RUN_ID}.csv"
    echo "  - data/sequences/lpmo_all_${RUN_ID}_raw.fasta"
    echo "  - data/metadata/m2_run_metadata_${RUN_ID}.json"
else
    echo "✗ MODULE 2.0 FAILED (exit code: $EXIT_CODE)"
    echo "======================================================================"
    echo ""
    echo "Troubleshooting:" 
    echo "  1. Check logs: logs/module_2_0_${SLURM_JOB_ID}.out/.err"
    echo "  2. Verify NCBI_EMAIL is set and valid"
    echo "  3. Confirm Python packages are installed (requests, biopython, bs4)"
    echo "  4. Confirm internet access to UniProt, NCBI, RCSB"
fi

echo "End time: $(date)"
echo "======================================================================"

exit $EXIT_CODE
