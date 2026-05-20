#!/bin/bash
#SBATCH --job-name=build_clust_fixture
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/build_primary_clustering_stage_fixture_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "[ERROR] Submit with sbatch, do not run directly on the login node." >&2
    echo "[INFO] Example: sbatch tests/run_tests_scripts/build_primary_clustering_stage_fixture.slurm" >&2
    exit 1
fi

REPO_ROOT="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse"
PYTHON_BIN="/cluster/work/projects/nn1003k/eirik/conda/analyse_env/bin/python"
SELECTED_GRID="${SELECTED_GRID:-$REPO_ROOT/tests/fixtures/clustering_stage_outputs/selected_primary_clustering_grid_results.tsv}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$REPO_ROOT/tests/fixtures/clustering_stage_outputs}"

cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT/src:${PYTHONPATH:-}"
export PATH="/cluster/work/projects/nn1003k/eirik/conda/analyse_env/bin:$PATH"

echo "[INFO] Selected grid: $SELECTED_GRID"
echo "[INFO] Output root:   $OUTPUT_ROOT"

"$PYTHON_BIN" "$REPO_ROOT/tests/run_tests_scripts/build_primary_clustering_stage_fixture.py" \
  --selected-grid-results "$SELECTED_GRID" \
  --output-root "$OUTPUT_ROOT"
