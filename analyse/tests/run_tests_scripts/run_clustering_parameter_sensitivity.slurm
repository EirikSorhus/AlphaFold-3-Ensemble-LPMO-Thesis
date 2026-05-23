#!/bin/bash
#SBATCH --job-name=clust_param_sens
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -euo pipefail

REPO_ROOT="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse"
PYTHON_BIN="/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin/python"
RUN_ROOT="${RUN_ROOT:-$REPO_ROOT/tests/tests_results/clustering_pilot_staged}"

if [[ $# -gt 0 ]]; then
    OUTPUT_ROOT=$(realpath -m "$1")
else
    timestamp=$(date +%Y%m%d_%H%M%S)
    OUTPUT_ROOT="$RUN_ROOT/clustering_parameter_sensitivity_runs/clustering_parameter_sensitivity_$timestamp"
fi

mkdir -p "$(dirname "$OUTPUT_ROOT")"

if [[ -e "$OUTPUT_ROOT" ]] && [[ -n "$(find "$OUTPUT_ROOT" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    echo "[ERROR] Refusing to write into non-empty output directory: $OUTPUT_ROOT" >&2
    exit 2
fi

cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT/src:${PYTHONPATH:-}"
export PATH="/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin:$PATH"

echo "[INFO] Input run root: $RUN_ROOT"
echo "[INFO] Output root:    $OUTPUT_ROOT"

"$PYTHON_BIN" "$REPO_ROOT/tests/run_tests_scripts/run_clustering_parameter_sensitivity.py" \
  --input-run-root "$RUN_ROOT" \
  --output-root "$OUTPUT_ROOT" \
  --hdbscan-min-cluster-sizes 3 5 10 \
  --hdbscan-min-samples none \
  --agglomerative-distance-thresholds 0.35 0.45 0.55 \
  --agglomerative-min-cluster-sizes 3 5 10
