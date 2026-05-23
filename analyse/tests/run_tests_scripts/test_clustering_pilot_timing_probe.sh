#!/bin/bash
#SBATCH --job-name=run_tests_clustering_pilot_timing_probe
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_clustering_pilot_timing_probe_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "[ERROR] This timing probe must be submitted with sbatch, not run directly on the login node."
    echo "[INFO] Example: sbatch tests/run_tests_scripts/test_clustering_pilot_timing_probe.sh"
    exit 1
fi

REPO_ROOT="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse"
PYTHON_BIN="/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin/python"
RUN_DIR="$REPO_ROOT/tests/tests_results/clustering_pilot_timing_probe_${SLURM_JOB_ID}"
N_JOBS_VALUES=(1 2 4)
MODES=(prepare_only pilot_like qc_full)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-dir)
            RUN_DIR=$(realpath -m "$2")
            shift 2
            ;;
        --n-jobs)
            IFS=',' read -r -a N_JOBS_VALUES <<< "$2"
            shift 2
            ;;
        --modes)
            IFS=',' read -r -a MODES <<< "$2"
            shift 2
            ;;
        --skip-normalization-profile)
            SKIP_PROFILE=1
            shift
            ;;
        --help)
            cat <<'EOF'
Usage:
    sbatch tests/run_tests_scripts/test_clustering_pilot_timing_probe.sh [options]

Options:
    --run-dir PATH                  Override output directory.
    --n-jobs CSV                    Worker counts, e.g. 1,2,4.
    --modes CSV                     Modes: prepare_only,pilot_like,qc_full.
    --skip-normalization-profile    Skip cProfile normalization pass.

Default:
    Runs the 12-pose timing probe with n_jobs=1,2,4 and all three modes.
EOF
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            exit 2
            ;;
    esac
done

cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT/src:${PYTHONPATH:-}"
export PATH="/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin:$PATH"
export XDG_CACHE_HOME="$RUN_DIR/cache"
mkdir -p "$RUN_DIR" "$XDG_CACHE_HOME"

echo "[INFO] Running focused timing-probe unit coverage first"
"$PYTHON_BIN" -m pytest tests/test_clustering_pilot_timing_probe.py -q

PROFILE_ARGS=()
if [[ "${SKIP_PROFILE:-0}" -eq 1 ]]; then
    PROFILE_ARGS+=(--skip-normalization-profile)
fi

echo "[INFO] Running 12-pose clustering-pilot timing probe"
"$PYTHON_BIN" tests/run_tests_scripts/run_clustering_pilot_timing_probe.py \
    --run-dir "$RUN_DIR" \
    --n-jobs "${N_JOBS_VALUES[@]}" \
    --modes "${MODES[@]}" \
    "${PROFILE_ARGS[@]}"

echo "[INFO] Timing probe complete"
echo "[INFO] Summary: $RUN_DIR/timing_summary.json"
echo "[INFO] Events:  $RUN_DIR/timing_events.tsv"
