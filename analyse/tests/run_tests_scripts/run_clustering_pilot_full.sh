#!/bin/bash
#SBATCH --job-name=run_tclustering_pilot_full
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=3
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_clustering_pilot_full_%j.log


set -euo pipefail

REPO_ROOT="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse"
SCRIPT_DIR="$REPO_ROOT/tests/run_tests_scripts"
PYTHON_BIN="/cluster/work/projects/nn1003k/eirik/conda/analyse_env/bin/python"
OVERVIEW_TSV="$REPO_ROOT/input_data/clustering_pilot_protein_overview.tsv"
DEFAULT_RUN_ROOT="$REPO_ROOT/tests/tests_results/clustering_pilot_full"

RUN_ROOT="$DEFAULT_RUN_ROOT"
PREPARE_ONLY=false
FORCE_STEPS=()
N_JOBS="${SLURM_CPUS_PER_TASK:-1}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-dir)
      RUN_ROOT=$(realpath -m "$2")
      shift 2
      ;;
    --prepare-only)
      PREPARE_ONLY=true
      shift
      ;;
    --force-step)
      FORCE_STEPS+=("$2")
      shift 2
      ;;
    --n-jobs)
      N_JOBS="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

if [[ "$PREPARE_ONLY" == false && -z "${SLURM_JOB_ID:-}" ]]; then
  echo "[ERROR] Actual pilot execution must be submitted with sbatch, not run directly on the login node." >&2
  echo "[INFO] Example: sbatch tests/run_tests_scripts/run_clustering_pilot_full.sh --run-dir tests/tests_results/clustering_pilot_full_run" >&2
  echo "[INFO] Use --prepare-only only when you intentionally want a non-sbatch readiness check." >&2
  exit 1
fi

cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT/src:${PYTHONPATH:-}"

LOG_DIR="$RUN_ROOT/logs"
MANIFEST_DIR="$RUN_ROOT/manifests"
SUMMARY_DIR="$RUN_ROOT/summaries"
DOMAIN_RUN_DIR="$RUN_ROOT/domain_only"
FULL_RUN_DIR="$RUN_ROOT/full_length"

mkdir -p "$LOG_DIR" "$MANIFEST_DIR" "$SUMMARY_DIR"

[[ -x "$PYTHON_BIN" ]] || { echo "Missing python interpreter: $PYTHON_BIN" >&2; exit 2; }
[[ -f "$OVERVIEW_TSV" ]] || { echo "Missing pilot overview TSV: $OVERVIEW_TSV" >&2; exit 2; }

FORCE_ARGS=()
for step in "${FORCE_STEPS[@]}"; do
  FORCE_ARGS+=(--force-step "$step")
done

run_step() {
  local step_name=$1
  shift
  local log_path="$LOG_DIR/${step_name}.log"
  echo "== $step_name =="
  if "$@" >"$log_path" 2>&1; then
    tail -n 20 "$log_path"
  else
    tail -n 40 "$log_path" >&2
    return 1
  fi
}

snapshot_summary() {
  local source_path=$1
  local target_name=$2
  cp "$source_path" "$SUMMARY_DIR/$target_name"
}

summary_step_status() {
  local summary_path=$1
  "$PYTHON_BIN" - "$summary_path" <<'PY'
import json
import sys
from pathlib import Path

summary = json.loads(Path(sys.argv[1]).read_text())
print(summary.get("run_step_status", "unknown"))
PY
}

LAST_PREPARE_STATUS=""

run_prepare_step() {
  local step_name=$1
  local run_dir=$2
  local manifest_path=$3
  local summary_name=$4
  local log_path="$LOG_DIR/${step_name}.log"
  local summary_path="$run_dir/clustering_pilot_real_case_summary.json"
  local exit_code

  echo "== $step_name =="
  set +e
  "$PYTHON_BIN" "$SCRIPT_DIR/run_clustering_pilot_real_case.py" \
    --run-dir "$run_dir" \
    --selection-manifest "$manifest_path" \
    --n-jobs "$N_JOBS" \
    "${FORCE_ARGS[@]}" >"$log_path" 2>&1
  exit_code=$?
  set -e

  if [[ -f "$log_path" ]]; then
    if [[ $exit_code -eq 0 ]]; then
      tail -n 20 "$log_path"
    else
      tail -n 40 "$log_path" >&2
    fi
  fi

  if [[ -f "$summary_path" ]]; then
    snapshot_summary "$summary_path" "$summary_name"
    LAST_PREPARE_STATUS=$(summary_step_status "$summary_path")
    echo "prepare status: $LAST_PREPARE_STATUS"
  else
    LAST_PREPARE_STATUS="missing_summary"
  fi

  if [[ $exit_code -ne 0 && "$LAST_PREPARE_STATUS" != "blocked_missing_selection_data" ]]; then
    return $exit_code
  fi

  return 0
}

run_step \
  01_build_manifests \
  "$PYTHON_BIN" "$SCRIPT_DIR/build_clustering_pilot_selection_manifests.py" \
  --overview-tsv "$OVERVIEW_TSV" \
  --output-dir "$MANIFEST_DIR"

cp "$OVERVIEW_TSV" "$MANIFEST_DIR/clustering_pilot_protein_overview.tsv"
snapshot_summary "$MANIFEST_DIR/clustering_pilot_manifest_summary.json" 01_manifest_generation.json

DOMAIN_MANIFEST="$MANIFEST_DIR/clustering_pilot_domain_only_selection.yaml"
FULL_MANIFEST="$MANIFEST_DIR/clustering_pilot_full_length_selection.yaml"
FULL_COUNT=$(
  "$PYTHON_BIN" - "$MANIFEST_DIR/clustering_pilot_manifest_summary.json" <<'PY'
import json
import sys
from pathlib import Path
summary = json.loads(Path(sys.argv[1]).read_text())
print(int(summary["n_full_length_selections"]))
PY
)

run_prepare_step \
  02_domain_prepare \
  "$DOMAIN_RUN_DIR" \
  "$DOMAIN_MANIFEST" \
  02_domain_prepare_summary.json
DOMAIN_PREPARE_STATUS="$LAST_PREPARE_STATUS"

if [[ "$PREPARE_ONLY" == false ]]; then
  [[ "$DOMAIN_PREPARE_STATUS" == "prepare_only" ]] || {
    echo "Domain-only prepare did not reach a runnable state: $DOMAIN_PREPARE_STATUS" >&2
    exit 1
  }
  run_step \
    03_domain_execute \
    "$PYTHON_BIN" "$SCRIPT_DIR/run_clustering_pilot_real_case.py" \
    --run-dir "$DOMAIN_RUN_DIR" \
    --selection-manifest "$DOMAIN_MANIFEST" \
    --execute \
    --n-jobs "$N_JOBS" \
    "${FORCE_ARGS[@]}"
  snapshot_summary "$DOMAIN_RUN_DIR/clustering_pilot_real_case_summary.json" 03_domain_execute_summary.json
fi

FULL_PREPARE_STATUS="not_requested"
if [[ "$FULL_COUNT" -gt 0 ]]; then
  run_prepare_step \
    04_full_length_prepare \
    "$FULL_RUN_DIR" \
    "$FULL_MANIFEST" \
    04_full_length_prepare_summary.json
  FULL_PREPARE_STATUS="$LAST_PREPARE_STATUS"

  if [[ "$PREPARE_ONLY" == false ]]; then
    [[ "$FULL_PREPARE_STATUS" == "prepare_only" ]] || {
      echo "Full-length prepare did not reach a runnable state: $FULL_PREPARE_STATUS" >&2
      exit 1
    }
    run_step \
      05_full_length_execute \
    "$PYTHON_BIN" "$SCRIPT_DIR/run_clustering_pilot_real_case.py" \
    --run-dir "$FULL_RUN_DIR" \
    --selection-manifest "$FULL_MANIFEST" \
    --execute \
    --n-jobs "$N_JOBS" \
    "${FORCE_ARGS[@]}"
    snapshot_summary "$FULL_RUN_DIR/clustering_pilot_real_case_summary.json" 05_full_length_execute_summary.json
  fi
fi

if [[ "$PREPARE_ONLY" == true ]]; then
  BLOCKED_PHASES=()
  [[ "$DOMAIN_PREPARE_STATUS" == "blocked_missing_selection_data" ]] && BLOCKED_PHASES+=("domain_only")
  [[ "$FULL_PREPARE_STATUS" == "blocked_missing_selection_data" ]] && BLOCKED_PHASES+=("full_length")

  if [[ ${#BLOCKED_PHASES[@]} -gt 0 ]]; then
    echo "Prepare-only audit found blocked phases: ${BLOCKED_PHASES[*]}" >&2
    exit 1
  fi
fi

echo "Pilot wrapper complete. Run root: $RUN_ROOT"
