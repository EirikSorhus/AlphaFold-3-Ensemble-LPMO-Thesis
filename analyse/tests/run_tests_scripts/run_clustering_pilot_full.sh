#!/bin/bash
#SBATCH --job-name=run_tclustering_pilot_full
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_clustering_pilot_full_%j.log


set -euo pipefail

REPO_ROOT="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse"
SCRIPT_DIR="$REPO_ROOT/tests/run_tests_scripts"
PYTHON_BIN="/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin/python"
DEFAULT_RUN_ROOT="$REPO_ROOT/tests/tests_results/clustering_pilot_full"

RUN_ROOT="$DEFAULT_RUN_ROOT"
PREPARE_ONLY=false
SKIP_SENSITIVITY=false
ACCOUNT="nn1003k"
PARTITION="small"
SHARD_SIZE=5
BALANCE_BY="estimated_pose_count"
ARRAY_LIMIT=""
EXECUTE_CPUS=8
EXECUTE_MEM="24G"
EXECUTE_TIME="24:00:00"
SUMMARY_CPUS=1
SUMMARY_MEM="4G"
SUMMARY_TIME="00:30:00"
SENSITIVITY_CPUS=1
SENSITIVITY_MEM="8G"
SENSITIVITY_TIME="04:00:00"
SENSITIVITY_OUTPUT_DIR=""

usage() {
  cat <<'EOF'
Usage:
    tests/run_tests_scripts/run_clustering_pilot_full.sh [options]

This wrapper submits the current staged clustering pilot and then chains the
existing clustering sensitivity job with the same parameter combinations used
in the earlier sensitivity run:
  - HDBSCAN min_cluster_size: 3 5 10
  - HDBSCAN min_samples: none
  - Agglomerative distance_threshold: 0.35 0.45 0.55
  - Agglomerative min_cluster_size: 3 5 10

Options:
    --run-dir PATH              Output root for pilot and sensitivity artifacts.
    --prepare-only             Do not submit jobs; print the staged and sensitivity commands.
    --skip-sensitivity         Submit only the staged pilot, not the parameter sensitivity phase.
    --account NAME             Slurm account (default: nn1003k).
    --partition NAME           Slurm partition (default: small).
    --shard-size N             Protein selections per staged shard task (default: 5).
    --balance-by MODE          Shard balancing: protein_count|estimated_pose_count (default: estimated_pose_count).
    --array-limit N            Optional staged array concurrency cap.
    --n-jobs N                 Alias for --execute-cpus.
    --execute-cpus N           CPUs per staged shard task (default: 8).
    --execute-mem VALUE        Memory per staged shard task (default: 24G).
    --execute-time VALUE       Walltime per staged shard task (default: 24:00:00).
    --summary-cpus N           CPUs for staged summary job (default: 1).
    --summary-mem VALUE        Memory for staged summary job (default: 4G).
    --summary-time VALUE       Walltime for staged summary job (default: 00:30:00).
    --sensitivity-cpus N       CPUs for parameter sensitivity job (default: 1).
    --sensitivity-mem VALUE    Memory for parameter sensitivity job (default: 8G).
    --sensitivity-time VALUE   Walltime for parameter sensitivity job (default: 04:00:00).
    --sensitivity-output-dir PATH
                               Optional explicit output directory for sensitivity outputs.
    --help                     Show this help.
EOF
}

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
    --skip-sensitivity)
      SKIP_SENSITIVITY=true
      shift
      ;;
    --account)
      ACCOUNT="$2"
      shift 2
      ;;
    --partition)
      PARTITION="$2"
      shift 2
      ;;
    --shard-size)
      SHARD_SIZE="$2"
      shift 2
      ;;
    --balance-by)
      BALANCE_BY="$2"
      shift 2
      ;;
    --array-limit)
      ARRAY_LIMIT="$2"
      shift 2
      ;;
    --n-jobs|--execute-cpus)
      EXECUTE_CPUS="$2"
      shift 2
      ;;
    --execute-mem)
      EXECUTE_MEM="$2"
      shift 2
      ;;
    --execute-time)
      EXECUTE_TIME="$2"
      shift 2
      ;;
    --summary-cpus)
      SUMMARY_CPUS="$2"
      shift 2
      ;;
    --summary-mem)
      SUMMARY_MEM="$2"
      shift 2
      ;;
    --summary-time)
      SUMMARY_TIME="$2"
      shift 2
      ;;
    --sensitivity-cpus)
      SENSITIVITY_CPUS="$2"
      shift 2
      ;;
    --sensitivity-mem)
      SENSITIVITY_MEM="$2"
      shift 2
      ;;
    --sensitivity-time)
      SENSITIVITY_TIME="$2"
      shift 2
      ;;
    --sensitivity-output-dir)
      SENSITIVITY_OUTPUT_DIR=$(realpath -m "$2")
      shift 2
      ;;
    --help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ "$SHARD_SIZE" -lt 1 ]]; then
  echo "[ERROR] --shard-size must be >= 1" >&2
  exit 2
fi

cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT/src:${PYTHONPATH:-}"

LOG_DIR="$RUN_ROOT/logs"
SUMMARY_DIR="$RUN_ROOT/summaries"
STAGED_METADATA_JSON="$SUMMARY_DIR/staged_submission_metadata.json"
ORCHESTRATOR_METADATA_JSON="$SUMMARY_DIR/full_wrapper_submission_metadata.json"

mkdir -p "$LOG_DIR" "$SUMMARY_DIR"

staged_cmd=(
  bash "$SCRIPT_DIR/submit_clustering_pilot_staged.sh"
  --run-dir "$RUN_ROOT"
  --account "$ACCOUNT"
  --partition "$PARTITION"
  --shard-size "$SHARD_SIZE"
  --balance-by "$BALANCE_BY"
  --execute-cpus "$EXECUTE_CPUS"
  --execute-mem "$EXECUTE_MEM"
  --execute-time "$EXECUTE_TIME"
  --summary-cpus "$SUMMARY_CPUS"
  --summary-mem "$SUMMARY_MEM"
  --summary-time "$SUMMARY_TIME"
  --submission-metadata-json "$STAGED_METADATA_JSON"
)

if [[ -n "$ARRAY_LIMIT" ]]; then
  staged_cmd+=(--array-limit "$ARRAY_LIMIT")
fi

if [[ "$PREPARE_ONLY" == true ]]; then
  staged_cmd+=(--dry-run)
fi

echo "== staged_pilot_submit =="
staged_output=$("${staged_cmd[@]}")
printf '%s\n' "$staged_output" | tee "$LOG_DIR/staged_pilot_submit.log"

sensitivity_cmd=(
  sbatch --parsable
  --job-name=clust_param_sens
  --account="$ACCOUNT"
  --partition="$PARTITION"
  --time="$SENSITIVITY_TIME"
  --mem="$SENSITIVITY_MEM"
  --cpus-per-task="$SENSITIVITY_CPUS"
  --output="$LOG_DIR/parameter_sensitivity_%j.log"
  --export="ALL,RUN_ROOT=$RUN_ROOT"
  "$SCRIPT_DIR/run_clustering_parameter_sensitivity.slurm"
)

if [[ -n "$SENSITIVITY_OUTPUT_DIR" ]]; then
  sensitivity_cmd+=("$SENSITIVITY_OUTPUT_DIR")
fi

if [[ "$PREPARE_ONLY" == true ]]; then
  if [[ "$SKIP_SENSITIVITY" == false ]]; then
    echo "== parameter_sensitivity_submit =="
    printf '[DRY-RUN]'
    printf ' %q' "${sensitivity_cmd[@]}"
    printf '\n'
  fi
  cat > "$ORCHESTRATOR_METADATA_JSON" <<EOF
{
  "run_root": "$RUN_ROOT",
  "balance_by": "$BALANCE_BY",
  "prepare_only": true,
  "skip_sensitivity": $([[ "$SKIP_SENSITIVITY" == true ]] && echo true || echo false),
  "staged_submission_metadata_json": "$STAGED_METADATA_JSON",
  "sensitivity_output_dir": "${SENSITIVITY_OUTPUT_DIR:-}"
}
EOF
  echo "Pilot orchestration dry-run complete. Run root: $RUN_ROOT"
  exit 0
fi

SUMMARY_JOB=$("$PYTHON_BIN" - "$STAGED_METADATA_JSON" <<'PY'
import json
import sys
from pathlib import Path

payload = json.loads(Path(sys.argv[1]).read_text())
print(payload.get("summary_job", ""))
PY
)

if [[ -z "$SUMMARY_JOB" ]]; then
  echo "[ERROR] Missing summary job id in $STAGED_METADATA_JSON" >&2
  exit 1
fi

SENSITIVITY_JOB=""
if [[ "$SKIP_SENSITIVITY" == false ]]; then
  sensitivity_cmd=(
    sbatch --parsable
    --dependency=afterok:"$SUMMARY_JOB"
    --job-name=clust_param_sens
    --account="$ACCOUNT"
    --partition="$PARTITION"
    --time="$SENSITIVITY_TIME"
    --mem="$SENSITIVITY_MEM"
    --cpus-per-task="$SENSITIVITY_CPUS"
    --output="$LOG_DIR/parameter_sensitivity_%j.log"
    --export="ALL,RUN_ROOT=$RUN_ROOT"
    "$SCRIPT_DIR/run_clustering_parameter_sensitivity.slurm"
  )

  if [[ -n "$SENSITIVITY_OUTPUT_DIR" ]]; then
    sensitivity_cmd+=("$SENSITIVITY_OUTPUT_DIR")
  fi

  echo "== parameter_sensitivity_submit =="
  SENSITIVITY_JOB=$("${sensitivity_cmd[@]}")
  printf 'Submitted parameter sensitivity job: %s\n' "$SENSITIVITY_JOB" | tee "$LOG_DIR/parameter_sensitivity_submit.log"
fi

cat > "$ORCHESTRATOR_METADATA_JSON" <<EOF
{
  "run_root": "$RUN_ROOT",
  "balance_by": "$BALANCE_BY",
  "prepare_only": false,
  "skip_sensitivity": $([[ "$SKIP_SENSITIVITY" == true ]] && echo true || echo false),
  "staged_submission_metadata_json": "$STAGED_METADATA_JSON",
  "summary_job": "$SUMMARY_JOB",
  "sensitivity_job": "$SENSITIVITY_JOB",
  "sensitivity_output_dir": "${SENSITIVITY_OUTPUT_DIR:-}"
}
EOF

echo "Pilot orchestration complete. Run root: $RUN_ROOT"
echo "Staged pilot summary job: $SUMMARY_JOB"
if [[ "$SKIP_SENSITIVITY" == false ]]; then
  echo "Parameter sensitivity job: $SENSITIVITY_JOB"
fi
