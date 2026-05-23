#!/usr/bin/env bash
#SBATCH --job-name=analysis_full_pipeline
#SBATCH --output=logs/analysis_full_pipeline_%j.out
#SBATCH --error=logs/analysis_full_pipeline_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

set -euo pipefail

# Usage:
#   sbatch scripts/run_analysis_full_pipeline.slurm.sh \
#     --config configs/production.analysis_core.example.yaml \
#     --output /path/to/output \
#     --del del_a \
#     --metadata /path/to/protein_metadata.tsv \
#     [--runtime-paths configs/runtime_paths.yaml] \
#     [--n-jobs 16] \
#     [--predictive-task all] \
#     [--random-state 42]

RUNTIME_PATHS="configs/runtime_paths.yaml"
CONFIG_PATH=""
OUTPUT_DIR=""
DEL_BRANCH="del_a"
N_JOBS="${SLURM_CPUS_PER_TASK:-1}"
PROTEIN_METADATA=""
PREDICTIVE_TASK="all"
RANDOM_STATE="42"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --runtime-paths)
      RUNTIME_PATHS="$2"
      shift 2
      ;;
    --config)
      CONFIG_PATH="$2"
      shift 2
      ;;
    --output)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --del)
      DEL_BRANCH="$2"
      shift 2
      ;;
    --n-jobs)
      N_JOBS="$2"
      shift 2
      ;;
    --metadata)
      PROTEIN_METADATA="$2"
      shift 2
      ;;
    --predictive-task)
      PREDICTIVE_TASK="$2"
      shift 2
      ;;
    --random-state)
      RANDOM_STATE="$2"
      shift 2
      ;;
    *)
      echo "[ERROR] Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "$CONFIG_PATH" || -z "$OUTPUT_DIR" || -z "$PROTEIN_METADATA" ]]; then
  echo "[ERROR] Missing required args. Need --config, --output, --metadata." >&2
  exit 1
fi

if [[ ! -f "$RUNTIME_PATHS" ]]; then
  echo "[ERROR] runtime_paths file not found: $RUNTIME_PATHS" >&2
  exit 1
fi

if [[ ! -f "$CONFIG_PATH" ]]; then
  echo "[ERROR] config file not found: $CONFIG_PATH" >&2
  exit 1
fi

if [[ ! -f "$PROTEIN_METADATA" ]]; then
  echo "[ERROR] protein metadata file not found: $PROTEIN_METADATA" >&2
  exit 1
fi

mkdir -p logs
mkdir -p "$OUTPUT_DIR"

export LPMO_PIPELINE_RUNTIME_PATHS_CONFIG="$RUNTIME_PATHS"

PYTHON_BIN="$(awk -F': ' '/^[[:space:]]*python_executable:/ {gsub(/"/,"",$2); print $2; exit}' "$RUNTIME_PATHS")"
if [[ -z "$PYTHON_BIN" ]]; then
  echo "[ERROR] Could not resolve runtime_settings.python_executable from $RUNTIME_PATHS" >&2
  exit 1
fi

if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "[ERROR] python_executable is not executable: $PYTHON_BIN" >&2
  exit 1
fi

CONDITION_TABLE="$OUTPUT_DIR/condition_table.tsv"
CLUSTER_TABLE="$OUTPUT_DIR/cluster_table.tsv"
PREDICTIVE_SUMMARY="$OUTPUT_DIR/10_predictive/predictive_summary.json"
CBM_SUMMARY="$OUTPUT_DIR/15_cbm_paired_analysis/cbm_paired_analysis_summary.json"

echo "[INFO] Runtime paths: $RUNTIME_PATHS"
echo "[INFO] Python: $PYTHON_BIN"
echo "[INFO] Output: $OUTPUT_DIR"
echo "[INFO] DEL branch: $DEL_BRANCH"
echo "[INFO] n_jobs: $N_JOBS"

# 1) Main production analysis (no tuning)
"$PYTHON_BIN" -m lpmo_pipeline.cli run \
  --mode production \
  --config "$CONFIG_PATH" \
  --output "$OUTPUT_DIR" \
  --del "$DEL_BRANCH" \
  --n-jobs "$N_JOBS"

if [[ ! -f "$CONDITION_TABLE" ]]; then
  echo "[ERROR] condition_table.tsv missing after run: $CONDITION_TABLE" >&2
  exit 1
fi

# 2) Predictive postprocess
"$PYTHON_BIN" -m lpmo_pipeline.cli predictive \
  --condition-table "$CONDITION_TABLE" \
  --protein-metadata "$PROTEIN_METADATA" \
  --output "$OUTPUT_DIR" \
  --task "$PREDICTIVE_TASK" \
  --n-folds 5 \
  --random-state "$RANDOM_STATE"

if [[ ! -f "$PREDICTIVE_SUMMARY" ]]; then
  echo "[ERROR] predictive summary missing: $PREDICTIVE_SUMMARY" >&2
  exit 1
fi

# 3) CBM paired postprocess
CBM_ARGS=(
  --condition-table "$CONDITION_TABLE"
  --output "$OUTPUT_DIR"
  --random-state "$RANDOM_STATE"
  --protein-metadata "$PROTEIN_METADATA"
)
if [[ -f "$CLUSTER_TABLE" ]]; then
  CBM_ARGS+=(--cluster-table "$CLUSTER_TABLE")
else
  echo "[WARN] cluster_table.tsv not found, running cbm-paired without --cluster-table"
fi

"$PYTHON_BIN" -m lpmo_pipeline.cli cbm-paired "${CBM_ARGS[@]}"

if [[ ! -f "$CBM_SUMMARY" ]]; then
  echo "[ERROR] CBM paired summary missing: $CBM_SUMMARY" >&2
  exit 1
fi

echo "[INFO] Completed full pipeline run"
echo "[INFO] condition_table: $CONDITION_TABLE"
echo "[INFO] predictive_summary: $PREDICTIVE_SUMMARY"
echo "[INFO] cbm_summary: $CBM_SUMMARY"
