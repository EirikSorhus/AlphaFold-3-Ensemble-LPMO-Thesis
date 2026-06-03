#!/usr/bin/env bash
#SBATCH --job-name=rerun_downstream_reuse
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=rerun_downstream_reuse_%j.out
#SBATCH --error=rerun_downstream_reuse_%j.err

# Submitter/wrapper for re-running analysis after prepare/QC.
#
# Usage from the analyse project root:
#   sbatch scripts/rerun_downstream_reuse_prepared.slurm.sh
#
# This short submitter job creates new shard configs and submits two Slurm arrays
# plus one dependent collect/postprocess job. Monitor with:
#   squeue -u "$USER"
#
# It does not delete, move, copy, or rewrite the old normalized/converted files.

set -euo pipefail

# ---------------------------------------------------------------------------
# User-adjustable paths and reuse policy
# ---------------------------------------------------------------------------
PROJECT_ROOT="${PROJECT_ROOT:-/cluster/work/projects/nn1003k/eirik/Masteroppgave_clean/analyse}"
OLD_RESULTS_DIR="${OLD_RESULTS_DIR:-$PROJECT_ROOT/results/2full_pipeline_array}"
NEW_RESULTS_DIR="${NEW_RESULTS_DIR:-$PROJECT_ROOT/results/2full_pipeline_array_rerun_downstream_${SLURM_JOB_ID:-manual}}"

REUSE_PREPARED_CASES="${REUSE_PREPARED_CASES:-true}"  # normalized.cif, for_posebusters.pdb, privateer_input.cif
REUSE_HARD_QC="${REUSE_HARD_QC:-true}"                # set false if QC thresholds or hard-QC geometry changed
REUSE_IFP_FOR_CONSTRUCTS="${REUSE_IFP_FOR_CONSTRUCTS:-domain_only}"
RECOMPUTE_IFP_FOR_CONSTRUCTS="${RECOMPUTE_IFP_FOR_CONSTRUCTS:-full_length}"
ALLOW_REUSE_FALLBACK_GENERATION="${ALLOW_REUSE_FALLBACK_GENERATION:-true}"
ALLOW_EXISTING_OUTPUT="${ALLOW_EXISTING_OUTPUT:-false}"
DRY_RUN="${DRY_RUN:-false}"

# Metadata/postprocess inputs. Defaults match the previous production run.
RUNTIME_PATHS="${RUNTIME_PATHS:-$PROJECT_ROOT/configs/runtime_paths.yaml}"
PROTEIN_METADATA="${PROTEIN_METADATA:-$PROJECT_ROOT/input_data/metadata_final_ec_fixed.tsv}"
CORE_FASTA="${CORE_FASTA:-$PROJECT_ROOT/input_data/lpmo_core_domain_2026-03-14_06-52-15_deduplicated.fasta}"
FAMILY_ALIGNMENT_DIR="${FAMILY_ALIGNMENT_DIR:-}"
MAFFT_EXECUTABLE="${MAFFT_EXECUTABLE:-mafft}"
PREDICTIVE_TASK="${PREDICTIVE_TASK:-all}"
RANDOM_STATE="${RANDOM_STATE:-42}"
RUN_GLOBAL_POSTPROCESS="${RUN_GLOBAL_POSTPROCESS:-true}"

# Slurm resources for the real work.
ACCOUNT="${ACCOUNT:-nn1003k}"
PARTITION="${PARTITION:-small}"
MAX_PARALLEL_SHARDS="${MAX_PARALLEL_SHARDS:-10}"
N_JOBS_PER_SHARD="${N_JOBS_PER_SHARD:-2}"
SHARD_CPUS="${SHARD_CPUS:-2}"
SHARD_MEM="${SHARD_MEM:-20G}"
SHARD_TIME="${SHARD_TIME:-05:00:00}"
SUMMARY_CPUS="${SUMMARY_CPUS:-2}"
SUMMARY_MEM="${SUMMARY_MEM:-20G}"
SUMMARY_TIME="${SUMMARY_TIME:-08:00:00}"

# ---------------------------------------------------------------------------
# Safety checks
# ---------------------------------------------------------------------------
if [[ "$REUSE_PREPARED_CASES" != "true" && "$REUSE_PREPARED_CASES" != "false" ]]; then
  echo "[ERROR] REUSE_PREPARED_CASES must be true or false" >&2
  exit 1
fi
if [[ "$REUSE_HARD_QC" != "true" && "$REUSE_HARD_QC" != "false" ]]; then
  echo "[ERROR] REUSE_HARD_QC must be true or false" >&2
  exit 1
fi
if [[ "$ALLOW_REUSE_FALLBACK_GENERATION" != "true" && "$ALLOW_REUSE_FALLBACK_GENERATION" != "false" ]]; then
  echo "[ERROR] ALLOW_REUSE_FALLBACK_GENERATION must be true or false" >&2
  exit 1
fi
if [[ "$RUN_GLOBAL_POSTPROCESS" != "true" && "$RUN_GLOBAL_POSTPROCESS" != "false" ]]; then
  echo "[ERROR] RUN_GLOBAL_POSTPROCESS must be true or false" >&2
  exit 1
fi
if [[ "$ALLOW_EXISTING_OUTPUT" != "true" && "$ALLOW_EXISTING_OUTPUT" != "false" ]]; then
  echo "[ERROR] ALLOW_EXISTING_OUTPUT must be true or false" >&2
  exit 1
fi

cd "$PROJECT_ROOT"

OLD_PLAN_JSON="$OLD_RESULTS_DIR/_plan/shard_plan.json"
for path in "$OLD_RESULTS_DIR" "$OLD_PLAN_JSON" "$RUNTIME_PATHS"; do
  if [[ ! -e "$path" ]]; then
    echo "[ERROR] Required path missing: $path" >&2
    exit 1
  fi
done

if [[ "$RUN_GLOBAL_POSTPROCESS" == "true" ]]; then
  for path in "$PROTEIN_METADATA" "$CORE_FASTA"; do
    if [[ ! -f "$path" ]]; then
      echo "[ERROR] Required postprocess input missing: $path" >&2
      exit 1
    fi
  done
fi

if [[ -d "$NEW_RESULTS_DIR" && -n "$(find "$NEW_RESULTS_DIR" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  if [[ "$ALLOW_EXISTING_OUTPUT" != "true" ]]; then
    echo "[ERROR] NEW_RESULTS_DIR already exists and is not empty: $NEW_RESULTS_DIR" >&2
    echo "[ERROR] Set ALLOW_EXISTING_OUTPUT=true only if you intend to add/overwrite downstream outputs there." >&2
    exit 1
  fi
fi

PYTHON_BIN="$(awk -F': ' '/^[[:space:]]*python_executable:/ {gsub(/"/,"",$2); print $2; exit}' "$RUNTIME_PATHS")"
if [[ -z "$PYTHON_BIN" || ! -x "$PYTHON_BIN" ]]; then
  echo "[ERROR] Could not resolve executable python from $RUNTIME_PATHS" >&2
  exit 1
fi

export PYTHONPATH="$PROJECT_ROOT/src:${PYTHONPATH:-}"
export LPMO_PIPELINE_RUNTIME_PATHS_CONFIG="$RUNTIME_PATHS"

mkdir -p "$NEW_RESULTS_DIR"
PLAN_DIR="$NEW_RESULTS_DIR/_plan"
CONFIG_OUT_DIR="$PLAN_DIR/generated_configs"
INDEX_DIR="$PLAN_DIR/shard_indexes"
JOB_SCRIPT_DIR="$PLAN_DIR/slurm_scripts"
LOG_DIR="$NEW_RESULTS_DIR/logs"
mkdir -p "$CONFIG_OUT_DIR" "$INDEX_DIR" "$JOB_SCRIPT_DIR" "$LOG_DIR"

NEW_PLAN_JSON="$PLAN_DIR/shard_plan.json"

echo "[INFO] Old results: $OLD_RESULTS_DIR"
echo "[INFO] New results: $NEW_RESULTS_DIR"
echo "[INFO] Reuse prepared cases: $REUSE_PREPARED_CASES"
echo "[INFO] Reuse hard QC: $REUSE_HARD_QC"
echo "[INFO] Domain-only: IFP reuse enabled when REUSE_IFP_FOR_CONSTRUCTS contains domain_only"
echo "[INFO] Full-length: IFP recompute enabled when RECOMPUTE_IFP_FOR_CONSTRUCTS contains full_length"
echo "[INFO] Fallback generation policy: $([[ "$ALLOW_REUSE_FALLBACK_GENERATION" == "true" ]] && echo 'allowed and logged' || echo 'disabled')"

# ---------------------------------------------------------------------------
# Build new shard configs from the old shard plan.
# ---------------------------------------------------------------------------
"$PYTHON_BIN" - <<PYEOF
import csv
import json
from pathlib import Path

import yaml

old_results_dir = Path("$OLD_RESULTS_DIR").resolve()
new_results_dir = Path("$NEW_RESULTS_DIR").resolve()
config_out_dir = Path("$CONFIG_OUT_DIR").resolve()
index_dir = Path("$INDEX_DIR").resolve()
reuse_prepared = "$REUSE_PREPARED_CASES" == "true"
reuse_hard_qc = "$REUSE_HARD_QC" == "true"
reuse_ifp_for_constructs = {value.strip() for value in "$REUSE_IFP_FOR_CONSTRUCTS".split(",") if value.strip()}
recompute_ifp_for_constructs = {value.strip() for value in "$RECOMPUTE_IFP_FOR_CONSTRUCTS".split(",") if value.strip()}
allow_reuse_fallback_generation = "$ALLOW_REUSE_FALLBACK_GENERATION" == "true"

old_plan = json.loads(Path("$OLD_PLAN_JSON").read_text())
new_plan = {"constructs": {}, "n_tasks": 0, "reuse": {
    "old_results_dir": str(old_results_dir),
    "reuse_prepared_cases": reuse_prepared,
    "reuse_hard_qc": reuse_hard_qc,
    "reuse_ifp_for_constructs": sorted(reuse_ifp_for_constructs),
    "recompute_ifp_for_constructs": sorted(recompute_ifp_for_constructs),
    "allow_reuse_fallback_generation": allow_reuse_fallback_generation,
}}

def _load_yaml(path: Path) -> dict:
    payload = yaml.safe_load(path.read_text())
    if not isinstance(payload, dict):
        raise TypeError(f"Config must be a mapping: {path}")
    return payload

def _validate_old_shard(old_output: Path) -> None:
    summary_path = old_output / "analysis_core_summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(f"missing old shard summary: {summary_path}")
    summary = json.loads(summary_path.read_text())
    cases = summary.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError(f"old shard summary has no cases: {summary_path}")
    if reuse_prepared:
        prep_errors = [str(case.get("pose_id")) for case in cases if case.get("status") == "prep_error"]
        if prep_errors:
            preview = ", ".join(prep_errors[:10])
            raise RuntimeError(
                f"{old_output} contains prep_error cases and cannot be fully reused: {preview}"
            )
        for case in cases:
            pose_id = str(case.get("pose_id") or "")
            for key in ("normalized_cif", "posebusters_pdb", "privateer_input_cif"):
                value = str(case.get(key) or "")
                if not value:
                    raise FileNotFoundError(f"{old_output}: {pose_id} missing {key} in summary")
                path = Path(value)
                if not path.is_file():
                    raise FileNotFoundError(f"{old_output}: {pose_id} missing {key}: {path}")
    if reuse_hard_qc and not (old_output / "qc_report.json").is_file():
        raise FileNotFoundError(f"missing old hard-QC report: {old_output / 'qc_report.json'}")
    if construct_type in reuse_ifp_for_constructs and not (old_output / "pose_ifp_table.tsv").is_file():
        raise FileNotFoundError(f"missing old IFP table for reuse: {old_output / 'pose_ifp_table.tsv'}")

for construct_type, construct in old_plan.get("constructs", {}).items():
    tasks = []
    cfg_dir = config_out_dir / construct_type
    cfg_dir.mkdir(parents=True, exist_ok=True)
    for old_task in construct.get("tasks", []):
        shard_id = old_task["shard_id"]
        old_output = Path(old_task["output_dir"]).resolve()
        _validate_old_shard(old_output)

        cfg = _load_yaml(Path(old_task["config_path"]))
        production = cfg.setdefault("production", {})
        if not isinstance(production, dict):
            raise TypeError(f"production must be a mapping in {old_task['config_path']}")
        production["construct_type"] = construct_type
        production["include_proteins"] = list(old_task.get("proteins") or [])
        production["run_id"] = f"{production.get('run_id') or 'analysis'}_reuse_{construct_type}_{shard_id}"
        if reuse_prepared:
            production["reuse_prepared_cases_from"] = str(old_output)
        else:
            production.pop("reuse_prepared_cases_from", None)
        if reuse_hard_qc:
            production["reuse_hard_qc_report_path"] = str(old_output / "qc_report.json")
        else:
            production.pop("reuse_hard_qc_report_path", None)
        production["reuse_analysis_artifacts_from"] = str(old_output)
        production["allow_reuse_fallback_generation"] = allow_reuse_fallback_generation
        if construct_type in reuse_ifp_for_constructs:
            production["reuse_ifp_table_path"] = str(old_output / "pose_ifp_table.tsv")
        else:
            production.pop("reuse_ifp_table_path", None)

        new_cfg_path = cfg_dir / f"{shard_id}.yaml"
        new_cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
        new_output = new_results_dir / construct_type / "shards" / shard_id
        task = {
            **old_task,
            "config_path": str(new_cfg_path),
            "output_dir": str(new_output),
            "reuse_prepared_cases_from": str(old_output) if reuse_prepared else "",
            "reuse_hard_qc_report_path": str(old_output / "qc_report.json") if reuse_hard_qc else "",
            "reuse_ifp_table_path": str(old_output / "pose_ifp_table.tsv") if construct_type in reuse_ifp_for_constructs else "",
            "reuse_analysis_artifacts_from": str(old_output),
            "allow_reuse_fallback_generation": allow_reuse_fallback_generation,
        }
        tasks.append(task)

    index_path = index_dir / f"{construct_type}_shard_index.tsv"
    with index_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["task_index", "shard_id", "del_branch", "config_path", "output_dir", "n_proteins"])
        for task in tasks:
            writer.writerow([
                task["task_index"],
                task["shard_id"],
                task["del_branch"],
                task["config_path"],
                task["output_dir"],
                task.get("n_proteins", ""),
            ])

    new_plan["constructs"][construct_type] = {
        **{key: value for key, value in construct.items() if key != "tasks"},
        "index_tsv": str(index_path),
        "tasks": tasks,
        "n_shards": len(tasks),
    }
    new_plan["n_tasks"] += len(tasks)

Path("$NEW_PLAN_JSON").write_text(json.dumps(new_plan, indent=2))
print(json.dumps({
    "n_tasks": new_plan["n_tasks"],
    "domain_only_shards": new_plan["constructs"].get("domain_only", {}).get("n_shards", 0),
    "full_length_shards": new_plan["constructs"].get("full_length", {}).get("n_shards", 0),
}, indent=2))
PYEOF

DOMAIN_INDEX_TSV="$INDEX_DIR/domain_only_shard_index.tsv"
FULL_INDEX_TSV="$INDEX_DIR/full_length_shard_index.tsv"

DOMAIN_SHARDS=$("$PYTHON_BIN" - <<PYEOF
import json
from pathlib import Path
plan = json.loads(Path("$NEW_PLAN_JSON").read_text())
print(plan["constructs"].get("domain_only", {}).get("n_shards", 0))
PYEOF
)
FULL_SHARDS=$("$PYTHON_BIN" - <<PYEOF
import json
from pathlib import Path
plan = json.loads(Path("$NEW_PLAN_JSON").read_text())
print(plan["constructs"].get("full_length", {}).get("n_shards", 0))
PYEOF
)

if [[ "$DOMAIN_SHARDS" -lt 1 && "$FULL_SHARDS" -lt 1 ]]; then
  echo "[ERROR] No shard tasks generated from old plan." >&2
  exit 1
fi

ARRAY_WORKER_SCRIPT="$JOB_SCRIPT_DIR/run_reuse_shard_array_task.sh"
cat > "$ARRAY_WORKER_SCRIPT" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

INDEX_TSV="$1"
REPO_ROOT="$2"
RUNTIME_PATHS="$3"
N_JOBS_PER_SHARD="$4"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT/src:${PYTHONPATH:-}"
export LPMO_PIPELINE_RUNTIME_PATHS_CONFIG="$RUNTIME_PATHS"

PYTHON_BIN="$(awk -F': ' '/^[[:space:]]*python_executable:/ {gsub(/"/,"",$2); print $2; exit}' "$RUNTIME_PATHS")"
if [[ -z "$PYTHON_BIN" || ! -x "$PYTHON_BIN" ]]; then
  echo "[ERROR] Could not resolve executable python from runtime paths: $RUNTIME_PATHS" >&2
  exit 2
fi

read -r SHARD_ID DEL_BRANCH CONFIG_PATH SHARD_OUTPUT N_PROTEINS < <(
  awk -F '\t' -v task="$TASK_ID" 'NR > 1 && $1 == task {print $2, $3, $4, $5, $6}' "$INDEX_TSV"
)
if [[ -z "${SHARD_ID:-}" || -z "${CONFIG_PATH:-}" || -z "${SHARD_OUTPUT:-}" ]]; then
  echo "[ERROR] No shard mapping for array task $TASK_ID in $INDEX_TSV" >&2
  exit 2
fi

mkdir -p "$SHARD_OUTPUT"
mkdir -p "$(dirname "$SHARD_OUTPUT")/../logs"
SHARD_LOG="$(dirname "$SHARD_OUTPUT")/../logs/${SHARD_ID}.log"

echo "[INFO] START reuse task=$TASK_ID shard=$SHARD_ID del=$DEL_BRANCH n_proteins=${N_PROTEINS:-unknown}" | tee -a "$SHARD_LOG"
"$PYTHON_BIN" -m lpmo_pipeline.cli run \
  --mode production \
  --config "$CONFIG_PATH" \
  --output "$SHARD_OUTPUT" \
  --del "$DEL_BRANCH" \
  --n-jobs "$N_JOBS_PER_SHARD" >>"$SHARD_LOG" 2>&1

if [[ ! -f "$SHARD_OUTPUT/cluster_table.tsv" ]]; then
  echo "[ERROR] Missing cluster_table.tsv for shard=$SHARD_ID" | tee -a "$SHARD_LOG" >&2
  exit 1
fi

echo "[INFO] DONE reuse task=$TASK_ID shard=$SHARD_ID" | tee -a "$SHARD_LOG"
EOF
chmod +x "$ARRAY_WORKER_SCRIPT"

array_spec() {
  local count="$1"
  printf '1-%s%%%s' "$count" "$MAX_PARALLEL_SHARDS"
}

submit_or_print_array() {
  local construct_type="$1"
  local count="$2"
  local index_tsv="$3"
  if [[ "$count" -lt 1 ]]; then
    echo ""
    return 0
  fi

  local cmd=(
    sbatch --parsable
    --array="$(array_spec "$count")"
    --job-name="reuse_${construct_type}_shards"
    --account="$ACCOUNT"
    --partition="$PARTITION"
    --time="$SHARD_TIME"
    --mem="$SHARD_MEM"
    --cpus-per-task="$SHARD_CPUS"
    --output="$LOG_DIR/${construct_type}_%A_%a.out"
    --error="$LOG_DIR/${construct_type}_%A_%a.err"
    "$ARRAY_WORKER_SCRIPT"
    "$index_tsv"
    "$PROJECT_ROOT"
    "$RUNTIME_PATHS"
    "$N_JOBS_PER_SHARD"
  )

  if [[ "$DRY_RUN" == "true" ]]; then
    {
      printf '[DRY-RUN]'
      printf ' %q' "${cmd[@]}"
      printf '\n'
    } >&2
    echo "dry_run_${construct_type}"
  else
    "${cmd[@]}"
  fi
}

COLLECT_SCRIPT="$JOB_SCRIPT_DIR/collect_and_postprocess.sh"
cat > "$COLLECT_SCRIPT" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$1"
RUNTIME_PATHS="$2"
OUTPUT_DIR="$3"
SHARD_PLAN_JSON="$4"
RUN_GLOBAL_POSTPROCESS="$5"
PROTEIN_METADATA="$6"
CORE_FASTA="$7"
PREDICTIVE_TASK="$8"
RANDOM_STATE="$9"
FAMILY_ALIGNMENT_DIR="${10}"
MAFFT_EXECUTABLE="${11}"

cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT/src:${PYTHONPATH:-}"
export LPMO_PIPELINE_RUNTIME_PATHS_CONFIG="$RUNTIME_PATHS"

PYTHON_BIN="$(awk -F': ' '/^[[:space:]]*python_executable:/ {gsub(/"/,"",$2); print $2; exit}' "$RUNTIME_PATHS")"
if [[ -z "$PYTHON_BIN" || ! -x "$PYTHON_BIN" ]]; then
  echo "[ERROR] Could not resolve executable python from runtime paths: $RUNTIME_PATHS" >&2
  exit 2
fi

MERGE_SUMMARY_JSON="$OUTPUT_DIR/merge_summary.json"

"$PYTHON_BIN" - <<PYEOF
import csv
import json
from pathlib import Path
from typing import Any

plan = json.loads(Path("$SHARD_PLAN_JSON").read_text())
output_dir = Path("$OUTPUT_DIR").resolve()
surfaces = [
    "condition_table.tsv",
    "cluster_table.tsv",
    "protein_condition_residue_scores.tsv",
    "protein_residue_regio_delta.tsv",
    "cluster_residue_signature.tsv",
    "cluster_ifp_signature.tsv",
    "pose_geometry.tsv",
    "pose_ifp_table.tsv",
    "pose_residue_contact_table.tsv",
    "cluster_assignments.tsv",
    "medoid_manifest.tsv",
    "condition_cluster_summary.tsv",
    "condition_patch_summary.tsv",
    "protein_patch_summary.tsv",
    "protein_summary_table.tsv",
    "crystal_anchor_table.tsv",
    "crystal_geometry_table.tsv",
]

def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})

summary: dict[str, Any] = {"constructs": {}, "combined": {}}
for construct_type, construct in sorted(plan.get("constructs", {}).items()):
    tasks = construct.get("tasks", [])
    merged_root = output_dir / construct_type / "merged"
    construct_summary = {}
    for surface in surfaces:
        rows: list[dict[str, Any]] = []
        source_files = []
        for task in tasks:
            src = Path(task["output_dir"]) / surface
            if not src.exists():
                continue
            source_files.append(str(src))
            for row in _read_tsv(src):
                rows.append({"_construct_type": construct_type, "_shard_id": task["shard_id"], **row})
        out_path = merged_root / surface
        _write_tsv(out_path, rows)
        construct_summary[surface] = {
            "merged_path": str(out_path) if out_path.exists() else "",
            "n_rows": len(rows),
            "n_source_files": len(source_files),
        }
    summary["constructs"][construct_type] = construct_summary

combined_root = output_dir / "combined"
for surface in surfaces:
    rows: list[dict[str, Any]] = []
    source_files = []
    for construct_type in ("domain_only", "full_length"):
        src = output_dir / construct_type / "merged" / surface
        if not src.exists():
            continue
        source_files.append(str(src))
        rows.extend(_read_tsv(src))
    out_path = combined_root / surface
    _write_tsv(out_path, rows)
    summary["combined"][surface] = {
        "merged_path": str(out_path) if out_path.exists() else "",
        "n_rows": len(rows),
        "n_source_files": len(source_files),
    }

Path("$MERGE_SUMMARY_JSON").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
PYEOF

COMBINED_CONDITION_TABLE="$OUTPUT_DIR/combined/condition_table.tsv"
COMBINED_CLUSTER_TABLE="$OUTPUT_DIR/combined/cluster_table.tsv"
COMBINED_CLUSTER_RESIDUE_SIGNATURE="$OUTPUT_DIR/combined/cluster_residue_signature.tsv"
COMBINED_STAGE16B_SCORES="$OUTPUT_DIR/combined/protein_condition_residue_scores.tsv"
COMBINED_STAGE16B_DELTA="$OUTPUT_DIR/combined/protein_residue_regio_delta.tsv"

if [[ ! -f "$COMBINED_CONDITION_TABLE" ]]; then
  echo "[ERROR] Missing merged combined condition table: $COMBINED_CONDITION_TABLE" >&2
  exit 1
fi

POSTPROCESS_OUTPUT="$OUTPUT_DIR/postprocess"
mkdir -p "$POSTPROCESS_OUTPUT"

if [[ "$RUN_GLOBAL_POSTPROCESS" == "true" ]]; then
  "$PYTHON_BIN" -m lpmo_pipeline.cli predictive \
    --condition-table "$COMBINED_CONDITION_TABLE" \
    --protein-metadata "$PROTEIN_METADATA" \
    --output "$POSTPROCESS_OUTPUT" \
    --task "$PREDICTIVE_TASK" \
    --n-folds 5 \
    --random-state "$RANDOM_STATE"

  CBM_ARGS=(
    --condition-table "$COMBINED_CONDITION_TABLE"
    --output "$POSTPROCESS_OUTPUT"
    --random-state "$RANDOM_STATE"
    --protein-metadata "$PROTEIN_METADATA"
  )
  if [[ -f "$COMBINED_CLUSTER_TABLE" ]]; then
    CBM_ARGS+=(--cluster-table "$COMBINED_CLUSTER_TABLE")
  fi
  if [[ -f "$COMBINED_CLUSTER_RESIDUE_SIGNATURE" ]]; then
    CBM_ARGS+=(--cluster-residue-signature-table "$COMBINED_CLUSTER_RESIDUE_SIGNATURE")
  fi
  "$PYTHON_BIN" -m lpmo_pipeline.cli cbm-paired "${CBM_ARGS[@]}"

  if [[ -f "$COMBINED_STAGE16B_SCORES" && -f "$COMBINED_STAGE16B_DELTA" ]]; then
    FAMILY_ARGS=(
      --protein-condition-residue-scores "$COMBINED_STAGE16B_SCORES"
      --protein-residue-regio-delta "$COMBINED_STAGE16B_DELTA"
      --protein-metadata "$PROTEIN_METADATA"
      --core-fasta "$CORE_FASTA"
      --output "$POSTPROCESS_OUTPUT"
      --families AA9 AA10
      --mafft-executable "$MAFFT_EXECUTABLE"
    )
    if [[ -n "$FAMILY_ALIGNMENT_DIR" ]]; then
      FAMILY_ARGS+=(--alignment-dir "$FAMILY_ALIGNMENT_DIR")
    fi
    "$PYTHON_BIN" -m lpmo_pipeline.cli family-enrichment "${FAMILY_ARGS[@]}"
  else
    echo "[WARN] Missing combined Stage 16b inputs; skipping family-enrichment"
  fi
fi

echo "[INFO] Completed collect/postprocess job"
echo "[INFO] Merge summary: $MERGE_SUMMARY_JSON"
EOF
chmod +x "$COLLECT_SCRIPT"

DOMAIN_ARRAY_JOB=$(submit_or_print_array "domain_only" "$DOMAIN_SHARDS" "$DOMAIN_INDEX_TSV")
FULL_ARRAY_JOB=$(submit_or_print_array "full_length" "$FULL_SHARDS" "$FULL_INDEX_TSV")

SUMMARY_DEPENDENCIES=()
[[ -n "$DOMAIN_ARRAY_JOB" ]] && SUMMARY_DEPENDENCIES+=("$DOMAIN_ARRAY_JOB")
[[ -n "$FULL_ARRAY_JOB" ]] && SUMMARY_DEPENDENCIES+=("$FULL_ARRAY_JOB")

SUMMARY_JOB=""
if [[ ${#SUMMARY_DEPENDENCIES[@]} -gt 0 ]]; then
  dependency_arg=$(IFS=:; echo "${SUMMARY_DEPENDENCIES[*]}")
  summary_cmd=(
    sbatch --parsable
    --dependency=afterok:"$dependency_arg"
    --job-name=reuse_collect_postprocess
    --account="$ACCOUNT"
    --partition="$PARTITION"
    --time="$SUMMARY_TIME"
    --mem="$SUMMARY_MEM"
    --cpus-per-task="$SUMMARY_CPUS"
    --output="$LOG_DIR/collect_%j.out"
    --error="$LOG_DIR/collect_%j.err"
    "$COLLECT_SCRIPT"
    "$PROJECT_ROOT"
    "$RUNTIME_PATHS"
    "$NEW_RESULTS_DIR"
    "$NEW_PLAN_JSON"
    "$RUN_GLOBAL_POSTPROCESS"
    "$PROTEIN_METADATA"
    "$CORE_FASTA"
    "$PREDICTIVE_TASK"
    "$RANDOM_STATE"
    "$FAMILY_ALIGNMENT_DIR"
    "$MAFFT_EXECUTABLE"
  )

  if [[ "$DRY_RUN" == "true" ]]; then
    {
      printf '[DRY-RUN]'
      printf ' %q' "${summary_cmd[@]}"
      printf '\n'
    } >&2
    SUMMARY_JOB="dry_run_collect"
  else
    SUMMARY_JOB=$("${summary_cmd[@]}")
  fi
fi

cat > "$NEW_RESULTS_DIR/submission_metadata.json" <<EOF
{
  "old_results_dir": "$(realpath -m "$OLD_RESULTS_DIR")",
  "new_results_dir": "$(realpath -m "$NEW_RESULTS_DIR")",
  "reuse_prepared_cases": $([[ "$REUSE_PREPARED_CASES" == "true" ]] && echo true || echo false),
  "reuse_hard_qc": $([[ "$REUSE_HARD_QC" == "true" ]] && echo true || echo false),
  "reuse_ifp_for_constructs": "$REUSE_IFP_FOR_CONSTRUCTS",
  "recompute_ifp_for_constructs": "$RECOMPUTE_IFP_FOR_CONSTRUCTS",
  "allow_reuse_fallback_generation": $([[ "$ALLOW_REUSE_FALLBACK_GENERATION" == "true" ]] && echo true || echo false),
  "runtime_paths": "$(realpath -m "$RUNTIME_PATHS")",
  "protein_metadata": "$(realpath -m "$PROTEIN_METADATA")",
  "core_fasta": "$(realpath -m "$CORE_FASTA")",
  "run_global_postprocess": $([[ "$RUN_GLOBAL_POSTPROCESS" == "true" ]] && echo true || echo false),
  "predictive_task": "$PREDICTIVE_TASK",
  "random_state": $RANDOM_STATE,
  "domain_shards": $DOMAIN_SHARDS,
  "full_shards": $FULL_SHARDS,
  "domain_array_job": "${DOMAIN_ARRAY_JOB:-}",
  "full_array_job": "${FULL_ARRAY_JOB:-}",
  "collect_job": "${SUMMARY_JOB:-}",
  "dry_run": $([[ "$DRY_RUN" == "true" ]] && echo true || echo false)
}
EOF

cat <<EOF
Submitted downstream rerun with prepare/QC reuse:
  old_results_dir:       $(realpath -m "$OLD_RESULTS_DIR")
  new_results_dir:       $(realpath -m "$NEW_RESULTS_DIR")
  reuse_prepared_cases:  $REUSE_PREPARED_CASES
  reuse_hard_qc:         $REUSE_HARD_QC
  reuse_ifp_for:         $REUSE_IFP_FOR_CONSTRUCTS
  recompute_ifp_for:     $RECOMPUTE_IFP_FOR_CONSTRUCTS
  fallback_generation:   $([[ "$ALLOW_REUSE_FALLBACK_GENERATION" == "true" ]] && echo 'allowed and logged' || echo 'disabled')
  domain_shards:         $DOMAIN_SHARDS
  full_shards:           $FULL_SHARDS
  domain_array_job:      ${DOMAIN_ARRAY_JOB:-none}
  full_array_job:        ${FULL_ARRAY_JOB:-none}
  collect_job:           ${SUMMARY_JOB:-none}
  generated_plan:        $NEW_PLAN_JSON

Monitor:
  squeue -u "$USER"
EOF
