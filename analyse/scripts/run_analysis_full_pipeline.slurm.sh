#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

RUNTIME_PATHS="configs/runtime_paths.yaml"
BASE_CONFIG=""
CORE_CONFIG=""
FULL_CONFIG=""
OUTPUT_DIR=""
PROTEIN_METADATA=""
CORE_FASTA="input_data/lpmo_core_domain_2026-03-14_06-52-15_deduplicated.fasta"
PROTEINS_PER_SHARD="5"
MAX_PARALLEL_SHARDS="10"
N_JOBS_PER_SHARD="2"
PREDICTIVE_TASK="all"
RANDOM_STATE="42"
SHARD_CPUS="2"
SHARD_MEM="32G"
SHARD_TIME="48:00:00"
SUMMARY_CPUS="2"
SUMMARY_MEM="12G"
SUMMARY_TIME="12:00:00"
RUN_GLOBAL_POSTPROCESS="true"
FAMILY_ALIGNMENT_DIR=""
MAFFT_EXECUTABLE="mafft"
ACCOUNT="nn1003k"
PARTITION="small"
SUBMISSION_METADATA_JSON=""
DRY_RUN="false"

usage() {
  cat <<'EOF'
Usage:
  scripts/run_analysis_full_pipeline.slurm.sh \
    --config configs/production.analysis_core.example.yaml \
    --output /path/to/output \
    [--core-config configs/core.yaml] \
    [--full-config configs/full.yaml] \
    [--runtime-paths configs/runtime_paths.yaml] \
    [--proteins-per-shard 5] \
    [--max-parallel-shards 10] \
    [--n-jobs-per-shard 2] \
    [--predictive-task all] \
    [--random-state 42] \
    [--shard-cpus 2] \
    [--shard-mem 32G] \
    [--shard-time 48:00:00] \
    [--summary-cpus 2] \
    [--summary-mem 12G] \
    [--summary-time 12:00:00] \
    [--run-global-postprocess true|false] \
    [--family-alignment-dir /path/to/alignments] \
    [--mafft-executable mafft] \
    [--account nn1003k] \
    [--partition small] \
    [--submission-metadata-json /path/to/submission.json] \
    [--dry-run]

Notes:
  - This script is a true Slurm-array submitter (not a single monolithic sbatch run).
  - It generates per-shard configs for domain_only and full_length and submits one array per construct.
  - The array workers run: lpmo_pipeline.cli run --mode production ...
  - A dependent summary job runs after arrays and performs merge + optional postprocess.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --runtime-paths)
      RUNTIME_PATHS="$2"
      shift 2
      ;;
    --config)
      BASE_CONFIG="$2"
      shift 2
      ;;
    --core-config)
      CORE_CONFIG="$2"
      shift 2
      ;;
    --full-config)
      FULL_CONFIG="$2"
      shift 2
      ;;
    --output)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --metadata)
      PROTEIN_METADATA="$2"
      shift 2
      ;;
    --core-fasta)
      CORE_FASTA="$2"
      shift 2
      ;;
    --proteins-per-shard)
      PROTEINS_PER_SHARD="$2"
      shift 2
      ;;
    --max-parallel-shards)
      MAX_PARALLEL_SHARDS="$2"
      shift 2
      ;;
    --n-jobs-per-shard)
      N_JOBS_PER_SHARD="$2"
      shift 2
      ;;
    --shard-cpus)
      SHARD_CPUS="$2"
      shift 2
      ;;
    --shard-mem)
      SHARD_MEM="$2"
      shift 2
      ;;
    --shard-time)
      SHARD_TIME="$2"
      shift 2
      ;;
    --account)
      ACCOUNT="$2"
      shift 2
      ;;
    --partition)
      PARTITION="$2"
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
    --run-global-postprocess)
      RUN_GLOBAL_POSTPROCESS="$2"
      shift 2
      ;;
    --family-alignment-dir)
      FAMILY_ALIGNMENT_DIR="$2"
      shift 2
      ;;
    --mafft-executable)
      MAFFT_EXECUTABLE="$2"
      shift 2
      ;;
    --submission-metadata-json)
      SUBMISSION_METADATA_JSON="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN="true"
      shift
      ;;
    --help)
      usage
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ -z "$OUTPUT_DIR" ]]; then
  echo "[ERROR] Missing required arg: --output" >&2
  exit 1
fi

if [[ -z "$BASE_CONFIG" && ( -z "$CORE_CONFIG" || -z "$FULL_CONFIG" ) ]]; then
  echo "[ERROR] Provide either --config OR both --core-config and --full-config." >&2
  exit 1
fi

if [[ -n "$BASE_CONFIG" ]]; then
  CORE_CONFIG="${CORE_CONFIG:-$BASE_CONFIG}"
  FULL_CONFIG="${FULL_CONFIG:-$BASE_CONFIG}"
fi

if [[ -z "$SHARD_CPUS" ]]; then
  SHARD_CPUS="$N_JOBS_PER_SHARD"
fi

if ! [[ "$PROTEINS_PER_SHARD" =~ ^[0-9]+$ ]] || [[ "$PROTEINS_PER_SHARD" -lt 1 ]]; then
  echo "[ERROR] --proteins-per-shard must be a positive integer" >&2
  exit 1
fi
if ! [[ "$MAX_PARALLEL_SHARDS" =~ ^[0-9]+$ ]] || [[ "$MAX_PARALLEL_SHARDS" -lt 1 ]]; then
  echo "[ERROR] --max-parallel-shards must be a positive integer" >&2
  exit 1
fi
if ! [[ "$N_JOBS_PER_SHARD" =~ ^[0-9]+$ ]] || [[ "$N_JOBS_PER_SHARD" -lt 1 ]]; then
  echo "[ERROR] --n-jobs-per-shard must be a positive integer" >&2
  exit 1
fi
if ! [[ "$SHARD_CPUS" =~ ^[0-9]+$ ]] || [[ "$SHARD_CPUS" -lt 1 ]]; then
  echo "[ERROR] --shard-cpus must be a positive integer" >&2
  exit 1
fi
if ! [[ "$SUMMARY_CPUS" =~ ^[0-9]+$ ]] || [[ "$SUMMARY_CPUS" -lt 1 ]]; then
  echo "[ERROR] --summary-cpus must be a positive integer" >&2
  exit 1
fi
if [[ "$RUN_GLOBAL_POSTPROCESS" != "true" && "$RUN_GLOBAL_POSTPROCESS" != "false" ]]; then
  echo "[ERROR] --run-global-postprocess must be true or false" >&2
  exit 1
fi

cd "$REPO_ROOT"

if [[ ! -f "$RUNTIME_PATHS" ]]; then
  echo "[ERROR] runtime_paths file not found: $RUNTIME_PATHS" >&2
  exit 1
fi
if [[ ! -f "$CORE_CONFIG" ]]; then
  echo "[ERROR] core config file not found: $CORE_CONFIG" >&2
  exit 1
fi
if [[ ! -f "$FULL_CONFIG" ]]; then
  echo "[ERROR] full config file not found: $FULL_CONFIG" >&2
  exit 1
fi
if [[ "$RUN_GLOBAL_POSTPROCESS" == "true" ]]; then
  if [[ -z "$PROTEIN_METADATA" || ! -f "$PROTEIN_METADATA" ]]; then
    echo "[ERROR] --metadata is required for postprocess and must point to an existing file" >&2
    exit 1
  fi
  if [[ -z "$CORE_FASTA" || ! -f "$CORE_FASTA" ]]; then
    echo "[ERROR] --core-fasta is required for postprocess and must point to an existing file" >&2
    exit 1
  fi
fi

export LPMO_PIPELINE_RUNTIME_PATHS_CONFIG="$RUNTIME_PATHS"

PYTHON_BIN="$(awk -F': ' '/^[[:space:]]*python_executable:/ {gsub(/"/,"",$2); print $2; exit}' "$RUNTIME_PATHS")"
if [[ -z "$PYTHON_BIN" ]]; then
  echo "[ERROR] Could not resolve python_executable from $RUNTIME_PATHS" >&2
  exit 1
fi
if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "[ERROR] python_executable is not executable: $PYTHON_BIN" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"
PLAN_DIR="$OUTPUT_DIR/_plan"
CONFIG_OUT_DIR="$PLAN_DIR/generated_configs"
INDEX_DIR="$PLAN_DIR/shard_indexes"
JOB_SCRIPT_DIR="$PLAN_DIR/slurm_scripts"
LOG_DIR="$OUTPUT_DIR/logs"
mkdir -p "$PLAN_DIR" "$CONFIG_OUT_DIR" "$INDEX_DIR" "$JOB_SCRIPT_DIR" "$LOG_DIR"

SHARD_PLAN_JSON="$PLAN_DIR/shard_plan.json"

"$PYTHON_BIN" - <<PYEOF
import copy
import csv
import json
import re
from pathlib import Path

import yaml

core_config_path = Path("$CORE_CONFIG").resolve()
full_config_path = Path("$FULL_CONFIG").resolve()
output_dir = Path("$OUTPUT_DIR").resolve()
config_out_dir = Path("$CONFIG_OUT_DIR").resolve()
index_dir = Path("$INDEX_DIR").resolve()
proteins_per_shard = int("$PROTEINS_PER_SHARD")


def _load(path: Path) -> dict:
    payload = yaml.safe_load(path.read_text())
    if not isinstance(payload, dict):
        raise TypeError(f"Config must be a mapping: {path}")
    return payload


def _dedupe(values):
    out = []
    seen = set()
    for value in values:
        s = str(value).strip()
        if not s or s in seen:
            continue
        seen.add(s)
        out.append(s)
    return out


def _as_list(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return list(value)
    return [value]


def _chunk(values, chunk_size):
    if not values:
        return []
    return [values[i:i + chunk_size] for i in range(0, len(values), chunk_size)]


def _discover_proteins(production: dict) -> list[str]:
    work_roots = production.get("work_roots") or {}
    construct_type = str(production.get("construct_type") or "").strip()

    work_root = ""
    if isinstance(work_roots, dict):
        work_root = str(work_roots.get(construct_type) or "").strip()
    if not work_root:
        work_root = str(production.get("work_root") or "").strip()
    if not work_root:
        return []

    root = Path(work_root)
    if not root.exists() or not root.is_dir():
        return []

    include_targets = _dedupe(_as_list(production.get("include_targets")))
    if include_targets:
        missing_targets = [t for t in include_targets if not (root / t).is_dir()]
        if missing_targets:
            print(
                f"[WARN] {construct_type}: configured include_targets not found under {root}: "
                f"{', '.join(missing_targets)}",
            )
        target_dirs = [root / t for t in include_targets if (root / t).is_dir()]
    else:
        target_dirs = sorted([p for p in root.iterdir() if p.is_dir() and p.name != "af3_msa"])

    if not target_dirs:
        return []

    proteins: list[str] = []
    seen: set[str] = set()
    for target_dir in target_dirs:
        target = target_dir.name
        input_dir = target_dir / "af3" / "input"
        if not input_dir.is_dir():
            continue

        pattern = re.compile(rf"^(?P<protein>.+)_{re.escape(target)}_data\\.json$")
        for json_path in sorted(input_dir.glob(f"*_{target}_data.json")):
            m = pattern.match(json_path.name)
            if not m:
                continue
            protein = m.group("protein").strip()
            if not protein or protein in seen:
                continue
            seen.add(protein)
            proteins.append(protein)
    return proteins


entries = [
    ("domain_only", "del_a", core_config_path),
    ("full_length", "del_b", full_config_path),
]

plan = {"constructs": {}, "n_tasks": 0}

for construct_type, del_branch, src_path in entries:
    base = _load(src_path)
    production = base.setdefault("production", {})
    if not isinstance(production, dict):
        raise TypeError(f"production must be a mapping in {src_path}")

    production["construct_type"] = construct_type

    configured_proteins = _dedupe(_as_list(production.get("include_proteins")))
    discovered_proteins = _discover_proteins(production)

    missing: list[str] = []
    if configured_proteins:
        if not discovered_proteins:
            raise RuntimeError(
                f"{construct_type}: include_proteins is configured, but no proteins were discovered "
                f"under this construct work_root. Refusing to submit unfiltered jobs because that "
                f"can run proteins that do not exist for this construct. Check production.work_roots, "
                f"production.work_root, production.include_targets, and the AF3 input layout."
            )

        discovered_set = set(discovered_proteins)
        proteins = [p for p in configured_proteins if p in discovered_set]
        missing = [p for p in configured_proteins if p not in discovered_set]
        protein_source = "include_proteins_filtered_by_auto_discovery"

        if missing:
            print(
                f"[WARN] {construct_type}: skipped {len(missing)} configured proteins "
                f"not found under this construct work_root",
            )
    else:
        proteins = discovered_proteins
        protein_source = "auto_discovery"

    if not proteins:
        print(f"[WARN] {construct_type}: no proteins selected; no shards will be submitted")

    shards = _chunk(proteins, proteins_per_shard)

    construct_cfg_dir = config_out_dir / construct_type
    construct_cfg_dir.mkdir(parents=True, exist_ok=True)

    run_id_base = str(production.get("run_id") or f"analysis_{construct_type}")
    tasks = []

    for idx, shard_proteins in enumerate(shards, start=1):
        shard_id = f"shard_{idx:03d}"
        shard_cfg = copy.deepcopy(base)
        shard_prod = shard_cfg.setdefault("production", {})
        shard_prod["construct_type"] = construct_type
        shard_prod["run_id"] = f"{run_id_base}_{construct_type}_{shard_id}"

        if proteins:
            shard_prod["include_proteins"] = shard_proteins

        shard_cfg_path = construct_cfg_dir / f"{shard_id}.yaml"
        shard_cfg_path.write_text(yaml.safe_dump(shard_cfg, sort_keys=False))

        shard_output_dir = output_dir / construct_type / "shards" / shard_id
        task = {
            "task_index": idx,
            "construct_type": construct_type,
            "del_branch": del_branch,
            "shard_id": shard_id,
            "config_path": str(shard_cfg_path),
            "output_dir": str(shard_output_dir),
            "n_proteins": len(shard_proteins),
            "proteins": shard_proteins,
        }
        tasks.append(task)

    index_path = index_dir / f"{construct_type}_shard_index.tsv"
    with index_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "task_index",
            "shard_id",
            "del_branch",
            "config_path",
            "output_dir",
            "n_proteins",
        ])
        for task in tasks:
            writer.writerow([
                task["task_index"],
                task["shard_id"],
                task["del_branch"],
                task["config_path"],
                task["output_dir"],
                task["n_proteins"],
            ])

    plan["constructs"][construct_type] = {
        "source_config": str(src_path),
        "index_tsv": str(index_path),
        "protein_source": protein_source,
        "n_configured_proteins": len(configured_proteins),
        "n_discovered_proteins": len(discovered_proteins),
        "n_selected_proteins": len(proteins),
        "n_skipped_configured_proteins": len(missing),
        "n_shards": len(tasks),
        "tasks": tasks,
    }
    plan["n_tasks"] += len(tasks)

Path("$SHARD_PLAN_JSON").write_text(json.dumps(plan, indent=2))
print(json.dumps({
    "n_tasks": plan["n_tasks"],
    "domain_only_shards": plan["constructs"].get("domain_only", {}).get("n_shards", 0),
    "full_length_shards": plan["constructs"].get("full_length", {}).get("n_shards", 0),
}, indent=2))
PYEOF

DOMAIN_INDEX_TSV="$INDEX_DIR/domain_only_shard_index.tsv"
FULL_INDEX_TSV="$INDEX_DIR/full_length_shard_index.tsv"

DOMAIN_SHARDS=$("$PYTHON_BIN" - <<PYEOF
import json
from pathlib import Path
plan = json.loads(Path("$SHARD_PLAN_JSON").read_text())
print(plan["constructs"].get("domain_only", {}).get("n_shards", 0))
PYEOF
)
FULL_SHARDS=$("$PYTHON_BIN" - <<PYEOF
import json
from pathlib import Path
plan = json.loads(Path("$SHARD_PLAN_JSON").read_text())
print(plan["constructs"].get("full_length", {}).get("n_shards", 0))
PYEOF
)

if [[ "$DOMAIN_SHARDS" -lt 1 && "$FULL_SHARDS" -lt 1 ]]; then
  echo "[ERROR] No shard tasks were generated." >&2
  exit 1
fi

ARRAY_WORKER_SCRIPT="$JOB_SCRIPT_DIR/run_shard_array_task.sh"
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

echo "[INFO] START task=$TASK_ID shard=$SHARD_ID del=$DEL_BRANCH n_proteins=${N_PROTEINS:-unknown}" | tee -a "$SHARD_LOG"
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

echo "[INFO] DONE task=$TASK_ID shard=$SHARD_ID" | tee -a "$SHARD_LOG"
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
    --job-name="analysis_${construct_type}_shards"
    --account="$ACCOUNT"
    --partition="$PARTITION"
    --time="$SHARD_TIME"
    --mem="$SHARD_MEM"
    --cpus-per-task="$SHARD_CPUS"
    --output="$LOG_DIR/${construct_type}_%A_%a.out"
    --error="$LOG_DIR/${construct_type}_%A_%a.err"
    "$ARRAY_WORKER_SCRIPT"
    "$index_tsv"
    "$REPO_ROOT"
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

DOMAIN_ARRAY_JOB=$(submit_or_print_array "domain_only" "$DOMAIN_SHARDS" "$DOMAIN_INDEX_TSV")
FULL_ARRAY_JOB=$(submit_or_print_array "full_length" "$FULL_SHARDS" "$FULL_INDEX_TSV")

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

per_construct_surfaces = [
    "condition_table.tsv",
    "cluster_table.tsv",
    "protein_condition_residue_scores.tsv",
    "protein_residue_regio_delta.tsv",
]
combined_surfaces = [
    "condition_table.tsv",
    "cluster_table.tsv",
    "protein_condition_residue_scores.tsv",
    "protein_residue_regio_delta.tsv",
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
        for key in row.keys():
            if key in seen:
                continue
            seen.add(key)
            fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


summary: dict[str, Any] = {
    "constructs": {},
    "combined": {},
}

construct_to_tasks: dict[str, list[dict[str, Any]]] = {}
for task in [t for c in plan["constructs"].values() for t in c.get("tasks", [])]:
    construct_to_tasks.setdefault(task["construct_type"], []).append(task)

for construct_type, tasks in sorted(construct_to_tasks.items()):
    merged_root = output_dir / construct_type / "merged"
    construct_summary = {}
    for surface in per_construct_surfaces:
        rows: list[dict[str, Any]] = []
        source_files: list[str] = []
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
for surface in combined_surfaces:
    rows: list[dict[str, Any]] = []
    source_files: list[str] = []
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
  else
    echo "[WARN] combined cluster_table.tsv missing; running cbm-paired without --cluster-table"
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

echo "[INFO] Completed collect job"
echo "[INFO] merge summary: $MERGE_SUMMARY_JSON"
EOF
chmod +x "$COLLECT_SCRIPT"

SUMMARY_DEPENDENCIES=()
[[ -n "$DOMAIN_ARRAY_JOB" ]] && SUMMARY_DEPENDENCIES+=("$DOMAIN_ARRAY_JOB")
[[ -n "$FULL_ARRAY_JOB" ]] && SUMMARY_DEPENDENCIES+=("$FULL_ARRAY_JOB")

SUMMARY_JOB=""
if [[ ${#SUMMARY_DEPENDENCIES[@]} -gt 0 ]]; then
  dependency_arg=$(IFS=:; echo "${SUMMARY_DEPENDENCIES[*]}")
  summary_cmd=(
    sbatch --parsable
    --dependency=afterok:"$dependency_arg"
    --job-name=analysis_full_collect
    --account="$ACCOUNT"
    --partition="$PARTITION"
    --time="$SUMMARY_TIME"
    --mem="$SUMMARY_MEM"
    --cpus-per-task="$SUMMARY_CPUS"
    --output="$LOG_DIR/collect_%j.out"
    --error="$LOG_DIR/collect_%j.err"
    "$COLLECT_SCRIPT"
    "$REPO_ROOT"
    "$RUNTIME_PATHS"
    "$OUTPUT_DIR"
    "$SHARD_PLAN_JSON"
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

PROTEIN_METADATA_ABS=""
if [[ -n "$PROTEIN_METADATA" ]]; then
  PROTEIN_METADATA_ABS="$(realpath -m "$PROTEIN_METADATA")"
fi

CORE_FASTA_ABS=""
if [[ -n "$CORE_FASTA" ]]; then
  CORE_FASTA_ABS="$(realpath -m "$CORE_FASTA")"
fi

FAMILY_ALIGNMENT_DIR_ABS=""
if [[ -n "$FAMILY_ALIGNMENT_DIR" ]]; then
  FAMILY_ALIGNMENT_DIR_ABS="$(realpath -m "$FAMILY_ALIGNMENT_DIR")"
fi

if [[ -n "$SUBMISSION_METADATA_JSON" ]]; then
  mkdir -p "$(dirname "$SUBMISSION_METADATA_JSON")"
  cat > "$SUBMISSION_METADATA_JSON" <<EOF
{
  "output_dir": "$(realpath -m "$OUTPUT_DIR")",
  "runtime_paths": "$(realpath -m "$RUNTIME_PATHS")",
  "core_config": "$(realpath -m "$CORE_CONFIG")",
  "full_config": "$(realpath -m "$FULL_CONFIG")",
  "proteins_per_shard": $PROTEINS_PER_SHARD,
  "max_parallel_shards": $MAX_PARALLEL_SHARDS,
  "n_jobs_per_shard": $N_JOBS_PER_SHARD,
  "shard_cpus": $SHARD_CPUS,
  "shard_mem": "$SHARD_MEM",
  "shard_time": "$SHARD_TIME",
  "summary_cpus": $SUMMARY_CPUS,
  "summary_mem": "$SUMMARY_MEM",
  "summary_time": "$SUMMARY_TIME",
  "run_global_postprocess": $([[ "$RUN_GLOBAL_POSTPROCESS" == "true" ]] && echo true || echo false),
  "predictive_task": "$PREDICTIVE_TASK",
  "random_state": $RANDOM_STATE,
  "protein_metadata": "$PROTEIN_METADATA_ABS",
  "core_fasta": "$CORE_FASTA_ABS",
  "family_alignment_dir": "$FAMILY_ALIGNMENT_DIR_ABS",
  "mafft_executable": "$MAFFT_EXECUTABLE",
  "domain_shards": $DOMAIN_SHARDS,
  "full_shards": $FULL_SHARDS,
  "domain_array_job": "${DOMAIN_ARRAY_JOB:-}",
  "full_array_job": "${FULL_ARRAY_JOB:-}",
  "collect_job": "${SUMMARY_JOB:-}",
  "dry_run": $([[ "$DRY_RUN" == "true" ]] && echo true || echo false)
}
EOF
fi

cat <<EOF
Submitted true Slurm arrays for sharded production run:
  output_dir:            $(realpath -m "$OUTPUT_DIR")
  proteins_per_shard:    $PROTEINS_PER_SHARD
  max_parallel_shards:   $MAX_PARALLEL_SHARDS
  n_jobs_per_shard:      $N_JOBS_PER_SHARD
  shard_cpus:            $SHARD_CPUS
  shard_mem:             $SHARD_MEM
  shard_time:            $SHARD_TIME
  summary_cpus:          $SUMMARY_CPUS
  summary_mem:           $SUMMARY_MEM
  summary_time:          $SUMMARY_TIME
  run_global_postprocess:$RUN_GLOBAL_POSTPROCESS
  domain_shards:         $DOMAIN_SHARDS
  full_shards:           $FULL_SHARDS
  domain_array_job:      ${DOMAIN_ARRAY_JOB:-none}
  full_array_job:        ${FULL_ARRAY_JOB:-none}
  collect_job:           ${SUMMARY_JOB:-none}

Generated plan:
  $SHARD_PLAN_JSON

Note:
  Collect job performs merge and optional predictive/cbm-paired/family-enrichment.
EOF
