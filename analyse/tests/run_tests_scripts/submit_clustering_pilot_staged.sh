#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
PYTHON_BIN="/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin/python"
OVERVIEW_TSV="$REPO_ROOT/input_data/clustering_pilot_protein_overview.tsv"
RUN_ROOT="$REPO_ROOT/tests/tests_results/clustering_pilot_staged"
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
RUN_GLOBAL_POSTPROCESS=true
PROTEIN_METADATA="$REPO_ROOT/input_data/metadata_final_ec_fixed.tsv"
CORE_FASTA="$REPO_ROOT/input_data/lpmo_core_domain_2026-03-14_06-52-15_deduplicated.fasta"
FAMILY_ALIGNMENT_DIR=""
MAFFT_EXECUTABLE="mafft"
SUBMISSION_METADATA_JSON=""
DRY_RUN=false

usage() {
    cat <<'EOF'
Usage:
    tests/run_tests_scripts/submit_clustering_pilot_staged.sh [options]

Options:
    --run-dir PATH          Output root for staged pilot artifacts.
    --account NAME          Slurm account (default: nn1003k).
    --partition NAME        Slurm partition (default: small).
    --shard-size N          Compatibility alias for --proteins-per-shard.
    --proteins-per-shard N  Protein selections per Slurm array task (default: 5).
    --balance-by MODE       Shard balancing: protein_count|estimated_pose_count (default: estimated_pose_count).
    --array-limit N         Optional Slurm array concurrency cap.
    --execute-cpus N        CPUs per shard task, passed to --n-jobs (default: 8).
    --execute-mem VALUE     Memory per shard task (default: 24G).
    --execute-time VALUE    Walltime per shard task (default: 24:00:00).
    --summary-cpus N        CPUs for aggregate summary job (default: 1).
    --summary-mem VALUE     Memory for aggregate summary job (default: 4G).
    --summary-time VALUE    Walltime for aggregate summary job (default: 00:30:00).
    --run-global-postprocess true|false
                            Run predictive/CBM/family postprocess after merge (default: true).
    --protein-metadata PATH Protein metadata TSV for global postprocess.
    --core-fasta PATH       Core-domain FASTA for global family enrichment.
    --family-alignment-dir PATH
                            Optional precomputed AA9/AA10 alignments.
    --mafft-executable NAME MAFFT executable name/path (default: mafft).
    --submission-metadata-json PATH
                            Optional JSON file describing submitted job ids.
    --dry-run               Write manifests/job scripts and print sbatch commands only.

This is the multi-node Slurm wrapper for the clustering pilot. It builds the
domain-only and full-length selection manifests, splits them into weighted
protein shards, submits one Slurm array per construct type, then collects and
merges production outputs after the arrays finish.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-dir)
            RUN_ROOT=$(realpath -m "$2")
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
        --shard-size)
            SHARD_SIZE="$2"
            shift 2
            ;;
        --proteins-per-shard)
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
        --execute-cpus|--prepare-cpus)
            EXECUTE_CPUS="$2"
            shift 2
            ;;
        --execute-mem|--prepare-mem)
            EXECUTE_MEM="$2"
            shift 2
            ;;
        --execute-time|--prepare-time)
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
        --run-global-postprocess)
            RUN_GLOBAL_POSTPROCESS="$2"
            shift 2
            ;;
        --protein-metadata)
            PROTEIN_METADATA=$(realpath -m "$2")
            shift 2
            ;;
        --core-fasta)
            CORE_FASTA=$(realpath -m "$2")
            shift 2
            ;;
        --family-alignment-dir)
            FAMILY_ALIGNMENT_DIR=$(realpath -m "$2")
            shift 2
            ;;
        --mafft-executable)
            MAFFT_EXECUTABLE="$2"
            shift 2
            ;;
        --submission-metadata-json)
            SUBMISSION_METADATA_JSON=$(realpath -m "$2")
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --help)
            usage
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ "$SHARD_SIZE" -lt 1 ]]; then
    echo "[ERROR] --shard-size must be >= 1" >&2
    exit 2
fi

if [[ "$BALANCE_BY" != "protein_count" && "$BALANCE_BY" != "estimated_pose_count" ]]; then
    echo "[ERROR] --balance-by must be protein_count or estimated_pose_count" >&2
    exit 2
fi

LOG_DIR="$RUN_ROOT/logs"
MANIFEST_DIR="$RUN_ROOT/manifests"
SHARD_DIR="$MANIFEST_DIR/shards"
JOB_SCRIPT_DIR="$RUN_ROOT/slurm_scripts"
SUMMARY_DIR="$RUN_ROOT/summaries"
MERGE_OUTPUT_ROOT="$RUN_ROOT/merged_production_output"
mkdir -p "$LOG_DIR" "$MANIFEST_DIR" "$SHARD_DIR" "$JOB_SCRIPT_DIR" "$SUMMARY_DIR"

cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT/src:${PYTHONPATH:-}"

[[ -x "$PYTHON_BIN" ]] || { echo "Missing python interpreter: $PYTHON_BIN" >&2; exit 2; }
[[ -f "$OVERVIEW_TSV" ]] || { echo "Missing pilot overview TSV: $OVERVIEW_TSV" >&2; exit 2; }

"$PYTHON_BIN" "$SCRIPT_DIR/build_clustering_pilot_selection_manifests.py" \
  --overview-tsv "$OVERVIEW_TSV" \
  --output-dir "$MANIFEST_DIR" > "$LOG_DIR/00_build_manifests_local.log" 2>&1
cp "$OVERVIEW_TSV" "$MANIFEST_DIR/clustering_pilot_protein_overview.tsv"
cp "$MANIFEST_DIR/clustering_pilot_manifest_summary.json" "$SUMMARY_DIR/01_manifest_generation.json"

DOMAIN_MANIFEST="$MANIFEST_DIR/clustering_pilot_domain_only_selection.yaml"
FULL_MANIFEST="$MANIFEST_DIR/clustering_pilot_full_length_selection.yaml"
DOMAIN_SHARD_DIR="$SHARD_DIR/domain_only"
FULL_SHARD_DIR="$SHARD_DIR/full_length"

"$PYTHON_BIN" "$SCRIPT_DIR/split_clustering_pilot_selection_manifest.py" \
  --input-manifest "$DOMAIN_MANIFEST" \
  --output-dir "$DOMAIN_SHARD_DIR" \
    --proteins-per-shard "$SHARD_SIZE" \
    --balance-by "$BALANCE_BY" \
    --output-root-template "$RUN_ROOT/{construct_type}_shards/{shard_id}/production_output" > "$LOG_DIR/00_split_domain_local.log" 2>&1

"$PYTHON_BIN" "$SCRIPT_DIR/split_clustering_pilot_selection_manifest.py" \
  --input-manifest "$FULL_MANIFEST" \
  --output-dir "$FULL_SHARD_DIR" \
    --proteins-per-shard "$SHARD_SIZE" \
    --balance-by "$BALANCE_BY" \
    --output-root-template "$RUN_ROOT/{construct_type}_shards/{shard_id}/production_output" > "$LOG_DIR/00_split_full_local.log" 2>&1

domain_shards=$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1]))["n_shards"])' "$DOMAIN_SHARD_DIR/shard_summary.json")
full_shards=$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1]))["n_shards"])' "$FULL_SHARD_DIR/shard_summary.json")

write_job_script() {
    local path=$1
    local body=$2
    printf '%s\n' "$body" > "$path"
    chmod +x "$path"
}

ARRAY_WORKER_SCRIPT="$JOB_SCRIPT_DIR/run_shard_array_task.sh"
write_job_script "$ARRAY_WORKER_SCRIPT" "#!/bin/bash
set -euo pipefail

INDEX_TSV=\"\$1\"
RUN_ROOT=\"\$2\"
CONSTRUCT_TYPE=\"\$3\"
PYTHON_BIN=\"$PYTHON_BIN\"
SCRIPT_DIR=\"$SCRIPT_DIR\"
REPO_ROOT=\"$REPO_ROOT\"
TASK_ID=\"\${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}\"

cd \"\$REPO_ROOT\"
export PYTHONPATH=\"\$REPO_ROOT/src:\${PYTHONPATH:-}\"
export PATH=\"/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin:\$PATH\"

read -r SHARD_ID MANIFEST_PATH < <(awk -F '\\t' -v task=\"\$TASK_ID\" 'NR > 1 && \$1 == task {print \$2, \$5}' \"\$INDEX_TSV\")
if [[ -z \"\${SHARD_ID:-}\" || -z \"\${MANIFEST_PATH:-}\" ]]; then
    echo \"[ERROR] No shard mapping for array task \$TASK_ID in \$INDEX_TSV\" >&2
    exit 2
fi

SHARD_RUN_DIR=\"\$RUN_ROOT/\${CONSTRUCT_TYPE}_shards/\$SHARD_ID\"
\"\$PYTHON_BIN\" \"\$SCRIPT_DIR/run_clustering_pilot_real_case.py\" \\
  --run-dir \"\$SHARD_RUN_DIR\" \\
  --run-id \"\${CONSTRUCT_TYPE}_\$SHARD_ID\" \\
  --selection-manifest \"\$MANIFEST_PATH\" \\
  --execute \\
  --n-jobs \"$EXECUTE_CPUS\"
"

COLLECT_SCRIPT="$JOB_SCRIPT_DIR/collect_shard_summaries.sh"
write_job_script "$COLLECT_SCRIPT" "#!/bin/bash
set -euo pipefail

cd \"$REPO_ROOT\"
export PYTHONPATH=\"$REPO_ROOT/src:\${PYTHONPATH:-}\"
collect_args=(
    --run-root \"$RUN_ROOT\"
    --output-json \"$SUMMARY_DIR/staged_array_summary.json\"
    --output-tsv \"$SUMMARY_DIR/staged_array_shards.tsv\"
    --merge-output-root \"$MERGE_OUTPUT_ROOT\"
    --protein-metadata \"$PROTEIN_METADATA\"
    --core-fasta \"$CORE_FASTA\"
    --mafft-executable \"$MAFFT_EXECUTABLE\"
)

if [[ -n \"$FAMILY_ALIGNMENT_DIR\" ]]; then
    collect_args+=(--family-alignment-dir \"$FAMILY_ALIGNMENT_DIR\")
fi
if [[ \"$RUN_GLOBAL_POSTPROCESS\" == \"true\" ]]; then
    collect_args+=(--run-global-postprocess)
fi

\"$PYTHON_BIN\" \"$SCRIPT_DIR/collect_clustering_pilot_shards.py\" \"\${collect_args[@]}\"
"

array_spec() {
    local count=$1
    if [[ -n "$ARRAY_LIMIT" ]]; then
        printf '1-%s%%%s' "$count" "$ARRAY_LIMIT"
    else
        printf '1-%s' "$count"
    fi
}

submit_or_print_array() {
    local construct_type=$1
    local count=$2
    local index_tsv=$3
    local output_pattern=$4

    if [[ "$count" -lt 1 ]]; then
        echo ""
        return 0
    fi

    local cmd=(
        sbatch --parsable
        --array="$(array_spec "$count")"
        --job-name="clust_${construct_type}_shards"
        --account="$ACCOUNT"
        --partition="$PARTITION"
        --time="$EXECUTE_TIME"
        --mem="$EXECUTE_MEM"
        --cpus-per-task="$EXECUTE_CPUS"
        --output="$output_pattern"
        "$ARRAY_WORKER_SCRIPT"
        "$index_tsv"
        "$RUN_ROOT"
        "$construct_type"
    )

    if [[ "$DRY_RUN" == true ]]; then
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

DOMAIN_ARRAY_JOB=$(submit_or_print_array \
    "domain_only" \
    "$domain_shards" \
    "$DOMAIN_SHARD_DIR/shard_index.tsv" \
    "$LOG_DIR/domain_only_shard_%A_%a.log")

FULL_ARRAY_JOB=$(submit_or_print_array \
    "full_length" \
    "$full_shards" \
    "$FULL_SHARD_DIR/shard_index.tsv" \
    "$LOG_DIR/full_length_shard_%A_%a.log")

SUMMARY_DEPENDENCIES=()
[[ -n "$DOMAIN_ARRAY_JOB" ]] && SUMMARY_DEPENDENCIES+=("$DOMAIN_ARRAY_JOB")
[[ -n "$FULL_ARRAY_JOB" ]] && SUMMARY_DEPENDENCIES+=("$FULL_ARRAY_JOB")

if [[ ${#SUMMARY_DEPENDENCIES[@]} -eq 0 ]]; then
    echo "[ERROR] No domain-only or full-length shards were generated." >&2
    exit 1
fi

dependency_arg=$(IFS=:; echo "${SUMMARY_DEPENDENCIES[*]}")
summary_cmd=(
    sbatch --parsable
    --dependency=afterany:"$dependency_arg"
    --job-name=clust_pilot_collect
    --account="$ACCOUNT"
    --partition="$PARTITION"
    --time="$SUMMARY_TIME"
    --mem="$SUMMARY_MEM"
    --cpus-per-task="$SUMMARY_CPUS"
    --output="$LOG_DIR/collect_%j.log"
    "$COLLECT_SCRIPT"
)

if [[ "$DRY_RUN" == true ]]; then
    {
        printf '[DRY-RUN]'
        printf ' %q' "${summary_cmd[@]}"
        printf '\n'
    } >&2
    SUMMARY_JOB="dry_run"
else
    SUMMARY_JOB=$("${summary_cmd[@]}")
fi

if [[ -n "$SUBMISSION_METADATA_JSON" ]]; then
        mkdir -p "$(dirname "$SUBMISSION_METADATA_JSON")"
        cat > "$SUBMISSION_METADATA_JSON" <<EOF
{
    "run_root": "$RUN_ROOT",
    "shard_size": $SHARD_SIZE,
    "balance_by": "$BALANCE_BY",
    "domain_shards": $domain_shards,
    "full_shards": $full_shards,
    "execute_cpus": $EXECUTE_CPUS,
    "run_global_postprocess": $([[ "$RUN_GLOBAL_POSTPROCESS" == true ]] && echo true || echo false),
    "domain_array_job": "${DOMAIN_ARRAY_JOB:-}",
    "full_array_job": "${FULL_ARRAY_JOB:-}",
    "summary_job": "$SUMMARY_JOB",
    "dry_run": $([[ "$DRY_RUN" == true ]] && echo true || echo false)
}
EOF
fi

cat <<EOF
Submitted staged clustering pilot arrays:
  run_root:          $RUN_ROOT
  shard_size:        $SHARD_SIZE
    balance_by:        $BALANCE_BY
  domain_shards:     $domain_shards
  full_shards:       $full_shards
  execute_cpus:      $EXECUTE_CPUS
  domain_array_job:  ${DOMAIN_ARRAY_JOB:-none}
  full_array_job:    ${FULL_ARRAY_JOB:-none}
  summary_job:       $SUMMARY_JOB

Production aggregate summary will be written to:
  $SUMMARY_DIR/staged_array_summary.json
  $SUMMARY_DIR/staged_array_shards.tsv
    $MERGE_OUTPUT_ROOT
EOF
