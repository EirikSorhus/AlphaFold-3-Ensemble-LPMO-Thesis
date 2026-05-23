#!/bin/bash
#SBATCH --job-name=run_tests_cbm_paired_real_validation
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_cbm_paired_real_validation_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "[ERROR] This script must be submitted with sbatch, not run directly on the login node."
    echo "[INFO] Example: sbatch tests/run_tests_scripts/test_cbm_paired_real_validation.sh"
    exit 1
fi

export PATH="/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin:$PATH"
export PIP_NO_CACHE_DIR=1

project_root="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse"
python_bin="/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin/python"
cd "$project_root"
export PYTHONPATH="$project_root/src:${PYTHONPATH:-}"

job_suffix="${SLURM_JOB_ID}"
run_id="cbm_paired_real_validation_${job_suffix}"
run_dir="$project_root/tests/tests_results/$run_id"
domain_only_root="$project_root/tests/tests_results/clustering_pilot_staged/domain_only_shards"
full_length_root="$project_root/tests/tests_results/clustering_pilot_staged/full_length_shards"
protein_metadata="$project_root/input_data/metadata_final_ec_fixed.tsv"
skip_unit_tests=0

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --run-dir)
            run_dir="$2"
            shift 2
            ;;
        --run-id)
            run_id="$2"
            shift 2
            ;;
        --domain-only-root)
            domain_only_root="$2"
            shift 2
            ;;
        --full-length-root)
            full_length_root="$2"
            shift 2
            ;;
        --protein-metadata)
            protein_metadata="$2"
            shift 2
            ;;
        --skip-unit-tests)
            skip_unit_tests=1
            shift
            ;;
        --help)
            cat <<'EOF'
Usage:
    sbatch tests/run_tests_scripts/test_cbm_paired_real_validation.sh [options]

Options:
    --run-dir PATH            Override the output directory.
    --run-id ID               Override the run identifier written to summary_paths.txt.
    --domain-only-root PATH   Staged domain-only shard root or one production output.
    --full-length-root PATH   Staged full-length shard root or one production output.
    --protein-metadata PATH   Protein metadata TSV used for family/CBM labels.
    --skip-unit-tests         Skip focused pytest coverage before the real-row harness.
    --help                    Show this help text.
EOF
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1"
            exit 1
            ;;
    esac
done

mkdir -p "$run_dir"

if [[ ! -d "$domain_only_root" ]]; then
    echo "[ERROR] Domain-only staged root does not exist: $domain_only_root"
    exit 1
fi

if [[ ! -d "$full_length_root" ]]; then
    echo "[ERROR] Full-length staged root does not exist: $full_length_root"
    exit 1
fi

if [[ ! -f "$protein_metadata" ]]; then
    echo "[ERROR] Protein metadata TSV does not exist: $protein_metadata"
    exit 1
fi

if [[ "$skip_unit_tests" -eq 0 ]]; then
    echo "[INFO] Running focused CBM pytest coverage"
    pytest tests/test_cbm_comparison.py -q
else
    echo "[INFO] Skipping focused CBM pytest coverage (--skip-unit-tests)"
fi

echo "[INFO] Running staged real-row CBM paired validation harness"
"$python_bin" tests/run_tests_scripts/run_cbm_paired_real_validation.py \
    --domain-only-root "$domain_only_root" \
    --full-length-root "$full_length_root" \
    --protein-metadata "$protein_metadata" \
    --output-dir "$run_dir"

validation_summary_path="$run_dir/cbm_paired_real_validation_summary.json"
summary_paths_file="$run_dir/summary_paths.txt"

if [[ ! -f "$validation_summary_path" ]]; then
    echo "[ERROR] CBM paired validation summary was not produced: $validation_summary_path"
    exit 1
fi

"$python_bin" - "$validation_summary_path" <<'PY'
from __future__ import annotations

import json
import pathlib
import sys

summary_path = pathlib.Path(sys.argv[1])
summary = json.loads(summary_path.read_text())
if not summary.get("validation_passed"):
    raise SystemExit(f"CBM paired validation failed: {summary.get('errors') or summary.get('error')}")

if int(summary.get("overlap_pair_key_count", 0)) <= 0:
    raise SystemExit("Expected at least one overlapping protein/substrate/dp pair between domain-only and full-length real rows")
if int(summary.get("n_paired_rows", 0)) <= 0:
    raise SystemExit("Expected cbm_paired_comparison_table.tsv to contain real paired rows")
if int(summary.get("n_nonempty_family_labels", 0)) != int(summary.get("n_paired_rows", 0)):
    raise SystemExit("Expected every paired CBM row to carry a non-empty family label")
if int(summary.get("n_nonempty_cbm_types", 0)) != int(summary.get("n_paired_rows", 0)):
    raise SystemExit("Expected every paired CBM row to carry a non-empty cbm_type label")

print("[INFO] Overlap pair keys:", summary.get("overlap_pair_key_count"))
print("[INFO] Paired rows:", summary.get("n_paired_rows"))
print("[INFO] CBM summary:", summary.get("cbm_paired_analysis_summary"))
PY

{
    echo "run_id=$run_id"
    echo "run_dir=$run_dir"
    echo "validation_summary=$validation_summary_path"
    echo "cbm_summary=$("$python_bin" - "$validation_summary_path" <<'PY'
from __future__ import annotations
import json
import pathlib
import sys
summary = json.loads(pathlib.Path(sys.argv[1]).read_text())
print(summary.get('cbm_paired_analysis_summary', ''))
PY
)"
} > "$summary_paths_file"

echo "[INFO] CBM paired real-row validation completed"
echo "[INFO] Summary paths: $summary_paths_file"
echo "[INFO] Validation summary: $validation_summary_path"