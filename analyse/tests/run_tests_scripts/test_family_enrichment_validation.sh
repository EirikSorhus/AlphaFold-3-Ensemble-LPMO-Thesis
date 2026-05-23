#!/bin/bash
#SBATCH --job-name=run_tests_family_enrichment_validation
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_family_enrichment_validation_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "[ERROR] This script must be submitted with sbatch, not run directly on the login node."
    echo "[INFO] Example: sbatch tests/run_tests_scripts/test_family_enrichment_validation.sh"
    exit 1
fi

export PATH="/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin:$PATH"
export PIP_NO_CACHE_DIR=1

project_root="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse"
python_bin="/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin/python"
cd "$project_root"
export PYTHONPATH="$project_root/src:${PYTHONPATH:-}"

out_dir="$project_root/tests/tests_results"
mkdir -p "$out_dir"

job_suffix="${SLURM_JOB_ID}"
run_id="family_enrichment_validation_${job_suffix}"
run_dir="$out_dir/$run_id"
production_output_root="$project_root/tests/tests_results/clustering_pilot_staged/domain_only_shards"
protein_metadata="$project_root/input_data/metadata_final_ec_fixed.tsv"
core_fasta="$project_root/input_data/lpmo_core_domain_2026-03-14_06-52-15_deduplicated.fasta"
alignment_dir=""
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
        --production-output-root)
            production_output_root="$2"
            shift 2
            ;;
        --protein-metadata)
            protein_metadata="$2"
            shift 2
            ;;
        --core-fasta)
            core_fasta="$2"
            shift 2
            ;;
        --alignment-dir)
            alignment_dir="$2"
            shift 2
            ;;
        --skip-unit-tests)
            skip_unit_tests=1
            shift
            ;;
        --help)
            cat <<'EOF'
Usage:
    sbatch tests/run_tests_scripts/test_family_enrichment_validation.sh [options]

Options:
    --run-dir PATH                 Override the output directory.
    --run-id ID                    Override the run identifier written to summary_paths.txt.
    --production-output-root PATH  Production output directory or staged shard root.
    --protein-metadata PATH        Protein metadata TSV used for family resolution.
    --core-fasta PATH              Deduplicated catalytic-core FASTA.
    --alignment-dir PATH           Optional directory with AA9.aligned.fasta / AA10.aligned.fasta.
    --skip-unit-tests              Skip focused family pytest coverage.
    --help                         Show this help text.
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

export XDG_CACHE_HOME="$run_dir/cache"
mkdir -p "$XDG_CACHE_HOME"

if [[ ! -d "$production_output_root" ]]; then
    echo "[ERROR] Production output root does not exist: $production_output_root"
    exit 1
fi

if [[ ! -f "$protein_metadata" ]]; then
    echo "[ERROR] Protein metadata TSV does not exist: $protein_metadata"
    exit 1
fi

if [[ ! -f "$core_fasta" ]]; then
    echo "[ERROR] Core FASTA does not exist: $core_fasta"
    exit 1
fi

if [[ -n "$alignment_dir" ]]; then
    if [[ ! -d "$alignment_dir" ]]; then
        echo "[ERROR] Alignment directory does not exist: $alignment_dir"
        exit 1
    fi
    echo "[INFO] Using precomputed family alignments from: $alignment_dir"
else
    if ! command -v mafft >/dev/null 2>&1; then
        echo "[ERROR] MAFFT is required for this validation when --alignment-dir is not provided."
        exit 1
    fi
    echo "[INFO] Using MAFFT from: $(command -v mafft)"
    echo "[INFO] MAFFT version: $(mafft --version 2>&1 | sed -n '1p')"
fi

if [[ "$skip_unit_tests" -eq 0 ]]; then
    echo "[INFO] Running focused family enrichment pytest suite"
    pytest \
        tests/test_residue_importance.py::test_compute_residue_importance_outputs_counts_only_positive_contact_conditions \
        tests/test_family_alignment.py \
        tests/test_family_residue_enrichment.py \
        tests/test_family_enrichment_postprocess.py \
        tests/test_cli_run.py::test_cmd_family_enrichment_calls_postprocess \
        -q
else
    echo "[INFO] Skipping focused family enrichment pytest suite (--skip-unit-tests)"
fi

python_args=(
    tests/run_tests_scripts/run_family_enrichment_validation.py
    --production-output-root "$production_output_root"
    --protein-metadata "$protein_metadata"
    --core-fasta "$core_fasta"
    --output-dir "$run_dir"
)

if [[ -n "$alignment_dir" ]]; then
    python_args+=(--alignment-dir "$alignment_dir")
fi

echo "[INFO] Running family enrichment validation harness on staged Stage 16b outputs"
"$python_bin" "${python_args[@]}"

validation_summary_path="$run_dir/family_enrichment_validation_summary.json"
summary_paths_file="$run_dir/summary_paths.txt"

if [[ ! -f "$validation_summary_path" ]]; then
    echo "[ERROR] Family enrichment validation summary was not produced: $validation_summary_path"
    exit 1
fi

"$python_bin" - "$validation_summary_path" <<'PY'
from __future__ import annotations

import json
import pathlib
import sys

validation_summary_path = pathlib.Path(sys.argv[1])
validation_summary = json.loads(validation_summary_path.read_text())
if validation_summary.get("validation_status") == "error":
    raise SystemExit(f"Validation harness reported error: {validation_summary.get('error', 'unknown error')}")

processed_families = validation_summary.get("processed_families") or {}
if not processed_families:
    raise SystemExit("Expected at least one processed family in real-data validation, but none were processed")

skipped_families = validation_summary.get("skipped_families") or {}
for family_label, payload in skipped_families.items():
    if payload.get("reason") == "mafft_not_available":
        raise SystemExit(f"Family {family_label} unexpectedly skipped because mafft_not_available")

family_summary_path = pathlib.Path(validation_summary["family_enrichment_summary_path"])
if not family_summary_path.exists():
    raise SystemExit(f"Family enrichment summary JSON was not produced: {family_summary_path}")

family_summary = json.loads(family_summary_path.read_text())
if int(family_summary.get("n_family_aligned_residue_rows", 0)) <= 0:
    raise SystemExit("Expected family_aligned_residue_table.tsv to contain data rows")
if int(family_summary.get("n_family_residue_enrichment_rows", 0)) <= 0:
    raise SystemExit("Expected family_residue_enrichment.tsv to contain data rows")

print("[INFO] Processed families:", ", ".join(sorted(processed_families)))
print("[INFO] Skipped families:", ", ".join(sorted(skipped_families)) if skipped_families else "none")
print("[INFO] Family aligned residue rows:", family_summary.get("n_family_aligned_residue_rows"))
print("[INFO] Family enrichment rows:", family_summary.get("n_family_residue_enrichment_rows"))
PY

{
    echo "run_id=$run_id"
    echo "run_dir=$run_dir"
    echo "production_output_root=$production_output_root"
    echo "validation_summary=$validation_summary_path"
    echo "family_summary=$run_dir/08_family_residue_enrichment/family_enrichment_summary.json"
    echo "family_aligned_residue_table=$run_dir/08_family_residue_enrichment/family_aligned_residue_table.tsv"
    echo "family_residue_enrichment_table=$run_dir/08_family_residue_enrichment/family_residue_enrichment.tsv"
    echo "family_alignment_manifest=$run_dir/08_family_residue_enrichment/family_alignment_manifest.tsv"
} > "$summary_paths_file"

echo "[INFO] Family enrichment validation run completed"
echo "[INFO] Summary paths: $summary_paths_file"
echo "[INFO] Validation summary: $validation_summary_path"
echo "[INFO] Output directory: $run_dir/08_family_residue_enrichment"