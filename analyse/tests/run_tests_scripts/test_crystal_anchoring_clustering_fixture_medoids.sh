#!/bin/bash
#SBATCH --job-name=run_tests_crystal_anchor_fixture_medoids
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_crystal_anchor_fixture_medoids_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "[ERROR] This script must be submitted with sbatch, not run directly on the login node."
    echo "[INFO] Example: sbatch tests/run_tests_scripts/test_crystal_anchoring_clustering_fixture_medoids.sh"
    exit 1
fi

export PATH="/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin:$PATH"
export PIP_NO_CACHE_DIR=1

project_root="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse"
cd "$project_root"
export PYTHONPATH="$project_root/src:${PYTHONPATH:-}"

out_dir="$project_root/tests/tests_results"
mkdir -p "$out_dir"

job_suffix="${SLURM_JOB_ID}"
run_id="crystal_anchoring_clustering_fixture_medoids_${job_suffix}"
run_dir="$out_dir/${run_id}"
fixture_root="$project_root/tests/fixtures/clustering_stage_outputs"
include_proteins=""
include_conditions=""
max_medoids=0
skip_unit_tests=0
allow_empty=0
fail_on_screen_error=0

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
        --fixture-root)
            fixture_root="$2"
            shift 2
            ;;
        --include-proteins)
            include_proteins="$2"
            shift 2
            ;;
        --include-conditions)
            include_conditions="$2"
            shift 2
            ;;
        --max-medoids)
            max_medoids="$2"
            shift 2
            ;;
        --skip-unit-tests)
            skip_unit_tests=1
            shift
            ;;
        --allow-empty)
            allow_empty=1
            shift
            ;;
        --fail-on-screen-error)
            fail_on_screen_error=1
            shift
            ;;
        --help)
            cat <<'EOF'
Usage:
    sbatch tests/run_tests_scripts/test_crystal_anchoring_clustering_fixture_medoids.sh [options]

Options:
    --run-dir PATH             Override output directory under tests/tests_results.
    --run-id ID                Override run identifier written to the summary.
    --fixture-root PATH        Override clustering-stage fixture root.
    --include-proteins CSV     Comma-separated protein IDs to include.
    --include-conditions CSV   Comma-separated condition IDs to include.
    --max-medoids N            Optional cap after filtering; 0 means all crystal-covered fixture medoids.
    --skip-unit-tests          Skip focused pytest checks before the harness.
    --allow-empty              Do not fail if no comparisons are produced.
    --fail-on-screen-error     Fail if any medoid screen raises an exception.
    --help                     Show this help text.

Default behavior screens all medoids in tests/fixtures/clustering_stage_outputs
whose protein has at least one crystal reference in input_data/pdb_structure_data.csv.
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

if [[ "$skip_unit_tests" -eq 0 ]]; then
    echo "[INFO] Running focused crystal anchoring unit tests"
    pytest tests/test_crystal_anchoring.py -q
else
    echo "[INFO] Skipping focused crystal anchoring unit tests (--skip-unit-tests)"
fi

echo "[INFO] Running crystal anchoring for clustering fixture medoids"
python_args=(
    tests/run_tests_scripts/run_crystal_anchoring_clustering_fixture_medoids.py
    --fixture-root "$fixture_root"
    --run-dir "$run_dir"
    --run-id "$run_id"
    --max-medoids "$max_medoids"
)

if [[ -n "$include_proteins" ]]; then
    python_args+=(--include-proteins "$include_proteins")
fi

if [[ -n "$include_conditions" ]]; then
    python_args+=(--include-conditions "$include_conditions")
fi

if [[ "$allow_empty" -eq 1 ]]; then
    python_args+=(--allow-empty)
fi

if [[ "$fail_on_screen_error" -eq 1 ]]; then
    python_args+=(--fail-on-screen-error)
fi

/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin/python "${python_args[@]}"

echo "[INFO] Crystal anchoring fixture-medoid run completed"
echo "[INFO] Results directory: $run_dir"
echo "[INFO] Summary: $run_dir/fixture_crystal_anchoring_summary.json"
echo "[INFO] Anchor table: $run_dir/fixture_crystal_anchor_table.tsv"
echo "[INFO] Geometry table: $run_dir/fixture_crystal_geometry_table.tsv"
echo "[INFO] Diagnostic table: $run_dir/fixture_crystal_ifp_diagnostic_summary.tsv"
