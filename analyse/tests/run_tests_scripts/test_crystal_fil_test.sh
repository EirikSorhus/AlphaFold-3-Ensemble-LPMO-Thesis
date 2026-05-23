#!/bin/bash
#SBATCH --job-name=run_tests_crystal_fil_test
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_crystal_fil_test_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "[ERROR] This script must be submitted with sbatch, not run directly on the login node."
    echo "[INFO] Example: sbatch tests/run_tests_scripts/test_crystal_fil_test.sh"
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
run_id="crystal_fil_test_${job_suffix}"
run_dir="$out_dir/${run_id}"
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
        --skip-unit-tests)
            skip_unit_tests=1
            shift
            ;;
        --help)
            cat <<'EOF'
Usage:
    sbatch tests/run_tests_scripts/test_crystal_fil_test.sh [options]

Options:
    --run-dir PATH       Override output directory.
    --run-id ID          Override run identifier written to the summary.
    --skip-unit-tests    Skip the focused crystal-anchoring pytest slice.
    --help               Show this help text.
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
    echo "[INFO] Running focused crystal-anchoring unit tests"
    pytest tests/test_crystal_anchoring.py -q -k 'filter_crystal_reference_records_matches_same_ligand_and_dp'
else
    echo "[INFO] Skipping focused crystal-anchoring unit tests (--skip-unit-tests)"
fi

echo "[INFO] Running crystal_fil_test harness"
/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin/python \
    tests/run_tests_scripts/run_crystal_fil_test.py \
    --run-dir "$run_dir" \
    --run-id "$run_id"

echo "[INFO] crystal_fil_test run completed"
echo "[INFO] Results directory: $run_dir"
echo "[INFO] Summary: $run_dir/crystal_fil_test_summary.json"
echo "[INFO] Matched pairs: $run_dir/matched_pairs.tsv"