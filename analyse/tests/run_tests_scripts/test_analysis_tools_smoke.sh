#!/bin/bash
#SBATCH --job-name=run_tests_analysis_tools_smoke
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_analysis_tools_smoke_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "[ERROR] This script must be submitted with sbatch, not run directly on the login node."
    echo "[INFO] Example: sbatch tests/run_tests_scripts/test_analysis_tools_smoke.sh"
    exit 1
fi

export PATH="/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin:$PATH"

project_root="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse"
cd "$project_root"
export PYTHONPATH="$project_root/src:${PYTHONPATH:-}"

out_dir="$project_root/tests/tests_results"
mkdir -p "$out_dir"

job_suffix="${SLURM_JOB_ID}"
run_id="analysis_tools_smoke_${job_suffix}"
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
    sbatch tests/run_tests_scripts/test_analysis_tools_smoke.sh [options]

Options:
    --run-dir PATH       Override output directory.
    --run-id ID          Override run identifier written to the summary.
    --skip-unit-tests    Skip pytest tests/test_manifest.py.
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

if [[ "$skip_unit_tests" -eq 0 ]]; then
    echo "[INFO] Running focused manifest/tool probe unit tests"
    pytest tests/test_manifest.py -q
else
    echo "[INFO] Skipping focused manifest/tool probe unit tests (--skip-unit-tests)"
fi

echo "[INFO] Running analysis tool smoke harness"
/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin/python \
    tests/run_tests_scripts/run_analysis_tools_smoke.py \
    --run-dir "$run_dir" \
    --run-id "$run_id"

echo "[INFO] Analysis tool smoke run completed"
echo "[INFO] Results directory: $run_dir"
echo "[INFO] Summary: $run_dir/analysis_tools_smoke_summary.json"