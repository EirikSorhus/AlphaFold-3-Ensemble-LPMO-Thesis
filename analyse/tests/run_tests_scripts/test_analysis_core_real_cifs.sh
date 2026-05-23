#!/bin/bash
#SBATCH --job-name=run_tests_analysis_core_real_cifs
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=00:35:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=3
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_analysis_core_real_cifs_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "[ERROR] This script must be submitted with sbatch, not run directly on the login node."
    echo "[INFO] Example: sbatch tests/run_tests_scripts/test_analysis_core_real_cifs.sh"
    exit 1
fi

export PATH="/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin:$PATH"
export PIP_NO_CACHE_DIR=1

project_root="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse"
cd "$project_root"
export PYTHONPATH="$project_root/src:${PYTHONPATH:-}"

out_dir="$project_root/tests/tests_results"
mkdir -p "$out_dir"

job_suffix="${SLURM_JOB_ID:-manual}"
run_dir="$out_dir/analysis_core_real_cifs_${job_suffix}"
mkdir -p "$run_dir"

export XDG_CACHE_HOME="$run_dir/cache"
mkdir -p "$XDG_CACHE_HOME"

cif_paths=(
    "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core/NAG4/af3/runs/408002/Q7SCE9_NAG4/seed-4_sample-1/Q7SCE9_NAG4_seed-4_sample-1_model.cif"
    "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core/STA6/af3/runs/408006/Q59930_STA6/seed-9_sample-0/Q59930_STA6_seed-9_sample-0_model.cif"
    "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core/STA4/af3/runs/408005/A0A0S2GKZ1_STA4/seed-9_sample-0/A0A0S2GKZ1_STA4_seed-9_sample-0_model.cif"
)

if [[ "$#" -gt 0 ]]; then
    cif_paths=("$@")
fi

for cif_path in "${cif_paths[@]}"; do
    if [[ ! -f "$cif_path" ]]; then
        echo "[ERROR] CIF does not exist: $cif_path"
        exit 1
    fi
done

echo "[INFO] Running focused production analysis-core tests and Stage 6 validation checks"
pytest tests/test_cluster_signatures.py tests/test_analysis_core_real_cifs_validation.py tests/test_analysis_orchestrator.py tests/test_cli_run.py -q

echo "[INFO] Running production analysis-core path on ${#cif_paths[@]} real CIFs"
/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin/python \
    tests/run_tests_scripts/run_analysis_core_real_cifs.py \
    --run-dir "$run_dir" \
    --run-id "analysis_core_real_cifs_${job_suffix}" \
    --del-branch del_a \
    "${cif_paths[@]}"

echo "[INFO] Analysis-core + clustering run completed"
echo "[INFO] Results directory: $run_dir"
echo "[INFO] Summary: $run_dir/analysis_core_real_cifs_summary.json"
echo "[INFO] Production output: $run_dir/production_output"
