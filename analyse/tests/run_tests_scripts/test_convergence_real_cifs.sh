#!/bin/bash
#SBATCH --job-name=run_tests_convergence_real_cifs
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=00:45:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_convergence_real_cifs_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "[ERROR] This script must be submitted with sbatch, not run directly on the login node."
    echo "[INFO] Example: sbatch tests/run_tests_scripts/test_convergence_real_cifs.sh"
    exit 1
fi

export PATH="/cluster/work/projects/nn1003k/eirik/conda/analyse_env/bin:$PATH"
export PIP_NO_CACHE_DIR=1

project_root="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse"
cd "$project_root"
export PYTHONPATH="$project_root/src:${PYTHONPATH:-}"

out_dir="$project_root/tests/tests_results"
mkdir -p "$out_dir"

job_suffix="${SLURM_JOB_ID}"
run_id="convergence_real_cifs_${job_suffix}"
run_dir="$out_dir/${run_id}"
mkdir -p "$run_dir"

export XDG_CACHE_HOME="$run_dir/cache"
mkdir -p "$XDG_CACHE_HOME"

declare -a cli_cif_paths=()

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
        --help)
            cat <<'EOF'
Usage:
    sbatch tests/run_tests_scripts/test_convergence_real_cifs.sh [options] [cif ...]

Options:
    --run-dir PATH       Override output directory.
    --run-id ID          Override run identifier written to the summary.
    --help               Show this help text.

Arguments:
    cif                  Optional explicit CIF paths. Defaults to the built-in multi-pose real-CIF set.
EOF
            exit 0
            ;;
        *)
            cli_cif_paths+=("$1")
            shift
            ;;
    esac
done

cif_paths=(
  "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/NAG4/af3/latest/Q7SCE9_NAG4/seed-1_sample-0/Q7SCE9_NAG4_seed-1_sample-0_model.cif"
  "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/NAG4/af3/latest/Q7SCE9_NAG4/seed-1_sample-1/Q7SCE9_NAG4_seed-1_sample-1_model.cif"
  "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/NAG4/af3/latest/Q7SCE9_NAG4/seed-4_sample-1/Q7SCE9_NAG4_seed-4_sample-1_model.cif"
  "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/STA6/af3/latest/Q59930_STA6/seed-1_sample-0/Q59930_STA6_seed-1_sample-0_model.cif"
  "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/STA6/af3/latest/Q59930_STA6/seed-1_sample-1/Q59930_STA6_seed-1_sample-1_model.cif"
  "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/STA6/af3/latest/Q59930_STA6/seed-2_sample-0/Q59930_STA6_seed-2_sample-0_model.cif"
)

if [[ "${#cli_cif_paths[@]}" -gt 0 ]]; then
    cif_paths=("${cli_cif_paths[@]}")
fi

for cif_path in "${cif_paths[@]}"; do
    if [[ ! -f "$cif_path" ]]; then
        echo "[ERROR] CIF does not exist: $cif_path"
        exit 1
    fi
done

echo "[INFO] Running focused convergence unit tests"
pytest tests/test_convergence_metrics.py -q

echo "[INFO] Running convergence real-case harness on ${#cif_paths[@]} CIFs"
/cluster/work/projects/nn1003k/eirik/conda/analyse_env/bin/python \
    tests/run_tests_scripts/run_convergence_real_cifs.py \
    --run-dir "$run_dir" \
    --run-id "$run_id" \
    "${cif_paths[@]}"

echo "[INFO] Convergence real-case run completed"
echo "[INFO] Results directory: $run_dir"
echo "[INFO] Summary: $run_dir/convergence_real_cifs_summary.json"
echo "[INFO] Pose convergence output: $run_dir/pose_convergence.tsv"
echo "[INFO] Condition summary output: $run_dir/condition_convergence_summary.tsv"