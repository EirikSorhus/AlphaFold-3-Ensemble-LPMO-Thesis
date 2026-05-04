#!/bin/bash
#SBATCH --job-name=run_tests_prolif_real_cifs
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_prolif_real_cifs_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "[ERROR] This script must be submitted with sbatch, not run directly on the login node."
    echo "[INFO] Example: sbatch tests/run_tests_scripts/test_prolif_real_cifs.sh"
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
run_id="prolif_real_cifs_${job_suffix}"
run_dir="$out_dir/${run_id}"

skip_unit_tests=0
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
                --skip-unit-tests)
                        skip_unit_tests=1
                        shift
                        ;;
                --help)
                        cat <<'EOF'
Usage:
    sbatch tests/run_tests_scripts/test_prolif_real_cifs.sh [options] [cif ...]

Options:
    --run-dir PATH       Override output directory.
    --run-id ID          Override run identifier written to the summary.
    --skip-unit-tests    Skip pytest tests/test_prolif_ifp.py.
    --help               Show this help text.

Arguments:
    cif                  Optional explicit CIF paths. Defaults to the standard three real cases.
EOF
                        exit 0
                        ;;
                *)
                        cli_cif_paths+=("$1")
                        shift
                        ;;
        esac
done

mkdir -p "$run_dir"

export XDG_CACHE_HOME="$run_dir/cache"
mkdir -p "$XDG_CACHE_HOME"

cif_paths=(
  "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/NAG4/af3/latest/Q7SCE9_NAG4/seed-4_sample-1/Q7SCE9_NAG4_seed-4_sample-1_model.cif"
  "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/STA6/af3/latest/Q59930_STA6/seed-2_sample-0/Q59930_STA6_seed-2_sample-0_model.cif"
  "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/STA4/af3/latest/A0A0S2GKZ1_STA4/seed-2_sample-2/A0A0S2GKZ1_STA4_seed-2_sample-2_model.cif"
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

if [[ "$skip_unit_tests" -eq 0 ]]; then
    echo "[INFO] Running focused ProLIF unit tests"
    pytest tests/test_prolif_ifp.py -q
else
    echo "[INFO] Skipping focused ProLIF unit tests (--skip-unit-tests)"
fi

echo "[INFO] Running ProLIF real-case harness on ${#cif_paths[@]} CIFs"
/cluster/work/projects/nn1003k/eirik/conda/analyse_env/bin/python \
    tests/run_tests_scripts/run_prolif_real_cifs.py \
    --run-dir "$run_dir" \
    --run-id "$run_id" \
    "${cif_paths[@]}"

echo "[INFO] ProLIF real-case run completed"
echo "[INFO] Results directory: $run_dir"
echo "[INFO] Summary: $run_dir/prolif_real_cifs_summary.json"
echo "[INFO] ProLIF output: $run_dir/prolif_output"