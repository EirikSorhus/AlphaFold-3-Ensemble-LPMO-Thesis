#!/bin/bash
#SBATCH --job-name=run_tests_geometry_real_cifs
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=00:20:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_geometry_real_cifs_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "[ERROR] This script must be submitted with sbatch, not run directly on the login node."
    echo "[INFO] Example: sbatch tests/run_tests_scripts/test_geometry_real_cifs.sh"
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
run_dir="$out_dir/geometry_real_cifs_${job_suffix}"
mkdir -p "$run_dir"

export XDG_CACHE_HOME="$run_dir/cache"
mkdir -p "$XDG_CACHE_HOME"

cif_paths=(
  "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/NAG4/af3/latest/Q7SCE9_NAG4/seed-4_sample-1/Q7SCE9_NAG4_seed-4_sample-1_model.cif"
  "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/STA6/af3/latest/Q59930_STA6/seed-2_sample-0/Q59930_STA6_seed-2_sample-0_model.cif"
  "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/STA4/af3/latest/A0A0S2GKZ1_STA4/seed-2_sample-2/A0A0S2GKZ1_STA4_seed-2_sample-2_model.cif"
)

if [[ "$#" -gt 0 ]]; then
    for extra_cif in "$@"; do
        cif_paths+=("$extra_cif")
    done
fi

for cif_path in "${cif_paths[@]}"; do
    if [[ ! -f "$cif_path" ]]; then
        echo "[ERROR] CIF does not exist: $cif_path"
        exit 1
    fi
done

echo "[INFO] Running focused downstream geometry tests"
pytest tests/test_mdanalysis_metrics.py -q

echo "[INFO] Running downstream geometry probe on ${#cif_paths[@]} CIFs"
/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin/python \
    tests/run_tests_scripts/run_geometry_real_cifs.py \
    --run-dir "$run_dir" \
    --run-id "geometry_real_cifs_${job_suffix}" \
    "${cif_paths[@]}"

echo "[INFO] Geometry test run completed"
echo "[INFO] Results directory: $run_dir"
echo "[INFO] Pose geometry table: $run_dir/pose_geometry.tsv"
echo "[INFO] Visual PDBs: $run_dir/<case>/geometry_debug.pdb"