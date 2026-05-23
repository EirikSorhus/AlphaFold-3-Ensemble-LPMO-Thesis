#!/bin/bash
#SBATCH --job-name=run_tests_contact_eligibility_real_cifs
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=3
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_contact_eligibility_real_cifs_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "[ERROR] This script must be submitted with sbatch, not run directly on the login node."
    echo "[INFO] Example: sbatch tests/run_tests_scripts/test_contact_eligibility_real_cifs.sh"
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
run_dir="$out_dir/contact_eligibility_real_cifs_${job_suffix}"
mkdir -p "$run_dir"

run_id="contact_eligibility_real_cifs_${job_suffix}"
del_branch="del_a"
work_root="/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core"
target="CEL8"
protein_id="A0A0S2GKZ1"
all_runs=0

export XDG_CACHE_HOME="$run_dir/cache"
mkdir -p "$XDG_CACHE_HOME"

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
        --del-branch)
            del_branch="$2"
            shift 2
            ;;
        --work-root)
            work_root="$2"
            shift 2
            ;;
        --target)
            target="$2"
            shift 2
            ;;
        --protein-id)
            protein_id="$2"
            shift 2
            ;;
        --all-runs)
            all_runs=1
            shift
            ;;
        --help)
            cat <<'EOF'
Usage:
    sbatch tests/run_tests_scripts/test_contact_eligibility_real_cifs.sh [options]

Options:
    --run-dir PATH       Override output directory.
    --run-id ID          Override run identifier written to the summary.
    --del-branch VALUE   Analysis branch passed to the production CLI (default: del_a).
    --work-root PATH     Existing structure_pipeline work root to analyze directly.
    --target TARGET      Ligand target to analyze (default: CEL8).
    --protein-id ID      UniProt protein id to analyze (default: A0A0S2GKZ1).
    --all-runs           Scan all numeric runs instead of only the latest run.
    --help               Show this help text.

Default behavior:
    Runs the full contact-eligibility + clustering production path on all poses
    for one protein-ligand condition from the live work_root, without staging.
EOF
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1"
            exit 1
            ;;
    esac
done

if [[ ! -d "$work_root" ]]; then
    echo "[ERROR] work_root does not exist: $work_root"
    exit 1
fi

echo "[INFO] Running focused contact-eligibility unit and integration tests"
pytest tests/test_prolif_ifp.py tests/test_clustering.py tests/test_analysis_orchestrator.py -q

echo "[INFO] Running contact-eligibility real-data harness on all poses for ${protein_id}/${target}"
/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin/python \
    tests/run_tests_scripts/run_contact_eligibility_real_cifs.py \
    --run-dir "$run_dir" \
    --run-id "$run_id" \
    --del-branch "$del_branch" \
    --work-root "$work_root" \
    --target "$target" \
    --protein-id "$protein_id" \
    $( [[ "$all_runs" -eq 1 ]] && printf '%s' '--all-runs' )

echo "[INFO] Contact-eligibility real-data run completed"
echo "[INFO] Results directory: $run_dir"
echo "[INFO] Summary: $run_dir/contact_eligibility_real_cifs_summary.json"
echo "[INFO] Production output: $run_dir/production_output"