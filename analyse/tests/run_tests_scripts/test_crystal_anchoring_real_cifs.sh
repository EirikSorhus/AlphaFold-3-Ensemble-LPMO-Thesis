#!/bin/bash
#SBATCH --job-name=run_tests_crystal_anchoring_real_cifs
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_crystal_anchoring_real_cifs_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "[ERROR] This script must be submitted with sbatch, not run directly on the login node."
    echo "[INFO] Example: sbatch tests/run_tests_scripts/test_crystal_anchoring_real_cifs.sh"
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
run_id="crystal_anchoring_real_cifs_${job_suffix}"
run_dir="$out_dir/${run_id}"

skip_unit_tests=0
representative_cif=""
protein_id="A0A0S2GKZ1"
ligand_id="CEL4"

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
        --representative-cif)
            representative_cif="$2"
            shift 2
            ;;
        --protein-id)
            protein_id="$2"
            shift 2
            ;;
        --ligand-id)
            ligand_id="$2"
            shift 2
            ;;
        --skip-unit-tests)
            skip_unit_tests=1
            shift
            ;;
        --help)
            cat <<'EOF'
Usage:
    sbatch tests/run_tests_scripts/test_crystal_anchoring_real_cifs.sh [options]

Options:
    --run-dir PATH             Override output directory.
    --run-id ID                Override run identifier written to the summary.
    --representative-cif PATH  Override the auto-selected best AF3 representative CIF.
    --protein-id ID            Protein identifier to resolve from the crystal index.
    --ligand-id ID             Ligand identifier written into the report metadata.
    --skip-unit-tests          Skip pytest tests/test_crystal_anchoring.py.
    --help                     Show this help text.
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

if [[ -n "$representative_cif" && ! -f "$representative_cif" ]]; then
    echo "[ERROR] Representative CIF does not exist: $representative_cif"
    exit 1
fi

if [[ "$skip_unit_tests" -eq 0 ]]; then
    echo "[INFO] Running focused crystal anchoring unit tests"
    pytest tests/test_crystal_anchoring.py -q
else
    echo "[INFO] Skipping focused crystal anchoring unit tests (--skip-unit-tests)"
fi

echo "[INFO] Running crystal anchoring real-case harness"
python_args=(
    tests/run_tests_scripts/run_crystal_anchoring_real_cifs.py
    --run-dir "$run_dir"
    --run-id "$run_id"
    --protein-id "$protein_id"
    --ligand-id "$ligand_id"
)

if [[ -n "$representative_cif" ]]; then
    python_args+=(--representative-cif "$representative_cif")
fi

/cluster/work/projects/nn1003k/eirik/conda/analyse_env/bin/python "${python_args[@]}"

echo "[INFO] Crystal anchoring real-case run completed"
echo "[INFO] Results directory: $run_dir"
echo "[INFO] Summary: $run_dir/crystal_anchoring_real_cifs_summary.json"
echo "[INFO] Crystal report: $run_dir/crystal_anchoring_output/crystal_reference_screen.json"