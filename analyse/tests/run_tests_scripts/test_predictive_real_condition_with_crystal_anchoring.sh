#!/bin/bash
#SBATCH --job-name=run_tests_predictive_real_condition
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=01:10:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=3
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_predictive_real_condition_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "[ERROR] This script must be submitted with sbatch, not run directly on the login node."
    echo "[INFO] Example: sbatch tests/run_tests_scripts/test_predictive_real_condition_with_crystal_anchoring.sh"
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
run_id="predictive_real_condition_${job_suffix}"
run_dir="$out_dir/${run_id}"
metadata_path="$project_root/input_data/metadata_final_ec_fixed.tsv"
task="all"
n_folds=5
random_state=42

cif_paths=(
    "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core/NAG4/af3/runs/408002/Q7SCE9_NAG4/seed-4_sample-1/Q7SCE9_NAG4_seed-4_sample-1_model.cif"
    "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core/STA6/af3/runs/408006/Q59930_STA6/seed-9_sample-0/Q59930_STA6_seed-9_sample-0_model.cif"
    "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core/STA4/af3/runs/408005/A0A0S2GKZ1_STA4/seed-9_sample-0/A0A0S2GKZ1_STA4_seed-9_sample-0_model.cif"
)

positional_cifs=()
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
        --metadata)
            metadata_path="$2"
            shift 2
            ;;
        --task)
            task="$2"
            shift 2
            ;;
        --n-folds)
            n_folds="$2"
            shift 2
            ;;
        --random-state)
            random_state="$2"
            shift 2
            ;;
        --help)
            cat <<'EOF'
Usage:
    sbatch tests/run_tests_scripts/test_predictive_real_condition_with_crystal_anchoring.sh [options] [CIF ...]

Options:
    --run-dir PATH         Override the top-level output directory.
    --run-id ID            Override the run identifier.
    --metadata PATH        Protein metadata TSV for predictive postprocess.
    --task NAME            Predictive task: all, c1_c4, or substrate.
    --n-folds N            Target number of grouped CV folds.
    --random-state N       Random seed for predictive CV.
    --help                 Show this help text.

Positional CIF arguments replace the default three real AF3 CIFs used for the
analysis-core staging run.
EOF
            exit 0
            ;;
        *)
            positional_cifs+=("$1")
            shift
            ;;
    esac
done

if [[ ${#positional_cifs[@]} -gt 0 ]]; then
    cif_paths=("${positional_cifs[@]}")
fi

mkdir -p "$run_dir"

export XDG_CACHE_HOME="$run_dir/cache"
mkdir -p "$XDG_CACHE_HOME"

if [[ ! -f "$metadata_path" ]]; then
    echo "[ERROR] Metadata TSV does not exist: $metadata_path"
    exit 1
fi

for cif_path in "${cif_paths[@]}"; do
    if [[ ! -f "$cif_path" ]]; then
        echo "[ERROR] CIF does not exist: $cif_path"
        exit 1
    fi
done

echo "[INFO] Running focused predictive + condition-summary + crystal-anchoring pytest suite"
pytest \
    tests/test_predictive_models.py \
    tests/test_predictive_postprocess.py \
    tests/test_cli_run.py \
    tests/test_condition_summary.py \
    tests/test_crystal_anchoring.py \
    tests/test_analysis_orchestrator.py::test_run_analysis_core_writes_qc_geometry_and_reports \
    -q

crystal_run_dir="$run_dir/crystal_anchoring_real_case"
analysis_run_dir="$run_dir/analysis_core_real_case"
summary_validation_dir="$run_dir/summary_validation"
production_output="$analysis_run_dir/production_output"
condition_table_path="$production_output/condition_table.tsv"
predictive_output_root="$production_output"
predictive_summary_path="$predictive_output_root/10_predictive/predictive_summary.json"
generated_condition_table_path="$summary_validation_dir/generated_condition_table.tsv"

echo "[INFO] Running crystal anchoring real-case harness with persistent artifacts"
/cluster/work/projects/nn1003k/eirik/conda/analyse_env/bin/python \
    tests/run_tests_scripts/run_crystal_anchoring_real_cifs.py \
    --run-dir "$crystal_run_dir" \
    --run-id "crystal_anchoring_real_cifs_${job_suffix}"

echo "[INFO] Running analysis-core real-CIF harness to produce a persistent production_output/condition_table.tsv"
/cluster/work/projects/nn1003k/eirik/conda/analyse_env/bin/python \
    tests/run_tests_scripts/run_analysis_core_real_cifs.py \
    --run-dir "$analysis_run_dir" \
    --run-id "analysis_core_real_cifs_${job_suffix}" \
    --del-branch del_a \
    "${cif_paths[@]}"

if [[ ! -f "$condition_table_path" ]]; then
    echo "[ERROR] Expected condition table was not produced: $condition_table_path"
    exit 1
fi

echo "[INFO] Running summary-table validation to leave a generated_condition_table.tsv artifact"
summary_validation_exit_code=0
set +e
/cluster/work/projects/nn1003k/eirik/conda/analyse_env/bin/python \
    tests/run_tests_scripts/run_summary_table_validation.py \
    --production-output "$production_output" \
    --output-dir "$summary_validation_dir"
summary_validation_exit_code=$?
set -e

if [[ ! -f "$generated_condition_table_path" ]]; then
    echo "[ERROR] Summary-table validation did not produce: $generated_condition_table_path"
    exit 1
fi

if [[ "$summary_validation_exit_code" -ne 0 ]]; then
    echo "[WARN] Summary-table validation returned exit code $summary_validation_exit_code; keeping generated tables because this tiny real-CIF run can legitimately lack retained clusters"
fi

echo "[INFO] Running predictive postprocess on the real condition_table.tsv"
/cluster/work/projects/nn1003k/eirik/conda/analyse_env/bin/python \
    -m lpmo_pipeline.cli predictive \
    --condition-table "$condition_table_path" \
    --protein-metadata "$metadata_path" \
    --output "$predictive_output_root" \
    --task "$task" \
    --n-folds "$n_folds" \
    --random-state "$random_state"

if [[ ! -f "$predictive_summary_path" ]]; then
    echo "[ERROR] Predictive postprocess did not produce: $predictive_summary_path"
    exit 1
fi

summary_paths_file="$run_dir/summary_paths.txt"
{
    echo "run_id=$run_id"
    echo "run_dir=$run_dir"
    echo "crystal_run_dir=$crystal_run_dir"
    echo "analysis_run_dir=$analysis_run_dir"
    echo "production_output=$production_output"
    echo "condition_table=$condition_table_path"
    echo "generated_condition_table=$generated_condition_table_path"
    echo "predictive_summary=$predictive_summary_path"
    echo "crystal_summary=$crystal_run_dir/crystal_anchoring_real_cifs_summary.json"
    echo "analysis_summary=$analysis_run_dir/analysis_core_real_cifs_summary.json"
    echo "summary_validation=$summary_validation_dir/summary_table_validation_summary.json"
    echo "summary_validation_exit_code=$summary_validation_exit_code"
} > "$summary_paths_file"

echo "[INFO] Combined predictive + crystal anchoring test run completed"
echo "[INFO] Summary paths: $summary_paths_file"
echo "[INFO] Condition table: $condition_table_path"
echo "[INFO] Generated condition table: $generated_condition_table_path"
echo "[INFO] Predictive summary: $predictive_summary_path"