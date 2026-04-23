#!/bin/bash
#SBATCH --job-name=run_tests_cif_to_pdb
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_cif_to_pdb_%j.log

set -euo pipefail
# Activate conda environment
export PATH="/cluster/work/projects/nn1003k/eirik/conda/analyse_env/bin:$PATH"

cd /cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse

out_dir="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/tests_results/"

mkdir -p "$out_dir"

job_suffix="${SLURM_JOB_ID:-manual}"
run_dir="$out_dir/cif_to_pdb_${job_suffix}"
mkdir -p "$run_dir"

export RUN_DIR="$run_dir"
python - <<'PY'
import os
from pathlib import Path
from lpmo_pipeline.io.cif_to_pdb import convert_cif_to_pdb

input_cif = Path("/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/tests_results/mapping_contracts_554656/normalize_for_mapping/normalized.cif")
output_dir = Path(os.environ["RUN_DIR"])
success, output_path = convert_cif_to_pdb(input_cif, output_dir)
print("success:", success)
print("output:", output_path)
PY
