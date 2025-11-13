#!/usr/bin/env bash
#SBATCH --job-name=AF3_7o1r_1_chain_pmm
#SBATCH --account=nn1003k
#SBATCH --partition=a100
#SBATCH --gpus=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-gpu=80G

set -euo pipefail
module purge
module --force swap StdEnv Zen2Env
module load CUDA/12.6.0

mkdir -p "${SUBMITDIR}/result"

AF3_DIR=/cluster/projects/nn1003k/prog/sif_alphafold3
AF3_IMAGE=${AF3_DIR}/alphafold3_v3.0.1.sif
AF3_INPUT_DIR=${SUBMITDIR}/input
AF3_OUTPUT_DIR=${SUBMITDIR}/result
AF3_MODEL_PARAMETERS_DIR=${AF3_DIR}/weights
AF3_DATABASES_DIR=/cluster/shared/databases/AlphaFold3/2025-04-15

singularity exec --nv \
  --pwd /app/alphafold \
  --bind "${AF3_INPUT_DIR}:/root/af_input" \
  --bind "${AF3_OUTPUT_DIR}:/root/af_output" \
  --bind "${AF3_MODEL_PARAMETERS_DIR}:/root/models" \
  --bind "${AF3_DATABASES_DIR}:/root/public_databases" \
  --bind "${AF3_DIR}:/root/image" \
  "${AF3_IMAGE}" \
  python run_alphafold.py \
    --json_path=/root/af_input/input_msa.json \
    --model_dir=/root/models \
    --db_dir=/root/public_databases \
    --output_dir=/root/af_output
