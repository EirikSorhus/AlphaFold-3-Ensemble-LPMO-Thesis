#!/usr/bin/env bash
#SBATCH --job-name=AF3_lpmo_test
#SBATCH --account=nn1003k
#SBATCH --partition=accel
#SBATCH --gpus=1
#SBATCH --time=00:30:00
#SBATCH --mem-per-gpu=20G

set -euo pipefail
module load NRIS/GPU

SUBMIT_DIR=${SUBMITDIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}
mkdir -p "${SUBMIT_DIR}/output"

AF3_DIR=/cluster/projects/nn1003k/prog/af3
AF3_IMAGE=${AF3_DIR}/af3_gpu_arm64.sif
AF3_INPUT_DIR=${SUBMIT_DIR}/input
AF3_OUTPUT_DIR=${SUBMIT_DIR}/output
AF3_MODEL_PARAMETERS_DIR=${AF3_DIR}/weights
AF3_DATABASES_SQUASHFS=/cluster/work/shared/alphafold_uncompressed.squashfs

apptainer exec --nv --cleanenv \
  --pwd /opt/alphafold3 \
  --bind "${AF3_INPUT_DIR}:/root/af_input" \
  --bind "${AF3_OUTPUT_DIR}:/root/af_output" \
  --bind "${AF3_MODEL_PARAMETERS_DIR}:/root/models" \
  --bind "${AF3_DATABASES_SQUASHFS}:/root/public_databases:image-src=/public_databases" \
  "${AF3_IMAGE}" \
  /opt/af3-venv/bin/python /opt/alphafold3/run_alphafold.py \
    --json_path=/root/af_input/msa_STA7.json \
    --model_dir=/root/models \
    --db_dir=/root/public_databases \
    --num_diffusion_samples=3 \
    --output_dir=/root/af_output

