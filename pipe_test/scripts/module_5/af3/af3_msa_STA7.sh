#!/usr/bin/env bash
#SBATCH --job-name=AF3_MSA
#SBATCH --partition=normal
#SBATCH --account=nn1003K
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=02:00:00

module load NRIS/CPU
set -euo pipefail


SUBMIT_DIR=${SUBMITDIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}

TMP_ON_SCRATCH="${SCRATCH}/af3_tmp_${SLURM_JOB_ID}"
mkdir -p "${SCRATCH}/input" "${SCRATCH}/result" "${TMP_ON_SCRATCH}"
export TMPDIR="${TMP_ON_SCRATCH}"
export APPTAINERENV_TMPDIR="${TMPDIR}"

AF3_DIR=/cluster/projects/nn1003k/prog/af3
AF3_IMAGE=${AF3_DIR}/af3_cpu_amd64.sif
AF3_INPUT_DIR=${SUBMIT_DIR}/input
AF3_OUTPUT_DIR=${SUBMIT_DIR}/result
AF3_MODEL_PARAMETERS_DIR=${AF3_DIR}/weights
AF3_DATABASES_SQUASHFS=/cluster/work/shared/alphafold_uncompressed.squashfs

# copy input to fast scratch
rsync -av --delete "${AF3_INPUT_DIR}/" "${SCRATCH}/input/"

# copy MSA json back when we exit (only if it exists)
trap 'shopt -s nullglob; files=(${SCRATCH}/*/*/*_data.json); if (( ${#files[@]} )); then cp -f "${files[0]}" "${AF3_INPUT_DIR}/msa_STA7.json"; fi' EXIT

apptainer exec --cleanenv \
  --env PATH=/hmmer/bin:/opt/af3-venv/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
  --env JAX_PLATFORM_NAME=cpu \
  --pwd /opt/alphafold3 \
  --bind "${SCRATCH}/input:/root/af_input" \
  --bind "${SCRATCH}/result:/root/af_output" \
  --bind "${TMP_ON_SCRATCH}:/tmp" \
  --bind "${AF3_MODEL_PARAMETERS_DIR}:/root/models" \
  --bind "${AF3_DATABASES_SQUASHFS}:/root/public_databases:image-src=/public_databases" \
  "${AF3_IMAGE}" \
  /opt/af3-venv/bin/python /opt/alphafold3/run_alphafold.py \
    --json_path=/root/af_input/af3_test_STA7.json \
    --model_dir=/root/models \
    --db_dir=/root/public_databases \
    --output_dir=/root/af_output \
    --run_inference=false \
    --jackhmmer_n_cpu=8 \
    --jackhmmer_binary_path=/hmmer/bin/jackhmmer \
    --nhmmer_binary_path=/hmmer/bin/nhmmer \
    --hmmalign_binary_path=/hmmer/bin/hmmalign \
    --hmmsearch_binary_path=/hmmer/bin/hmmsearch \
    --hmmbuild_binary_path=/hmmer/bin/hmmbuild

