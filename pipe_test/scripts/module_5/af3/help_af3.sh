#!/usr/bin/env bash
#SBATCH --job-name=AF3_HELPFULL
#SBATCH --partition=normal
#SBATCH --account=nn1003K
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:05:00

module load NRIS/CPU
set -euo pipefail

SUBMIT_DIR=${SUBMITDIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}

# Skriv help-hit
OUTPUT_FILE="${SUBMIT_DIR}/alphafold3_help.txt"

# Midlertidig scratch (greit å ha, i tilfelle containeren forventer /tmp)
TMP_ON_SCRATCH="${SCRATCH}/af3_tmp_${SLURM_JOB_ID}"
mkdir -p "${TMP_ON_SCRATCH}"
export TMPDIR="${TMP_ON_SCRATCH}"
export APPTAINERENV_TMPDIR="${TMPDIR}"

AF3_DIR=/cluster/projects/nn1003k/prog/af3
AF3_IMAGE=${AF3_DIR}/af3_cpu_amd64.sif
AF3_MODEL_PARAMETERS_DIR=${AF3_DIR}/weights
AF3_DATABASES_SQUASHFS=/cluster/work/shared/alphafold_uncompressed.squashfs

echo "Writing AlphaFold3 help to: ${OUTPUT_FILE}"

apptainer exec --cleanenv \
  --env PATH=/hmmer/bin:/opt/af3-venv/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
  --env JAX_PLATFORM_NAME=cpu \
  --pwd /opt/alphafold3 \
  --bind "${TMP_ON_SCRATCH}:/tmp" \
  --bind "${AF3_MODEL_PARAMETERS_DIR}:/root/models" \
  --bind "${AF3_DATABASES_SQUASHFS}:/root/public_databases:image-src=/public_databases" \
  "${AF3_IMAGE}" \
  /opt/af3-venv/bin/python /opt/alphafold3/run_alphafold.py --helpfull \
  > "${OUTPUT_FILE}" 2>&1

echo "Done."