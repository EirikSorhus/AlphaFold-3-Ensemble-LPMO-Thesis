#!/usr/bin/env bash
#SBATCH --job-name=AF3_MSA_7o1r_1_chain_pmm
#SBATCH --partition=bigmem
#SBATCH --account=NN1003K
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=30G
#SBATCH --time=24:00:00

set -euo pipefail
TMP_ON_SCRATCH="${SCRATCH}/af3_tmp_${SLURM_JOB_ID}"
mkdir -p "${SCRATCH}/input" "${SCRATCH}/result" "${TMP_ON_SCRATCH}"
export TMPDIR="${TMP_ON_SCRATCH}"
export SINGULARITYENV_TMPDIR="${TMPDIR}"

AF3_DIR=/cluster/projects/nn1003k/prog/sif_alphafold3
AF3_IMAGE=${AF3_DIR}/alphafold3_v3.0.1.sif
AF3_MODEL_PARAMETERS_DIR=${AF3_DIR}/weights
AF3_DATABASES_DIR=/cluster/shared/databases/AlphaFold3/2025-04-15

# copy input to fast scratch
rsync -av --delete "${SUBMITDIR}/input/" "${SCRATCH}/input/"

# copy MSA json back when we exit
trap 'cp -r ${SCRATCH}/*/*/*_data.json ${SUBMITDIR}/input/input_msa.json' EXIT

singularity exec \
  --pwd /app/alphafold \
  --bind "${SCRATCH}/input:/root/af_input" \
  --bind "${SCRATCH}/result:/root/af_output" \
  --bind "${TMP_ON_SCRATCH}:/tmp" \
  --bind "${AF3_MODEL_PARAMETERS_DIR}:/root/models" \
  --bind "${AF3_DATABASES_DIR}:/root/public_databases" \
  "${AF3_IMAGE}" \
  python run_alphafold.py \
      --json_path=/root/af_input/alphafold_input.json \
      --model_dir=/root/models \
      --db_dir=/root/public_databases \
      --output_dir=/root/af_output \
      --norun_inference \
      --jackhmmer_n_cpu=8
