#!/usr/bin/env bash
#SBATCH --job-name=rf3_fold
#SBATCH --account=nn1003k
#SBATCH --partition=accel
#SBATCH --gpus=1
#SBATCH --time=00:20:00
#SBATCH --mem-per-gpu=30G
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/pipe_test/scripts/module_5/rf3/slurm_%j.out
#SBATCH --error=/cluster/work/projects/nn1003k/eirik/Masteroppgave/pipe_test/scripts/module_5/rf3/slurm_%j.err

set -euo pipefail
module load NRIS/GPU

require_file() {
  local path="$1"
  local label="$2"
  if [[ ! -f "${path}" ]]; then
    echo "ERROR: Missing ${label}: ${path}" >&2
    exit 1
  fi
}

FOUNDRY_ROOT=/cluster/projects/nn1003k/prog/foundry
SIF=${FOUNDRY_ROOT}/foundry.sif

CKPT_DIR=${FOUNDRY_ROOT}/checkpoints
CKPT_FILE=${CKPT_DIR}/rf3_foundry_01_24_latest_remapped.ckpt

OUTBASE=/cluster/work/projects/nn1003k/eirik/Masteroppgave/pipe_test/scripts/module_5/rf3
OUT=${OUTBASE}/${SLURM_JOB_ID}

TMPBASE=/cluster/work/projects/nn1003k/eirik/tmp
TMPDIR=${TMPBASE}/foundry_${USER}

INPUT_TEMPLATE=${FOUNDRY_ROOT}/models/rf3/tests/data/5vht_from_json.json
MSA_PATH=/cluster/work/projects/nn1003k/eirik/Masteroppgave/pipe_test/scripts/module_5/rf3/input/UniProtIDs_Q7S111.a3m

INPUT_JSON="${OUTBASE}/input/rf3_test.json"

mkdir -p "${OUT}" "${OUTBASE}" "${TMPDIR}" #Fjernet "${CKPT_DIR}" fordi jeg ikke har skrivetilgang i mappen
require_file "${SIF}" "Apptainer image"
require_file "${CKPT_FILE}" "RF3 checkpoint"
require_file "${INPUT_TEMPLATE}" "input template JSON"
require_file "${MSA_PATH}" "MSA file"
require_file "${INPUT_JSON}" "RF3 input JSON"

set -a
[[ -f "${FOUNDRY_ROOT}/.env" ]] && source "${FOUNDRY_ROOT}/.env"
set +a

export APPTAINERENV_FOUNDRY_CHECKPOINT_DIRS="${CKPT_DIR}"
export APPTAINERENV_XDG_CACHE_HOME="${OUT}/cache"
export APPTAINERENV_TMPDIR="${TMPDIR}"
# Prefer source tree configs over installed package paths
export APPTAINERENV_PYTHONPATH="${FOUNDRY_ROOT}/src:${FOUNDRY_ROOT}/models/rf3/src:${FOUNDRY_ROOT}/models/mpnn/src:${FOUNDRY_ROOT}/models/rfd3/src:${PYTHONPATH:-}"
export APPTAINERENV_PROJECT_ROOT="${FOUNDRY_ROOT}"
export APPTAINERENV_HYDRA_FULL_ERROR=1
mkdir -p "${OUT}/cache"

TRITON_VERSION=3.5.0
TRITON_SITE=${TMPDIR}/triton_site
mkdir -p "${TRITON_SITE}"


echo "=== Sanity: show checkpoint + env inside container ==="
apptainer exec --nv \
  --bind "${FOUNDRY_ROOT}:${FOUNDRY_ROOT}" \
  --bind "${TMPBASE}:${TMPBASE}" \
  --bind "/cluster/work/projects/nn1003k:/cluster/work/projects/nn1003k" \
  "${SIF}" bash -lc "
    echo FOUNDRY_CHECKPOINT_DIRS=\$FOUNDRY_CHECKPOINT_DIRS
    ls -lh '${CKPT_FILE}'
    python -c \"import torch; print(torch.__version__); print(torch.cuda.is_available())\"
  "

echo "=== Run RF3 fold ==="
apptainer exec --nv \
  --bind "${FOUNDRY_ROOT}:${FOUNDRY_ROOT}" \
  --bind "${TMPBASE}:${TMPBASE}" \
  --bind "/cluster/work/projects/nn1003k:/cluster/work/projects/nn1003k" \
  "${SIF}" bash -lc "
    set -e
    export PYTHONPATH=\"${TRITON_SITE}:\${PYTHONPATH:-}\"
    if ! python - <<'PY'
from importlib.metadata import version, PackageNotFoundError
import sys
try:
    v = version('triton')
    sys.exit(0 if v.startswith('3.5.') else 1)
except PackageNotFoundError:
    sys.exit(1)
PY
    then
      echo 'Installing triton==${TRITON_VERSION} into ${TRITON_SITE}...'
      pip install --no-cache-dir --target \"${TRITON_SITE}\" \"triton==${TRITON_VERSION}\"
      export PYTHONPATH=\"${TRITON_SITE}:\${PYTHONPATH:-}\"
    fi
    python \"${FOUNDRY_ROOT}/models/rf3/src/rf3/cli.py\" fold \
      inputs=\"${INPUT_JSON}\" \
      ckpt_path=\"${CKPT_FILE}\" \
      out_dir=\"${OUT}\" \
      skip_existing=False \
      dump_predictions=True \
      dump_trajectories=False \
      n_recycles=10 \
      diffusion_batch_size=1 \
      num_steps=50
  "

echo "Done. Outputs in: ${OUT}"

