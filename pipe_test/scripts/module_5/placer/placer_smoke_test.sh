#!/usr/bin/env bash
#SBATCH --job-name=placer_smoke
#SBATCH --account=nn1003k
#SBATCH --partition=accel
#SBATCH --gpus=1
#SBATCH --time=00:02:00
#SBATCH --mem-per-gpu=20G
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/pipe_test/scripts/module_5/placer/slurm_%j.out
#SBATCH --error=/cluster/work/projects/nn1003k/eirik/Masteroppgave/pipe_test/scripts/module_5/placer/slurm_%j.err

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

# In sbatch jobs, BASH_SOURCE can point to Slurm spool copy under /var/spool.
# Use SLURM_SUBMIT_DIR so outputs are written to a user-writable project path.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PLACER_DIR=/cluster/projects/nn1003k/prog/placer/build

DEFAULT_INPUT_CIF=/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/CEL6/af3/latest/B6EQJ6_CEL6/B6EQJ6_CEL6_model.cif
INPUT_CIF="${1:-${DEFAULT_INPUT_CIF}}"
NSAMPLES="${NSAMPLES:-5}"

if [[ -n "${PLACER_SIF:-}" ]]; then
  SIF="${PLACER_SIF}"
elif [[ -f "${PLACER_DIR}/placer.sif" ]]; then
  SIF="${PLACER_DIR}/placer.sif"
else
  echo "ERROR: Missing PLACER runtime image: ${PLACER_DIR}/placer.sif" >&2
  echo "'placer-base.sif' is a base image and does not include the PLACER CLI." >&2
  echo "Build failed logs exist in ${PLACER_DIR}/slurm-*.out (see openbabel header errors)." >&2
  echo "Use the original build process there to produce placer.sif, or set PLACER_SIF to another existing runtime image." >&2
  exit 1
fi

OUTBASE="${SCRIPT_DIR}/output"
RUN_TAG="${SLURM_JOB_ID:-manual_$(date +%Y%m%d_%H%M%S)}"
OUT_DIR="${OUTBASE}/placer_smoke_${RUN_TAG}"
HELP_TXT="${OUT_DIR}/placer_help.txt"
RUN_LOG="${OUT_DIR}/placer_run.log"

mkdir -p "${OUT_DIR}"

require_file "${SIF}" "Apptainer image"
require_file "${INPUT_CIF}" "input CIF"

echo "=== PLACER smoke test ==="
echo "SIF:       ${SIF}"
echo "Input CIF: ${INPUT_CIF}"
echo "Output:    ${OUT_DIR}"
echo "nsamples:  ${NSAMPLES}"

echo ""
echo "[1/3] Collecting PLACER help..."
if ! apptainer run --nv --cleanenv "${SIF}" --help >"${HELP_TXT}" 2>&1; then
  apptainer exec --nv --cleanenv "${SIF}" bash -lc '
    set -e
    if command -v placer >/dev/null 2>&1; then
      placer --help
    elif [[ -f /opt/PLACER/run_PLACER.py ]]; then
      python /opt/PLACER/run_PLACER.py --help
    else
      echo "ERROR: Could not find placer CLI or /opt/PLACER/run_PLACER.py in container." >&2
      exit 1
    fi
  ' >"${HELP_TXT}" 2>&1
fi

echo "[2/3] Running PLACER on test CIF..."
apptainer exec --nv --cleanenv \
  --bind "/cluster/work/projects/nn1003k:/cluster/work/projects/nn1003k" \
  --bind "${OUT_DIR}:/out" \
  "${SIF}" bash -lc "
    set -euo pipefail
    if command -v placer >/dev/null 2>&1; then
      placer --ifile '${INPUT_CIF}' --odir /out --nsamples '${NSAMPLES}' --suffix smoke
    elif [[ -f /opt/PLACER/run_PLACER.py ]]; then
      python /opt/PLACER/run_PLACER.py --ifile '${INPUT_CIF}' --odir /out --nsamples '${NSAMPLES}' --suffix smoke
    else
      echo 'ERROR: Could not find placer CLI or /opt/PLACER/run_PLACER.py in container.' >&2
      exit 1
    fi
  " >"${RUN_LOG}" 2>&1

echo "[3/3] Validating output..."
model_count=$(find "${OUT_DIR}" -type f \( -name "*.pdb" -o -name "*.cif" \) | wc -l)
score_count=$(find "${OUT_DIR}" -type f \( -name "*.csv" -o -name "*score*" -o -name "*rank*" \) | wc -l)

if [[ "${model_count}" -lt 1 ]]; then
  echo "ERROR: No model files (.pdb/.cif) were produced by PLACER." >&2
  echo "Check run log: ${RUN_LOG}" >&2
  exit 1
fi

if [[ "${score_count}" -lt 1 ]]; then
  echo "WARNING: No obvious score/rank files found. PLACER run produced models, but inspect logs/output format manually."
fi

echo "Smoke test PASSED."
echo "Help file: ${HELP_TXT}"
echo "Run log:   ${RUN_LOG}"
echo "Models:    ${model_count}"
echo "Scores:    ${score_count}"

echo ""
echo "Ligand note: External MOL/SDF is usually NOT required for this smoke test because ligand coordinates are in CIF."
echo "Use --ligand_file (or CCD hints) only if atom typing/connectivity issues appear."
