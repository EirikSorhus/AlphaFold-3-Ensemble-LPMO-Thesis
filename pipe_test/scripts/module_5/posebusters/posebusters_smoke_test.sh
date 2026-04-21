#!/usr/bin/env bash
#SBATCH --job-name=posebusters_smoke
#SBATCH --account=nn1003k
#SBATCH --partition=normal
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/pipe_test/scripts/module_5/posebusters/slurm_%j.out
#SBATCH --error=/cluster/work/projects/nn1003k/eirik/Masteroppgave/pipe_test/scripts/module_5/posebusters/slurm_%j.err

set -euo pipefail

# Usage:
#   sbatch posebusters_smoke_test.sh [input_structure]
#   bash posebusters_smoke_test.sh [input_structure]
# Optional environment overrides:
#   POSEBUSTERS_SIF=/path/to/posebusters.sif
#   PB_OUTFMT=csv|long|short   (default: csv)
#   PB_FULL_REPORT=1|0         (default: 1)

require_file() {
  local path="$1"
  local label="$2"
  if [[ ! -f "${path}" ]]; then
    echo "ERROR: Missing ${label}: ${path}" >&2
    exit 1
  fi
}

resolve_sif() {
  local primary=/cluster/projects/nn1003k/prog/posebusters/posebusters.sif
  local secondary=/cluster/projects/nn1003k/prog/posebusters/build/posebusters.sif

  if [[ -n "${POSEBUSTERS_SIF:-}" ]]; then
    echo "${POSEBUSTERS_SIF}"
    return 0
  fi

  if [[ -f "${primary}" ]]; then
    echo "${primary}"
    return 0
  fi

  if [[ -f "${secondary}" ]]; then
    echo "${secondary}"
    return 0
  fi

  echo "ERROR: Could not find PoseBusters runtime image." >&2
  echo "Tried:" >&2
  echo "  1) POSEBUSTERS_SIF (env var)" >&2
  echo "  2) ${primary}" >&2
  echo "  3) ${secondary}" >&2
  return 1
}

# In sbatch jobs, BASH_SOURCE can point to Slurm spool copy under /var/spool.
# Use SLURM_SUBMIT_DIR so outputs are written to a user-writable project path.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

DEFAULT_INPUT=/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/CEL6/af3/latest/B6EQJ6_CEL6/B6EQJ6_CEL6_model.cif
INPUT_PATH="${1:-${DEFAULT_INPUT}}"
PB_OUTFMT="${PB_OUTFMT:-csv}"
PB_FULL_REPORT="${PB_FULL_REPORT:-1}"

SIF="$(resolve_sif)"
OUTBASE="${SCRIPT_DIR}/output"
RUN_TAG="${SLURM_JOB_ID:-manual_$(date +%Y%m%d_%H%M%S)}"
OUT_DIR="${OUTBASE}/posebusters_smoke_${RUN_TAG}"
HELP_TXT="${OUT_DIR}/posebusters_help.txt"
RUN_LOG="${OUT_DIR}/posebusters_run.log"
RESULT_CSV="${OUT_DIR}/posebusters_results.csv"
INPUT_FOR_BUST="${INPUT_PATH}"

mkdir -p "${OUT_DIR}"

require_file "${SIF}" "Apptainer image"
require_file "${INPUT_PATH}" "input structure"

if [[ "${PB_OUTFMT}" != "short" && "${PB_OUTFMT}" != "long" && "${PB_OUTFMT}" != "csv" ]]; then
  echo "ERROR: PB_OUTFMT must be one of: short, long, csv (got '${PB_OUTFMT}')." >&2
  exit 1
fi

echo "=== PoseBusters smoke test ==="
echo "SIF:         ${SIF}"
echo "Input file:  ${INPUT_PATH}"
echo "Output dir:  ${OUT_DIR}"
echo "Output fmt:  ${PB_OUTFMT}"
echo "Full report: ${PB_FULL_REPORT}"

ext="${INPUT_PATH##*.}"
ext="${ext,,}"
if [[ "${ext}" != "pdb" && "${ext}" != "sdf" && "${ext}" != "mol" && "${ext}" != "mol2" ]]; then
  echo "ERROR: Unsupported input format '${ext}'." >&2
  echo "PoseBusters in this environment supports: .sdf, .mol, .mol2, .pdb" >&2
  echo "Current input: ${INPUT_PATH}" >&2
  if [[ "${ext}" == "cif" || "${ext}" == "mmcif" ]]; then
    echo "Hint: submit this smoke test again when a .pdb export is available." >&2
  fi
  exit 1
fi

echo ""
echo "[1/3] Collecting PoseBusters help..."
MODE=""
if apptainer exec --cleanenv "${SIF}" bash -lc 'command -v bust >/dev/null 2>&1'; then
  MODE="exec"
  apptainer exec --cleanenv "${SIF}" bust --help >"${HELP_TXT}" 2>&1
elif apptainer run --cleanenv "${SIF}" --help >"${HELP_TXT}" 2>&1; then
  MODE="run"
else
  echo "ERROR: Could not find a working PoseBusters CLI via apptainer exec or run." >&2
  exit 1
fi

echo "[2/3] Running PoseBusters..."
FULL_REPORT_FLAG=()
if [[ "${PB_FULL_REPORT}" == "1" ]]; then
  FULL_REPORT_FLAG+=(--full-report)
fi

if [[ "${MODE}" == "exec" ]]; then
  apptainer exec --cleanenv \
    --bind "/cluster/work/projects/nn1003k:/cluster/work/projects/nn1003k" \
    "${SIF}" \
    bust "${INPUT_FOR_BUST}" --outfmt "${PB_OUTFMT}" --output "${RESULT_CSV}" "${FULL_REPORT_FLAG[@]}" \
    >"${RUN_LOG}" 2>&1
else
  apptainer run --cleanenv \
    --bind "/cluster/work/projects/nn1003k:/cluster/work/projects/nn1003k" \
    "${SIF}" \
    "${INPUT_FOR_BUST}" --outfmt "${PB_OUTFMT}" --output "${RESULT_CSV}" "${FULL_REPORT_FLAG[@]}" \
    >"${RUN_LOG}" 2>&1
fi

echo "[3/3] Validating output..."
require_file "${RUN_LOG}" "run log"
require_file "${RESULT_CSV}" "PoseBusters output"

if [[ ! -s "${RESULT_CSV}" ]]; then
  echo "ERROR: Output file exists but is empty: ${RESULT_CSV}" >&2
  exit 1
fi

line_count=$(wc -l <"${RESULT_CSV}")
if [[ "${line_count}" -lt 2 ]]; then
  echo "ERROR: Expected CSV header + at least one data line in ${RESULT_CSV}, got ${line_count} line(s)." >&2
  echo "Run log: ${RUN_LOG}" >&2
  exit 1
fi

header=$(head -n 1 "${RESULT_CSV}")
if [[ "${header}" != *","* ]]; then
  echo "WARNING: CSV header does not contain commas. Check format manually: ${RESULT_CSV}" >&2
fi

echo "Smoke test PASSED."
echo "Mode:      ${MODE}"
echo "Input run: ${INPUT_FOR_BUST}"
echo "Help file: ${HELP_TXT}"
echo "Run log:   ${RUN_LOG}"
echo "Result:    ${RESULT_CSV}"
echo "Lines:     ${line_count}"
