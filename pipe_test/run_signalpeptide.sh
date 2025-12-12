#!/bin/bash
#SBATCH --job-name=lpmo_signalp              # Navn på jobben i SLURM-køen
#SBATCH --account=nn1003k                    # ENDRE: Ditt prosjekt-ID
#SBATCH --time=01:00:00                      # Maks kjøretid (1 time)
#SBATCH --mem=20G                            # Minne per node
#SBATCH --cpus-per-task=8                    # Antall CPU-kjerner
#SBATCH --ntasks=1
#SBATCH --output=logs/lpmo_signalp_%j.out    # Output-fil (%j = job ID)
#SBATCH --error=logs/lpmo_signalp_%j.err     # Error-fil

################################################################################
# SLURM BATCH SCRIPT FOR MODULE SIGNALPEPTIDE: SIGNAL PEPTIDE PREDICTION
################################################################################

# ============================================================================
# 1. SLURM JOB INFORMASJON
# ============================================================================
echo "======================================================================"
echo "LPMO Signal Peptide Prediction and Trimming (SignalPeptide)"
echo "======================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo ""

# ============================================================================
# 2. MILJØKONFIGURASJON
# ============================================================================
set -euo pipefail

module purge
module load NRIS/CPU
module load hpc-container-wrapper
module load Python/3.11.5-GCCcore-13.2.0

# Pek PATH til lpmo_pipe_env
export PATH=/cluster/projects/nn1003k/eirik/conda/lpmo_pipe_env/bin:$PATH

# Aktiver Python .venv
source /cluster/projects/nn1003k/eirik/Masteroppgave/.venv/bin/activate

# Fail-fast hvis venv ikke faktisk er aktiv
: "${VIRTUAL_ENV:?ERROR: VIRTUAL_ENV is not set (venv not activated?)}"

# CPU / threading control
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"

# Debug / verifiser at riktig python brukes + hard fail hvis for gammel python
echo "VIRTUAL_ENV: $VIRTUAL_ENV"
echo "python: $(which python)"
echo "python --version: $(python --version)"
python - <<'PY'
import sys
print("sys.executable:", sys.executable)
print("sys.version:", sys.version)
assert sys.version_info >= (3,7), "ERROR: Python < 3.7 (need >= 3.7 for __future__ annotations)"
PY
echo ""

# ============================================================================
# 3. ARBEIDSMAPPE
# ============================================================================
cd /cluster/projects/nn1003k/eirik/Masteroppgave/pipe_test || {
    echo "[ERROR] Could not cd to project directory"
    exit 1
}

mkdir -p data/sequences
mkdir -p data/metadata
mkdir -p data/signalpeptide
mkdir -p logs

echo "Working directory: $(pwd)"
echo ""

# ============================================================================
# 4. SIGNALP6 KONFIGURASJON
# ============================================================================
export SIGNALP6_PATH="${SIGNALP6_PATH:-/cluster/home/eisorhus/.local/bin/signalp6}" # endre for annen bruker

if [ ! -f "$SIGNALP6_PATH" ]; then
    echo "[ERROR] SignalP6 not found at: $SIGNALP6_PATH"
    echo "Please set SIGNALP6_PATH environment variable or edit this script"
    exit 1
fi

echo "SignalP6 path: $SIGNALP6_PATH"
echo "SignalP6 version:"
"$SIGNALP6_PATH" --version 2>/dev/null || echo "  (version check unavailable)"
echo ""

# ============================================================================
# 5. GENERATE RUN ID (same format as Module 2)
# ============================================================================
RUN_ID=$(date +"%Y%m%d_%H%M%S")
echo "Generated Run ID: $RUN_ID"
echo ""

# ============================================================================
# 6. PIPELINE INNSTILLINGER
# ============================================================================
SIGNALP6_MODE="slow"
SIGNALP6_FORMAT="txt"
MIN_SP_LENGTH="10"
MAX_SP_LENGTH="70"
TREAT_MISSING_AS_NO_SP="true"

INPUT_FASTA=""
INPUT_METADATA=""
OUTPUT_FASTA="data/sequences/lpmo_mature_${RUN_ID}.fasta"
OUTPUT_METADATA="data/metadata/m2_sequence_3d_metadata_with_sp_${RUN_ID}.csv"

CONFIG_FILE=""

# ============================================================================
# 7. KOMMANDOLINJE ARGUMENTER
# ============================================================================
while [[ $# -gt 0 ]]; do
    case $1 in
        --mode)
            SIGNALP6_MODE="$2"; shift 2 ;;
        --format)
            SIGNALP6_FORMAT="$2"; shift 2 ;;
        --min-sp)
            MIN_SP_LENGTH="$2"; shift 2 ;;
        --max-sp)
            MAX_SP_LENGTH="$2"; shift 2 ;;
        --fasta)
            INPUT_FASTA="$2"; shift 2 ;;
        --metadata)
            INPUT_METADATA="$2"; shift 2 ;;
        --output-fasta)
            OUTPUT_FASTA="$2"; shift 2 ;;
        --output-metadata)
            OUTPUT_METADATA="$2"; shift 2 ;;
        --config)
            CONFIG_FILE="$2"; shift 2 ;;
        *)
            echo "[ERROR] Unknown option: $1"
            exit 1
            ;;
    esac
done

echo "Pipeline configuration:"
echo "  Mode: $SIGNALP6_MODE"
echo "  Format: $SIGNALP6_FORMAT"
echo "  Min SP length: $MIN_SP_LENGTH"
echo "  Max SP length: $MAX_SP_LENGTH"
echo "  Treat missing as no SP: $TREAT_MISSING_AS_NO_SP"
if [ -n "$INPUT_FASTA" ]; then echo "  Input FASTA: $INPUT_FASTA"; fi
if [ -n "$INPUT_METADATA" ]; then echo "  Input metadata: $INPUT_METADATA"; fi
echo "  Output FASTA: $OUTPUT_FASTA"
echo "  Output metadata: $OUTPUT_METADATA"
if [ -n "$CONFIG_FILE" ]; then echo "  Config file: $CONFIG_FILE"; fi
echo ""

# ============================================================================
# 8. KJØR SIGNALPEPTIDE MODULEN
# ============================================================================
echo "Starting Signal Peptide Prediction Pipeline..."
echo "Start time: $(date)"
echo "======================================================================"
echo ""

# Bruk eksplisitt python fra venv og bygg kommandoen som array (ingen eval)
PYTHON="$VIRTUAL_ENV/bin/python"
CMD=("$PYTHON" -m scripts.module_signalpeptide.signalp_pipeline)

if [ -n "$CONFIG_FILE" ]; then
    CMD+=("--config" "$CONFIG_FILE")
fi

if [ -n "$INPUT_FASTA" ]; then
    CMD+=("--fasta" "$INPUT_FASTA")
fi

if [ -n "$INPUT_METADATA" ]; then
    CMD+=("--metadata" "$INPUT_METADATA")
fi

CMD+=("--output-fasta" "$OUTPUT_FASTA" "--output-metadata" "$OUTPUT_METADATA")

echo "Executing: ${CMD[*]}"
echo ""

# Kjør (set -e gjør at jobben stopper ved feil automatisk)
"${CMD[@]}"

echo ""
echo "======================================================================"
echo "Signal Peptide Prediction COMPLETED SUCCESSFULLY"
echo "End time: $(date)"
echo ""
echo "Output files:"
echo "  - Mature FASTA: $OUTPUT_FASTA"
echo "  - Updated metadata: $OUTPUT_METADATA"
echo "  - Parsed predictions: data/signalpeptide/signalp6_parsed.tsv"
echo "  - Run metadata: data/metadata/signalpeptide_run_metadata.json"
echo "======================================================================"
