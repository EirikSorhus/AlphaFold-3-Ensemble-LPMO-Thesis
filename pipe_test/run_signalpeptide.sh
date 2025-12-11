#!/bin/bash
#SBATCH --job-name=lpmo_signalp        # Navn på jobben i SLURM-køen
#SBATCH --account=nn1003k              # ENDRE: Ditt prosjekt-ID
#SBATCH --time=01:00:00                # Maks kjøretid (1 time)
#SBATCH --mem=20                       # Minne per node. min 14, men ber om mer for sikkerhets skyld. Vet ikke nøyaktig hvor mye signalp6 bruker.
#SBATCH --cpus-per-task=8              # Antall CPU-kjerner
#SBATCH --ntasks=1
#SBATCH --output=logs/lpmo_signalp_%j.out   # Output-fil (%j = job ID)
#SBATCH --error=logs/lpmo_signalp_%j.err    # Error-fil

################################################################################
# SLURM BATCH SCRIPT FOR MODULE SIGNALPEPTIDE: SIGNAL PEPTIDE PREDICTION
################################################################################
#
# Dette scriptet kjører module_signalpeptide som bruker SignalP6 til å
# predikere og trimme signalpeptider fra LPMO-sekvenser produsert av Module 2.
#
# Steg:
#   1. Leser M2 output FASTA (lpmo_all_*_raw.fasta) og metadata
#   2. Grupperer sekvenser etter kingdom (Eukaryota vs other)
#   3. Kjører SignalP6 per gruppe med riktig --organism flag
#   4. Parser SignalP6 output
#   5. Genererer trimmet FASTA (lpmo_mature.fasta)
#   6. Oppdaterer metadata med signal peptide-informasjon
#
# VIKTIG: Les gjennom og tilpass følgende seksjoner før du kjører:
#   1. SLURM-parametere (linjene over)
#   2. Miljøkonfigurasjon (Python-versjon)
#   3. SignalP6 binær-path (PÅKREVD)
#   4. Pipeline-alternativer (valgfritt)
#
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
# Setter pipefail slik at scriptet stopper ved feil
set -euo pipefail

# Last nødvendig Python-miljø
# ENDRE til din Python-versjon/modul hvis nødvendig:
module purge
module load NRIS/CPU
module load hpc-container-wrapper

# Pek PATH til lpmo_pipe_env
export PATH=/cluster/projects/nn1003k/eirik/conda/lpmo_pipe_env/bin:$PATH

# Aktiver Python .venv
source /cluster/projects/nn1003k/eirik/Masteroppgave/.venv/bin/activate

echo "Python version: $(python --version)"
echo "Python location: $(which python)"
echo ""

# ============================================================================
# 3. ARBEIDSMAPPE
# ============================================================================
# Sørg for at du er i riktig mappe (project_root)
# ENDRE til din prosjektmappe:
cd /cluster/projects/nn1003k/eirik/Masteroppgave/pipe_test || {
    echo "[ERROR] Could not cd to project directory"
    exit 1
}
# Opprett nødvendige mapper hvis de ikke finnes
mkdir -p data/sequences
mkdir -p data/metadata
mkdir -p data/signalpeptide
mkdir -p logs

echo "Working directory: $(pwd)"
echo ""

# ============================================================================
# 4. SIGNALP6 KONFIGURASJON
# ============================================================================
# SignalP6 binær-path (PÅKREVD)
# ENDRE til hvor SignalP6 er installert på ditt system:
export SIGNALP6_PATH="${SIGNALP6_PATH:-/cluster/home/eisorhus/.local/bin/signalp6}"

# Sjekk at SignalP6 finnes
if [ ! -f "$SIGNALP6_PATH" ]; then
    echo "[ERROR] SignalP6 not found at: $SIGNALP6_PATH"
    echo "Please set SIGNALP6_PATH environment variable or edit this script"
    exit 1
fi

echo "SignalP6 path: $SIGNALP6_PATH"
echo "SignalP6 version:"
$SIGNALP6_PATH --version 2>/dev/null || echo "  (version check unavailable)"
echo ""

# ============================================================================
# 5. PIPELINE INNSTILLINGER
# ============================================================================
# Disse parametrene kan overstyres via kommandolinje eller config file

# Standard innstillinger:
SIGNALP6_MODE="slow"           # "fast" eller "slow" Jeg bruker "slow" for bedre nøyaktighet, men "fast" kan brukes for raskere kjøring.
SIGNALP6_FORMAT="txt"          # "txt" eller annet format SignalP6 støtter
MIN_SP_LENGTH="10"             # Minimum signalpeptid-lengde
MAX_SP_LENGTH="70"             # Maksimum signalpeptid-lengde
TREAT_MISSING_AS_NO_SP="true"  # Behandle sekvenser uten prediksjon som "no signal peptide"

# Input/output filer (auto-detect eller override):
INPUT_FASTA=""                 # Tomt = auto-detect latest lpmo_all_*_raw.fasta
INPUT_METADATA=""              # Tomt = auto-detect latest m2_sequence_3d_metadata_*.csv
OUTPUT_FASTA="data/sequences/lpmo_mature.fasta"
OUTPUT_METADATA="data/metadata/m2_sequence_3d_metadata_with_sp.csv"

# Optional: Config file override
CONFIG_FILE=""                 # Sett til path hvis du bruker config_signalp.yaml

# ============================================================================
# 6. KOMMANDOLINJE ARGUMENTER (valgfritt)
# ============================================================================
# Eksempler:
#   ./run_signalpeptide.sh --mode slow --min-sp 15
#   ./run_signalpeptide.sh --fasta data/sequences/lpmo_all_20250101_120000_raw.fasta
#   ./run_signalpeptide.sh --config config_signalp.yaml

while [[ $# -gt 0 ]]; do
    case $1 in
        --mode)
            SIGNALP6_MODE="$2"
            shift 2
            ;;
        --format)
            SIGNALP6_FORMAT="$2"
            shift 2
            ;;
        --min-sp)
            MIN_SP_LENGTH="$2"
            shift 2
            ;;
        --max-sp)
            MAX_SP_LENGTH="$2"
            shift 2
            ;;
        --fasta)
            INPUT_FASTA="$2"
            shift 2
            ;;
        --metadata)
            INPUT_METADATA="$2"
            shift 2
            ;;
        --output-fasta)
            OUTPUT_FASTA="$2"
            shift 2
            ;;
        --output-metadata)
            OUTPUT_METADATA="$2"
            shift 2
            ;;
        --config)
            CONFIG_FILE="$2"
            shift 2
            ;;
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
if [ -n "$INPUT_FASTA" ]; then
    echo "  Input FASTA: $INPUT_FASTA"
fi
if [ -n "$INPUT_METADATA" ]; then
    echo "  Input metadata: $INPUT_METADATA"
fi
echo "  Output FASTA: $OUTPUT_FASTA"
echo "  Output metadata: $OUTPUT_METADATA"
if [ -n "$CONFIG_FILE" ]; then
    echo "  Config file: $CONFIG_FILE"
fi
echo ""

# ============================================================================
# 7. KJØR SIGNALPEPTIDE MODULEN
# ============================================================================
echo "Starting Signal Peptide Prediction Pipeline..."
echo "Start time: $(date)"
echo "======================================================================"
echo ""

# Bygg kommandoen
CMD="python -m scripts.module_signalpeptide.signalp_pipeline"

# Legg til konfig-fil hvis spesifisert
if [ -n "$CONFIG_FILE" ]; then
    CMD="$CMD --config $CONFIG_FILE"
fi

# Legg til input-filer hvis spesifisert
if [ -n "$INPUT_FASTA" ]; then
    CMD="$CMD --fasta $INPUT_FASTA"
fi

if [ -n "$INPUT_METADATA" ]; then
    CMD="$CMD --metadata $INPUT_METADATA"
fi

# Legg til output-filer
CMD="$CMD --output-fasta $OUTPUT_FASTA --output-metadata $OUTPUT_METADATA"

echo "Executing: $CMD"
echo ""

# Kjør kommandoen med feilhåndtering
if eval "$CMD"; then
    EXIT_CODE=$?
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
    exit $EXIT_CODE
else
    EXIT_CODE=$?
    echo ""
    echo "======================================================================"
    echo "Signal Peptide Prediction FAILED (exit code: $EXIT_CODE)"
    echo "End time: $(date)"
    echo "======================================================================"
    exit $EXIT_CODE
fi
