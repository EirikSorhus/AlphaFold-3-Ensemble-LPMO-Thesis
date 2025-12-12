#!/bin/bash
#SBATCH --job-name=lpmo_m3             # Navn på jobben i SLURM-køen
#SBATCH --account=nn1003k              # ENDRE: Ditt prosjekt-ID
#SBATCH --time=02:00:00                # Maks kjøretid (2 timer)
#SBATCH --mem=8G                       # Minne per node
#SBATCH --cpus-per-task=4              # Antall CPU-kjerner for hmmscan
#SBATCH --ntasks=1
#SBATCH --output=logs/lpmo_m3_%j.out   # Output-fil (%j = job ID)
#SBATCH --error=logs/lpmo_m3_%j.err    # Error-fil

################################################################################
# SLURM BATCH SCRIPT FOR MODULE 3: DOMAIN ANNOTATION
################################################################################
#
# Dette scriptet kjører Module 3 som annoterer LPMO- og CBM-domener
# ved hjelp av hmmscan mot dbCAN HMM-databaser.
#
# VIKTIG: Les gjennom og tilpass følgende seksjoner før du kjører:
#   1. SLURM-parametere (linjene over)
#   2. Miljøkonfigurasjon (Python-versjon, HMMER)
#   3. HMM-database-stier
#   4. Input-filer fra Module 2
#
################################################################################

# ============================================================================
# 1. SLURM JOB INFORMASJON
# ============================================================================
echo "======================================================================"
echo "LPMO Domain Annotation Pipeline (Module 3)"
echo "======================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo ""

# ============================================================================
# 2. MILJØKONFIGURASJON
# ============================================================================
# Setter pipefail
set -euo pipefail

# Last nødvendig Python-miljø og HMMER
# ENDRE til din Python-versjon/modul:
module purge
module load NRIS/CPU
module load hpc-container-wrapper
module load Python/3.11.5-GCCcore-13.2.0

# Pek PATH til lpmo_pipe_env (conda-containerize laget dette bin/-området)
export PATH=/cluster/projects/nn1003k/eirik/conda/lpmo_pipe_env/bin:$PATH

# Aktiver Python .venv
source /cluster/projects/nn1003k/eirik/Masteroppgave/.venv/bin/activate

# Last HMMER-modul (hmmscan trengs for domene-søk)
# ENDRE til din HMMER-modul hvis annerledes:
# module load HMMER/3.3.2

echo "Python version: $(python --version)"
echo "Python location: $(which python)"
echo "hmmscan version: $(hmmscan -h | head -n1)"
echo "hmmscan location: $(which hmmscan)"
echo ""

# ============================================================================
# 3. ARBEIDSMAPPE
# ============================================================================
# Sørg for at du er i riktig mappe (project_root)
# ENDRE til din prosjektmappe:
cd /cluster/projects/nn1003k/eirik/Masteroppgave/pipe_test

# Alternativ: Bruk SLURM submit directory
# cd $SLURM_SUBMIT_DIR

# Opprett nødvendige mapper hvis de ikke finnes
mkdir -p data/domains/dbcan
mkdir -p data/metadata
mkdir -p data/sequences

echo "Project directory: $(pwd)"
echo ""

# ============================================================================
# 4. VERIFISER HMM-DATABASER
# ============================================================================
echo "Verifying HMM databases..."
echo "----------------------------------------------------------------------"

DBCAN_HMM="data/domains/dbcan/dbCAN-HMMdb-V14.hmm"
DBCAN_SUB_HMM="data/domains/dbcan/dbCAN_sub.hmm"

if [ ! -f "$DBCAN_HMM" ]; then
    echo "ERROR: dbCAN.hmm not found at $DBCAN_HMM"
    echo "Please download and hmmpress dbCAN.hmm"
    exit 1
fi

if [ ! -f "${DBCAN_HMM}.h3i" ]; then
    echo "ERROR: dbCAN.hmm index files not found (.h3i, .h3f, .h3m, .h3p)"
    echo "Run: hmmpress $DBCAN_HMM"
    exit 1
fi

echo "✓ dbCAN.hmm found and indexed"

if [ -f "$DBCAN_SUB_HMM" ]; then
    if [ -f "${DBCAN_SUB_HMM}.h3i" ]; then
        echo "✓ dbCAN_sub.hmm found and indexed"
    else
        echo "WARNING: dbCAN_sub.hmm found but not indexed"
        echo "Run: hmmpress $DBCAN_SUB_HMM"
        echo "Continuing without subfamily HMMs..."
    fi
else
    echo "WARNING: dbCAN_sub.hmm not found (optional)"
    echo "Module 3 will run without subfamily HMMs"
fi

echo ""

# ============================================================================
# 5. VERIFISER INPUT FRA MODULE 2
# ============================================================================
echo "Verifying Module 2 outputs..."
echo "----------------------------------------------------------------------"

# Finn nyeste FASTA-fil fra Module 2
LATEST_FASTA=$(ls -t data/sequences/lpmo_all_*_raw.fasta 2>/dev/null | head -n1)
if [ -z "$LATEST_FASTA" ]; then
    echo "ERROR: No Module 2 FASTA files found (lpmo_all_*_raw.fasta)"
    echo "Please run Module 2 first"
    exit 1
fi
echo "✓ Input FASTA: $LATEST_FASTA"

# Finn nyeste metadata-fil fra Module 2
LATEST_METADATA=$(ls -t data/metadata/m2_sequence_3d_metadata_*.csv 2>/dev/null | head -n1)
if [ -z "$LATEST_METADATA" ]; then
    echo "WARNING: No Module 2 metadata files found (m2_sequence_3d_metadata_*.csv)"
    echo "Module 3 will run without metadata (family info may be limited)"
fi
echo "✓ Input metadata: $LATEST_METADATA"

echo ""

# ============================================================================
# 6. TEST PYTHON ENVIRONMENT
# ============================================================================
echo "Testing Python environment..."
python -c "from scripts.module_3.config_m3 import Module3Config; print('[TEST] Module 3 imports OK')" || { echo "[ERROR] Python/imports failed"; exit 1; }
echo ""

# ============================================================================
# 7. KJØR MODULE 3
# ============================================================================
echo "======================================================================"
echo "Starting Module 3 execution..."
echo "======================================================================"

# Kjør med timeout (30 minutter - hmmscan kan ta tid)
timeout 1800 python -m scripts.module_3.run_module3 \
    --cpu $SLURM_CPUS_PER_TASK \
    --fasta "$LATEST_FASTA" \
    --metadata "$LATEST_METADATA"

PIPELINE_EXIT=$?

if [ $PIPELINE_EXIT -eq 124 ]; then
    echo ""
    echo "======================================================================"
    echo "ERROR: Module 3 timed out after 30 minutes"
    echo "This may indicate hmmscan is stuck or input file is very large"
    echo "======================================================================"
    exit 1
fi

# ============================================================================
# 8. SJEKK EXIT CODE OG OUTPUT
# ============================================================================
echo ""
echo "======================================================================"
if [ $PIPELINE_EXIT -eq 0 ]; then
    echo "✓ MODULE 3 COMPLETED SUCCESSFULLY"
    echo "======================================================================"
    echo ""
    echo "Output files:"
    echo "  Domain annotations:"
    echo "    - data/domains/m3_domains_parsed.csv"
    echo "    - data/domains/dbcan/m3_dbcan_raw.tbl"
    echo "    - data/domains/dbcan/m3_dbcan_sub_raw.tbl (if subfamilies used)"
    echo ""
    echo "  Sequence metadata:"
    echo "    - data/metadata/m3_domain_metadata.csv"
    echo "    - data/metadata/m3_run_metadata.json"
    echo ""
    echo "  FASTA variants:"
    echo "    - data/sequences/lpmo_full_length.fasta"
    echo "    - data/sequences/lpmo_catalytic_domain.fasta"
    echo ""
    echo "Next steps:"
    echo "  1. Check m3_domain_metadata.csv for domain statistics"
    echo "  2. Verify FASTA files contain expected sequences"
    echo "  3. Use lpmo_full_length.fasta and lpmo_catalytic_domain.fasta for structure prediction"
else
    echo "✗ MODULE 3 FAILED (exit code: $PIPELINE_EXIT)"
    echo "======================================================================"
    echo ""
    echo "Troubleshooting:"
    echo "  1. Check error messages in the output above"
    echo "  2. Verify hmmscan is working: hmmscan -h"
    echo "  3. Check HMM database files are properly indexed"
    echo "  4. Verify input FASTA from Module 2 is valid"
    echo "  5. Check Python dependencies (biopython, pyyaml)"
    echo ""
    echo "For more help, see: scripts/module_3/README.md (if exists)"
fi

echo "End time: $(date)"
echo "======================================================================"

exit $PIPELINE_EXIT
