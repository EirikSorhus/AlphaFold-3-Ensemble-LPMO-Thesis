#!/bin/bash
#SBATCH --job-name=lpmo_fetch          # Navn på jobben i SLURM-køen
#SBATCH --account=nn1003k              # ENDRE: Ditt prosjekt-ID
#SBATCH --time=00:30:00                # Maks kjøretid (30 minutter)
#SBATCH --mem=1G                       # Minne per node
#SBATCH --cpus-per-task=1              # Antall CPU-kjerner
#SBATCH --ntasks=1
#SBATCH --output=logs/lpmo_m2_%j.out      # Output-fil (%j = job ID)
#SBATCH --error=logs/lpmo_m2_%j.err       # Error-fil

################################################################################
# SLURM BATCH SCRIPT FOR MODULE 2: LPMO SEQUENCE FETCHING
################################################################################
#
# Dette scriptet kjører den forbedrede versjonen av Module 2 som henter
# LPMO-proteinsekvenser fra UniProt, CAZy og NCBI.
#
# VIKTIG: Les gjennom og tilpass følgende seksjoner før du kjører:
#   1. SLURM-parametere (linjene over)
#   2. Miljøkonfigurasjon (Python-versjon)
#   3. NCBI email (PÅKREVD)
#   4. Pipeline-alternativer (valgfritt)
#
################################################################################

# ============================================================================
# 1. SLURM JOB INFORMASJON
# ============================================================================
echo "======================================================================"
echo "LPMO Sequence Fetching Pipeline"
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

# Last nødvendig Python-miljø
# ENDRE til din Python-versjon/modul:
module purge
module load NRIS/CPU
module load hpc-container-wrapper
module load Python/3.11.5-GCCcore-13.2.0

# Pek PATH til lpmo_pipe_env (conda-containerize laget dette bin/-området)
export PATH=/cluster/projects/nn1003k/eirik/conda/lpmo_pipe_env/bin:$PATH

# Aktiver Python .venv
source /cluster/projects/nn1003k/eirik/Masteroppgave/.venv/bin/activate

# Alternativ: Bruk conda environment
# module load Anaconda3
# source activate myenv

echo "Python version: $(python --version)"
echo "Python location: $(which python)"
echo "hmmscan location: $(which hmmscan)"
echo ""

# ============================================================================
# 3. NCBI EMAIL (PÅKREVD)
# ============================================================================
# NCBI krever en gyldig e-postadresse for API-tilgang.
# 
# VELG ÉN AV DISSE METODENE:

# Metode A: Sett direkte i script (ENDRE EMAIL!)
export NCBI_EMAIL="eirik.sorhus@nmbu.no"

# Metode B: Hent fra en fil (mer sikkert for delte scripts)
# export NCBI_EMAIL=$(cat ~/.ncbi_email)

# Metode C: Bruk CLI-argument (se nedenfor under "Pipeline Options")
# (da trenger du ikke export her)

echo "NCBI Email: $NCBI_EMAIL"
echo ""

# ============================================================================
# 4. ARBEIDSMAPPE
# ============================================================================
# Sørg for at du er i riktig mappe (project_root)
# ENDRE til din prosjektmappe:
cd /cluster/projects/nn1003k/eirik/Masteroppgave/pipe_test

# Alternativ: Bruk SLURM submit directory
# cd $SLURM_SUBMIT_DIR

# Opprett data-mappe hvis den ikke finnes
mkdir -p data

echo "Project directory: $(pwd)"
echo ""

# ============================================================================
# 5. PIPELINE OPTIONS
# ============================================================================
# Her kan du velge hvordan du vil kjøre pipelinen.
# Les kommentarene for å forstå hvert alternativ.

# Test at Python fungerer og imports er OK
echo "Testing Python environment..."
python -c "from scripts.module_2.config_improved import check_dependencies; print('[TEST] Imports OK'); check_dependencies()" 2>&1 || { echo "[ERROR] Python/imports failed"; exit 1; }
echo ""

# ----------------------------------------------------------------------------
# ALTERNATIV 1: STANDARD KJØRING (anbefalt - alle LPMO-familier)
# ----------------------------------------------------------------------------
# Kjører alle tre steg: UniProt → CAZy/NCBI → Merge
# Bruker standard familier: AA9, AA10, AA11, AA13
# Bruker NCBI_EMAIL fra miljøvariabel satt over



# Kjører koden
echo "======================================================================" 
echo "Starting pipeline execution with 10-minute timeout per step..."
echo "======================================================================"

# Kjør med timeout
timeout 600 python -m scripts.module_2.run_module2_improved --families AA9
PIPELINE_EXIT=$?

if [ $PIPELINE_EXIT -eq 124 ]; then
    echo ""
    echo "======================================================================" 
    echo "ERROR: Pipeline timed out after 10 minutes"
    echo "This usually means UniProt or CAZy request is hanging"
    echo "======================================================================"
    exit 1
fi



# ----------------------------------------------------------------------------
# ALTERNATIV 1B: VELG SPESIFIKKE FAMILIER
# ----------------------------------------------------------------------------
# Spesifiser hvilke CAZy-familier du vil hente med --families flagget
# Du kan velge én eller flere familier

# Eksempel: Kun AA9 og AA10
# python -m scripts.module_2.run_module2_improved --families AA9 AA10

# Eksempel: Kun AA13
# python -m scripts.module_2.run_module2_improved --families AA13

# Eksempel: Alle standard LPMO-familier (eksplisitt)
# python -m scripts.module_2.run_module2_improved --families AA9 AA10 AA11 AA13

# ----------------------------------------------------------------------------
# ALTERNATIV 2: MED CLI EMAIL (istedenfor miljøvariabel)
# ----------------------------------------------------------------------------
# Hvis du ikke vil sette NCBI_EMAIL som miljøvariabel, kan du sende den
# direkte som argument:

# python -m scripts.module_2.run_module2_improved --email din.epost@institusjon.no

# ----------------------------------------------------------------------------
# ALTERNATIV 3: HOPP OVER UNIPROT (hvis du allerede har filen)
# ----------------------------------------------------------------------------
# Hvis du allerede har kjørt UniProt-steget og har filen
# 'data/uniprot_LPMO_raw.fasta', kan du hoppe over dette steget:

# python -m scripts.module_2.run_module2_improved --skip-uniprot

# ----------------------------------------------------------------------------
# ALTERNATIV 4: HOPP OVER CAZY/NCBI (hvis du kun vil ha UniProt-data)
# ----------------------------------------------------------------------------
# Hvis du kun vil hente fra UniProt og hoppe over CAZy/NCBI:

# python -m scripts.module_2.run_module2_improved --skip-cazy

# ----------------------------------------------------------------------------
# ALTERNATIV 5: HOPP OVER MERGE (kun hent data, ikke slå sammen)
# ----------------------------------------------------------------------------
# Hvis du vil hente data men ikke merge ennå:

# python -m scripts.module_2.run_module2_improved --skip-merge

# ----------------------------------------------------------------------------
# ALTERNATIV 6: KOMBINER FLERE FLAGG
# ----------------------------------------------------------------------------
# Du kan kombinere flere alternativer:

# Eksempel: Hopp over UniProt og bruk spesifikk email
# python -m scripts.module_2.run_module2_improved \
#     --skip-uniprot \
#     --email din.epost@institusjon.no

# Eksempel: Kun kjør merge-steget (krever at de to andre filene finnes)
# python -m scripts.module_2.run_module2_improved \
#     --skip-uniprot \
#     --skip-cazy

# ----------------------------------------------------------------------------
# ALTERNATIV 7: VELG DEDUPLISERINGSSTRATEGI
# ----------------------------------------------------------------------------
# Standard er å deduplisere basert på sekvens.
# Du kan velge andre strategier:

# Kun sekvens (standard) - fjerner identiske sekvenser uavhengig av ID
# python -m scripts.module_2.run_module2_improved --deduplicate-by sequence

# Kun ID - fjerner entries med samme ID
# python -m scripts.module_2.run_module2_improved --deduplicate-by id

# Både ID og sekvens må matche
# python -m scripts.module_2.run_module2_improved --deduplicate-by both

# Ingen deduplisering
# python -m scripts.module_2.run_module2_improved --deduplicate-by none

# ============================================================================
# 6. SJEKK EXIT CODE
# ============================================================================
# Pipeline returnerer 0 ved suksess, 1 ved feil
EXIT_CODE=$?

echo ""
echo "======================================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ PIPELINE COMPLETED SUCCESSFULLY"
    echo "======================================================================"
    echo ""
    echo "Output files:"
    echo "  - data/uniprot_LPMO_raw.fasta"
    echo "  - data/cazy_LPMO_raw.fasta"
    echo "  - data/LPMO_all_raw.fasta  (merged)"
    echo ""
    echo "Next steps:"
    echo "  1. Verify output files exist and contain sequences"
    echo "  2. Check sequence counts in the output above"
    echo "  3. Proceed to next module in your pipeline"
else
    echo "✗ PIPELINE FAILED (exit code: $EXIT_CODE)"
    echo "======================================================================"
    echo ""
    echo "Troubleshooting:"
    echo "  1. Check error messages in the output above"
    echo "  2. Verify NCBI_EMAIL is set correctly"
    echo "  3. Check network connectivity to UniProt/CAZy/NCBI"
    echo "  4. Verify all Python packages are installed"
    echo ""
    echo "For more help, see: scripts/module_2/README_IMPROVED.md"
fi

echo "End time: $(date)"
echo "======================================================================"

exit $EXIT_CODE
