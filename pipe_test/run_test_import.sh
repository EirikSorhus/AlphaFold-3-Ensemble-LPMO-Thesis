#!/bin/bash
#SBATCH --job-name=test_import              # Navn på jobben i SLURM-køen
#SBATCH --account=nn1003k                    # ENDRE: Ditt prosjekt-ID
#SBATCH --time=01:00:00                      # Maks kjøretid (1 time)
#SBATCH --mem=1G                            # Minne per node
#SBATCH --cpus-per-task=1                    # Antall CPU-kjerner
#SBATCH --ntasks=1
#SBATCH --output=logs/test_import_%j.out    # Output-fil (%j = job ID)
#SBATCH --error=logs/test_import_%j.err     # Error-fil

################################################################################
# SLURM BATCH SCRIPT FOR MODULE SIGNALPEPTIDE: test import
################################################################################

# ============================================================================
# 1. SLURM JOB INFORMASJON
# ============================================================================
echo "======================================================================"
echo "Test Import Script"
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
module load Python/3.11.5-GCCcore-13.2.0


# Aktiver Python .venv
source /cluster/projects/nn1003k/eirik/Masteroppgave/.venv/bin/activate

python get_cazy_sequences.py -f AA9 AA13

