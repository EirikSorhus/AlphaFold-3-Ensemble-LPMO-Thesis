#!/bin/bash
#SBATCH --job-name=cazy_uniprot_metadata      # Jobbnavn i SLURM-køen
#SBATCH --account=nn1003k                     # Prosjekt-ID
#SBATCH --time=00:10:00                       # Maks kjøretid
#SBATCH --mem=500M                            # Minne
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=logs/cazy_uniprot_%j.out     # STDOUT
#SBATCH --error=logs/cazy_uniprot_%j.err      # STDERR

###############################################################################
# CAZy → UniProt metadata pipeline (module_2)
###############################################################################

set -euo pipefail

echo "======================================================================"
echo "CAZy → UniProt metadata pipeline"
echo "======================================================================"
echo "Job ID:        $SLURM_JOB_ID"
echo "Node:          $SLURM_NODELIST"
echo "Start time:    $(date)"
echo "Working dir:   $(pwd)"
echo "======================================================================"
echo ""

# ---------------------------------------------------------------------------
# 1) Miljø
# ---------------------------------------------------------------------------
module purge
module load NRIS/CPU
module load Python/3.11.5-GCCcore-13.2.0

# Aktiver prosjektets venv
source /cluster/projects/nn1003k/eirik/Masteroppgave/.venv/bin/activate

echo "Python executable: $(which python)"
python --version
echo ""

# ---------------------------------------------------------------------------
# 2) Forbered mapper
# ---------------------------------------------------------------------------
mkdir -p logs
mkdir -p data/run
mkdir -p data/metadata

# ---------------------------------------------------------------------------
# 3) Kjør pipeline
# ---------------------------------------------------------------------------
echo "Running cazy_to_uniprot_driver.py ..."
echo ""

python scripts/module_2/cazy_to_uniprot_driver.py \
  --families AA9 AA13 \
  --format tsv \
  --contact-email eirik.sorhus@nmbu.no

# ---------------------------------------------------------------------------
# 4) Ferdig
# ---------------------------------------------------------------------------
echo ""
echo "======================================================================"
echo "Pipeline finished"
echo "End time: $(date)"
echo "======================================================================"
