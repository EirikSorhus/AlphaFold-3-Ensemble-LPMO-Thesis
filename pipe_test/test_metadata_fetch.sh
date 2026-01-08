#!/bin/bash
#SBATCH --job-name=metadata_fetch             # Navn på jobben i SLURM-køen
#SBATCH --account=nn1003k              # ENDRE: Ditt prosjekt-ID
#SBATCH --time=00:10:00                # Maks kjøretid (10 minutter)
#SBATCH --mem=500M                       # Minne per node
#SBATCH --cpus-per-task=1              # Antall CPU-kjerner
#SBATCH --ntasks=1
#SBATCH --output=logs/metadata_fetch_%j.out   # Output-fil (%j = job ID)
#SBATCH --error=logs/metadata_fetch_%j.err    # Error-fil
set -euo pipefail


module purge
module load NRIS/CPU
module load Python/3.11.5-GCCcore-13.2.0

echo "== Test UniProt metadata pipeline =="

# Aktiver Python .venv
source /cluster/projects/nn1003k/eirik/Masteroppgave/.venv/bin/activate

# Sørg for at output-mapper finnes
mkdir -p data/metadata
mkdir -p logs

# Lag en liten testliste med UniProt IDs (LPMOer)
# !!! Dette er ikke uniprot ID dette er NCBI og vil ikke fungere !!!
cat > lpmo_test_ids.txt << EOF
WQF89147.1
WQF76824.1
WDK15276.1
CAQ16217.1
WDK18612.1
WYZ42323.1
UQC83797.1
EOF

echo "Input IDs written to lpmo_test_ids.txt"

# Kjør Python-koden (antatt at den ligger i scripts/module_2)
python uniprot_metadata.py \
  --input lpmo_test_ids.txt \
  --output data/metadata/metadata_test_lpmo.tsv

echo "== Ferdig =="
echo "Output: data/metadata/metadata_test_lpmo.tsv"
