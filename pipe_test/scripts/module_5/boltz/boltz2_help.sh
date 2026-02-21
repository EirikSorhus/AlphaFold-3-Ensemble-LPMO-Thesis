#!/usr/bin/env bash
#SBATCH --job-name=boltz2_help
#SBATCH --account=nn1003k
#SBATCH --partition=accel
#SBATCH --gpus=1
#SBATCH --time=00:10:00
#SBATCH --mem-per-gpu=10G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail
module load NRIS/GPU

IMG=/cluster/projects/nn1003k/prog/boltz/boltz2_alt.sif
WEIGHTS_HOST=/cluster/projects/nn1003k/prog/boltz/weights

OUTBASE=$SUBMITDIR
OUT_DIR=$OUTBASE/boltz2_help_${SLURM_JOB_ID}

[[ -f "$IMG" ]] || { echo "ERROR: image not found: $IMG"; exit 2; }

mkdir -p "$OUT_DIR"

echo "=== Writing Boltz 2 help documentation to text files ==="

apptainer exec --nv \
  --bind "$SUBMITDIR:/work" \
  --bind "$WEIGHTS_HOST:/weights" \
  --env HF_HOME=/weights/huggingface \
  "$IMG" bash -lc "
    boltz --help > /work/boltz2_help_${SLURM_JOB_ID}/boltz_help.txt
    boltz predict --help > /work/boltz2_help_${SLURM_JOB_ID}/boltz_predict_help.txt
  "

echo "Done."
echo "Help files written to:"
echo "  ${OUT_DIR}/boltz_help.txt"
echo "  ${OUT_DIR}/boltz_predict_help.txt"
