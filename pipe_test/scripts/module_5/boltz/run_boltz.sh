#!/usr/bin/env bash
#SBATCH --job-name=boltz2_run
#SBATCH --account=nn1003k
#SBATCH --partition=accel
#SBATCH --gpus=1
#SBATCH --time=00:20:00
#SBATCH --mem-per-gpu=80G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail
module load NRIS/GPU

IMG=/cluster/projects/nn1003k/prog/boltz/boltz2_alt.sif
WEIGHTS_HOST=/cluster/projects/nn1003k/prog/boltz/weights

INPUT="${SUBMITDIR}/input/boltz_test_pkl_NAG6.yaml"
OUT_DIR="${SUBMITDIR}/result"

[[ -f "$IMG" ]]   || { echo "ERROR: image not found: $IMG"; exit 2; }
[[ -f "$INPUT" ]] || { echo "ERROR: input not found: $INPUT"; exit 2; }

mkdir -p "$OUT_DIR"

apptainer exec --nv \
  --bind "${SUBMITDIR}:/work" \
  --bind "${WEIGHTS_HOST}:/weights" \
  --env HF_HOME=/weights/huggingface \
  "$IMG" boltz predict "/work/input/$(basename "$INPUT")" \
    --out_dir /work/result \
    --cache /weights/boltz \
    --override

echo "Boltz run completed. Results in: $OUT_DIR"
