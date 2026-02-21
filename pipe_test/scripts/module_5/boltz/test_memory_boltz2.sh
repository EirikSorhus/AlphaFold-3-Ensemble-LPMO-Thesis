#!/usr/bin/env bash
#SBATCH --job-name=boltz2_profile
#SBATCH --account=nn1003k
#SBATCH --partition=accel
#SBATCH --gpus=1
#SBATCH --time=00:20:00
#SBATCH --mem-per-gpu=30G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail
module load NRIS/GPU

# ---- Paths (host) ----
IMG=/cluster/projects/nn1003k/prog/boltz/boltz2_alt.sif
WEIGHTS_HOST=/cluster/projects/nn1003k/prog/boltz/weights

INPUT="${SUBMITDIR}/input/boltz_eksempel_amylase.yaml"
OUT_DIR="${SUBMITDIR}/result"

[[ -f "$IMG" ]]   || { echo "ERROR: image not found: $IMG"; exit 2; }
[[ -f "$INPUT" ]] || { echo "ERROR: input not found: $INPUT"; exit 2; }

mkdir -p "$OUT_DIR"

# ---- Logs ----
GPU_LOG="${OUT_DIR}/gpu_usage_${SLURM_JOB_ID}.log"
RAM_LOG="${OUT_DIR}/ram_job_${SLURM_JOB_ID}.log"

# ---- cgroup path (Olivia / Slurm). SLURM_UID is not always set; use id -u. ----
USER_UID="$(id -u)"
CGROUP_DIR="/sys/fs/cgroup/slurm/uid_${USER_UID}/job_${SLURM_JOB_ID}"

# Fallback: try to locate the job cgroup if path differs on the node
if [[ ! -d "$CGROUP_DIR" ]]; then
  CGROUP_DIR="$(find /sys/fs/cgroup -maxdepth 6 -type d -name "job_${SLURM_JOB_ID}" 2>/dev/null | head -n 1 || true)"
fi

cleanup() {
  kill "${GPU_MONITOR_PID:-}" "${RAM_MONITOR_PID:-}" 2>/dev/null || true
}
trap cleanup EXIT

echo "JobID: ${SLURM_JOB_ID}"
echo "Node(s): ${SLURM_JOB_NODELIST}"
echo "Output dir: ${OUT_DIR}"
echo "GPU log: ${GPU_LOG}"
echo "RAM log: ${RAM_LOG}"
echo "cgroup dir: ${CGROUP_DIR:-NA}"
echo

# ---- Start GPU monitor (5s) ----
( while sleep 5; do
    echo "$(date '+%F %T'), $(nvidia-smi --query-gpu=memory.used,utilization.gpu --format=csv,noheader,nounits)" \
      >> "$GPU_LOG"
  done ) &
GPU_MONITOR_PID=$!

# ---- Start job-level RAM monitor from cgroup (5s) ----
( while sleep 5; do
    if [[ -n "${CGROUP_DIR:-}" && -r "${CGROUP_DIR}/memory.current" ]]; then
      CUR="$(cat "${CGROUP_DIR}/memory.current")"
      if [[ -r "${CGROUP_DIR}/memory.peak" ]]; then
        PEAK="$(cat "${CGROUP_DIR}/memory.peak")"
      else
        PEAK="NA"
      fi
      echo "$(date '+%F %T'), current=${CUR}, peak=${PEAK}" >> "$RAM_LOG"
    else
      echo "$(date '+%F %T'), cgroup_memory_unavailable" >> "$RAM_LOG"
    fi
  done ) &
RAM_MONITOR_PID=$!

# ---- Run Boltz ----
apptainer exec --nv \
  --bind "${SUBMITDIR}:/work" \
  --bind "${WEIGHTS_HOST}:/weights" \
  --env HF_HOME=/weights/huggingface \
  "$IMG" boltz predict "/work/input/$(basename "$INPUT")" \
    --out_dir /work/result \
    --cache /weights/boltz \
    --override

# ---- Summary (best-effort) ----
echo
echo "=== Profiling summary (best-effort) ==="

if [[ -s "$GPU_LOG" ]]; then
  PEAK_GPU_MIB="$(awk -F',' 'NF>=2 {gsub(/ /,"",$2); print $2}' "$GPU_LOG" | sort -n | tail -1 || true)"
  echo "Peak GPU memory.used: ${PEAK_GPU_MIB:-NA} MiB"
else
  echo "Peak GPU memory.used: NA (no gpu log)"
fi

PEAK_RAM_BYTES="NA"
if [[ -n "${CGROUP_DIR:-}" && -r "${CGROUP_DIR}/memory.peak" ]]; then
  PEAK_RAM_BYTES="$(cat "${CGROUP_DIR}/memory.peak" || true)"
elif [[ -s "$RAM_LOG" ]]; then
  PEAK_RAM_BYTES="$(awk -F'peak=' 'NF==2 {print $2}' "$RAM_LOG" | grep -v NA | sort -n | tail -1 || true)"
fi

if [[ "$PEAK_RAM_BYTES" =~ ^[0-9]+$ ]]; then
  export PEAK_RAM_BYTES
  PEAK_RAM_GIB="$(python3 - <<'PY'
import os
b = int(os.environ["PEAK_RAM_BYTES"])
print(f"{b/1024/1024/1024:.2f}")
PY
)"
  echo "Peak job RAM (cgroup): ${PEAK_RAM_GIB} GiB (${PEAK_RAM_BYTES} bytes)"
else
  echo "Peak job RAM (cgroup): NA"
fi

echo "Logs written to:"
echo "  $GPU_LOG"
echo "  $RAM_LOG"
