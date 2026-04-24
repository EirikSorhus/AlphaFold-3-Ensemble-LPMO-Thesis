#!/bin/bash
#SBATCH --job-name=run_tests_privateer_real_cifs
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=00:15:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_privateer_real_cifs_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "[ERROR] This script must be submitted with sbatch, not run directly on the login node."
    echo "[INFO] Example: sbatch tests/run_tests_scripts/test_privateer_real_cifs.sh"
    exit 1
fi

export PATH="/cluster/work/projects/nn1003k/eirik/conda/analyse_env/bin:$PATH"

project_root="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse"
cd "$project_root"
export PYTHONPATH="$project_root/src:${PYTHONPATH:-}"

out_dir="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/tests_results"
mkdir -p "$out_dir"

job_suffix="${SLURM_JOB_ID:-manual}"
run_dir="$out_dir/privateer_real_cifs_${job_suffix}"
mkdir -p "$run_dir"

cif_paths=(
  "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/NAG4/af3/latest/Q7SCE9_NAG4/seed-4_sample-1/Q7SCE9_NAG4_seed-4_sample-1_model.cif"
  "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/STA6/af3/latest/Q59930_STA6/seed-2_sample-0/Q59930_STA6_seed-2_sample-0_model.cif"
  "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/STA4/af3/latest/A0A0S2GKZ1_STA4/seed-2_sample-2/A0A0S2GKZ1_STA4_seed-2_sample-2_model.cif"
)

if [[ "$#" -gt 0 ]]; then
    for extra_cif in "$@"; do
        cif_paths+=("$extra_cif")
    done
fi

for cif_path in "${cif_paths[@]}"; do
  if [[ ! -f "$cif_path" ]]; then
    echo "[ERROR] CIF does not exist: $cif_path"
    exit 1
  fi
done

echo "[INFO] Running targeted Privateer tests"
pytest tests/test_privateer_runner.py tests/test_qc_gates.py::TestPrivateerGate -q

echo "[INFO] Probing Privateer CLI metadata"
if ! command -v privateer >/dev/null 2>&1; then
  echo "[ERROR] privateer is not available in analyse_env PATH"
  exit 1
fi

privateer_bin="$(command -v privateer)"
echo "$privateer_bin" > "$run_dir/privateer_path.txt"
privateer --version > "$run_dir/privateer_version.txt" 2>&1 || true
privateer --help > "$run_dir/privateer_help.txt" 2>&1 || true

cat > "$run_dir/privateer_command_templates.txt" <<'EOF'
privateer --json {input}
privateer -json {input}
privateer validate --input {input} --json
privateer validate --json {input}
privateer -pdbin {input} -mode glycan_validation -json
EOF

export RUN_DIR="$run_dir"
export EXTRA_CIFS="${*:-}"

python - <<'PY'
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path


RUN_DIR = Path(os.environ["RUN_DIR"])

cif_paths = [
    Path("/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/NAG4/af3/latest/Q7SCE9_NAG4/seed-4_sample-1/Q7SCE9_NAG4_seed-4_sample-1_model.cif"),
    Path("/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/STA6/af3/latest/Q59930_STA6/seed-2_sample-0/Q59930_STA6_seed-2_sample-0_model.cif"),
    Path("/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/STA4/af3/latest/A0A0S2GKZ1_STA4/seed-2_sample-2/A0A0S2GKZ1_STA4_seed-2_sample-2_model.cif"),
]

extra_cifs = os.environ.get("EXTRA_CIFS", "").strip()
if extra_cifs:
    cif_paths.extend(Path(part) for part in shlex.split(extra_cifs))

templates = [
    "privateer --json {input}",
    "privateer -json {input}",
    "privateer validate --input {input} --json",
    "privateer validate --json {input}",
    "privateer -pdbin {input} -mode glycan_validation -json",
]


def json_summary(text: str) -> dict[str, object]:
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        return {
            "parsed": False,
            "error": f"{exc.__class__.__name__}: {exc}",
        }

    if isinstance(payload, dict):
        return {
            "parsed": True,
            "top_level_type": "dict",
            "top_level_keys": sorted(payload.keys()),
        }
    if isinstance(payload, list):
        return {
            "parsed": True,
            "top_level_type": "list",
            "length": len(payload),
            "first_item_type": type(payload[0]).__name__ if payload else None,
        }
    return {
        "parsed": True,
        "top_level_type": type(payload).__name__,
    }


results: list[dict[str, object]] = []

for index, cif_path in enumerate(cif_paths, start=1):
    cif_dir = RUN_DIR / f"cif_{index}"
    cif_dir.mkdir(parents=True, exist_ok=True)

    cif_result = {
        "index": index,
        "cif_path": str(cif_path),
        "attempts": [],
    }

    for attempt_index, template in enumerate(templates, start=1):
        rendered = template.format(input=str(cif_path))
        cmd = shlex.split(rendered)

        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
        )

        stdout_path = cif_dir / f"attempt_{attempt_index}_stdout.txt"
        stderr_path = cif_dir / f"attempt_{attempt_index}_stderr.txt"
        stdout_path.write_text(proc.stdout)
        stderr_path.write_text(proc.stderr)

        summary = json_summary(proc.stdout)
        attempt = {
            "template": template,
            "command": cmd,
            "returncode": proc.returncode,
            "stdout_path": str(stdout_path),
            "stderr_path": str(stderr_path),
            "stdout_nonempty": bool(proc.stdout.strip()),
            "stderr_nonempty": bool(proc.stderr.strip()),
            "json_summary": summary,
        }
        cif_result["attempts"].append(attempt)

        if summary.get("parsed"):
            cif_result["first_json_success_template"] = template
            break

    results.append(cif_result)

summary = {
    "privateer_path": (RUN_DIR / "privateer_path.txt").read_text().strip(),
    "version_file": str(RUN_DIR / "privateer_version.txt"),
    "help_file": str(RUN_DIR / "privateer_help.txt"),
    "results": results,
}

summary_path = RUN_DIR / "privateer_probe_summary.json"
summary_path.write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))

if not any(
    attempt.get("json_summary", {}).get("parsed")
    for result in results
    for attempt in result["attempts"]
):
    print("[ERROR] None of the attempted Privateer command templates produced parseable JSON.")
    sys.exit(1)
PY

echo "[INFO] Privateer real-CIF probe completed"
echo "[INFO] Results directory: $run_dir"