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

echo "[INFO] Probing Privateer SIF metadata"
privateer_sif="/cluster/projects/nn1003k/prog/privateer/privateer.sif"
if [[ ! -f "$privateer_sif" ]]; then
    echo "[ERROR] Privateer SIF not found: $privateer_sif"
    exit 1
fi

echo "$privateer_sif" > "$run_dir/privateer_path.txt"
apptainer run --cleanenv "$privateer_sif" -list > "$run_dir/privateer_list.txt" 2>&1 || true
apptainer run --cleanenv "$privateer_sif" help > "$run_dir/privateer_help.txt" 2>&1 || true

cat > "$run_dir/privateer_command_templates.txt" <<'EOF'
apptainer run --cleanenv /cluster/projects/nn1003k/prog/privateer/privateer.sif -pdbin {input}
apptainer run --cleanenv /cluster/projects/nn1003k/prog/privateer/privateer.sif -pdbin {input} -mode ccp4i2
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
    "apptainer run --cleanenv /cluster/projects/nn1003k/prog/privateer/privateer.sif -pdbin {input}",
    "apptainer run --cleanenv /cluster/projects/nn1003k/prog/privateer/privateer.sif -pdbin {input} -mode ccp4i2",
]

subprocess_env = {
    "PATH": "/usr/bin:/bin:/usr/local/bin",
    "HOME": os.environ.get("HOME", "/tmp"),
    "USER": os.environ.get("USER", "unknown"),
    "LOGNAME": os.environ.get("LOGNAME", os.environ.get("USER", "unknown")),
    "LANG": os.environ.get("LANG", "C.UTF-8"),
    "LC_ALL": os.environ.get("LC_ALL", "C.UTF-8"),
    "TERM": os.environ.get("TERM", "xterm"),
}


def output_summary(text: str) -> dict[str, object]:
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        return {
            "parsed_as_json": False,
            "error": f"{exc.__class__.__name__}: {exc}",
        }

    if isinstance(payload, dict):
        return {
            "parsed_as_json": True,
            "top_level_type": "dict",
            "top_level_keys": sorted(payload.keys()),
        }
    if isinstance(payload, list):
        return {
            "parsed_as_json": True,
            "top_level_type": "list",
            "length": len(payload),
            "first_item_type": type(payload[0]).__name__ if payload else None,
        }
    return {
        "parsed_as_json": True,
        "top_level_type": type(payload).__name__,
    }


def classify_failure(returncode: int, stderr_text: str) -> dict[str, object]:
    stderr_lower = stderr_text.lower()
    category = "unknown"
    if "mmdbfile: read_file error" in stderr_lower:
        category = "mmdb_read_error"
    elif "exec: privateer: not found" in stderr_lower:
        category = "privateer_not_found_in_container"
    elif "unrecognised" in stderr_lower:
        category = "unsupported_cli_flag"
    elif "jsondecodeerror" in stderr_lower:
        category = "json_parse_error"

    return {
        "category": category,
        "returncode": returncode,
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
            cwd=str(cif_dir),
            env=subprocess_env,
        )

        stdout_path = cif_dir / f"attempt_{attempt_index}_stdout.txt"
        stderr_path = cif_dir / f"attempt_{attempt_index}_stderr.txt"
        stdout_path.write_text(proc.stdout)
        stderr_path.write_text(proc.stderr)

        generated_files = sorted(
            str(path.relative_to(cif_dir))
            for path in cif_dir.rglob("*")
            if path.is_file()
            and path.name not in {stdout_path.name, stderr_path.name}
        )
        summary = output_summary(proc.stdout)
        failure_summary = classify_failure(proc.returncode, proc.stderr)
        attempt = {
            "template": template,
            "command": cmd,
            "returncode": proc.returncode,
            "stdout_path": str(stdout_path),
            "stderr_path": str(stderr_path),
            "stdout_nonempty": bool(proc.stdout.strip()),
            "stderr_nonempty": bool(proc.stderr.strip()),
            "stdout_summary": summary,
            "failure_summary": failure_summary,
            "generated_files": generated_files,
        }
        cif_result["attempts"].append(attempt)

    results.append(cif_result)

summary = {
    "privateer_path": (RUN_DIR / "privateer_path.txt").read_text().strip(),
    "list_file": str(RUN_DIR / "privateer_list.txt"),
    "help_file": str(RUN_DIR / "privateer_help.txt"),
    "results": results,
}

any_success = any(
    attempt.get("returncode") == 0
    for result in results
    for attempt in result["attempts"]
)

failure_counts: dict[str, int] = {}
for result in results:
    for attempt in result["attempts"]:
        category = str(attempt.get("failure_summary", {}).get("category", "unknown"))
        failure_counts[category] = failure_counts.get(category, 0) + 1
summary["failure_counts"] = failure_counts

if not any_success:
    summary["warning"] = (
        "No attempted Privateer command completed successfully. "
        "See per-attempt stderr/stdout artifacts for container and input-format diagnostics."
    )

summary_path = RUN_DIR / "privateer_probe_summary.json"
summary_path.write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))

if not any_success:
    print("[WARN] None of the attempted Privateer command templates completed successfully.")
PY

echo "[INFO] Privateer real-CIF probe completed"
echo "[INFO] Results directory: $run_dir"