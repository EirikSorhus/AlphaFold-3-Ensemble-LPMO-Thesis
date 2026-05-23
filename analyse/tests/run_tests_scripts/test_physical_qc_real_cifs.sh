#!/bin/bash
#SBATCH --job-name=run_tests_physical_qc_real_cifs
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_physical_qc_real_cifs_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "[ERROR] This script must be submitted with sbatch, not run directly on the login node."
    echo "[INFO] Example: sbatch tests/run_tests_scripts/test_physical_qc_real_cifs.sh"
    exit 1
fi

export PATH="/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin:$PATH"

# Avoid writing pip/http cache files into repo-adjacent .cache paths.
export PIP_NO_CACHE_DIR=1

project_root="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse"
cd "$project_root"
export PYTHONPATH="$project_root/src:${PYTHONPATH:-}"

out_dir="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/tests_results"
mkdir -p "$out_dir"

job_suffix="${SLURM_JOB_ID:-manual}"
run_dir="$out_dir/physical_qc_real_cifs_${job_suffix}"
mkdir -p "$run_dir"

# Route generic user-cache writes to run-specific output, not workspace root.
export XDG_CACHE_HOME="$run_dir/cache"
mkdir -p "$XDG_CACHE_HOME"

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

echo "[INFO] Running targeted QC unit tests"
pytest \
    tests/test_qc_gates.py::TestCuHisGate \
    tests/test_qc_gates.py::TestCuHisGateAtomSelection \
    tests/test_hard_qc_orchestrator.py::TestHardQCOrchestrator::test_pre_qc_failure_skips_pb_and_geometry \
    -q

echo "[INFO] Running physical QC probe on real CIFs"
export RUN_DIR="$run_dir"
export EXTRA_CIFS="${*:-}"

python - <<'PY'
import json
import os
import shlex
import traceback
from pathlib import Path

import gemmi

from lpmo_pipeline.io.normalize_mmcif import NormalizeMMCIFRunner
from lpmo_pipeline.qc.active_site_proximity import (
    ACTIVE_SITE_PROXIMITY_MAX_A,
    check_active_site_proximity,
)
from lpmo_pipeline.qc.custom_geometry_checks import (
    CU_HIS_MAX,
    CU_HIS_MIN,
    check_geometry,
)
from lpmo_pipeline.qc.qc_report import compute_verdict


RUN_DIR = Path(os.environ["RUN_DIR"])

cif_paths = [
    Path("/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/NAG4/af3/latest/Q7SCE9_NAG4/seed-4_sample-1/Q7SCE9_NAG4_seed-4_sample-1_model.cif"),
    Path("/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/STA6/af3/latest/Q59930_STA6/seed-2_sample-0/Q59930_STA6_seed-2_sample-0_model.cif"),
    Path("/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/STA4/af3/latest/A0A0S2GKZ1_STA4/seed-2_sample-2/A0A0S2GKZ1_STA4_seed-2_sample-2_model.cif"),
]

extra_cifs = os.environ.get("EXTRA_CIFS", "").strip()
if extra_cifs:
    cif_paths.extend(Path(part) for part in shlex.split(extra_cifs))


def finite_or_none(value):
    if isinstance(value, float) and (value == float("inf") or value == float("-inf")):
        return None
    return value


results = []
for index, cif_path in enumerate(cif_paths, start=1):
    case_dir = RUN_DIR / f"cif_{index}"
    case_dir.mkdir(parents=True, exist_ok=True)

    pose_id = cif_path.stem
    case = {
        "index": index,
        "cif_path": str(cif_path),
        "pose_id": pose_id,
        "status": "ok",
        "error": None,
    }

    try:
        normalize_out = case_dir / "normalize"
        normalize_out.mkdir(parents=True, exist_ok=True)

        ok_norm, normalized_path = NormalizeMMCIFRunner(cif_path, normalize_out).run()
        case["normalize"] = {
            "ok": bool(ok_norm),
            "normalized_path": str(normalized_path) if normalized_path else None,
            "report_path": str(normalize_out / "normalize_report.json"),
        }

        if not ok_norm or not normalized_path:
            raise RuntimeError("Normalization failed; cannot run QC")

        normalized_path = Path(normalized_path)
        if not normalized_path.exists():
            raise FileNotFoundError(f"normalized.cif not found: {normalized_path}")

        structure = gemmi.read_structure(str(normalized_path))

        prox = check_active_site_proximity(structure=structure, pose_id=pose_id)
        geom = check_geometry(structure=structure, pose_id=pose_id)

        verdict = compute_verdict(
            pose_id=pose_id,
            pb_result=None,
            priv_result=None,
            geom_result=geom,
            proximity_result=prox,
        )

        cu_his_measurements = [
            {
                "his_chain": m.his_chain,
                "his_resnum": m.his_resnum,
                "his_atom": m.his_atom,
                "distance_angstrom": m.distance_angstrom,
                "in_range": m.in_range,
            }
            for m in geom.cu_his_measurements
        ]

        case["thresholds"] = {
            "active_site_proximity_max_a": ACTIVE_SITE_PROXIMITY_MAX_A,
            "cu_his_min_a": CU_HIS_MIN,
            "cu_his_max_a": CU_HIS_MAX,
            "cu_his_gate_selected_atoms": [
                "His1:N",
                "His1:ND1",
                "His!=1:(ND1|NE2) nearest within cutoff",
            ],
        }
        case["active_site_proximity"] = {
            "passed": prox.passed,
            "min_cu_ligand_distance": finite_or_none(prox.min_cu_ligand_distance),
            "min_cu_c1": finite_or_none(prox.min_cu_c1),
            "min_cu_c4": finite_or_none(prox.min_cu_c4),
            "nearest_ligand_atom": prox.nearest_ligand_atom,
            "nearest_c1_atom": prox.nearest_c1_atom,
            "nearest_c4_atom": prox.nearest_c4_atom,
            "failure_reasons": list(prox.failure_reasons),
            "warnings": list(prox.warnings),
        }
        case["geometry"] = {
            "passed": geom.passed,
            "cu_found": geom.cu_found,
            "cu_position": list(geom.cu_position),
            "cu_his_all_in_range": geom.cu_his_all_in_range,
            "cu_his_measurements": cu_his_measurements,
            "cu_his_n_measurements": len(cu_his_measurements),
            "cu_his_expected_selected_count": 3,
            "min_cu_c1": finite_or_none(geom.min_cu_c1),
            "min_cu_c4": finite_or_none(geom.min_cu_c4),
            "his_brace_angle": finite_or_none(geom.his_brace_angle),
            "failure_reasons": list(geom.failure_reasons),
        }
        case["qc_verdict"] = {
            "status": verdict.status,
            "drop_reasons": list(verdict.drop_reasons),
            "warnings": list(verdict.warnings),
            "metrics": verdict.metrics,
            "cu_geometry": verdict.cu_geometry,
        }

        case_report_path = case_dir / "physical_qc_case_report.json"
        case["case_report_path"] = str(case_report_path)
        case_report_path.write_text(json.dumps(case, indent=2))

    except Exception as exc:  # pragma: no cover
        case["status"] = "error"
        case["error"] = f"{exc.__class__.__name__}: {exc}"
        case["traceback"] = traceback.format_exc()

    results.append(case)

summary = {
    "run_dir": str(RUN_DIR),
    "n_inputs": len(cif_paths),
    "n_ok": sum(1 for r in results if r.get("status") == "ok"),
    "n_error": sum(1 for r in results if r.get("status") == "error"),
    "results": results,
}

summary_path = RUN_DIR / "physical_qc_probe_summary.json"
summary_path.write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))

if summary["n_ok"] == 0:
    raise SystemExit("[ERROR] No CIF completed physical QC probe successfully")
PY

echo "[INFO] Physical QC probe completed"
echo "[INFO] Results directory: $run_dir"
