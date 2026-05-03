#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import logging
import re
import traceback
from pathlib import Path
from typing import Any

from jsonschema import validate

from lpmo_pipeline.io.cif_to_pdb import convert_cif_to_pdb
from lpmo_pipeline.io.gemmi_compat import gemmi
from lpmo_pipeline.io.normalize_mmcif import NormalizeMMCIFRunner
from lpmo_pipeline.qc.hard_qc_orchestrator import HardQCInput, run_hard_qc
from lpmo_pipeline.qc.privateer_runner import get_privateer_version, prepare_privateer_input
from lpmo_pipeline.qc.qc_report import write_qc_report


DEFAULT_CIF_PATHS: tuple[Path, ...] = (
    Path(
        "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/"
        "NAG4/af3/latest/Q7SCE9_NAG4/seed-4_sample-1/Q7SCE9_NAG4_seed-4_sample-1_model.cif"
    ),
    Path(
        "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/"
        "STA6/af3/latest/Q59930_STA6/seed-2_sample-0/Q59930_STA6_seed-2_sample-0_model.cif"
    ),
    Path(
        "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/"
        "STA4/af3/latest/A0A0S2GKZ1_STA4/seed-2_sample-2/A0A0S2GKZ1_STA4_seed-2_sample-2_model.cif"
    ),
)

SCHEMA_PATH = Path(__file__).resolve().parents[2] / "schemas" / "qc_report_schema.json"
LOGGER = logging.getLogger("hard_qc_real_cifs")


def _slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run full hard QC on a small set of real AF3 CIFs.",
    )
    parser.add_argument(
        "--run-dir",
        required=True,
        help="Output directory for the prepared artifacts and qc_report.json.",
    )
    parser.add_argument(
        "--run-id",
        default="",
        help="Optional QC run identifier written into qc_report.json.",
    )
    parser.add_argument(
        "cifs",
        nargs="*",
        help="Optional explicit CIF paths. If omitted, the three default real-CIF cases are used.",
    )
    return parser.parse_args()


def _load_schema() -> dict[str, Any]:
    return json.loads(SCHEMA_PATH.read_text())


def _safe_privateer_version() -> str:
    try:
        return get_privateer_version()
    except Exception as exc:
        return f"error:{exc}"


def _prepare_case(
    cif_path: Path,
    *,
    index: int,
    run_dir: Path,
) -> tuple[dict[str, Any], HardQCInput | None]:
    pose_id = cif_path.stem
    case_label = f"{index:02d}_{_slug(pose_id)}"
    case_dir = run_dir / case_label
    case_dir.mkdir(parents=True, exist_ok=True)

    case: dict[str, Any] = {
        "index": index,
        "pose_id": pose_id,
        "cif_path": str(cif_path),
        "case_dir": str(case_dir),
        "status": "preparing",
    }

    try:
        normalize_dir = case_dir / "normalize"
        normalize_ok, normalized_path = NormalizeMMCIFRunner(cif_path, normalize_dir).run()
        case["normalize_report_path"] = str(normalize_dir / "normalize_report.json")
        if not normalize_ok or normalized_path is None:
            raise RuntimeError("Normalization failed")

        normalized_path = Path(normalized_path).resolve()
        if not normalized_path.exists():
            raise FileNotFoundError(f"normalized.cif not found: {normalized_path}")

        posebusters_dir = case_dir / "posebusters_input"
        posebusters_ok, posebusters_pdb = convert_cif_to_pdb(normalized_path, posebusters_dir)
        case["cif_to_pdb_report_path"] = str(posebusters_dir / "cif_to_pdb_report.json")
        if not posebusters_ok or posebusters_pdb is None:
            raise RuntimeError("cif_to_pdb failed")

        posebusters_pdb = Path(posebusters_pdb).resolve()
        if not posebusters_pdb.exists():
            raise FileNotFoundError(f"for_posebusters.pdb not found: {posebusters_pdb}")

        privateer_input_path = prepare_privateer_input(
            normalized_path,
            case_dir / "privateer_input.cif",
        )
        privateer_output_dir = privateer_input_path.parent / "privateer_output" / _slug(pose_id)

        structure = gemmi.read_structure(str(normalized_path))

        case.update(
            {
                "status": "prepared",
                "normalized_cif": str(normalized_path),
                "posebusters_pdb": str(posebusters_pdb),
                "privateer_input_cif": str(privateer_input_path),
                "privateer_output_dir": str(privateer_output_dir),
            }
        )
        qc_input = HardQCInput(
            pose_id=pose_id,
            mol_pred_path=posebusters_pdb,
            structure=structure,
            privateer_cif_path=privateer_input_path,
        )
        return case, qc_input
    except Exception as exc:
        case.update(
            {
                "status": "prep_error",
                "error": f"{exc.__class__.__name__}: {exc}",
                "traceback": traceback.format_exc(),
            }
        )
        return case, None


def _attach_qc_results(summary: dict[str, Any], run_dir: Path) -> None:
    verdicts = {
        verdict["pose_id"]: verdict
        for verdict in summary.get("qc_report_payload", {}).get("verdicts", [])
    }
    poses = {
        pose["pose_id"]: pose
        for pose in summary.get("qc_report_payload", {}).get("poses", [])
    }

    for case in summary["cases"]:
        pose_id = case.get("pose_id")
        if case.get("status") != "prepared" or pose_id not in verdicts:
            continue

        verdict = verdicts[pose_id]
        pose_payload = poses.get(pose_id, {})
        case["qc_verdict"] = verdict
        case["qc_pose_payload"] = pose_payload

        privateer_output_dir = Path(case["privateer_output_dir"])
        artifact_paths = {
            "validation_data": privateer_output_dir / "validation_data-privateer",
            "stdout": privateer_output_dir / "privateer_stdout.txt",
            "stderr": privateer_output_dir / "privateer_stderr.txt",
        }
        case["privateer_artifacts"] = {
            key: str(path)
            for key, path in artifact_paths.items()
            if path.exists()
        }

    summary["results_root"] = str(run_dir)


def main() -> int:
    args = _parse_args()
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

    run_dir = Path(args.run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    cif_paths = [Path(path).resolve() for path in args.cifs] if args.cifs else list(DEFAULT_CIF_PATHS)
    run_id = args.run_id or run_dir.name

    summary: dict[str, Any] = {
        "run_dir": str(run_dir),
        "run_id": run_id,
        "schema_path": str(SCHEMA_PATH),
        "requested_cifs": [str(path) for path in cif_paths],
        "n_requested": len(cif_paths),
        "privateer_version": _safe_privateer_version(),
        "cases": [],
        "n_prepared": 0,
        "n_prep_errors": 0,
        "schema_valid": False,
    }

    LOGGER.info("Preparing %d CIF inputs for hard QC", len(cif_paths))
    prepared_inputs: list[HardQCInput] = []
    for index, cif_path in enumerate(cif_paths, start=1):
        case, qc_input = _prepare_case(cif_path, index=index, run_dir=run_dir)
        summary["cases"].append(case)
        if qc_input is not None:
            prepared_inputs.append(qc_input)

    summary["n_prepared"] = len(prepared_inputs)
    summary["n_prep_errors"] = sum(1 for case in summary["cases"] if case["status"] == "prep_error")

    if prepared_inputs:
        LOGGER.info("Running hard QC on %d prepared poses", len(prepared_inputs))
        report = run_hard_qc(prepared_inputs, run_id=run_id)
        qc_report_path = run_dir / "qc_report.json"
        write_qc_report(report, qc_report_path)

        qc_report_payload = json.loads(qc_report_path.read_text())
        validate(instance=qc_report_payload, schema=_load_schema())

        summary.update(
            {
                "qc_report_path": str(qc_report_path),
                "qc_report_payload": qc_report_payload,
                "schema_valid": True,
                "qc_counts": {
                    "total": report.total,
                    "passed": report.passed,
                    "flagged": report.flagged,
                    "dropped": report.dropped,
                },
            }
        )
    else:
        summary["qc_report_error"] = "No CIFs were prepared successfully for hard QC"

    _attach_qc_results(summary, run_dir)

    summary_path = run_dir / "hard_qc_real_cifs_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    if summary["n_prepared"] != summary["n_requested"]:
        return 1
    if not summary["schema_valid"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())