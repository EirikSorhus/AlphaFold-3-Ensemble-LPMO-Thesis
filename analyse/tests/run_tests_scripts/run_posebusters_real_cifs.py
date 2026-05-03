#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import logging
import re
import traceback
from pathlib import Path
from typing import Any

from lpmo_pipeline.io.cif_to_pdb import convert_cif_to_pdb
from lpmo_pipeline.io.normalize_mmcif import NormalizeMMCIFRunner
from lpmo_pipeline.qc.posebusters_runner import run_posebusters_single


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

LOGGER = logging.getLogger("posebusters_real_cifs")


def _slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare real AF3 CIFs and run PoseBusters only.",
    )
    parser.add_argument(
        "--run-dir",
        required=True,
        help="Output directory for prepared artifacts and the PoseBusters summary.",
    )
    parser.add_argument(
        "--run-id",
        default="",
        help="Optional run identifier written into the summary.",
    )
    parser.add_argument(
        "cifs",
        nargs="*",
        help="Optional explicit CIF paths. If omitted, the three default real-CIF cases are used.",
    )
    return parser.parse_args()


def _summarize_posebusters_chains(pdb_path: Path) -> dict[str, list[str]]:
    chains: dict[str, set[str]] = {}
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        chain_id = line[21].strip() or "_"
        resname = line[17:20].strip()
        chains.setdefault(chain_id, set()).add(resname)
    return {
        chain_id: sorted(resnames)
        for chain_id, resnames in sorted(chains.items())
    }


def _run_case(cif_path: Path, *, index: int, run_dir: Path) -> dict[str, Any]:
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
        posebusters_dir = case_dir / "posebusters_input"
        posebusters_ok, posebusters_pdb = convert_cif_to_pdb(normalized_path, posebusters_dir)
        case["cif_to_pdb_report_path"] = str(posebusters_dir / "cif_to_pdb_report.json")
        if not posebusters_ok or posebusters_pdb is None:
            raise RuntimeError("cif_to_pdb failed")

        posebusters_pdb = Path(posebusters_pdb).resolve()
        chain_summary = _summarize_posebusters_chains(posebusters_pdb)
        auto_split_mode = "A" in chain_summary and any(chain != "A" for chain in chain_summary)
        pb_result = run_posebusters_single(posebusters_pdb, pose_id)

        case.update(
            {
                "status": "tested",
                "normalized_cif": str(normalized_path),
                "posebusters_pdb": str(posebusters_pdb),
                "posebusters_chain_summary": chain_summary,
                "posebusters_input_contract": (
                    "dock_auto_split_from_combined_pdb" if auto_split_mode else "direct_input"
                ),
                "posebusters_result": {
                    "passed": pb_result.passed,
                    "critical_errors": pb_result.critical_errors,
                    "warnings": pb_result.warnings,
                    "details": pb_result.details,
                },
            }
        )
        return case
    except Exception as exc:
        case.update(
            {
                "status": "error",
                "error": f"{exc.__class__.__name__}: {exc}",
                "traceback": traceback.format_exc(),
            }
        )
        return case


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
        "requested_cifs": [str(path) for path in cif_paths],
        "n_requested": len(cif_paths),
        "cases": [],
    }

    LOGGER.info("Preparing %d CIF inputs for PoseBusters-only validation", len(cif_paths))
    for index, cif_path in enumerate(cif_paths, start=1):
        summary["cases"].append(_run_case(cif_path, index=index, run_dir=run_dir))

    summary["n_tested"] = sum(1 for case in summary["cases"] if case.get("status") == "tested")
    summary["n_errors"] = sum(1 for case in summary["cases"] if case.get("status") == "error")
    summary["n_passed"] = sum(
        1
        for case in summary["cases"]
        if case.get("posebusters_result", {}).get("passed") is True
    )
    summary["n_failed"] = sum(
        1
        for case in summary["cases"]
        if case.get("posebusters_result", {}).get("passed") is False
    )

    summary_path = run_dir / "posebusters_real_cifs_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    return 1 if summary["n_errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())