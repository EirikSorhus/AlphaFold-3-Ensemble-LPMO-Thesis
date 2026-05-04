#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from typing import Any

from lpmo_pipeline.qc.posebusters_runner import POSEBUSTERS_SIF_CANDIDATES
from lpmo_pipeline.utils.manifest import ToolVersionFetcher


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run lightweight smoke checks for the analysis toolchain.",
    )
    parser.add_argument(
        "--run-dir",
        required=True,
        help="Output directory for the smoke summary JSON.",
    )
    parser.add_argument(
        "--run-id",
        default="",
        help="Optional run identifier written into the summary output.",
    )
    return parser.parse_args()


def _selected_posebusters_sif() -> Path | None:
    return next((path for path in POSEBUSTERS_SIF_CANDIDATES if path.exists()), None)


def main() -> int:
    args = _parse_args()
    run_dir = Path(args.run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    run_id = args.run_id or run_dir.name

    tool_versions = ToolVersionFetcher.get_versions()
    dependency_statuses = ToolVersionFetcher.get_analysis_dependency_statuses()

    apptainer_path = shutil.which("apptainer")
    selected_posebusters_sif = _selected_posebusters_sif()

    failures: list[str] = []
    if any(status != "ok" for status in dependency_statuses.values()):
        failures.append("One or more Python dependency probes failed")
    if apptainer_path is None:
        failures.append("apptainer is not on PATH")
    if selected_posebusters_sif is None:
        failures.append("No PoseBusters SIF candidate was found")
    if tool_versions.get("privateer") in {None, "", "not_found"}:
        failures.append("Privateer version probe failed")

    summary: dict[str, Any] = {
        "run_dir": str(run_dir),
        "run_id": run_id,
        "tool_versions": tool_versions,
        "dependency_statuses": dependency_statuses,
        "external_runtime": {
            "apptainer_path": apptainer_path,
            "posebusters_python_api_available": tool_versions.get("posebusters") != "not_installed",
            "posebusters_sif_candidates": [
                {"path": str(path), "exists": path.exists()}
                for path in POSEBUSTERS_SIF_CANDIDATES
            ],
            "selected_posebusters_sif": str(selected_posebusters_sif) if selected_posebusters_sif else None,
            "posebusters_runtime_ready": apptainer_path is not None and selected_posebusters_sif is not None,
            "privateer_runtime_ready": tool_versions.get("privateer") not in {None, "", "not_found"},
        },
        "success": not failures,
        "failures": failures,
    }

    summary_path = run_dir / "analysis_tools_smoke_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())