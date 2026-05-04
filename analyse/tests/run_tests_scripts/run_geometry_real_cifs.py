#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import logging
import re
import traceback
from pathlib import Path
from typing import Any

from lpmo_pipeline.analysis.mdanalysis_metrics import (
    compute_pose_metrics_from_structure,
    write_geometry_debug_pdb,
    write_geometry_metrics,
    write_pose_geometry_tsv,
)
from lpmo_pipeline.io.gemmi_compat import gemmi
from lpmo_pipeline.io.normalize_mmcif import NormalizeMMCIFRunner


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

LOGGER = logging.getLogger("geometry_real_cifs")


def _slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def _parse_ids(pose_id: str) -> tuple[str, str, str]:
    parts = pose_id.split("_")
    if len(parts) >= 2:
        return parts[0], parts[1], "af3"
    return "", "", "af3"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run downstream geometry analysis on a small set of real AF3 CIFs.",
    )
    parser.add_argument(
        "--run-dir",
        required=True,
        help="Output directory for prepared artifacts, pose_geometry.tsv, and debug PDBs.",
    )
    parser.add_argument(
        "--run-id",
        default="",
        help="Optional run identifier written into the JSON summary.",
    )
    parser.add_argument(
        "cifs",
        nargs="*",
        help="Optional explicit CIF paths. Defaults to the three standard real-CIF cases.",
    )
    return parser.parse_args()


def _prepare_case(
    cif_path: Path,
    *,
    index: int,
    run_dir: Path,
) -> tuple[dict[str, Any], Any | None, Any | None]:
    pose_id = cif_path.stem
    protein_id, ligand_id, model = _parse_ids(pose_id)
    case_label = f"{index:02d}_{_slug(pose_id)}"
    case_dir = run_dir / case_label
    case_dir.mkdir(parents=True, exist_ok=True)

    case: dict[str, Any] = {
        "index": index,
        "pose_id": pose_id,
        "protein_id": protein_id,
        "ligand_id": ligand_id,
        "model": model,
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

        structure = gemmi.read_structure(str(normalized_path))
        metrics = compute_pose_metrics_from_structure(
            structure=structure,
            pose_id=pose_id,
            protein_id=protein_id,
            ligand_id=ligand_id,
            model=model,
        )

        metrics_path = case_dir / "geometry_metrics.json"
        write_geometry_metrics(metrics, metrics_path)

        debug_pdb_path = case_dir / "geometry_debug.pdb"
        write_geometry_debug_pdb(structure, metrics, debug_pdb_path)

        case.update(
            {
                "status": "prepared",
                "normalized_cif": str(normalized_path),
                "geometry_metrics_json": str(metrics_path),
                "geometry_debug_pdb": str(debug_pdb_path),
                "geometry_row": metrics.to_row(),
                "debug_coordinates": {
                    "original_cu_coordinates": metrics.original_cu_coordinates,
                    "repositioned_cu_coordinates": metrics.repositioned_cu_coordinates,
                    "virtual_oxyl_coordinates": metrics.virtual_oxyl_coordinates,
                    "virtual_h_c1_coordinates": metrics.virtual_h_c1_coordinates,
                    "virtual_h_c4_coordinates": metrics.virtual_h_c4_coordinates,
                },
            }
        )
        return case, metrics, structure
    except Exception as exc:
        case.update(
            {
                "status": "prep_error",
                "error": f"{exc.__class__.__name__}: {exc}",
                "traceback": traceback.format_exc(),
            }
        )
        return case, None, None


def main() -> int:
    args = _parse_args()
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

    run_dir = Path(args.run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    run_id = args.run_id or run_dir.name

    cif_paths = [Path(path).resolve() for path in args.cifs] if args.cifs else list(DEFAULT_CIF_PATHS)
    summary: dict[str, Any] = {
        "run_dir": str(run_dir),
        "run_id": run_id,
        "requested_cifs": [str(path) for path in cif_paths],
        "n_requested": len(cif_paths),
        "cases": [],
    }

    metrics_rows = []
    LOGGER.info("Preparing %d CIF inputs for downstream geometry test", len(cif_paths))
    for index, cif_path in enumerate(cif_paths, start=1):
        case, metrics, _ = _prepare_case(cif_path, index=index, run_dir=run_dir)
        summary["cases"].append(case)
        if metrics is not None:
            metrics_rows.append(metrics)

    summary["n_prepared"] = len(metrics_rows)
    summary["n_prep_errors"] = sum(1 for case in summary["cases"] if case["status"] == "prep_error")

    if metrics_rows:
        pose_geometry_path = run_dir / "pose_geometry.tsv"
        write_pose_geometry_tsv(metrics_rows, pose_geometry_path)
        summary["pose_geometry_tsv"] = str(pose_geometry_path)

    summary_path = run_dir / "geometry_real_cifs_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    if summary["n_prepared"] == 0:
        return 1
    if summary["n_prepared"] != summary["n_requested"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())