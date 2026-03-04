# src/lpmo_pipeline/report/build_metrics_csv.py
"""
Responsibility: Build flat metrics.csv from all per-pose data.
Input:  Per-pose geometry metrics, IFP data, cluster labels, QC verdicts
Output: metrics.csv conforming to schemas/metrics_csv_schema.json

Columns: run_id, protein_id, ligand_id, model, seed, pose_id, cluster_id,
         qc_status, Cu_C1, Cu_C4, his_brace_angle, ifp_similarity_crystal,
         n_contacts, activity_class, substrate_type, dp, cbm_present, ...

NOTE: Dropped poses are EXCLUDED from metrics.csv (logged in failures only).
      Flagged poses are INCLUDED with qc_status="flagged".
"""
from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

# Column order (must match metrics_csv_schema.json)
METRICS_COLUMNS = [
    "run_id",
    "protein_id",
    "ligand_id",
    "model",
    "seed",
    "pose_id",
    "cluster_id",
    "qc_status",
    "posebusters_passed",
    "privateer_passed",
    "geometry_passed",
    "min_cu_c1",
    "min_cu_c4",
    "his_brace_angle_deg",
    "core_rmsd_vs_reference",
    "pocket_rmsd_vs_crystal",
    "ifp_similarity_crystal",
    "n_ifp_contacts",
    "activity_class",
    "substrate_type",
    "dp",
    "del_variant",
    "cbm_present",
    "cbm_ligand_min_dist",
    "n_warnings",
]


def build_metrics_csv(
    pose_records: list[dict[str, Any]],
    output_path: Path,
    exclude_dropped: bool = True,
) -> int:
    """Build flat metrics CSV from per-pose records.

    Args:
        pose_records: List of per-pose dicts. Each must have at minimum:
            {pose_id, protein_id, ligand_id, model, qc_status, ...}
        output_path: Path to write CSV.
        exclude_dropped: If True, exclude QC-dropped poses (default).

    Returns:
        Number of rows written.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    filtered = pose_records
    if exclude_dropped:
        filtered = [r for r in pose_records if r.get("qc_status") != "dropped"]
        n_dropped = len(pose_records) - len(filtered)
        if n_dropped > 0:
            logger.info("Excluded %d dropped poses from metrics CSV", n_dropped)

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=METRICS_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for record in filtered:
            writer.writerow(record)

    logger.info("Wrote metrics.csv (%d rows) to %s", len(filtered), output_path)
    return len(filtered)


def merge_pose_record(
    pose_id: str,
    run_id: str,
    protein_id: str,
    ligand_id: str,
    model: str,
    seed: int,
    qc_verdict: dict[str, Any],
    geometry: dict[str, Any],
    cluster_label: int = -1,
    crystal_sim: float | None = None,
    n_ifp_contacts: int = 0,
    activity_class: str = "",
    substrate_type: str = "",
    dp: int = 0,
    del_variant: str = "a",
    cbm_present: bool = False,
    cbm_ligand_min_dist: float | None = None,
) -> dict[str, Any]:
    """Merge data sources into a single flat record for metrics CSV.

    Assembles one row from QC, geometry, cluster, and annotation sources.
    """
    return {
        "run_id": run_id,
        "protein_id": protein_id,
        "ligand_id": ligand_id,
        "model": model,
        "seed": seed,
        "pose_id": pose_id,
        "cluster_id": cluster_label,
        "qc_status": qc_verdict.get("status", "unknown"),
        "posebusters_passed": qc_verdict.get("posebusters_passed", True),
        "privateer_passed": qc_verdict.get("privateer_passed", True),
        "geometry_passed": qc_verdict.get("geometry_passed", True),
        "min_cu_c1": geometry.get("min_cu_c1"),
        "min_cu_c4": geometry.get("min_cu_c4"),
        "his_brace_angle_deg": geometry.get("his_brace_angle_deg"),
        "core_rmsd_vs_reference": geometry.get("core_rmsd_vs_reference"),
        "pocket_rmsd_vs_crystal": geometry.get("pocket_rmsd_vs_crystal"),
        "ifp_similarity_crystal": crystal_sim,
        "n_ifp_contacts": n_ifp_contacts,
        "activity_class": activity_class,
        "substrate_type": substrate_type,
        "dp": dp,
        "del_variant": del_variant,
        "cbm_present": cbm_present,
        "cbm_ligand_min_dist": cbm_ligand_min_dist,
        "n_warnings": len(qc_verdict.get("warnings", [])),
    }
