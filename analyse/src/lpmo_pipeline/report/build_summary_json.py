# src/lpmo_pipeline/report/build_summary_json.py
"""
Responsibility: Aggregate all pipeline outputs into summary.json.
Input:  QC reports, cluster results, geometry stats, tuning results, crystal anchoring
Output: summary.json conforming to schemas/summary_json_schema.json

STEP 13a in masterplan.
"""
from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def build_summary(
    run_id: str,
    mode: str,
    del_variant: str,
    qc_reports: list[dict[str, Any]],
    cluster_results: list[dict[str, Any]],
    geometry_stats: list[dict[str, Any]],
    crystal_reports: list[dict[str, Any]] | None = None,
    tuning_summary: dict[str, Any] | None = None,
    activity_features: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Build the unified summary.json.

    Args:
        run_id: Pipeline run identifier.
        mode: "tune" | "production".
        del_variant: "a" | "b".
        qc_reports: List of per-run QC report dicts.
        cluster_results: List of per-protein×ligand cluster result dicts.
        geometry_stats: List of per-pose geometry metric dicts.
        crystal_reports: Optional crystal anchoring results.
        tuning_summary: Optional tuning summary (if mode=tune).
        activity_features: Optional activity feature dicts.

    Returns:
        summary dict ready for JSON serialization.
    """
    # --- Dataset stats ---
    n_proteins = len(set(
        r.get("protein_id", "") for g in geometry_stats for r in [g]
    ))
    n_ligands = len(set(
        r.get("ligand_id", "") for g in geometry_stats for r in [g]
    ))
    n_models = len(set(
        r.get("model", "") for g in geometry_stats for r in [g]
    ))
    n_total_poses = len(geometry_stats)

    # --- QC stats ---
    total_qc = sum(q.get("total", 0) for q in qc_reports)
    total_passed = sum(q.get("passed", 0) for q in qc_reports)
    total_flagged = sum(q.get("flagged", 0) for q in qc_reports)
    total_dropped = sum(q.get("dropped", 0) for q in qc_reports)

    # --- Cluster stats ---
    all_n_clusters = [c.get("n_clusters", 0) for c in cluster_results]
    all_outlier_rates = [c.get("outlier_rate", 0.0) for c in cluster_results]

    import numpy as np
    mean_clusters = float(np.mean(all_n_clusters)) if all_n_clusters else 0.0
    mean_outlier = float(np.mean(all_outlier_rates)) if all_outlier_rates else 0.0

    # --- Geometry stats ---
    cu_c1_vals = [
        g.get("min_cu_c1", float("inf"))
        for g in geometry_stats
        if g.get("min_cu_c1", float("inf")) < float("inf")
    ]
    cu_c4_vals = [
        g.get("min_cu_c4", float("inf"))
        for g in geometry_stats
        if g.get("min_cu_c4", float("inf")) < float("inf")
    ]

    summary: dict[str, Any] = {
        "run_id": run_id,
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "mode": mode,
        "del_variant": del_variant,
        "dataset_stats": {
            "n_proteins": n_proteins,
            "n_ligands": n_ligands,
            "n_models": n_models,
            "n_total_poses": n_total_poses,
        },
        "qc_stats": {
            "total_poses_qc": total_qc,
            "passed": total_passed,
            "flagged": total_flagged,
            "dropped": total_dropped,
            "pass_rate": total_passed / total_qc if total_qc > 0 else 0.0,
        },
        "cluster_stats": {
            "mean_n_clusters": mean_clusters,
            "mean_outlier_rate": mean_outlier,
            "n_protein_ligand_combinations": len(cluster_results),
        },
        "geometry_stats": {
            "cu_c1_range": [float(min(cu_c1_vals)), float(max(cu_c1_vals))] if cu_c1_vals else None,
            "cu_c4_range": [float(min(cu_c4_vals)), float(max(cu_c4_vals))] if cu_c4_vals else None,
            "n_measurements": len(cu_c1_vals),
        },
    }

    if tuning_summary:
        summary["tuning_results"] = tuning_summary

    if crystal_reports:
        tanimotos = [
            c.get("best_tanimoto", 0.0) for c in crystal_reports
        ]
        summary["crystal_stats"] = {
            "n_comparisons": len(crystal_reports),
            "mean_best_tanimoto": float(np.mean(tanimotos)) if tanimotos else 0.0,
        }

    if activity_features:
        summary["activity_stats"] = {
            "n_proteins_annotated": len(activity_features),
        }

    logger.info("Built summary.json: %d proteins, %d poses, QC pass=%.1f%%",
                n_proteins, n_total_poses,
                summary["qc_stats"]["pass_rate"] * 100)
    return summary


def write_summary_json(summary: dict[str, Any], output_path: Path) -> None:
    """Write summary to JSON file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(summary, f, indent=2)
    logger.info("Wrote summary.json to %s", output_path)
