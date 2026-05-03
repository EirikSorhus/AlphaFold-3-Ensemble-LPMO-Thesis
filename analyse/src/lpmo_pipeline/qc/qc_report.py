# src/lpmo_pipeline/qc/qc_report.py
"""
Responsibility: Aggregate QC results (pre-QC proximity + PoseBusters + Privateer + geometry)
                into a unified per-pose QC report.
Input:  ActiveSiteProximityResult, PoseBustersBatchResult, PrivateerResult, GeometryResult (per pose)
Output: qc_report.json conforming to schemas/qc_report_schema.json

Implements failure policy:
  - Hard-fail → pose marked "dropped" with reason
  - Soft-flag → pose marked "flagged" with warning list
  - Pass      → pose marked "passed"
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from lpmo_pipeline.qc.active_site_proximity import ActiveSiteProximityResult
from lpmo_pipeline.qc.posebusters_runner import PoseBustersSingleResult
from lpmo_pipeline.qc.privateer_runner import PrivateerResult
from lpmo_pipeline.qc.custom_geometry_checks import GeometryResult

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Per-pose QC verdict
# ---------------------------------------------------------------------------
@dataclass
class PoseQCVerdict:
    """Unified QC verdict for a single pose."""

    pose_id: str
    status: str  # "passed" | "flagged" | "dropped"
    proximity_passed: bool = True
    posebusters_passed: bool = True
    privateer_passed: bool = True
    geometry_passed: bool = True
    drop_reasons: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    metrics: dict[str, Any] = field(default_factory=dict)
    posebusters: dict[str, Any] = field(default_factory=dict)
    privateer: dict[str, Any] = field(default_factory=dict)
    cu_geometry: dict[str, Any] = field(default_factory=dict)


def compute_verdict(
    pose_id: str,
    pb_result: PoseBustersSingleResult | None,
    priv_result: PrivateerResult | None,
    geom_result: GeometryResult | None,
    proximity_result: ActiveSiteProximityResult | None = None,
) -> PoseQCVerdict:
    """Compute unified QC verdict for a single pose.

    Failure policy:
            - Pre-QC active-site proximity fail -> drop
      - Any critical PoseBusters error → drop
      - Privateer recognition < 100% or anomer fail → drop
        - Cu–His outside hard gate window → drop
        - Cu–His outside preferred QC window → flag
      - Soft PoseBusters warnings → flag (keep)
      - Crystal similarity low → flag (keep)
    """
    verdict = PoseQCVerdict(pose_id=pose_id, status="passed")

    # --- Pre-QC active-site proximity ---
    if proximity_result is not None:
        verdict.proximity_passed = proximity_result.passed
        verdict.metrics["min_cu_ligand_distance"] = proximity_result.min_cu_ligand_distance
        verdict.metrics["min_cu_c1"] = proximity_result.min_cu_c1
        verdict.metrics["min_cu_c4"] = proximity_result.min_cu_c4
        if proximity_result.nearest_ligand_atom:
            verdict.metrics["nearest_ligand_atom"] = proximity_result.nearest_ligand_atom
        if not proximity_result.passed:
            verdict.drop_reasons.extend(
                [f"active_site_proximity:{r}" for r in proximity_result.failure_reasons]
            )
        if proximity_result.warnings:
            verdict.warnings.extend(
                [f"active_site_proximity:{w}" for w in proximity_result.warnings]
            )

    # --- PoseBusters ---
    if pb_result is not None:
        verdict.posebusters_passed = pb_result.passed
        error_counts: dict[str, int] = {}
        for err in pb_result.critical_errors:
            error_counts[err] = error_counts.get(err, 0) + 1
        verdict.posebusters = {
            "passed": pb_result.passed,
            "errors": error_counts,
            "critical": len(pb_result.critical_errors) > 0,
            "warnings": list(pb_result.warnings),
        }
        if not pb_result.passed:
            verdict.drop_reasons.extend(
                [f"posebusters_critical:{e}" for e in pb_result.critical_errors]
            )
        if pb_result.warnings:
            verdict.warnings.extend(
                [f"posebusters_soft:{w}" for w in pb_result.warnings]
            )

    # --- Privateer ---
    if priv_result is not None:
        verdict.privateer_passed = priv_result.all_pass
        privateer_errors: list[str] = []
        if priv_result.runner_error:
            privateer_errors.append(f"runner_error:{priv_result.runner_error}")
        if priv_result.recognition_rate < 1.0:
            privateer_errors.append("recognition_rate_below_100")
        if priv_result.anomer_pass < priv_result.total_sugars:
            privateer_errors.append("anomer_failure")
        if priv_result.ring_pucker_pass < priv_result.total_sugars:
            privateer_errors.append("ring_pucker_failure")
        if priv_result.linkage_pass < priv_result.total_sugars:
            privateer_errors.append("linkage_failure")

        verdict.privateer = {
            "recognized_sugars": priv_result.recognized,
            "total_sugars": priv_result.total_sugars,
            "recognition_rate": priv_result.recognition_rate,
            "anomer_ok": priv_result.anomer_pass == priv_result.total_sugars,
            "ring_pucker_ok": priv_result.ring_pucker_pass == priv_result.total_sugars,
            "linkage_ok": priv_result.linkage_pass == priv_result.total_sugars,
            "errors": privateer_errors,
        }
        if not priv_result.all_pass:
            if priv_result.runner_error:
                verdict.drop_reasons.append(f"privateer_runner_error:{priv_result.runner_error}")
            elif priv_result.recognition_rate < 1.0:
                unrecognized = [
                    r.resname for r in priv_result.residues if not r.sugar_recognized
                ]
                verdict.drop_reasons.append(
                    f"privateer_unrecognized:{unrecognized}"
                )
            bad_anomers = [
                f"{r.chain}:{r.resname}{r.resnum}"
                for r in priv_result.residues if not r.anomer_ok
            ]
            if bad_anomers:
                verdict.drop_reasons.append(f"privateer_anomer:{bad_anomers}")

    # --- Geometry ---
    if geom_result is not None:
        verdict.geometry_passed = geom_result.passed
        if not geom_result.passed:
            verdict.drop_reasons.extend(geom_result.failure_reasons)
        if geom_result.warnings:
            verdict.warnings.extend(
                [f"geometry_soft:{warning}" for warning in geom_result.warnings]
            )
        cu_his_distances = [m.distance_angstrom for m in geom_result.cu_his_measurements]
        verdict.cu_geometry = {
            "cu_his_distances_a": cu_his_distances,
            "all_in_range": geom_result.cu_his_all_in_range,
            "all_in_soft_range": geom_result.cu_his_all_in_soft_range,
            "cu_c1_dist_a": None if geom_result.min_cu_c1 == float("inf") else geom_result.min_cu_c1,
            "cu_c4_dist_a": None if geom_result.min_cu_c4 == float("inf") else geom_result.min_cu_c4,
        }
        # Preserve pre-QC Cu-substrate metrics when present; geometry keeps its
        # own explicit distance keys so both measurement surfaces remain usable.
        verdict.metrics["cu_found"] = geom_result.cu_found
        verdict.metrics.setdefault("min_cu_c1", geom_result.min_cu_c1)
        verdict.metrics.setdefault("min_cu_c4", geom_result.min_cu_c4)
        verdict.metrics["geometry_min_cu_c1"] = geom_result.min_cu_c1
        verdict.metrics["geometry_min_cu_c4"] = geom_result.min_cu_c4
        verdict.metrics["cu_his_all_in_range"] = geom_result.cu_his_all_in_range
        verdict.metrics["cu_his_all_in_soft_range"] = geom_result.cu_his_all_in_soft_range
        if geom_result.his_brace_angle is not None:
            verdict.metrics["his_brace_angle"] = geom_result.his_brace_angle

    # --- Final status ---
    if verdict.drop_reasons:
        verdict.status = "dropped"
    elif verdict.warnings:
        verdict.status = "flagged"
    else:
        verdict.status = "passed"

    logger.info("QC verdict for %s: %s", pose_id, verdict.status)
    return verdict


# ---------------------------------------------------------------------------
# Batch report
# ---------------------------------------------------------------------------
@dataclass
class QCReport:
    """Full QC report for a run (all poses)."""

    run_id: str = ""
    verdicts: list[PoseQCVerdict] = field(default_factory=list)
    total: int = 0
    passed: int = 0
    flagged: int = 0
    dropped: int = 0


def build_qc_report(
    run_id: str,
    verdicts: list[PoseQCVerdict],
) -> QCReport:
    """Build aggregated QC report."""
    return QCReport(
        run_id=run_id,
        verdicts=verdicts,
        total=len(verdicts),
        passed=sum(1 for v in verdicts if v.status == "passed"),
        flagged=sum(1 for v in verdicts if v.status == "flagged"),
        dropped=sum(1 for v in verdicts if v.status == "dropped"),
    )


def write_qc_report(report: QCReport, output_path: Path) -> None:
    """Write QC report to JSON conforming to qc_report_schema.json."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).isoformat()
    data = {
        "timestamp": timestamp,
        "run_id": report.run_id,
        "poses": [
            {
                "pose_id": v.pose_id,
                "overall_status": (
                    "pass" if v.status == "passed"
                    else "soft_flag" if v.status == "flagged"
                    else "hard_fail"
                ),
                "posebusters": v.posebusters or {
                    "passed": v.posebusters_passed,
                    "errors": {},
                    "critical": not v.posebusters_passed,
                    "warnings": [],
                },
                "privateer": v.privateer or {
                    "recognized_sugars": 0,
                    "total_sugars": 0,
                    "recognition_rate": 0.0,
                    "anomer_ok": True,
                    "ring_pucker_ok": True,
                    "linkage_ok": True,
                    "errors": [],
                },
                "cu_geometry": v.cu_geometry or {
                    "cu_his_distances_a": [],
                    "all_in_range": v.geometry_passed,
                    "cu_c1_dist_a": None,
                    "cu_c4_dist_a": None,
                },
            }
            for v in report.verdicts
        ],
        "total": report.total,
        "passed": report.passed,
        "flagged": report.flagged,
        "dropped": report.dropped,
        "verdicts": [
            {
                "pose_id": v.pose_id,
                "status": v.status,
                "proximity_passed": v.proximity_passed,
                "posebusters_passed": v.posebusters_passed,
                "privateer_passed": v.privateer_passed,
                "geometry_passed": v.geometry_passed,
                "drop_reasons": v.drop_reasons,
                "warnings": v.warnings,
                "metrics": v.metrics,
            }
            for v in report.verdicts
        ],
    }
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    logger.info("Wrote QC report (%d poses) to %s", report.total, output_path)
