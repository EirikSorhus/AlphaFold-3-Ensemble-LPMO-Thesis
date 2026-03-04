# src/lpmo_pipeline/qc/qc_report.py
"""
Responsibility: Aggregate QC results (PoseBusters + Privateer + geometry)
                into a unified per-pose QC report.
Input:  PoseBustersBatchResult, PrivateerResult, GeometryResult (per pose)
Output: qc_report.json conforming to schemas/qc_report_schema.json

Implements failure policy:
  - Hard-fail → pose marked "dropped" with reason
  - Soft-flag → pose marked "flagged" with warning list
  - Pass      → pose marked "passed"
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

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
    posebusters_passed: bool = True
    privateer_passed: bool = True
    geometry_passed: bool = True
    drop_reasons: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    metrics: dict[str, Any] = field(default_factory=dict)


def compute_verdict(
    pose_id: str,
    pb_result: PoseBustersSingleResult | None,
    priv_result: PrivateerResult | None,
    geom_result: GeometryResult | None,
) -> PoseQCVerdict:
    """Compute unified QC verdict for a single pose.

    Failure policy:
      - Any critical PoseBusters error → drop
      - Privateer recognition < 100% or anomer fail → drop
      - Cu–His outside 1.9–2.6 Å → drop
      - Soft PoseBusters warnings → flag (keep)
      - Crystal similarity low → flag (keep)
    """
    verdict = PoseQCVerdict(pose_id=pose_id)

    # --- PoseBusters ---
    if pb_result is not None:
        verdict.posebusters_passed = pb_result.passed
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
        if not priv_result.all_pass:
            if priv_result.recognition_rate < 1.0:
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
        # Store metrics regardless
        verdict.metrics["cu_found"] = geom_result.cu_found
        verdict.metrics["min_cu_c1"] = geom_result.min_cu_c1
        verdict.metrics["min_cu_c4"] = geom_result.min_cu_c4
        verdict.metrics["cu_his_all_in_range"] = geom_result.cu_his_all_in_range
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
    data = {
        "run_id": report.run_id,
        "total": report.total,
        "passed": report.passed,
        "flagged": report.flagged,
        "dropped": report.dropped,
        "verdicts": [
            {
                "pose_id": v.pose_id,
                "status": v.status,
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
