# src/lpmo_pipeline/placer/rank_ensemble.py
"""
Responsibility: Rank PLACER ensemble poses using multi-criteria scoring.
Input:  PlacerResult + QC results (PoseBusters, Privateer, geometry)
Output: Ranked list with composite score

Rank criteria (descending priority):
  1. PLACER_score (lower = better)
  2. Privateer quality (all sugars pass)
  3. Cu geometry sanity (His-brace in range)
  4. Clash score (lower = better, from PoseBusters)
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

from lpmo_pipeline.placer.run_placer import PlacerPose, PlacerResult
from lpmo_pipeline.qc.qc_report import PoseQCVerdict

logger = logging.getLogger(__name__)


@dataclass
class RankedPose:
    """Pose with composite ranking information."""

    pose_id: str
    rank: int
    placer_score: float
    privateer_ok: bool
    cu_geom_ok: bool
    clash_count: int
    composite_score: float
    qc_status: str  # "passed" | "flagged" | "dropped"
    pdb_path: str = ""


def rank_ensemble(
    placer_result: PlacerResult,
    qc_verdicts: dict[str, PoseQCVerdict],
    weights: dict[str, float] | None = None,
) -> list[RankedPose]:
    """Re-rank PLACER poses incorporating QC results.

    Dropped poses are pushed to the bottom.
    Flagged poses are kept but penalized slightly.

    Args:
        placer_result: PLACER output with initial scores.
        qc_verdicts: {pose_id: PoseQCVerdict} from qc_report.
        weights: Optional {criterion: weight} for composite scoring.
                 Defaults: placer=0.4, privateer=0.2, cu_geom=0.2, clash=0.2

    Returns:
        Sorted list of RankedPose (best first).
    """
    if weights is None:
        weights = {
            "placer": 0.4,
            "privateer": 0.2,
            "cu_geom": 0.2,
            "clash": 0.2,
        }

    ranked: list[RankedPose] = []

    for pose in placer_result.poses:
        verdict = qc_verdicts.get(pose.pose_id)
        privateer_ok = verdict.privateer_passed if verdict else True
        cu_geom_ok = verdict.geometry_passed if verdict else True
        qc_status = verdict.status if verdict else "unknown"
        clash_count = pose.clash_count

        # Composite score (lower = better)
        # Normalize PLACER score to 0–1 range (placeholder: divide by max_score)
        placer_norm = pose.placer_score  # Already lower=better
        priv_penalty = 0.0 if privateer_ok else 100.0
        cu_penalty = 0.0 if cu_geom_ok else 100.0
        clash_norm = float(clash_count)

        composite = (
            weights["placer"] * placer_norm
            + weights["privateer"] * priv_penalty
            + weights["cu_geom"] * cu_penalty
            + weights["clash"] * clash_norm
        )

        # Dropped poses get infinite score (pushed to bottom)
        if qc_status == "dropped":
            composite = float("inf")

        ranked.append(RankedPose(
            pose_id=pose.pose_id,
            rank=0,  # Assigned after sorting
            placer_score=pose.placer_score,
            privateer_ok=privateer_ok,
            cu_geom_ok=cu_geom_ok,
            clash_count=clash_count,
            composite_score=composite,
            qc_status=qc_status,
            pdb_path=str(pose.pdb_path),
        ))

    # Sort by composite score (ascending)
    ranked.sort(key=lambda r: r.composite_score)
    for i, r in enumerate(ranked, 1):
        r.rank = i

    n_usable = sum(1 for r in ranked if r.qc_status != "dropped")
    logger.info(
        "Ensemble ranking: %d total, %d usable, best composite=%.3f (pose %s)",
        len(ranked), n_usable,
        ranked[0].composite_score if ranked else float("nan"),
        ranked[0].pose_id if ranked else "N/A",
    )
    return ranked
