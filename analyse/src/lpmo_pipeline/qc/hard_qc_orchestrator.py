# src/lpmo_pipeline/qc/hard_qc_orchestrator.py
"""
Responsibility: Orchestrate the hard QC sequence for a batch of poses.

Sequence per pose:
  1. Run PoseBusters (if available; graceful fallback on failure)
  2. Run Cu-His + substrate geometry checks
  3. Skip Privateer (disabled until download issues resolved)
  4. Aggregate into unified QC verdict via compute_verdict()

INVARIANT — "Ikke-slett regel":
  All numeric metrics are preserved in the verdict even when the pose is
  dropped by a hard gate.  Nothing computed is discarded.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from lpmo_pipeline.qc.custom_geometry_checks import GeometryResult, check_geometry
from lpmo_pipeline.qc.posebusters_runner import PoseBustersSingleResult, run_posebusters_single
from lpmo_pipeline.qc.qc_report import (
    PoseQCVerdict,
    QCReport,
    build_qc_report,
    compute_verdict,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Input model
# ---------------------------------------------------------------------------
@dataclass
class HardQCInput:
    """Input for a single pose entering hard QC.

    Attributes:
        pose_id: Unique identifier for this pose.
        mol_pred_path: Path to the PDB file (for PoseBusters).
        structure: Parsed structure object (gemmi) for geometry checks.
        protein_path: Optional protein PDB for PoseBusters dock/redock mode.
        reference_path: Optional reference ligand for PoseBusters redock RMSD.
        cu_chain: Expected Cu chain (default ``"E"`` per normalisation).
        glycan_chains: Expected glycan chains (default ``["B", "C", "D"]``).
    """

    pose_id: str
    mol_pred_path: Path
    structure: Any  # gemmi.Structure
    protein_path: Path | None = None
    reference_path: Path | None = None
    cu_chain: str = "E"
    glycan_chains: list[str] | None = None


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------
def run_hard_qc(
    poses: list[HardQCInput],
    run_id: str = "",
) -> QCReport:
    """Run the full hard-QC sequence on a list of poses.

    For each pose the following checks are executed:
      1. **PoseBusters** — chemical/stereochemical validation.
         If the PoseBusters backend raises an unexpected exception the pose
         is still processed (``pb_result`` set to ``None``).
      2. **Cu-His + substrate geometry** — LPMO-specific active-site check.
      3. **Privateer** — *currently disabled* (``priv_result=None``).

    All numeric metrics are retained regardless of pass/fail status
    ("Ikke-slett regel").

    Args:
        poses: List of :class:`HardQCInput` objects.
        run_id: Identifier for this QC run (used in the report).

    Returns:
        :class:`QCReport` with per-pose verdicts and aggregate counts.
    """
    verdicts: list[PoseQCVerdict] = []

    for pose in poses:
        logger.info("Hard QC: processing pose %s", pose.pose_id)

        # --- Step 1: PoseBusters ---
        pb_result: PoseBustersSingleResult | None = None
        try:
            pb_result = run_posebusters_single(
                pdb_path=pose.mol_pred_path,
                pose_id=pose.pose_id,
                protein_path=pose.protein_path,
                reference_path=pose.reference_path,
            )
        except Exception:
            logger.exception(
                "PoseBusters raised an unexpected error for pose %s; "
                "continuing without PB result",
                pose.pose_id,
            )

        # --- Step 2: Geometry (Cu-His + Cu-substrate) ---
        geom_result: GeometryResult | None = None
        try:
            geom_result = check_geometry(
                structure=pose.structure,
                pose_id=pose.pose_id,
                cu_chain=pose.cu_chain,
                glycan_chains=pose.glycan_chains,
            )
        except Exception:
            logger.exception(
                "Geometry check raised an unexpected error for pose %s",
                pose.pose_id,
            )

        # --- Step 3: Privateer (disabled) ---
        priv_result = None

        # --- Step 4: Aggregate verdict ---
        verdict = compute_verdict(
            pose_id=pose.pose_id,
            pb_result=pb_result,
            priv_result=priv_result,
            geom_result=geom_result,
        )
        verdicts.append(verdict)

        logger.info(
            "Hard QC verdict for %s: %s (drop_reasons=%s)",
            pose.pose_id,
            verdict.status,
            verdict.drop_reasons or "none",
        )

    report = build_qc_report(run_id=run_id, verdicts=verdicts)
    logger.info(
        "Hard QC complete: %d total, %d passed, %d flagged, %d dropped",
        report.total,
        report.passed,
        report.flagged,
        report.dropped,
    )
    return report