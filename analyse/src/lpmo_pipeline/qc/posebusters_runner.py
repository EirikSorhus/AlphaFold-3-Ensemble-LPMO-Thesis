# src/lpmo_pipeline/qc/posebusters_runner.py
"""
Responsibility: Run PoseBusters on PLACER-refined poses.
Input:  PDB file (for_posebusters.pdb) from protonation step
Output: PoseBustersResult with per-test pass/fail + error types

GATE: no_critical_posebusters_errors = true  (hard-fail → drop pose)
SOFT: minor warnings → keep pose, mark flag
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Result models
# ---------------------------------------------------------------------------
@dataclass
class PoseBustersSingleResult:
    """Result for a single pose."""

    pose_id: str
    passed: bool
    critical_errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    details: dict[str, Any] = field(default_factory=dict)


@dataclass
class PoseBustersBatchResult:
    """Aggregated result for a batch of poses."""

    results: list[PoseBustersSingleResult] = field(default_factory=list)
    total: int = 0
    passed: int = 0
    failed: int = 0
    pass_rate: float = 0.0
    error_distribution: dict[str, int] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Critical vs soft error classification
# ---------------------------------------------------------------------------
CRITICAL_ERROR_TYPES: set[str] = {
    "steric_clash",
    "chem_valence_violation",
    "bond_length_outlier",
    "ring_not_planar",
    "chirality_error",
}

SOFT_WARNING_TYPES: set[str] = {
    "minor_angle_deviation",
    "unusual_torsion",
    "slight_bump",
}


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------
def run_posebusters_single(
    pdb_path: Path,
    pose_id: str,
) -> PoseBustersSingleResult:
    """Run PoseBusters on a single PDB pose file.

    PSEUDOCODE — wraps the posebusters Python API.

    Args:
        pdb_path: Path to PDB file with CONECT records.
        pose_id: Identifier for this pose.

    Returns:
        PoseBustersSingleResult.
    """
    logger.info("Running PoseBusters on pose %s: %s", pose_id, pdb_path)

    # --- Step 1: Import and run PoseBusters ---
    # from posebusters import PoseBusters
    # pb = PoseBusters(config="redock")   # or "mol" depending on mode
    # results_df = pb.bust(pdb_path, ...)
    # PSEUDOCODE: parse results_df into structured output

    # --- Step 2: Classify errors ---
    all_tests: dict[str, Any] = {}  # {test_name: {passed: bool, value: float}}
    # PSEUDOCODE: iterate over results_df columns
    # for col in results_df.columns:
    #     all_tests[col] = {"passed": bool(results_df[col].iloc[0]), ...}

    critical: list[str] = []
    warnings: list[str] = []
    for test_name, test_result in all_tests.items():
        if not test_result.get("passed", True):
            if test_name in CRITICAL_ERROR_TYPES:
                critical.append(test_name)
            elif test_name in SOFT_WARNING_TYPES:
                warnings.append(test_name)
            else:
                # Unknown test failure → treat as critical (conservative)
                critical.append(test_name)

    passed = len(critical) == 0

    return PoseBustersSingleResult(
        pose_id=pose_id,
        passed=passed,
        critical_errors=critical,
        warnings=warnings,
        details=all_tests,
    )


def run_posebusters_batch(
    pdb_dir: Path,
    pose_ids: list[str],
) -> PoseBustersBatchResult:
    """Run PoseBusters on all poses in a directory.

    Args:
        pdb_dir: Directory containing per-pose PDB files.
        pose_ids: List of pose IDs (filenames without extension).

    Returns:
        PoseBustersBatchResult with per-pose results and aggregates.
    """
    results: list[PoseBustersSingleResult] = []

    for pid in pose_ids:
        pdb_path = pdb_dir / f"{pid}.pdb"
        if not pdb_path.exists():
            logger.warning("PDB not found for pose %s: %s", pid, pdb_path)
            results.append(PoseBustersSingleResult(
                pose_id=pid, passed=False,
                critical_errors=["file_not_found"],
            ))
            continue
        result = run_posebusters_single(pdb_path, pid)
        results.append(result)

    # Aggregate
    total = len(results)
    passed = sum(1 for r in results if r.passed)
    failed = total - passed

    # Error distribution
    error_dist: dict[str, int] = {}
    for r in results:
        for err in r.critical_errors + r.warnings:
            error_dist[err] = error_dist.get(err, 0) + 1

    batch = PoseBustersBatchResult(
        results=results,
        total=total,
        passed=passed,
        failed=failed,
        pass_rate=passed / total if total > 0 else 0.0,
        error_distribution=error_dist,
    )

    logger.info(
        "PoseBusters batch: %d/%d passed (%.1f%%). Critical errors: %s",
        passed, total, batch.pass_rate * 100, error_dist,
    )
    return batch


def write_posebusters_report(
    batch: PoseBustersBatchResult,
    output_path: Path,
) -> None:
    """Write PoseBusters results to JSON."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "total": batch.total,
        "passed": batch.passed,
        "failed": batch.failed,
        "pass_rate": batch.pass_rate,
        "error_distribution": batch.error_distribution,
        "per_pose": [
            {
                "pose_id": r.pose_id,
                "passed": r.passed,
                "critical_errors": r.critical_errors,
                "warnings": r.warnings,
            }
            for r in batch.results
        ],
    }
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    logger.info("Wrote PoseBusters report to %s", output_path)
