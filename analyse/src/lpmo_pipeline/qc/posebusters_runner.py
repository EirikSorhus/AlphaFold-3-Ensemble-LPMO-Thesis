# src/lpmo_pipeline/qc/posebusters_runner.py
"""
Responsibility: Run PoseBusters on PLACER-refined poses.
Input:  PDB file (for_posebusters.pdb) from protonation step
Output: PoseBustersResult with per-test pass/fail + error types

GATE: no_critical_posebusters_errors = true  (hard-fail → drop pose)
SOFT: minor warnings → keep pose, mark flag
"""
from __future__ import annotations

import csv
import io
import json
import logging
import subprocess
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
    "sanitization",
    "all_atoms_connected",
    "no_radicals",
    "internal_steric_clash",
    "bond_lengths",
    "bond_angles",
    "tetrahedral_chirality",
    "volume_overlap_with_protein",
}

SOFT_WARNING_TYPES: set[str] = {
    "aromatic_ring_flatness",
    "double_bond_flatness",
    "double_bond_stereochemistry",
    "non-aromatic_ring_non-flatness",
    "internal_energy",
    "inchi_convertible",
}

# Columns from PoseBusters that track file loading, not test results
_LOADING_COLUMNS: set[str] = {"mol_pred_loaded", "mol_true_loaded", "mol_cond_loaded"}

# SIF container path for fallback execution
POSEBUSTERS_SIF: Path = Path(
    "/cluster/projects/nn1003k/prog/posebusters/build/posebusters.sif"
)


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------
def run_posebusters_single(
    pdb_path: Path,
    pose_id: str,
    protein_path: Path | None = None,
    reference_path: Path | None = None,
) -> PoseBustersSingleResult:
    """Run PoseBusters on a single PDB pose file.

    Uses the PoseBusters Python API as primary backend.  Falls back to the
    SIF container at ``POSEBUSTERS_SIF`` if the library is not importable.

    Args:
        pdb_path: Path to PDB file with CONECT records.
        pose_id: Identifier for this pose.
        protein_path: Optional protein PDB for dock/redock mode.
        reference_path: Optional true ligand for redock RMSD.

    Returns:
        PoseBustersSingleResult.
    """
    logger.info("Running PoseBusters on pose %s: %s", pose_id, pdb_path)

    try:
        all_tests = _run_pb_python_api(pdb_path, protein_path, reference_path)
    except _PBImportError:
        logger.info("PoseBusters Python API unavailable, trying SIF fallback")
        all_tests = _run_pb_sif(pdb_path, protein_path, reference_path)

    return _classify_results(pose_id, all_tests)


# ---------------------------------------------------------------------------
# PoseBusters backends
# ---------------------------------------------------------------------------
class _PBImportError(Exception):
    """PoseBusters Python package not importable."""


def _run_pb_python_api(
    pdb_path: Path,
    protein_path: Path | None,
    reference_path: Path | None,
) -> dict[str, bool]:
    """Run PoseBusters via Python API.  Returns ``{test_name: passed}``."""
    try:
        from posebusters import PoseBusters  # type: ignore[import-untyped]
    except ImportError as exc:
        raise _PBImportError("posebusters not installed") from exc

    if protein_path and reference_path:
        config = "redock"
    elif protein_path:
        config = "dock"
    else:
        config = "mol"

    pb = PoseBusters(config=config, max_workers=0)

    try:
        results_df = pb.bust(
            mol_pred=pdb_path,
            mol_cond=protein_path,
            mol_true=reference_path,
        )
    except Exception:
        logger.exception("PoseBusters bust() raised for %s", pdb_path)
        return {}

    if results_df.empty:
        logger.warning("PoseBusters returned empty DataFrame for %s", pdb_path)
        return {}

    row = results_df.iloc[0]
    return {
        col: bool(row[col])
        for col in results_df.columns
        if col not in _LOADING_COLUMNS
    }


def _run_pb_sif(
    pdb_path: Path,
    protein_path: Path | None,
    reference_path: Path | None,
) -> dict[str, bool]:
    """Run PoseBusters via SIF container.  Returns ``{test_name: passed}``."""
    if not POSEBUSTERS_SIF.exists():
        logger.error("PoseBusters SIF not found: %s", POSEBUSTERS_SIF)
        return {}

    cmd: list[str] = [
        "apptainer", "exec", str(POSEBUSTERS_SIF),
        "bust", str(pdb_path),
    ]
    if protein_path:
        cmd.extend(["-p", str(protein_path)])
    if reference_path:
        cmd.extend(["-l", str(reference_path)])
    cmd.extend(["--outfmt", "csv"])

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, check=True, timeout=300,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired,
            FileNotFoundError) as exc:
        logger.error("PoseBusters SIF execution failed: %s", exc)
        return {}

    return _parse_csv_output(proc.stdout)


def _parse_csv_output(csv_text: str) -> dict[str, bool]:
    """Parse PoseBusters CSV stdout into ``{test_name: passed}`` dict."""
    skip = {"file", "molecule", "position"}
    reader = csv.DictReader(io.StringIO(csv_text))
    for row in reader:
        return {
            k: v.strip().lower() == "true"
            for k, v in row.items()
            if k not in skip and k not in _LOADING_COLUMNS and v is not None
        }
    return {}


def _classify_results(
    pose_id: str,
    all_tests: dict[str, bool],
) -> PoseBustersSingleResult:
    """Classify PoseBusters test results into critical errors / soft warnings."""
    if not all_tests:
        return PoseBustersSingleResult(
            pose_id=pose_id,
            passed=False,
            critical_errors=["posebusters_no_results"],
        )

    critical: list[str] = []
    warnings: list[str] = []
    for test_name, passed in all_tests.items():
        if passed:
            continue
        if test_name in CRITICAL_ERROR_TYPES:
            critical.append(test_name)
        elif test_name in SOFT_WARNING_TYPES:
            warnings.append(test_name)
        else:
            # Unknown failure → conservative: treat as critical
            critical.append(test_name)

    return PoseBustersSingleResult(
        pose_id=pose_id,
        passed=len(critical) == 0,
        critical_errors=critical,
        warnings=warnings,
        details={k: {"passed": v} for k, v in all_tests.items()},
    )


def run_posebusters_batch(
    pdb_dir: Path,
    pose_ids: list[str],
    protein_path: Path | None = None,
) -> PoseBustersBatchResult:
    """Run PoseBusters on all poses in a directory.

    Args:
        pdb_dir: Directory containing per-pose PDB files.
        pose_ids: List of pose IDs (filenames without extension).
        protein_path: Optional protein PDB for dock/redock checks.

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
        result = run_posebusters_single(pdb_path, pid, protein_path=protein_path)
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
