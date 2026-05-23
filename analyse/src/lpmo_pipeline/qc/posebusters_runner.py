# src/lpmo_pipeline/qc/posebusters_runner.py
"""
Responsibility: Run PoseBusters on AF3 poses after normalization.
Input:  AF3 PoseBusters PDB export. Combined protein+glycan exports are
    auto-split into ligand-only ``mol_pred`` and protein ``mol_cond``
    so PoseBusters is invoked with the documented ``dock``/``redock``
    contract instead of treating the whole complex as a standalone ligand.
Output: PoseBustersResult with per-test pass/fail + error types

Current runtime contract:
        - combined AF3 exports are auto-split into ligand-only ``mol_pred`` and
            protein ``mol_cond`` inputs, then PoseBusters runs in built-in ``dock``
            mode
        - the pipeline does not override PoseBusters dock thresholds
        - in the installed PoseBusters build, the intermolecular-distance module
            uses ``max_distance=5.0`` and ``search_distance=6.0``
        - ``protein-ligand_maximum_distance`` is the far-away check, while
            ``minimum_distance_to_protein`` is the renamed ``no_clashes`` result

GATE: no_critical_posebusters_errors = true  (hard-fail → drop pose)
SOFT: minor warnings → keep pose, mark flag
"""
from __future__ import annotations

import csv
from contextlib import contextmanager
import io
import json
import logging
import math
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator

from lpmo_pipeline.config import load_defaults_config, load_runtime_paths_config

logger = logging.getLogger(__name__)

_RUNTIME_PATHS = load_runtime_paths_config()
_DEFAULTS_CONFIG = load_defaults_config()
APPTAINER_EXECUTABLE = _RUNTIME_PATHS.external_tools.apptainer_executable


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
    "protein-ligand_maximum_distance",
    "minimum_distance_to_protein",
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

# SIF container candidates (first existing path is used).
# Here, "fallback" means path selection between known equivalent SIF files,
# not skipping PoseBusters as a QC stage.
POSEBUSTERS_SIF_CANDIDATES = _RUNTIME_PATHS.external_tools.posebusters_sif_candidates

_PROTEIN_CHAIN_ID = str(((_DEFAULTS_CONFIG.get("chain_schema") or {}).get("protein")) or "A")


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

    Args:
        pdb_path: Path to AF3 PoseBusters PDB export. Combined protein+ligand
            files are auto-split into ligand/protein inputs before PoseBusters
            is invoked.
        pose_id: Identifier for this pose.
        protein_path: Optional protein PDB for dock/redock mode.
        reference_path: Optional true ligand for redock RMSD.

    Returns:
        PoseBustersSingleResult.
    """
    logger.info("Running PoseBusters on pose %s: %s", pose_id, pdb_path)

    all_tests = _run_pb_sif(pdb_path, protein_path, reference_path)

    return _classify_results(pose_id, all_tests)


# ---------------------------------------------------------------------------
# PoseBusters backend (SIF-only)
# ---------------------------------------------------------------------------
def _run_pb_sif(
    pdb_path: Path,
    protein_path: Path | None,
    reference_path: Path | None,
) -> dict[str, bool]:
    """Run PoseBusters via SIF container.  Returns ``{test_name: passed}``."""
    with _prepare_posebusters_inputs(
        pdb_path,
        protein_path=protein_path,
        reference_path=reference_path,
    ) as (mol_pred_path, mol_cond_path, mol_true_path):
        api_results = _run_pb_api(mol_pred_path, mol_cond_path, mol_true_path)
        if api_results is not None:
            return api_results

        sif_path = next((p for p in POSEBUSTERS_SIF_CANDIDATES if p.exists()), None)
        if sif_path is None:
            logger.error("PoseBusters SIF not found in candidates: %s", POSEBUSTERS_SIF_CANDIDATES)
            return {}

        cmd: list[str] = [
            APPTAINER_EXECUTABLE, "exec", str(sif_path),
            "bust", str(mol_pred_path),
        ]
        if mol_cond_path:
            cmd.extend(["-p", str(mol_cond_path)])
        if mol_true_path:
            cmd.extend(["-l", str(mol_true_path)])
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


def _run_pb_api(
    pdb_path: Path,
    protein_path: Path | None,
    reference_path: Path | None,
) -> dict[str, bool] | None:
    """Run PoseBusters through the installed Python API when available."""
    try:
        from posebusters import PoseBusters
    except ImportError:
        return None

    config = _select_pb_api_config(protein_path=protein_path, reference_path=reference_path)
    kwargs: dict[str, str] = {}
    if protein_path is not None:
        kwargs["mol_cond"] = str(protein_path)
    if protein_path is not None and reference_path is not None:
        kwargs["mol_true"] = str(reference_path)

    logger.info(
        "PoseBusters Python API config=%s mol_pred=%s mol_cond=%s mol_true=%s",
        config,
        pdb_path,
        protein_path,
        reference_path,
    )

    try:
        results_table = PoseBusters(config=config, max_workers=0).bust(
            str(pdb_path),
            full_report=False,
            **kwargs,
        )
    except Exception as exc:
        logger.warning("PoseBusters Python API failed for %s: %s", pdb_path, exc)
        return None

    return _parse_result_table(results_table)


def _select_pb_api_config(
    *,
    protein_path: Path | None,
    reference_path: Path | None,
) -> str:
    if protein_path is not None and reference_path is not None:
        return "redock"
    if protein_path is not None:
        return "dock"
    return "mol"


@contextmanager
def _prepare_posebusters_inputs(
    pdb_path: Path,
    *,
    protein_path: Path | None,
    reference_path: Path | None,
) -> Iterator[tuple[Path, Path | None, Path | None]]:
    """Yield PoseBusters inputs that follow the documented ligand/protein contract.

    AF3 exports enter this module as combined protein+glycan PDB files. PoseBusters,
    however, expects ``mol_pred`` to be the ligand and ``mol_cond`` to be the
    conditioning protein for ``dock``/``redock`` mode. When we detect a combined
    export, split it on the fly into ligand-only and protein-only inputs.
    """

    split_lines = _split_combined_pose_pdb_lines(pdb_path)
    if split_lines is None:
        yield pdb_path, protein_path, reference_path
        return

    ligand_lines, protein_lines = split_lines
    with tempfile.TemporaryDirectory(prefix="posebusters_input_") as tmp_dir_name:
        tmp_dir = Path(tmp_dir_name)
        ligand_path = tmp_dir / f"{pdb_path.stem}_ligand.pdb"
        auto_protein_path = tmp_dir / f"{pdb_path.stem}_protein.pdb"
        ligand_path.write_text("\n".join(ligand_lines) + "\n")
        auto_protein_path.write_text("\n".join(protein_lines) + "\n")

        effective_protein_path = protein_path or auto_protein_path
        logger.info(
            "PoseBusters auto-split combined input %s into ligand=%s protein=%s",
            pdb_path,
            ligand_path,
            effective_protein_path,
        )
        yield ligand_path, effective_protein_path, reference_path


def _split_combined_pose_pdb_lines(pdb_path: Path) -> tuple[list[str], list[str]] | None:
    """Split a combined AF3 PDB export into ligand-only and protein-only records.

    Returns ``None`` when the input does not look like the combined AF3 export we
    generate for hard QC.
    """

    try:
        lines = pdb_path.read_text().splitlines()
    except OSError as exc:
        logger.warning("Failed to read PoseBusters input %s for auto-split: %s", pdb_path, exc)
        return None

    protein_lines: list[str] = []
    ligand_lines: list[str] = []
    ligand_serials: set[int] = set()
    has_protein = False
    has_ligand = False

    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue

        serial = _parse_pdb_serial(line[6:11])
        if serial is None:
            continue

        chain_id = line[21].strip()
        if chain_id == _PROTEIN_CHAIN_ID:
            protein_lines.append(line)
            has_protein = True
        else:
            ligand_lines.append(line)
            ligand_serials.add(serial)
            has_ligand = True

    if not (has_protein and has_ligand):
        return None

    for line in lines:
        if not line.startswith("CONECT"):
            continue
        serials = _parse_pdb_conect_serials(line)
        kept_serials = [serial for serial in serials if serial in ligand_serials]
        if len(kept_serials) >= 2:
            ligand_lines.append(_format_pdb_conect_line(kept_serials))

    protein_lines.append("END")
    ligand_lines.append("END")
    return ligand_lines, protein_lines


def _parse_pdb_serial(field: str) -> int | None:
    try:
        return int(field.strip())
    except ValueError:
        return None


def _parse_pdb_conect_serials(line: str) -> list[int]:
    serials: list[int] = []
    for token in line[6:].split():
        serial = _parse_pdb_serial(token)
        if serial is not None:
            serials.append(serial)
    return serials


def _format_pdb_conect_line(serials: list[int]) -> str:
    return "CONECT" + "".join(f"{serial:5d}" for serial in serials)


def _parse_result_table(results_table: Any) -> dict[str, bool]:
    """Convert the first PoseBusters result row into ``{test_name: passed}``."""
    try:
        row = results_table.iloc[0]
    except Exception:
        return {}

    return {
        str(key): _coerce_pb_value(value)
        for key, value in row.items()
        if key not in _LOADING_COLUMNS
    }


def _coerce_pb_value(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() == "true"
    try:
        if math.isnan(value):
            return False
    except TypeError:
        pass
    return bool(value)


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
