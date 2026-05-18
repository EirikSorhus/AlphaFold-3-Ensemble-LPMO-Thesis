# src/lpmo_pipeline/qc/privateer_runner.py
"""
Responsibility: Run Privateer validation on glycan residues.
Input:  privateer_input.cif (must contain ONLY valid CCD monosaccharides)
Output: PrivateerResult with per-residue sugar_identity, anomer, ring_pucker, linkage

HARD RULES:
  - privateer_recognized_sugars must be 100%
  - Anomer errors → hard fail (drop pose)
  - Input MUST use valid CCD monosaccharide codes (enforced in STEP 3)
"""
from __future__ import annotations

import json
import logging
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from lpmo_pipeline.config import load_runtime_paths_config
from lpmo_pipeline.utils.exceptions import PrivateerCCDError

logger = logging.getLogger(__name__)


# Privateer must run from SIF, never from local env/PATH.
_RUNTIME_PATHS = load_runtime_paths_config()
APPTAINER_EXECUTABLE = _RUNTIME_PATHS.external_tools.apptainer_executable
PRIVATEER_SIF_CANDIDATES = _RUNTIME_PATHS.external_tools.privateer_sif_candidates
PRIVATEER_DEFAULT_MODE = _RUNTIME_PATHS.runtime_settings.privateer_default_mode

_CCD_MONOSACCHARIDE_CODES: frozenset[str] = frozenset(
    {
        "NAG", "BGC", "GLC", "MAN", "BMA", "GAL", "FUC", "XYS",
        "ARA", "FRU", "NDG", "SIA", "NAN", "STA",
    }
)
_GLYCAN_SENTINEL_ATOMS: frozenset[str] = frozenset(
    {"C1", "C2", "C3", "C4", "C5", "C6", "O3", "O4", "O5", "O6", "N2", "C7", "C8", "O7"}
)


# ---------------------------------------------------------------------------
# Result models
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class PrivateerBatchInput:
    """Batch input for Privateer execution."""

    cif_path: Path
    pose_id: str
    output_dir: Path | None = None
    dry_run: bool = False


@dataclass(frozen=True)
class PrivateerInvocation:
    """Resolved Privateer container invocation."""

    command: list[str]
    bind_mounts: list[str]
    work_dir: str


@dataclass(frozen=True)
class PrivateerSummary:
    """Counts parsed from Privateer ccp4i2 summary text."""

    wrong_anomer: int = 0
    wrong_configuration: int = 0
    unphysical_puckering: int = 0
    high_energy_conformations: int = 0
    issues_detected: int = 0
    affected_sugars: int = 0
    total_sugars_reported: int = 0


@dataclass
class PrivateerResidueResult:
    """Privateer result for a single glycan residue."""

    chain: str
    resname: str
    resnum: int
    sugar_recognized: bool
    anomer_ok: bool
    ring_pucker_ok: bool
    ring_pucker_conformation: str = ""  # e.g. "4C1", "1C4"
    linkage_ok: bool = True
    diagnostics: dict[str, Any] = field(default_factory=dict)


@dataclass
class PrivateerResult:
    """Aggregated Privateer result for a structure."""

    pose_id: str = ""
    residues: list[PrivateerResidueResult] = field(default_factory=list)
    total_sugars: int = 0
    recognized: int = 0
    recognition_rate: float = 0.0
    anomer_pass: int = 0
    ring_pucker_pass: int = 0
    linkage_pass: int = 0
    all_pass: bool = False
    runner_error: str = ""
    command: list[str] = field(default_factory=list)
    validation_data_path: str = ""
    stdout_path: str = ""
    stderr_path: str = ""
    work_dir: str = ""
    returncode: int = -1
    dry_run: bool = False


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------
def run_privateer(
    cif_path: Path,
    pose_id: str = "",
    output_dir: Path | None = None,
    structure_factor_cif: Path | None = None,
    structure_factor_mtz: Path | None = None,
    colin_fo: str = "",
    codein: str = "",
    showgeom: bool = False,
    radiusin: float | None = None,
    dry_run: bool = False,
    probe_container: bool = True,
    debug_output: bool = False,
) -> PrivateerResult:
    """Run Privateer on a CIF file and parse results.

    Args:
        cif_path: Path to privateer_input.cif (with CCD monosaccharides).
        pose_id: Pose identifier for logging context.

    Returns:
        PrivateerResult with per-residue and aggregate metrics.
    """
    cif_path = Path(cif_path)
    resolved_pose_id = pose_id or cif_path.stem

    if not cif_path.exists():
        return PrivateerResult(
            pose_id=resolved_pose_id,
            runner_error=f"Privateer input does not exist: {cif_path}",
        )

    if (
        output_dir is not None
        or dry_run
        or structure_factor_cif is not None
        or structure_factor_mtz is not None
    ):
        resolved_output_dir = Path(output_dir) if output_dir is not None else _default_privateer_output_dir(
            cif_path, resolved_pose_id
        )
        invocation = build_privateer_invocation(
            cif_path=cif_path,
            output_dir=resolved_output_dir,
            structure_factor_cif=structure_factor_cif,
            structure_factor_mtz=structure_factor_mtz,
            colin_fo=colin_fo,
            codein=codein,
            showgeom=showgeom,
            radiusin=radiusin,
        )
        validation_data_path = resolved_output_dir / "validation_data-privateer"
        if dry_run:
            return PrivateerResult(
                pose_id=resolved_pose_id,
                command=invocation.command,
                validation_data_path=str(validation_data_path),
                work_dir=invocation.work_dir,
                dry_run=True,
            )
        if probe_container:
            _probe_privateer_container(invocation)
        return _run_privateer_ccp4i2(
            cif_path=cif_path,
            pose_id=resolved_pose_id,
            output_dir=resolved_output_dir,
            invocation=invocation,
            debug_output=debug_output,
        )

    logger.info("Running Privateer on %s (pose=%s)", cif_path, pose_id)

    try:
        raw = _run_privateer_cli(cif_path)
        residue_results = _parse_privateer_output(raw)
        result = _aggregate_privateer_results(residue_results)
        result.pose_id = pose_id
    except Exception as exc:
        logger.error("Privateer execution failed for %s: %s", cif_path, exc)
        # Conservative behavior: failed run is treated as failed gate.
        return PrivateerResult(
            pose_id=pose_id,
            residues=[],
            total_sugars=0,
            recognized=0,
            recognition_rate=0.0,
            anomer_pass=0,
            ring_pucker_pass=0,
            linkage_pass=0,
            all_pass=False,
            runner_error=str(exc),
        )

    return _finalize_privateer_result(result, residue_results, cif_path=cif_path, pose_id=pose_id)


def run_privateer_batch(
    inputs: list[PrivateerBatchInput],
    *,
    max_workers: int | None = None,
) -> list[PrivateerResult]:
    """Run Privateer on a batch of prepared CIF inputs."""
    _ = max_workers
    results: list[PrivateerResult] = []
    for item in inputs:
        results.append(
            run_privateer(
                item.cif_path,
                pose_id=item.pose_id,
                output_dir=item.output_dir or _default_privateer_output_dir(item.cif_path, item.pose_id),
                dry_run=item.dry_run,
            )
        )
    return results


def prepare_privateer_input(input_cif: Path, output_cif: Path | None = None) -> Path:
    """Prepare Privateer input CIF.

    The normalized mmCIF is already in the form expected by the current
    Privateer container path, so preparation is a controlled copy into the
    case-local working location.
    """
    input_cif = Path(input_cif)
    if not input_cif.exists():
        raise FileNotFoundError(input_cif)

    text = input_cif.read_text()
    _validate_privateer_ccd_codes(text)

    output_cif = Path(output_cif) if output_cif is not None else input_cif.with_name("privateer_input.cif")
    output_cif.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(input_cif, output_cif)
    return output_cif.resolve()


def get_privateer_version() -> str:
    """Return a stable identifier for the configured Privateer backend."""
    sif_path = next((p for p in PRIVATEER_SIF_CANDIDATES if p.exists()), None)
    if sif_path is None:
        raise RuntimeError(
            f"Privateer SIF not found in candidates: {PRIVATEER_SIF_CANDIDATES}"
        )
    try:
        proc = subprocess.run(
            [APPTAINER_EXECUTABLE, "run", "--cleanenv", str(sif_path), "--version"],
            capture_output=True,
            text=True,
            check=False,
            timeout=30,
            env=_privateer_subprocess_env(),
        )
    except Exception:
        return sif_path.name

    text = f"{proc.stdout}\n{proc.stderr}"
    matched = re.search(r"Version\s+([A-Za-z0-9_.-]+)", text, flags=re.IGNORECASE)
    if matched:
        return matched.group(1).rstrip(".,;:")
    return sif_path.name


def _default_privateer_output_dir(cif_path: Path, pose_id: str) -> Path:
    return cif_path.parent / "privateer_output" / _slug(pose_id)


def _slug(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in value)


def build_privateer_invocation(
    *,
    cif_path: Path,
    output_dir: Path,
    structure_factor_cif: Path | None = None,
    structure_factor_mtz: Path | None = None,
    colin_fo: str = "",
    codein: str = "",
    showgeom: bool = False,
    radiusin: float | None = None,
    mode: str = PRIVATEER_DEFAULT_MODE,
) -> PrivateerInvocation:
    sif_path = next((p for p in PRIVATEER_SIF_CANDIDATES if p.exists()), None)
    if sif_path is None:
        raise RuntimeError(
            f"Privateer SIF not found in candidates: {PRIVATEER_SIF_CANDIDATES}"
        )

    case_dir = cif_path.resolve().parent
    output_dir = output_dir.resolve()

    bind_mounts = [f"{case_dir}:{case_dir}", f"{output_dir}:{output_dir}"]
    if structure_factor_cif is not None:
        sf_cif_dir = Path(structure_factor_cif).resolve().parent
        mount = f"{sf_cif_dir}:{sf_cif_dir}"
        if mount not in bind_mounts:
            bind_mounts.append(mount)
    if structure_factor_mtz is not None:
        sf_mtz_dir = Path(structure_factor_mtz).resolve().parent
        mount = f"{sf_mtz_dir}:{sf_mtz_dir}"
        if mount not in bind_mounts:
            bind_mounts.append(mount)

    command = [APPTAINER_EXECUTABLE, "run", "--cleanenv"]
    for bind in bind_mounts:
        command.extend(["--bind", bind])
    command.extend([str(sif_path), "-pdbin", str(cif_path.resolve())])
    if structure_factor_cif is not None:
        command.extend(["-cifin", str(Path(structure_factor_cif).resolve())])
    if structure_factor_mtz is not None:
        command.extend(["-mtzin", str(Path(structure_factor_mtz).resolve())])
    if colin_fo:
        command.extend(["-colin-fo", colin_fo])
    if codein:
        command.extend(["-codein", codein])
    if showgeom:
        command.append("-showgeom")
    if radiusin is not None:
        command.extend(["-radiusin", str(radiusin)])
    command.extend(["-mode", mode])
    return PrivateerInvocation(command=command, bind_mounts=bind_mounts, work_dir=str(output_dir))


def _build_privateer_ccp4i2_command(cif_path: Path, output_dir: Path) -> list[str]:
    return build_privateer_invocation(cif_path=cif_path, output_dir=output_dir).command


def _validate_privateer_ccd_codes(cif_text: str) -> None:
    residue_names: set[str] = set()
    for raw_line in cif_text.splitlines():
        line = raw_line.strip()
        if not line.startswith(("ATOM", "HETATM")):
            continue
        fields = line.split()
        if len(fields) < 6:
            continue
        atom_name = fields[3].upper()
        comp_id = fields[5].upper()
        if atom_name in _GLYCAN_SENTINEL_ATOMS:
            residue_names.add(comp_id)

    invalid = sorted(name for name in residue_names if name not in _CCD_MONOSACCHARIDE_CODES)
    if invalid:
        raise PrivateerCCDError(f"non-CCD glycan residues: {invalid}")


def _probe_privateer_container(_invocation: PrivateerInvocation) -> None:
    """Hook for optional container probing before the full run."""
    return None


def _run_privateer_ccp4i2(
    *,
    cif_path: Path,
    pose_id: str,
    output_dir: Path,
    invocation: PrivateerInvocation,
    debug_output: bool = False,
) -> PrivateerResult:
    output_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = output_dir / "privateer_stdout.txt"
    stderr_path = output_dir / "privateer_stderr.txt"
    validation_data_path = output_dir / "validation_data-privateer"

    logger.info(
        "Running Privateer on %s (pose=%s) with %s",
        cif_path,
        pose_id,
        " ".join(invocation.command),
    )

    try:
        proc = subprocess.run(
            invocation.command,
            capture_output=True,
            text=True,
            check=False,
            timeout=300,
            cwd=output_dir,
            env=_privateer_subprocess_env(),
        )
    except Exception as exc:
        return PrivateerResult(
            pose_id=pose_id,
            runner_error=str(exc),
            command=invocation.command,
            validation_data_path=str(validation_data_path),
            work_dir=invocation.work_dir,
        )

    residue_results: list[PrivateerResidueResult] = []
    runner_error = ""
    summary = _parse_privateer_summary(proc.stdout)
    if validation_data_path.exists() and validation_data_path.stat().st_size > 0:
        try:
            rows = _parse_validation_rows(validation_data_path.read_text())
            residue_results = [_row_to_residue_result(row, summary) for row in rows]
        except Exception as exc:
            runner_error = f"validation_data_parse_error:{exc}"
    elif proc.stdout.strip().lstrip().startswith(("{", "[")):
        try:
            residue_results = _parse_privateer_output(proc.stdout)
        except Exception as exc:
            runner_error = str(exc)
    else:
        runner_error = (
            f"Privateer did not produce validation_data-privateer or JSON stdout "
            f"(returncode={proc.returncode})"
        )

    if residue_results:
        result = _aggregate_privateer_results(
            residue_results=residue_results,
            pose_id=pose_id,
            summary=summary,
            invocation=invocation,
            returncode=proc.returncode,
        )
        result.pose_id = pose_id
        result.command = invocation.command
        result.validation_data_path = str(validation_data_path)
        result.work_dir = invocation.work_dir
        result.returncode = proc.returncode
        result.runner_error = runner_error
        if debug_output:
            stdout_path.write_text(proc.stdout or "")
            stderr_path.write_text(proc.stderr or "")
            result.stdout_path = str(stdout_path)
            result.stderr_path = str(stderr_path)
        return _finalize_privateer_result(result, residue_results, cif_path=cif_path, pose_id=pose_id)

    stdout_path.write_text(proc.stdout or "")
    stderr_path.write_text(proc.stderr or "")
    if not runner_error:
        runner_error = f"Privateer failed with return code {proc.returncode}"
    elif f"returncode={proc.returncode}" in runner_error:
        runner_error = f"Privateer failed with return code {proc.returncode}"

    return PrivateerResult(
        pose_id=pose_id,
        runner_error=runner_error,
        command=invocation.command,
        validation_data_path=str(validation_data_path),
        stdout_path=str(stdout_path),
        stderr_path=str(stderr_path),
        work_dir=invocation.work_dir,
        returncode=proc.returncode,
    )


def _finalize_privateer_result(
    result: PrivateerResult,
    residue_results: list[PrivateerResidueResult],
    *,
    cif_path: Path,
    pose_id: str,
) -> PrivateerResult:
    # Gate check
    if result.recognition_rate < 1.0:
        unrecognized = [
            r.resname for r in residue_results if not r.sugar_recognized
        ]
        logger.error(
            "GATE FAIL: privateer_recognized_sugars = %.1f%% (must be 100%%). "
            "Unrecognized: %s",
            result.recognition_rate * 100,
            unrecognized,
        )

    if result.anomer_pass < result.total_sugars:
        bad_anomers = [
            f"{r.chain}:{r.resname}{r.resnum}"
            for r in residue_results
            if not r.anomer_ok
        ]
        logger.error(
            "GATE FAIL: Privateer anomer errors on: %s", bad_anomers
        )

    logger.info(
        "Privateer: %d/%d recognized, %d/%d anomer OK, %d/%d ring OK",
        result.recognized,
        result.total_sugars,
        result.anomer_pass,
        result.total_sugars,
        result.ring_pucker_pass,
        result.total_sugars,
    )
    return result


def _run_privateer_cli(cif_path: Path) -> str:
    """Run Privateer from SIF and return stdout JSON."""
    sif_path = next((p for p in PRIVATEER_SIF_CANDIDATES if p.exists()), None)
    if sif_path is None:
        raise RuntimeError(
            f"Privateer SIF not found in candidates: {PRIVATEER_SIF_CANDIDATES}"
        )

    cmd = [
        APPTAINER_EXECUTABLE,
        "run",
        "--cleanenv",
        str(sif_path),
        "-pdbin",
        str(cif_path),
    ]

    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True,
        timeout=300,
        env=_privateer_subprocess_env(),
    )
    if not proc.stdout.strip():
        raise RuntimeError("Privateer produced empty stdout")
    if not proc.stdout.lstrip().startswith(("{", "[")):
        raise RuntimeError(
            "Installed Privateer build did not emit JSON. "
            "This container exposes text/HTML output with -pdbin and no supported JSON flag."
        )
    return proc.stdout


def _privateer_subprocess_env() -> dict[str, str]:
    """Return a minimal host environment for Apptainer subprocesses.

    The analyse conda environment can interfere with Apptainer's view of the
    container filesystem when invoked from Python. Keeping only the basic host
    variables makes the container execution match the working shell behavior.
    """
    return {
        "PATH": "/usr/bin:/bin:/usr/local/bin",
        "HOME": os.environ.get("HOME", "/tmp"),
        "USER": os.environ.get("USER", "unknown"),
        "LOGNAME": os.environ.get("LOGNAME", os.environ.get("USER", "unknown")),
        "LANG": os.environ.get("LANG", "C.UTF-8"),
        "LC_ALL": os.environ.get("LC_ALL", "C.UTF-8"),
        "TERM": os.environ.get("TERM", "xterm"),
    }


def _parse_privateer_output(raw_output: str) -> list[PrivateerResidueResult]:
    """Parse Privateer output into per-residue records.

    Supported format (primary): JSON with residue records under one of
    ``residues``, ``glycans``, ``results``, or as a top-level list.
    """
    payload: Any
    try:
        payload = json.loads(raw_output)
    except json.JSONDecodeError as exc:
        raise RuntimeError("Privateer output is not valid JSON") from exc

    records: list[dict[str, Any]] = []
    if isinstance(payload, list):
        records = [r for r in payload if isinstance(r, dict)]
    elif isinstance(payload, dict):
        for key in ("residues", "glycans", "results"):
            section = payload.get(key)
            if isinstance(section, list):
                records = [r for r in section if isinstance(r, dict)]
                break

    residue_results: list[PrivateerResidueResult] = []
    for rec in records:
        chain = str(rec.get("chain") or rec.get("chain_id") or "")
        resname = str(rec.get("resname") or rec.get("residue_name") or "")
        resnum = int(rec.get("resnum") or rec.get("residue_number") or 0)

        recognized = _as_bool(
            rec.get("sugar_recognized", rec.get("is_recognized", False))
        )
        anomer_ok = _as_bool(
            rec.get("anomer_ok", rec.get("anomer_correct", False))
        )
        ring_ok = _as_bool(
            rec.get("ring_pucker_ok", rec.get("puckering_ok", False))
        )
        linkage_ok = _as_bool(rec.get("linkage_ok", True))

        residue_results.append(
            PrivateerResidueResult(
                chain=chain,
                resname=resname,
                resnum=resnum,
                sugar_recognized=recognized,
                anomer_ok=anomer_ok,
                ring_pucker_ok=ring_ok,
                ring_pucker_conformation=str(rec.get("ring_pucker_conformation", "")),
                linkage_ok=linkage_ok,
                diagnostics=rec,
            )
        )

    return residue_results


def _parse_privateer_validation_data(raw_output: str) -> list[PrivateerResidueResult]:
    """Parse Privateer's ccp4i2 validation-data text artifact."""
    summary = PrivateerSummary()
    return [_row_to_residue_result(row, summary) for row in _parse_validation_rows(raw_output)]


def _parse_validation_rows(raw_output: str) -> list[dict[str, str]]:
    """Parse tabular validation_data-privateer rows into dicts."""
    rows: list[dict[str, str]] = []
    for raw_line in raw_output.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        parts = [part.strip() for part in raw_line.split("\t") if part.strip()]
        if len(parts) < 10:
            parts = line.split()
        if len(parts) < 10:
            continue

        rows.append(
            {
                "record": parts[0],
                "sugar_label": parts[1],
                "score": parts[2],
                "torsion": parts[3],
                "distance": parts[4],
                "sugar_type": parts[5],
                "pucker": parts[6],
                "rscc": parts[7],
                "linkage": parts[8],
                "ok": parts[9],
            }
        )
    return rows


def _parse_privateer_summary(raw_output: str) -> PrivateerSummary:
    """Parse the ccp4i2 textual summary counts from Privateer stdout."""

    def _match(pattern: str) -> int:
        matched = re.search(pattern, raw_output, flags=re.IGNORECASE)
        if not matched:
            return 0
        return int(matched.group(1))

    issues_match = re.search(
        r"Privateer has identified\s+(\d+)\s+issues,\s+with\s+(\d+)\s+of\s+(\d+)\s+sugars\s+affected",
        raw_output,
        flags=re.IGNORECASE,
    )

    return PrivateerSummary(
        wrong_anomer=_match(r"Wrong anomer:\s*(\d+)"),
        wrong_configuration=_match(r"Wrong configuration:\s*(\d+)"),
        unphysical_puckering=_match(r"Unphysical puckering amplitude:\s*(\d+)"),
        high_energy_conformations=_match(r"In higher-energy conformations:\s*(\d+)"),
        issues_detected=int(issues_match.group(1)) if issues_match else 0,
        affected_sugars=int(issues_match.group(2)) if issues_match else 0,
        total_sugars_reported=int(issues_match.group(3)) if issues_match else 0,
    )


def _row_to_residue_result(
    row: dict[str, str],
    summary: PrivateerSummary,
) -> PrivateerResidueResult:
    residue_id = row.get("sugar_label", "")
    resname = ""
    chain = ""
    resnum = 0
    try:
        resname, chain, resnum_txt = residue_id.split("-", 2)
        resnum = int(resnum_txt)
    except ValueError:
        pass

    ok = row.get("ok", "").strip().lower() == "yes"
    anomer_ok = ok or summary.wrong_anomer == 0
    ring_pucker_ok = ok or (
        summary.unphysical_puckering == 0 and summary.high_energy_conformations == 0
    )
    if not ok and summary.wrong_anomer == 0 and summary.unphysical_puckering == 0 and summary.high_energy_conformations == 0:
        anomer_ok = False
        ring_pucker_ok = False

    return PrivateerResidueResult(
        chain=chain,
        resname=resname,
        resnum=resnum,
        sugar_recognized=True,
        anomer_ok=anomer_ok,
        ring_pucker_ok=ring_pucker_ok and row.get("pucker", "") not in {"", "?", "unknown"},
        ring_pucker_conformation=row.get("pucker", ""),
        linkage_ok=True,
        diagnostics={**row, "source": "validation_data-privateer"},
    )


def _aggregate_privateer_results(
    residue_results: list[PrivateerResidueResult],
    pose_id: str = "",
    summary: PrivateerSummary | None = None,
    invocation: PrivateerInvocation | None = None,
    returncode: int = -1,
) -> PrivateerResult:
    """Aggregate per-residue records into the QC gate summary."""
    _ = summary
    total = len(residue_results)
    recognized = sum(1 for r in residue_results if r.sugar_recognized)
    anomer_pass = sum(1 for r in residue_results if r.anomer_ok)
    ring_pass = sum(1 for r in residue_results if r.ring_pucker_ok)
    linkage_pass = sum(1 for r in residue_results if r.linkage_ok)

    return PrivateerResult(
        pose_id=pose_id,
        residues=residue_results,
        total_sugars=total,
        recognized=recognized,
        recognition_rate=recognized / total if total > 0 else 0.0,
        anomer_pass=anomer_pass,
        ring_pucker_pass=ring_pass,
        linkage_pass=linkage_pass,
        all_pass=(
            recognized == total
            and anomer_pass == total
            and ring_pass == total
            and linkage_pass == total
        ),
        command=invocation.command if invocation is not None else [],
        work_dir=invocation.work_dir if invocation is not None else "",
        returncode=returncode,
    )


def _as_bool(value: Any) -> bool:
    """Coerce common JSON/string truthy values to bool."""
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "y", "ok", "pass"}
    return False


def write_privateer_report(
    result: PrivateerResult,
    output_path: Path,
) -> None:
    """Write Privateer results to JSON."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "total_sugars": result.total_sugars,
        "recognized": result.recognized,
        "recognition_rate": result.recognition_rate,
        "anomer_pass": result.anomer_pass,
        "ring_pucker_pass": result.ring_pucker_pass,
        "linkage_pass": result.linkage_pass,
        "all_pass": result.all_pass,
        "per_residue": [
            {
                "chain": r.chain,
                "resname": r.resname,
                "resnum": r.resnum,
                "sugar_recognized": r.sugar_recognized,
                "anomer_ok": r.anomer_ok,
                "ring_pucker_ok": r.ring_pucker_ok,
                "ring_pucker_conformation": r.ring_pucker_conformation,
                "linkage_ok": r.linkage_ok,
            }
            for r in result.residues
        ],
    }
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    logger.info("Wrote Privateer report to %s", output_path)
