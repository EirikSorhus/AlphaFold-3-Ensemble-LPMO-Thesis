# src/lpmo_pipeline/qc/privateer_runner.py
"""Run Privateer from its SIF image and parse validation artifacts.

This module intentionally treats the Privateer container as the source of truth.
The installed SIF on this cluster does not expose a JSON CLI. The reliable path is
``apptainer run --cleanenv ... -mode ccp4i2`` with explicit bind mounts, followed
by parsing the ``validation_data-privateer`` artifact written in the working
directory.
"""
from __future__ import annotations

import csv
import json
import logging
import os
import re
import shlex
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from lpmo_pipeline.io.ccd_lookup import (
    REJECTED_OLIGOMER_ALIASES,
    validate_all_glycan_residues,
    validate_comp_id,
)
from lpmo_pipeline.utils.exceptions import PrivateerCCDError

logger = logging.getLogger(__name__)


# Privateer must run from SIF, never from local env/PATH.
PRIVATEER_SIF_CANDIDATES: tuple[Path, ...] = (
    Path("/cluster/projects/nn1003k/prog/privateer/privateer.sif"),
)
PRIVATEER_RUNTIME = "apptainer"
PRIVATEER_DEFAULT_MODE = "ccp4i2"
PRIVATEER_VALIDATION_FILENAME = "validation_data-privateer"
PRIVATEER_SUPPORTED_MODES = frozenset({"normal", "ccp4i2"})

_PROTEIN_RESIDUE_NAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "MSE", "PHE", "PRO", "SER", "THR", "TRP", "TYR",
    "VAL",
}
_METAL_OR_SOLVENT_COMP_IDS = {"CU", "ZN", "FE", "MN", "CO", "NI", "MG", "CA", "HOH", "WAT"}


@dataclass(frozen=True)
class PrivateerInvocation:
    """Concrete execution plan for one Privateer run."""

    runtime: str
    sif_path: Path
    input_path: Path
    output_dir: Path
    mode: str
    bind_mounts: list[str]
    command: list[str]
    stdout_path: Path
    stderr_path: Path
    validation_data_path: Path

    @property
    def quoted_command(self) -> str:
        return shlex.join(self.command)


@dataclass(frozen=True)
class PrivateerBatchInput:
    """One queued Privateer run, optionally for batch execution."""

    cif_path: Path
    pose_id: str = ""
    output_dir: Path | None = None
    structure_factor_cif: Path | None = None
    structure_factor_mtz: Path | None = None
    colin_fo: str | None = None
    codein: str | None = None
    showgeom: bool = False
    radiusin: float | None = None
    mtzout: Path | None = None
    mode: str = PRIVATEER_DEFAULT_MODE
    dry_run: bool = False


@dataclass(frozen=True)
class _PrivateerSummary:
    wrong_anomer: int = 0
    wrong_configuration: int = 0
    unphysical_puckering: int = 0
    high_energy_conformations: int = 0
    issues_detected: int = 0
    affected_sugars: int = 0
    total_sugars_reported: int = 0


# ---------------------------------------------------------------------------
# Result models
# ---------------------------------------------------------------------------
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
    sugar_label: str = ""
    detected_type: str = ""
    quality_score: float | None = None
    phi: float | None = None
    theta: float | None = None
    context: str = ""
    privateer_ok: bool = True
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
    wrong_anomer_count: int = 0
    wrong_configuration_count: int = 0
    unphysical_puckering_count: int = 0
    high_energy_conformation_count: int = 0
    issues_detected: int = 0
    affected_sugars: int = 0
    all_pass: bool = False
    runner_error: str = ""
    runtime: str = ""
    sif_path: str = ""
    bind_mounts: list[str] = field(default_factory=list)
    command: list[str] = field(default_factory=list)
    work_dir: str = ""
    stdout_path: str = ""
    stderr_path: str = ""
    validation_data_path: str = ""
    returncode: int | None = None
    dry_run: bool = False


# ---------------------------------------------------------------------------
# Prep helpers
# ---------------------------------------------------------------------------
def prepare_privateer_input(
    input_path: Path,
    output_path: Path | None = None,
) -> Path:
    """Write a validated ``privateer_input.cif`` artifact.

    The current AF3-normalized inputs already use monomer CCD codes in the cases
    we care about, so the implemented prep step validates the glycan comp_ids and
    rewrites the structure as mmCIF under the canonical artifact name.
    """
    input_path = Path(input_path).resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Privateer input source does not exist: {input_path}")

    output_path = Path(output_path).resolve() if output_path else input_path.with_name("privateer_input.cif")
    comp_ids = [resname for _, resname, _ in _extract_nonpolymer_residues(input_path)]

    all_valid, results = validate_all_glycan_residues(comp_ids)
    if not all_valid:
        invalid = [result.comp_id for result in results if not result.is_valid_ccd_mono]
        raise PrivateerCCDError(
            f"Privateer prep rejected non-CCD glycan residues: {invalid}",
            step="privateer_prep",
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if input_path != output_path:
        shutil.copyfile(input_path, output_path)
    logger.info("Prepared Privateer input %s -> %s", input_path, output_path)
    return output_path


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------
def run_privateer(
    cif_path: Path,
    pose_id: str = "",
    *,
    output_dir: Path | None = None,
    structure_factor_cif: Path | None = None,
    structure_factor_mtz: Path | None = None,
    colin_fo: str | None = None,
    codein: str | None = None,
    showgeom: bool = False,
    radiusin: float | None = None,
    mtzout: Path | None = None,
    mode: str = PRIVATEER_DEFAULT_MODE,
    dry_run: bool = False,
    probe_container: bool = True,
    timeout: int = 300,
    runtime: str = PRIVATEER_RUNTIME,
) -> PrivateerResult:
    """Run Privateer on one input structure.

    Args:
        cif_path: Path to privateer_input.cif (with CCD monosaccharides).
        pose_id: Pose identifier for logging context.

    Returns:
        PrivateerResult with per-residue and aggregate metrics.
    """
    cif_path = Path(cif_path).resolve()
    resolved_output_dir = _resolve_privateer_output_dir(cif_path, pose_id, output_dir)

    try:
        invocation = build_privateer_invocation(
            cif_path=cif_path,
            output_dir=resolved_output_dir,
            structure_factor_cif=structure_factor_cif,
            structure_factor_mtz=structure_factor_mtz,
            colin_fo=colin_fo,
            codein=codein,
            showgeom=showgeom,
            radiusin=radiusin,
            mtzout=mtzout,
            mode=mode,
            runtime=runtime,
        )
    except Exception as exc:
        return _failed_privateer_result(pose_id, str(exc))

    logger.info(
        "Running Privateer on %s (pose=%s) with %s",
        cif_path,
        pose_id,
        invocation.quoted_command,
    )

    if dry_run:
        return PrivateerResult(
            pose_id=pose_id,
            all_pass=False,
            runtime=invocation.runtime,
            sif_path=str(invocation.sif_path),
            bind_mounts=list(invocation.bind_mounts),
            command=list(invocation.command),
            work_dir=str(invocation.output_dir),
            stdout_path=str(invocation.stdout_path),
            stderr_path=str(invocation.stderr_path),
            validation_data_path=str(invocation.validation_data_path),
            dry_run=True,
        )

    try:
        _preflight_privateer_host(
            invocation,
            structure_factor_cif=structure_factor_cif,
            structure_factor_mtz=structure_factor_mtz,
            mtzout=mtzout,
        )
        if probe_container:
            _probe_privateer_container(
                invocation,
                structure_factor_cif=structure_factor_cif,
                structure_factor_mtz=structure_factor_mtz,
            )

        stdout_text, stderr_text, returncode = _run_privateer_cli(invocation, timeout=timeout)
        validation_text = _read_validation_text(invocation.validation_data_path)
        summary = _parse_privateer_summary(stdout_text)
        expected_residues = _extract_expected_glycan_residues(cif_path)
        residue_results = _build_residue_results(
            validation_text=validation_text,
            expected_residues=expected_residues,
            summary=summary,
        )
        result = _aggregate_privateer_results(
            pose_id=pose_id,
            residue_results=residue_results,
            summary=summary,
            invocation=invocation,
            returncode=returncode,
        )
    except Exception as exc:
        logger.error("Privateer execution failed for %s: %s", cif_path, exc)
        return _failed_privateer_result(
            pose_id,
            str(exc),
            invocation=invocation,
        )

    if result.recognition_rate < 1.0:
        unrecognized = [f"{r.chain}:{r.resname}{r.resnum}" for r in residue_results if not r.sugar_recognized]
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


def run_privateer_batch(
    inputs: list[PrivateerBatchInput],
    *,
    max_workers: int | None = None,
    probe_container: bool = True,
) -> list[PrivateerResult]:
    """Run multiple Privateer jobs in parallel.

    The workload is subprocess-bound, so a thread pool is sufficient here.
    """
    if not inputs:
        return []

    worker_count = max_workers or min(len(inputs), os.cpu_count() or 1)
    with ThreadPoolExecutor(max_workers=max(1, worker_count)) as executor:
        futures = [
            executor.submit(
                run_privateer,
                item.cif_path,
                item.pose_id,
                output_dir=item.output_dir,
                structure_factor_cif=item.structure_factor_cif,
                structure_factor_mtz=item.structure_factor_mtz,
                colin_fo=item.colin_fo,
                codein=item.codein,
                showgeom=item.showgeom,
                radiusin=item.radiusin,
                mtzout=item.mtzout,
                mode=item.mode,
                dry_run=item.dry_run,
                probe_container=probe_container,
            )
            for item in inputs
        ]
        return [future.result() for future in futures]


def build_privateer_invocation(
    cif_path: Path,
    *,
    output_dir: Path,
    structure_factor_cif: Path | None = None,
    structure_factor_mtz: Path | None = None,
    colin_fo: str | None = None,
    codein: str | None = None,
    showgeom: bool = False,
    radiusin: float | None = None,
    mtzout: Path | None = None,
    mode: str = PRIVATEER_DEFAULT_MODE,
    runtime: str = PRIVATEER_RUNTIME,
) -> PrivateerInvocation:
    """Build the concrete Apptainer command and bind plan for Privateer."""
    if mode not in PRIVATEER_SUPPORTED_MODES:
        raise ValueError(f"Unsupported Privateer mode '{mode}'")
    if structure_factor_cif is not None and structure_factor_mtz is not None:
        raise ValueError("Provide either structure_factor_cif or structure_factor_mtz, not both")

    cif_path = Path(cif_path).resolve()
    output_dir = Path(output_dir).resolve()
    sif_path = _find_privateer_sif()

    bind_paths = [cif_path.parent, output_dir]
    if structure_factor_cif is not None:
        bind_paths.append(Path(structure_factor_cif).resolve().parent)
    if structure_factor_mtz is not None:
        bind_paths.append(Path(structure_factor_mtz).resolve().parent)
    if mtzout is not None:
        bind_paths.append(Path(mtzout).resolve().parent)

    bind_mounts = _build_bind_mounts(bind_paths)
    command = [runtime, "run", "--cleanenv"]
    for bind in bind_mounts:
        command.extend(["--bind", bind])
    command.extend([str(sif_path), "-pdbin", str(cif_path), "-mode", mode])

    if structure_factor_cif is not None:
        command.extend(["-cifin", str(Path(structure_factor_cif).resolve())])
    if structure_factor_mtz is not None:
        command.extend(["-mtzin", str(Path(structure_factor_mtz).resolve())])
    if colin_fo is not None:
        command.extend(["-colin-fo", colin_fo])
    if codein is not None:
        command.extend(["-codein", codein])
    if showgeom:
        command.append("-showgeom")
    if radiusin is not None:
        command.extend(["-radiusin", str(radiusin)])
    if mtzout is not None:
        command.extend(["-mtzout", str(Path(mtzout).resolve())])

    return PrivateerInvocation(
        runtime=runtime,
        sif_path=sif_path,
        input_path=cif_path,
        output_dir=output_dir,
        mode=mode,
        bind_mounts=bind_mounts,
        command=command,
        stdout_path=output_dir / "privateer_stdout.txt",
        stderr_path=output_dir / "privateer_stderr.txt",
        validation_data_path=output_dir / PRIVATEER_VALIDATION_FILENAME,
    )


def get_privateer_version(
    *,
    runtime: str = PRIVATEER_RUNTIME,
    timeout: int = 30,
) -> str:
    """Query the Privateer version from the SIF image itself."""
    sif_path = _find_privateer_sif()
    proc = subprocess.run(
        [runtime, "run", "--cleanenv", str(sif_path), "-list"],
        capture_output=True,
        text=True,
        check=False,
        timeout=timeout,
        env=_privateer_subprocess_env(),
    )
    text = f"{proc.stdout}\n{proc.stderr}"
    patterns = (
        r"Privateer Version\s+([A-Za-z0-9_.-]+)",
        r"CCP4\s+[0-9.]+:\s+Privateer\s+version\s+([A-Za-z0-9_.-]+)",
    )
    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            return match.group(1).rstrip(".")
    return "unknown"


def _run_privateer_cli(
    invocation: PrivateerInvocation,
    *,
    timeout: int,
) -> tuple[str, str, int]:
    """Run Privateer from SIF and persist stdout/stderr to the output directory."""
    invocation.output_dir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        invocation.command,
        capture_output=True,
        text=True,
        check=False,
        timeout=timeout,
        cwd=invocation.output_dir,
        env=_privateer_subprocess_env(),
    )
    invocation.stdout_path.write_text(proc.stdout)
    invocation.stderr_path.write_text(proc.stderr)

    if proc.returncode != 0:
        raise RuntimeError(
            f"Privateer failed with return code {proc.returncode}: {_summarize_privateer_failure(proc.stderr)}"
        )
    return proc.stdout, proc.stderr, proc.returncode


def _find_privateer_sif() -> Path:
    sif_path = next((p for p in PRIVATEER_SIF_CANDIDATES if p.exists()), None)
    if sif_path is None:
        raise RuntimeError(
            f"Privateer SIF not found in candidates: {PRIVATEER_SIF_CANDIDATES}"
        )
    return sif_path


def _resolve_privateer_output_dir(
    cif_path: Path,
    pose_id: str,
    output_dir: Path | None,
) -> Path:
    if output_dir is not None:
        return Path(output_dir).resolve()
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", pose_id or cif_path.stem)
    return cif_path.parent / "privateer_output" / label


def _build_bind_mounts(paths: list[Path]) -> list[str]:
    mounts: list[str] = []
    seen: set[Path] = set()
    for path in paths:
        resolved = Path(path).resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        mounts.append(f"{resolved}:{resolved}")
    return mounts


def _preflight_privateer_host(
    invocation: PrivateerInvocation,
    *,
    structure_factor_cif: Path | None,
    structure_factor_mtz: Path | None,
    mtzout: Path | None,
) -> None:
    runtime_path = shutil.which(invocation.runtime)
    if runtime_path is None:
        raise RuntimeError(f"Container runtime not found on PATH: {invocation.runtime}")

    if not invocation.sif_path.exists():
        raise RuntimeError(f"Privateer SIF does not exist: {invocation.sif_path}")
    if not invocation.input_path.exists():
        raise FileNotFoundError(f"Privateer input does not exist: {invocation.input_path}")
    if structure_factor_cif is not None and not Path(structure_factor_cif).exists():
        raise FileNotFoundError(f"Privateer structure-factor CIF does not exist: {structure_factor_cif}")
    if structure_factor_mtz is not None and not Path(structure_factor_mtz).exists():
        raise FileNotFoundError(f"Privateer MTZ does not exist: {structure_factor_mtz}")

    invocation.output_dir.mkdir(parents=True, exist_ok=True)
    if not os.access(invocation.output_dir, os.W_OK | os.X_OK):
        raise RuntimeError(f"Privateer output directory is not writable: {invocation.output_dir}")

    if mtzout is not None:
        Path(mtzout).resolve().parent.mkdir(parents=True, exist_ok=True)

    logger.debug(
        "Privateer preflight runtime=%s sif=%s binds=%s cwd=%s input=%s stdout=%s stderr=%s",
        runtime_path,
        invocation.sif_path,
        invocation.bind_mounts,
        invocation.output_dir,
        invocation.input_path,
        invocation.stdout_path,
        invocation.stderr_path,
    )


def _probe_privateer_container(
    invocation: PrivateerInvocation,
    *,
    structure_factor_cif: Path | None,
    structure_factor_mtz: Path | None,
) -> None:
    probe_paths = [str(invocation.input_path)]
    if structure_factor_cif is not None:
        probe_paths.append(str(Path(structure_factor_cif).resolve()))
    if structure_factor_mtz is not None:
        probe_paths.append(str(Path(structure_factor_mtz).resolve()))

    probe_command = [invocation.runtime, "exec"]
    for bind in invocation.bind_mounts:
        probe_command.extend(["--bind", bind])
    probe_command.extend(
        [
            str(invocation.sif_path),
            "/bin/sh",
            "-lc",
            (
                'command -v privateer >/dev/null && '
                'for p in "$@"; do test -r "$p" || exit 12; done && '
                'test -d "$PRIVATEER_OUTDIR" && test -w "$PRIVATEER_OUTDIR"'
            ),
            "_",
            *probe_paths,
        ]
    )

    probe_env = dict(_privateer_subprocess_env())
    probe_env["PRIVATEER_OUTDIR"] = str(invocation.output_dir)
    proc = subprocess.run(
        probe_command,
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
        env=probe_env,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            "Privateer container preflight failed: "
            f"{_summarize_privateer_failure(proc.stderr or proc.stdout)}"
        )


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


def _read_validation_text(validation_path: Path) -> str:
    if not validation_path.exists():
        raise RuntimeError(
            f"Privateer did not write {PRIVATEER_VALIDATION_FILENAME} in {validation_path.parent}"
        )
    validation_text = validation_path.read_text().strip()
    if not validation_text:
        raise RuntimeError(f"Privateer wrote an empty {PRIVATEER_VALIDATION_FILENAME} file")
    return validation_text


def _parse_privateer_summary(stdout_text: str) -> _PrivateerSummary:
    patterns = {
        "wrong_anomer": r"Wrong anomer:\s*(\d+)",
        "wrong_configuration": r"Wrong configuration:\s*(\d+)",
        "unphysical_puckering": r"Unphysical puckering amplitude:\s*(\d+)",
        "high_energy_conformations": r"In higher-energy conformations:\s*(\d+)",
        "issues_detected": r"Privateer has identified\s*(\d+)\s*issues",
        "affected_sugars": r"with\s*(\d+)\s*of\s*(\d+)\s*sugars affected",
    }

    values: dict[str, int] = {
        "wrong_anomer": 0,
        "wrong_configuration": 0,
        "unphysical_puckering": 0,
        "high_energy_conformations": 0,
        "issues_detected": 0,
        "affected_sugars": 0,
        "total_sugars_reported": 0,
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, stdout_text)
        if not match:
            continue
        if key == "affected_sugars":
            values["affected_sugars"] = int(match.group(1))
            values["total_sugars_reported"] = int(match.group(2))
        else:
            values[key] = int(match.group(1))

    return _PrivateerSummary(**values)


def _build_residue_results(
    *,
    validation_text: str,
    expected_residues: list[tuple[str, str, int]],
    summary: _PrivateerSummary,
) -> list[PrivateerResidueResult]:
    rows = _parse_validation_rows(validation_text)
    parsed_rows = {(_row_chain(row), _row_resname(row), _row_resnum(row)): row for row in rows}

    residue_results: list[PrivateerResidueResult] = []
    for chain, resname, resnum in expected_residues:
        row = parsed_rows.pop((chain, resname, resnum), None)
        if row is None:
            residue_results.append(
                PrivateerResidueResult(
                    chain=chain,
                    resname=resname,
                    resnum=resnum,
                    sugar_recognized=False,
                    anomer_ok=False,
                    ring_pucker_ok=False,
                    linkage_ok=False,
                    sugar_label=f"{resname}-{chain}-{resnum}",
                    privateer_ok=False,
                    diagnostics={"reason": "missing_from_privateer_validation_data"},
                )
            )
            continue

        residue_results.append(_row_to_residue_result(row, summary))

    for row in parsed_rows.values():
        residue_results.append(_row_to_residue_result(row, summary))

    return residue_results


def _parse_validation_rows(validation_text: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    reader = csv.reader(validation_text.splitlines(), delimiter="\t")
    for raw_fields in reader:
        fields = [field.strip() for field in raw_fields if field.strip()]
        if not fields:
            continue
        if len(fields) < 10:
            raise RuntimeError(f"Unexpected Privateer validation row format: {raw_fields}")
        rows.append(
            {
                "pdb_code": fields[0],
                "sugar_label": fields[1],
                "quality_score": fields[2],
                "phi": fields[3],
                "theta": fields[4],
                "detected_type": fields[5],
                "conformation": fields[6],
                "bfactor": fields[7],
                "context": fields[8],
                "ok": fields[9].lower(),
            }
        )
    return rows


def _row_to_residue_result(
    row: dict[str, str],
    summary: _PrivateerSummary,
) -> PrivateerResidueResult:
    privateer_ok = row["ok"] == "yes"
    has_anomer_issue = (summary.wrong_anomer > 0 or summary.wrong_configuration > 0) and not privateer_ok
    has_ring_issue = (
        summary.unphysical_puckering > 0 or summary.high_energy_conformations > 0
    ) and not privateer_ok
    linkage_ok = privateer_ok or (has_anomer_issue or has_ring_issue)

    return PrivateerResidueResult(
        chain=_row_chain(row),
        resname=_row_resname(row),
        resnum=_row_resnum(row),
        sugar_recognized=True,
        anomer_ok=not has_anomer_issue,
        ring_pucker_ok=not has_ring_issue,
        ring_pucker_conformation=row["conformation"],
        linkage_ok=linkage_ok,
        sugar_label=row["sugar_label"],
        detected_type=row["detected_type"],
        quality_score=_as_float(row["quality_score"]),
        phi=_as_float(row["phi"]),
        theta=_as_float(row["theta"]),
        context=row["context"],
        privateer_ok=privateer_ok,
        diagnostics=dict(row),
    )


def _row_chain(row: dict[str, str]) -> str:
    return row["sugar_label"].split("-")[1]


def _row_resname(row: dict[str, str]) -> str:
    return row["sugar_label"].split("-")[0]


def _row_resnum(row: dict[str, str]) -> int:
    return int(row["sugar_label"].split("-")[2])


def _extract_expected_glycan_residues(cif_path: Path) -> list[tuple[str, str, int]]:
    return [
        residue
        for residue in _extract_nonpolymer_residues(cif_path)
        if _is_candidate_glycan_residue(residue[1])
    ]


def _extract_nonpolymer_residues(cif_path: Path) -> list[tuple[str, str, int]]:
    text = cif_path.read_text()
    if _looks_like_pdb(text):
        return _extract_nonpolymer_residues_from_pdb_text(text)
    return _extract_nonpolymer_residues_from_mmcif_text(text)


def _looks_like_pdb(text: str) -> bool:
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped == "#":
            continue
        if stripped.startswith(("data_", "loop_", "_")):
            return False
        return line.startswith(("ATOM  ", "HETATM"))
    return False


def _extract_nonpolymer_residues_from_pdb_text(text: str) -> list[tuple[str, str, int]]:
    residues: dict[tuple[str, str, int], None] = {}
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        resname = line[17:20].strip().upper()
        if _should_skip_residue_name(resname):
            continue
        chain = line[21].strip() or "?"
        try:
            resnum = int(line[22:26].strip())
        except ValueError:
            continue
        residues[(chain, resname, resnum)] = None
    return list(residues)


def _extract_nonpolymer_residues_from_mmcif_text(text: str) -> list[tuple[str, str, int]]:
    residues: dict[tuple[str, str, int], None] = {}
    lines = text.splitlines()
    index = 0
    while index < len(lines):
        if lines[index].strip() != "loop_":
            index += 1
            continue

        index += 1
        tags: list[str] = []
        while index < len(lines) and lines[index].strip().startswith("_"):
            tags.append(lines[index].strip())
            index += 1

        if not tags or not all(tag.startswith("_atom_site.") for tag in tags):
            continue

        comp_idx = _find_atom_site_tag(tags, "_atom_site.auth_comp_id", "_atom_site.label_comp_id")
        chain_idx = _find_atom_site_tag(tags, "_atom_site.auth_asym_id", "_atom_site.label_asym_id")
        resnum_idx = _find_atom_site_tag(tags, "_atom_site.auth_seq_id", "_atom_site.label_seq_id")

        if comp_idx is None or chain_idx is None or resnum_idx is None:
            continue

        while index < len(lines):
            stripped = lines[index].strip()
            if not stripped:
                index += 1
                continue
            if stripped == "#":
                index += 1
                break
            if stripped == "loop_" or stripped.startswith("_"):
                break

            fields = shlex.split(lines[index], posix=True)
            index += 1
            if len(fields) < len(tags):
                continue

            resname = fields[comp_idx].strip().upper()
            if _should_skip_residue_name(resname):
                continue
            chain = fields[chain_idx].strip() or "?"
            try:
                resnum = int(float(fields[resnum_idx]))
            except ValueError:
                continue
            residues[(chain, resname, resnum)] = None
    return list(residues)


def _find_atom_site_tag(tags: list[str], *candidates: str) -> int | None:
    for candidate in candidates:
        if candidate in tags:
            return tags.index(candidate)
    return None


def _is_candidate_glycan_residue(resname: str) -> bool:
    if _should_skip_residue_name(resname):
        return False
    lookup = validate_comp_id(resname)
    return lookup.is_valid_ccd_mono or resname in REJECTED_OLIGOMER_ALIASES


def _should_skip_residue_name(resname: str) -> bool:
    return resname in _PROTEIN_RESIDUE_NAMES or resname in _METAL_OR_SOLVENT_COMP_IDS


def _aggregate_privateer_results(
    *,
    pose_id: str,
    residue_results: list[PrivateerResidueResult],
    summary: _PrivateerSummary,
    invocation: PrivateerInvocation,
    returncode: int,
) -> PrivateerResult:
    """Aggregate per-residue records into the QC gate summary."""
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
        wrong_anomer_count=summary.wrong_anomer,
        wrong_configuration_count=summary.wrong_configuration,
        unphysical_puckering_count=summary.unphysical_puckering,
        high_energy_conformation_count=summary.high_energy_conformations,
        issues_detected=summary.issues_detected,
        affected_sugars=summary.affected_sugars,
        all_pass=(
            total == 0
            or (
                recognized == total
                and anomer_pass == total
                and ring_pass == total
                and linkage_pass == total
            )
        ),
        runtime=invocation.runtime,
        sif_path=str(invocation.sif_path),
        bind_mounts=list(invocation.bind_mounts),
        command=list(invocation.command),
        work_dir=str(invocation.output_dir),
        stdout_path=str(invocation.stdout_path),
        stderr_path=str(invocation.stderr_path),
        validation_data_path=str(invocation.validation_data_path),
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


def _as_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _failed_privateer_result(
    pose_id: str,
    runner_error: str,
    *,
    invocation: PrivateerInvocation | None = None,
) -> PrivateerResult:
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
        runner_error=runner_error,
        runtime=invocation.runtime if invocation else "",
        sif_path=str(invocation.sif_path) if invocation else "",
        bind_mounts=list(invocation.bind_mounts) if invocation else [],
        command=list(invocation.command) if invocation else [],
        work_dir=str(invocation.output_dir) if invocation else "",
        stdout_path=str(invocation.stdout_path) if invocation else "",
        stderr_path=str(invocation.stderr_path) if invocation else "",
        validation_data_path=str(invocation.validation_data_path) if invocation else "",
        returncode=None,
    )


def _summarize_privateer_failure(text: str) -> str:
    cleaned = " ".join(line.strip() for line in text.splitlines() if line.strip())
    return cleaned[:400] if cleaned else "no stderr captured"


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
        "wrong_anomer_count": result.wrong_anomer_count,
        "wrong_configuration_count": result.wrong_configuration_count,
        "unphysical_puckering_count": result.unphysical_puckering_count,
        "high_energy_conformation_count": result.high_energy_conformation_count,
        "issues_detected": result.issues_detected,
        "affected_sugars": result.affected_sugars,
        "all_pass": result.all_pass,
        "runner_error": result.runner_error,
        "runtime": result.runtime,
        "sif_path": result.sif_path,
        "bind_mounts": result.bind_mounts,
        "command": result.command,
        "work_dir": result.work_dir,
        "stdout_path": result.stdout_path,
        "stderr_path": result.stderr_path,
        "validation_data_path": result.validation_data_path,
        "returncode": result.returncode,
        "dry_run": result.dry_run,
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
                "sugar_label": r.sugar_label,
                "detected_type": r.detected_type,
                "quality_score": r.quality_score,
                "phi": r.phi,
                "theta": r.theta,
                "context": r.context,
                "privateer_ok": r.privateer_ok,
            }
            for r in result.residues
        ],
    }
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    logger.info("Wrote Privateer report to %s", output_path)
