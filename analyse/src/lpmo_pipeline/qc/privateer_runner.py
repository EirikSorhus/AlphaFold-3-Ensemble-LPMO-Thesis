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
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


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
    diagnostics: dict[str, Any] = field(default_factory=dict)


@dataclass
class PrivateerResult:
    """Aggregated Privateer result for a structure."""

    residues: list[PrivateerResidueResult] = field(default_factory=list)
    total_sugars: int = 0
    recognized: int = 0
    recognition_rate: float = 0.0
    anomer_pass: int = 0
    ring_pucker_pass: int = 0
    linkage_pass: int = 0
    all_pass: bool = False
    runner_error: str = ""


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------
def run_privateer(
    cif_path: Path,
    pose_id: str = "",
) -> PrivateerResult:
    """Run Privateer on a CIF file and parse results.

    PSEUDOCODE — calls Privateer CLI or Python API.

    Args:
        cif_path: Path to privateer_input.cif (with CCD monosaccharides).
        pose_id: Pose identifier for logging context.

    Returns:
        PrivateerResult with per-residue and aggregate metrics.
    """
    logger.info("Running Privateer on %s (pose=%s)", cif_path, pose_id)

    try:
        raw = _run_privateer_cli(cif_path)
        residue_results = _parse_privateer_output(raw)
        result = _aggregate_privateer_results(residue_results)
    except Exception as exc:
        logger.error("Privateer execution failed for %s: %s", cif_path, exc)
        # Conservative behavior: failed run is treated as failed gate.
        return PrivateerResult(
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

    if anomer_pass < total:
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
    """Run Privateer CLI and return stdout.

    Environment override:
      - ``LPMO_PRIVATEER_CMD`` can define a full command template containing
        ``{input}``, for example:
        ``privateer validate --input {input} --json``.
    """
    cmd_template = os.environ.get("LPMO_PRIVATEER_CMD", "privateer --json {input}")
    rendered = cmd_template.replace("{input}", str(cif_path))
    cmd = rendered.split()

    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True,
        timeout=300,
    )
    if not proc.stdout.strip():
        raise RuntimeError("Privateer produced empty stdout")
    return proc.stdout


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


def _aggregate_privateer_results(
    residue_results: list[PrivateerResidueResult],
) -> PrivateerResult:
    """Aggregate per-residue records into the QC gate summary."""
    total = len(residue_results)
    recognized = sum(1 for r in residue_results if r.sugar_recognized)
    anomer_pass = sum(1 for r in residue_results if r.anomer_ok)
    ring_pass = sum(1 for r in residue_results if r.ring_pucker_ok)
    linkage_pass = sum(1 for r in residue_results if r.linkage_ok)

    return PrivateerResult(
        residues=residue_results,
        total_sugars=total,
        recognized=recognized,
        recognition_rate=recognized / total if total > 0 else 0.0,
        anomer_pass=anomer_pass,
        ring_pucker_pass=ring_pass,
        linkage_pass=linkage_pass,
        all_pass=(recognized == total and anomer_pass == total),
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
