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

    # --- Step 1: Call Privateer ---
    # Option A: CLI
    # cmd = ["privateer", "-pdbin", str(cif_path), "-mode", "glycan_validation"]
    # proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    # raw_output = proc.stdout
    #
    # Option B: Python API (privateer module)
    # import privateer
    # results = privateer.validate_glycans(str(cif_path))

    # --- Step 2: Parse Privateer output ---
    # PSEUDOCODE: Privateer outputs per-sugar diagnostics
    # For each sugar residue in the output:
    residue_results: list[PrivateerResidueResult] = []

    # PSEUDOCODE: iterate parsed output
    # for sugar_entry in privateer_parsed_output:
    #     rr = PrivateerResidueResult(
    #         chain=sugar_entry["chain"],
    #         resname=sugar_entry["resname"],
    #         resnum=sugar_entry["resnum"],
    #         sugar_recognized=sugar_entry["is_recognized"],
    #         anomer_ok=sugar_entry["anomer_correct"],
    #         ring_pucker_ok=sugar_entry["ring_pucker_ok"],
    #         ring_pucker_conformation=sugar_entry.get("conformation", ""),
    #         linkage_ok=sugar_entry.get("linkage_ok", True),
    #         diagnostics=sugar_entry,
    #     )
    #     residue_results.append(rr)

    # --- Step 3: Aggregate ---
    total = len(residue_results)
    recognized = sum(1 for r in residue_results if r.sugar_recognized)
    anomer_pass = sum(1 for r in residue_results if r.anomer_ok)
    ring_pass = sum(1 for r in residue_results if r.ring_pucker_ok)
    linkage_pass = sum(1 for r in residue_results if r.linkage_ok)

    result = PrivateerResult(
        residues=residue_results,
        total_sugars=total,
        recognized=recognized,
        recognition_rate=recognized / total if total > 0 else 0.0,
        anomer_pass=anomer_pass,
        ring_pucker_pass=ring_pass,
        linkage_pass=linkage_pass,
        all_pass=(recognized == total and anomer_pass == total),
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
        recognized, total, anomer_pass, total, ring_pass, total,
    )
    return result


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
