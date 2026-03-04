# src/lpmo_pipeline/qc/custom_geometry_checks.py
"""
Responsibility: LPMO-specific active-site geometry validation.
Input:  Structure (gemmi or MDAnalysis Universe)
Output: GeometryResult with Cu–His distances, Cu–C1/C4 distances, angles

GATES:
  - Cu–His distance must be 1.9–2.6 Å (hard fail)
  - Cu–C1, Cu–C4 measured and flagged if > 7 Å (soft flag)

INVARIANT: Atom names not assumed consistent across models.
           Use element + residue context to identify Cu and His-brace atoms.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Thresholds (from configs/thresholds.yaml)
# ---------------------------------------------------------------------------
CU_HIS_MIN: float = 1.9   # Å
CU_HIS_MAX: float = 2.6   # Å
CU_SUBSTRATE_FLAG: float = 7.0  # Å — flag if Cu–C1/C4 > this


@dataclass
class CuHisMeasurement:
    """A single Cu–His distance measurement."""

    his_chain: str
    his_resnum: int
    his_atom: str       # e.g. NE2 or ND1
    cu_chain: str
    cu_resnum: int
    distance_angstrom: float
    in_range: bool = False


@dataclass
class CuSubstrateMeasurement:
    """Cu distance to a specific substrate carbon."""

    carbon_type: str    # "C1" or "C4"
    chain: str
    resnum: int
    atom_name: str
    distance_angstrom: float
    flagged: bool = False  # True if > CU_SUBSTRATE_FLAG


@dataclass
class GeometryResult:
    """Full geometry validation result."""

    pose_id: str = ""
    cu_found: bool = False
    cu_position: tuple[float, float, float] = (0.0, 0.0, 0.0)

    cu_his_measurements: list[CuHisMeasurement] = field(default_factory=list)
    cu_his_all_in_range: bool = False

    cu_c1_measurements: list[CuSubstrateMeasurement] = field(default_factory=list)
    cu_c4_measurements: list[CuSubstrateMeasurement] = field(default_factory=list)
    min_cu_c1: float = float("inf")
    min_cu_c4: float = float("inf")

    his_brace_angle: float | None = None  # Angle between His-N–Cu–His-N plane
    active_site_orientation: dict[str, Any] = field(default_factory=dict)

    passed: bool = False
    failure_reasons: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------
def check_geometry(
    structure: Any,  # gemmi.Structure or MDAnalysis.Universe
    pose_id: str = "",
    cu_chain: str = "E",
    glycan_chains: list[str] | None = None,
) -> GeometryResult:
    """Run all LPMO-specific geometry checks.

    Args:
        structure: Parsed structure (gemmi preferred).
        pose_id: Identifier for logging.
        cu_chain: Expected chain for Cu (default "E" per normalization).
        glycan_chains: Expected glycan chains (default ["B", "C", "D"]).

    Returns:
        GeometryResult with all measurements and pass/fail.
    """
    if glycan_chains is None:
        glycan_chains = ["B", "C", "D"]

    result = GeometryResult(pose_id=pose_id)

    # --- Step 1: Find Cu atom ---
    cu_atom = _find_cu_atom(structure, cu_chain)
    if cu_atom is None:
        result.failure_reasons.append("Cu atom not found in chain " + cu_chain)
        logger.error("GEOMETRY: Cu not found in chain %s for pose %s", cu_chain, pose_id)
        return result

    result.cu_found = True
    cu_pos = np.array([cu_atom.pos.x, cu_atom.pos.y, cu_atom.pos.z])
    result.cu_position = tuple(cu_pos.tolist())

    # --- Step 2: Cu–His distances ---
    his_residues = _find_his_brace_residues(structure, protein_chain="A")
    for his_res in his_residues:
        for atom_name in ["NE2", "ND1"]:
            his_atom = _get_atom_by_name_or_element(his_res, atom_name, "N")
            if his_atom is None:
                continue
            his_pos = np.array([his_atom.pos.x, his_atom.pos.y, his_atom.pos.z])
            dist = float(np.linalg.norm(cu_pos - his_pos))
            in_range = CU_HIS_MIN <= dist <= CU_HIS_MAX

            result.cu_his_measurements.append(CuHisMeasurement(
                his_chain="A",
                his_resnum=his_res.seqid.num,
                his_atom=atom_name,
                cu_chain=cu_chain,
                cu_resnum=cu_atom.residue_seqid if hasattr(cu_atom, 'residue_seqid') else 1,
                distance_angstrom=dist,
                in_range=in_range,
            ))

    result.cu_his_all_in_range = all(m.in_range for m in result.cu_his_measurements)
    if not result.cu_his_all_in_range:
        bad = [
            f"His{m.his_resnum}:{m.his_atom}={m.distance_angstrom:.2f}Å"
            for m in result.cu_his_measurements
            if not m.in_range
        ]
        result.failure_reasons.append(f"Cu-His distance out of range: {bad}")
        logger.error("GATE FAIL: Cu-His distance out of 1.9–2.6 Å: %s (pose %s)", bad, pose_id)

    # --- Step 3: Cu–C1/C4 distances ---
    for gly_chain_name in glycan_chains:
        gly_chain = _get_chain(structure, gly_chain_name)
        if gly_chain is None:
            continue

        for residue in gly_chain:
            # Find C1 atom (element C, bonded context ~ ring C1)
            c1 = _get_atom_by_name_or_element(residue, "C1", "C")
            if c1 is not None:
                c1_pos = np.array([c1.pos.x, c1.pos.y, c1.pos.z])
                d_c1 = float(np.linalg.norm(cu_pos - c1_pos))
                result.cu_c1_measurements.append(CuSubstrateMeasurement(
                    carbon_type="C1",
                    chain=gly_chain_name,
                    resnum=residue.seqid.num,
                    atom_name=c1.name,
                    distance_angstrom=d_c1,
                    flagged=d_c1 > CU_SUBSTRATE_FLAG,
                ))
                result.min_cu_c1 = min(result.min_cu_c1, d_c1)

            # Find C4 atom
            c4 = _get_atom_by_name_or_element(residue, "C4", "C")
            if c4 is not None:
                c4_pos = np.array([c4.pos.x, c4.pos.y, c4.pos.z])
                d_c4 = float(np.linalg.norm(cu_pos - c4_pos))
                result.cu_c4_measurements.append(CuSubstrateMeasurement(
                    carbon_type="C4",
                    chain=gly_chain_name,
                    resnum=residue.seqid.num,
                    atom_name=c4.name,
                    distance_angstrom=d_c4,
                    flagged=d_c4 > CU_SUBSTRATE_FLAG,
                ))
                result.min_cu_c4 = min(result.min_cu_c4, d_c4)

    # --- Step 4: His-brace angle (optional) ---
    if len(result.cu_his_measurements) >= 2:
        result.his_brace_angle = _compute_his_brace_angle(
            result.cu_his_measurements, cu_pos
        )

    # --- Step 5: Final verdict ---
    result.passed = (
        result.cu_found
        and result.cu_his_all_in_range
        and len(result.failure_reasons) == 0
    )

    logger.info(
        "Geometry check pose=%s: passed=%s, Cu-His OK=%s, "
        "min Cu-C1=%.2f, min Cu-C4=%.2f",
        pose_id, result.passed, result.cu_his_all_in_range,
        result.min_cu_c1, result.min_cu_c4,
    )
    return result


# ---------------------------------------------------------------------------
# Helpers (PSEUDOCODE — adapt to actual gemmi API)
# ---------------------------------------------------------------------------
def _find_cu_atom(structure: Any, chain_name: str) -> Any | None:
    """Find Cu atom in the specified chain."""
    model = structure[0]
    for chain in model:
        if chain.name != chain_name:
            continue
        for residue in chain:
            for atom in residue:
                if atom.element.name == "Cu" or atom.element.name == "CU":
                    return atom
    return None


def _find_his_brace_residues(structure: Any, protein_chain: str = "A") -> list[Any]:
    """Find histidine residues that form the His-brace (Cu-coordinating).

    In LPMOs, the His-brace typically consists of:
    - His1 (N-terminal His, coordinates Cu via amino-N and ND1)
    - A second His (coordinates Cu via NE2)

    Returns list of gemmi.Residue objects for His residues.
    """
    model = structure[0]
    his_residues = []
    for chain in model:
        if chain.name != protein_chain:
            continue
        for residue in chain:
            if residue.name == "HIS":
                his_residues.append(residue)
    return his_residues


def _get_chain(structure: Any, chain_name: str) -> Any | None:
    """Get a chain by name from the first model."""
    model = structure[0]
    for chain in model:
        if chain.name == chain_name:
            return chain
    return None


def _get_atom_by_name_or_element(
    residue: Any, preferred_name: str, element: str,
) -> Any | None:
    """Get atom by preferred name; fallback to element match.

    INVARIANT: Atom names not assumed consistent across models.
    """
    # Try exact name first
    for atom in residue:
        if atom.name.strip() == preferred_name:
            return atom
    # Fallback: first atom of matching element
    for atom in residue:
        if atom.element.name == element:
            return atom
    return None


def _compute_his_brace_angle(
    measurements: list[CuHisMeasurement],
    cu_pos: np.ndarray,
) -> float:
    """Compute angle between the two His-N–Cu vectors (His-brace angle).

    Returns angle in degrees.
    """
    # PSEUDOCODE: get the two His-N positions, compute CuN1-Cu-CuN2 angle
    # For now, return placeholder
    return 0.0
