# src/lpmo_pipeline/qc/active_site_proximity.py
"""
Responsibility: Pre-QC active-site proximity gate before PoseBusters/Privateer.
Input:  Structure object (gemmi-like)
Output: ActiveSiteProximityResult

Hard gate:
  - min Cu-ligand heavy-atom distance <= ACTIVE_SITE_PROXIMITY_MAX_A

Important:
  - Keep all computed distances as metadata, even when gate fails.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from lpmo_pipeline.config import load_defaults_config

logger = logging.getLogger(__name__)

_DEFAULTS_CONFIG = load_defaults_config()
_CHAIN_SCHEMA = _DEFAULTS_CONFIG.get("chain_schema") or {}
DEFAULT_CU_CHAIN = str(_CHAIN_SCHEMA.get("metal") or "E")
DEFAULT_GLYCAN_CHAINS = tuple(str(chain) for chain in (_CHAIN_SCHEMA.get("glycans") or ["B", "C", "D"]))


ACTIVE_SITE_PROXIMITY_MAX_A: float = 10.0
CU_C_SOFT_FLAG_A: float = 7.0


def _atom_position(atom: Any) -> np.ndarray:
    """Return atom coordinates for native gemmi and compat fallback atoms."""
    pos = getattr(atom, "pos", None)
    if pos is not None and all(hasattr(pos, axis) for axis in ("x", "y", "z")):
        return np.array([float(pos.x), float(pos.y), float(pos.z)], dtype=float)
    if all(hasattr(atom, axis) for axis in ("x", "y", "z")):
        return np.array([float(atom.x), float(atom.y), float(atom.z)], dtype=float)
    raise AttributeError("Atom has no compatible coordinate attributes")


@dataclass
class ActiveSiteProximityResult:
    """Result of active-site proximity pre-check for one pose."""

    pose_id: str = ""
    cu_found: bool = False
    min_cu_ligand_distance: float = float("inf")
    min_cu_c1: float = float("inf")
    min_cu_c4: float = float("inf")
    nearest_ligand_atom: str = ""
    nearest_c1_atom: str = ""
    nearest_c4_atom: str = ""
    passed: bool = False
    failure_reasons: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def check_active_site_proximity(
    structure: Any,
    pose_id: str = "",
    cu_chain: str = DEFAULT_CU_CHAIN,
    glycan_chains: list[str] | None = None,
    hard_cutoff_a: float = ACTIVE_SITE_PROXIMITY_MAX_A,
    cu_c_soft_flag_a: float = CU_C_SOFT_FLAG_A,
) -> ActiveSiteProximityResult:
    """Run pre-QC active-site proximity check.

    Args:
        structure: Parsed structure (gemmi preferred).
        pose_id: Pose identifier for logs.
        cu_chain: Expected Cu chain (default E).
        glycan_chains: Ligand chains to evaluate (default B/C/D).
        hard_cutoff_a: Hard gate threshold for nearest Cu-ligand distance.
        cu_c_soft_flag_a: Soft warning threshold for nearest Cu-C1/C4 distance.

    Returns:
        ActiveSiteProximityResult with metrics and pass/fail.
    """
    if glycan_chains is None:
        glycan_chains = list(DEFAULT_GLYCAN_CHAINS)

    result = ActiveSiteProximityResult(pose_id=pose_id)

    cu_atom = _find_cu_atom(structure, chain_name=cu_chain)
    if cu_atom is None:
        result.failure_reasons.append(f"Cu atom not found in chain {cu_chain}")
        logger.error("PRE-QC FAIL pose=%s: Cu atom not found", pose_id)
        return result

    result.cu_found = True
    cu_pos = _atom_position(cu_atom)

    for gly_chain_name in glycan_chains:
        gly_chain = _get_chain(structure, gly_chain_name)
        if gly_chain is None:
            continue

        for residue in gly_chain:
            for atom in residue:
                if atom.element.name == "H":
                    continue
                atom_pos = _atom_position(atom)
                dist = float(np.linalg.norm(cu_pos - atom_pos))
                if dist < result.min_cu_ligand_distance:
                    result.min_cu_ligand_distance = dist
                    result.nearest_ligand_atom = (
                        f"{gly_chain_name}:{residue.name}{residue.seqid.num}:{atom.name.strip()}"
                    )

            c1_atom = _get_atom_by_name_or_element(residue, preferred_name="C1", element="C")
            if c1_atom is not None:
                c1_pos = _atom_position(c1_atom)
                c1_dist = float(np.linalg.norm(cu_pos - c1_pos))
                if c1_dist < result.min_cu_c1:
                    result.min_cu_c1 = c1_dist
                    result.nearest_c1_atom = (
                        f"{gly_chain_name}:{residue.name}{residue.seqid.num}:{c1_atom.name.strip()}"
                    )

            c4_atom = _get_atom_by_name_or_element(residue, preferred_name="C4", element="C")
            if c4_atom is not None:
                c4_pos = _atom_position(c4_atom)
                c4_dist = float(np.linalg.norm(cu_pos - c4_pos))
                if c4_dist < result.min_cu_c4:
                    result.min_cu_c4 = c4_dist
                    result.nearest_c4_atom = (
                        f"{gly_chain_name}:{residue.name}{residue.seqid.num}:{c4_atom.name.strip()}"
                    )

    if result.min_cu_ligand_distance == float("inf"):
        result.failure_reasons.append("No ligand heavy atoms found in glycan chains")
        logger.error("PRE-QC FAIL pose=%s: No ligand atoms found", pose_id)
        return result

    if result.min_cu_ligand_distance > hard_cutoff_a:
        result.failure_reasons.append(
            "Ligand too far from active site: "
            f"min Cu-ligand={result.min_cu_ligand_distance:.2f} A > {hard_cutoff_a:.2f} A"
        )

    if result.min_cu_c1 != float("inf") and result.min_cu_c1 > cu_c_soft_flag_a:
        result.warnings.append(
            f"Cu-C1 distance flagged: {result.min_cu_c1:.2f} A > {cu_c_soft_flag_a:.2f} A"
        )

    if result.min_cu_c4 != float("inf") and result.min_cu_c4 > cu_c_soft_flag_a:
        result.warnings.append(
            f"Cu-C4 distance flagged: {result.min_cu_c4:.2f} A > {cu_c_soft_flag_a:.2f} A"
        )

    result.passed = len(result.failure_reasons) == 0
    logger.info(
        "PRE-QC pose=%s passed=%s min_cu_ligand=%.2f min_cu_c1=%.2f min_cu_c4=%.2f",
        pose_id,
        result.passed,
        result.min_cu_ligand_distance,
        result.min_cu_c1,
        result.min_cu_c4,
    )
    return result


def _find_cu_atom(structure: Any, chain_name: str) -> Any | None:
    """Find Cu atom in the expected chain."""
    model = structure[0]
    for chain in model:
        if chain.name != chain_name:
            continue
        for residue in chain:
            for atom in residue:
                if atom.element.name in {"CU", "Cu"}:
                    return atom
    return None


def _get_chain(structure: Any, chain_name: str) -> Any | None:
    """Get chain by name from first model."""
    model = structure[0]
    for chain in model:
        if chain.name == chain_name:
            return chain
    return None


def _get_atom_by_name_or_element(residue: Any, preferred_name: str, element: str) -> Any | None:
    """Get atom by preferred name; fallback to element match."""
    for atom in residue:
        if atom.name.strip() == preferred_name:
            return atom
    for atom in residue:
        if atom.element.name == element:
            return atom
    return None
