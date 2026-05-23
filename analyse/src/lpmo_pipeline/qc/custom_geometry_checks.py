# src/lpmo_pipeline/qc/custom_geometry_checks.py
"""
Responsibility: LPMO-specific active-site geometry validation.
Input:  Structure (gemmi or MDAnalysis Universe)
Output: GeometryResult with Cu–His distances, Cu–C1/C4 distances, angles

GATES:
  - Cu-His distance must be 1.5-3.0 A (hard fail)
  - Cu-His distance outside 1.8-2.6 A is retained as a soft warning
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

from lpmo_pipeline.config import load_defaults_config

logger = logging.getLogger(__name__)

_DEFAULTS_CONFIG = load_defaults_config()
_CHAIN_SCHEMA = _DEFAULTS_CONFIG.get("chain_schema") or {}
DEFAULT_PROTEIN_CHAIN = str(_CHAIN_SCHEMA.get("protein") or "A")
DEFAULT_CU_CHAIN = str(_CHAIN_SCHEMA.get("metal") or "E")
DEFAULT_GLYCAN_CHAINS = tuple(str(chain) for chain in (_CHAIN_SCHEMA.get("glycans") or ["B", "C", "D"]))


# ---------------------------------------------------------------------------
# Thresholds (from configs/thresholds.yaml)
# ---------------------------------------------------------------------------
CU_HIS_MIN: float = 1.5   # Å
CU_HIS_MAX: float = 3.0   # Å
CU_HIS_SOFT_MIN: float = 1.8   # Å
CU_HIS_SOFT_MAX: float = 2.6   # Å
CU_SUBSTRATE_FLAG: float = 7.0  # Å — flag if Cu–C1/C4 > this
HIS_BRACE_MAX_SEARCH_A: float = 3.0  # Å


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
    his_atom_position: tuple[float, float, float] = (0.0, 0.0, 0.0)


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
    cu_his_all_in_soft_range: bool = False

    cu_c1_measurements: list[CuSubstrateMeasurement] = field(default_factory=list)
    cu_c4_measurements: list[CuSubstrateMeasurement] = field(default_factory=list)
    min_cu_c1: float = float("inf")
    min_cu_c4: float = float("inf")

    his_brace_angle: float | None = None  # Angle between His-N–Cu–His-N plane
    active_site_orientation: dict[str, Any] = field(default_factory=dict)

    passed: bool = False
    failure_reasons: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------
def check_geometry(
    structure: Any,  # gemmi.Structure or MDAnalysis.Universe
    pose_id: str = "",
    cu_chain: str = DEFAULT_CU_CHAIN,
    glycan_chains: list[str] | None = None,
    hard_cu_his_min_a: float = CU_HIS_MIN,
    hard_cu_his_max_a: float = CU_HIS_MAX,
    soft_cu_his_min_a: float = CU_HIS_SOFT_MIN,
    soft_cu_his_max_a: float = CU_HIS_SOFT_MAX,
    cu_c_soft_flag_a: float = CU_SUBSTRATE_FLAG,
    his_brace_max_search_a: float = HIS_BRACE_MAX_SEARCH_A,
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
        glycan_chains = list(DEFAULT_GLYCAN_CHAINS)

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

    # --- Step 2: Cu–His distances (locked brace rules) ---
    his_residues = _find_his_brace_residues(
        structure,
        protein_chain=str(_CHAIN_SCHEMA.get("protein") or "A"),
    )
    candidate_n_atoms: list[dict[str, Any]] = []
    for his_res in his_residues:
        resnum = his_res.seqid.num
        if resnum == 1:
            atom_names = ["N", "ND1", "NE2"]
        else:
            atom_names = ["ND1", "NE2"]

        for atom_name in atom_names:
            atom = _get_atom_by_exact_name(his_res, atom_name)
            if atom is None:
                continue
            pos = np.array([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float)
            candidate_n_atoms.append(
                {
                    "his_chain": "A",
                    "his_resnum": resnum,
                    "atom_name": atom_name,
                    "position": pos,
                }
            )

    selected_n_atoms, selection_errors = _select_cu_his_gate_nitrogens(
        candidate_n_atoms=candidate_n_atoms,
        cu_pos=cu_pos,
        max_search_a=his_brace_max_search_a,
    )

    for atom_rec in selected_n_atoms:
        dist = float(atom_rec["distance_angstrom"])
        in_range = hard_cu_his_min_a <= dist <= hard_cu_his_max_a
        his_pos = np.asarray(atom_rec["position"], dtype=float)

        result.cu_his_measurements.append(
            CuHisMeasurement(
                his_chain=str(atom_rec["his_chain"]),
                his_resnum=int(atom_rec["his_resnum"]),
                his_atom=str(atom_rec["atom_name"]),
                cu_chain=cu_chain,
                cu_resnum=cu_atom.residue_seqid if hasattr(cu_atom, "residue_seqid") else 1,
                distance_angstrom=dist,
                in_range=in_range,
                his_atom_position=tuple(his_pos.tolist()),
            )
        )

    result.cu_his_all_in_range = (
        len(result.cu_his_measurements) == 3
        and all(m.in_range for m in result.cu_his_measurements)
    )
    result.cu_his_all_in_soft_range = (
        len(result.cu_his_measurements) == 3
        and all(
            soft_cu_his_min_a <= m.distance_angstrom <= soft_cu_his_max_a
            for m in result.cu_his_measurements
        )
    )
    if not result.cu_his_all_in_range:
        if selection_errors:
            result.failure_reasons.extend(selection_errors)
        bad = [
            f"His{m.his_resnum}:{m.his_atom}={m.distance_angstrom:.2f}Å"
            for m in result.cu_his_measurements
            if not m.in_range
        ]
        if bad:
            result.failure_reasons.append(
                "Cu-His distance out of range: "
                f"{bad} not in [{hard_cu_his_min_a:.2f}, {hard_cu_his_max_a:.2f}] A"
            )
            logger.error(
                "GATE FAIL: Cu-His distance outside hard range %.2f-%.2f A: %s (pose %s)",
                hard_cu_his_min_a,
                hard_cu_his_max_a,
                bad,
                pose_id,
            )
    elif not result.cu_his_all_in_soft_range:
        soft_bad = [
            f"His{m.his_resnum}:{m.his_atom}={m.distance_angstrom:.2f}Å"
            for m in result.cu_his_measurements
            if not (soft_cu_his_min_a <= m.distance_angstrom <= soft_cu_his_max_a)
        ]
        if soft_bad:
            warning = (
                "Cu-His distance outside preferred QC range: "
                f"{soft_bad} not in [{soft_cu_his_min_a:.2f}, {soft_cu_his_max_a:.2f}] A"
            )
            result.warnings.append(warning)
            logger.warning("SOFT QC: %s (pose %s)", warning, pose_id)

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
                    flagged=d_c1 > cu_c_soft_flag_a,
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
                    flagged=d_c4 > cu_c_soft_flag_a,
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
        "Geometry check pose=%s: passed=%s, Cu-His hard=%s, Cu-His soft=%s, "
        "min Cu-C1=%.2f, min Cu-C4=%.2f",
        pose_id,
        result.passed,
        result.cu_his_all_in_range,
        result.cu_his_all_in_soft_range,
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


def _find_his_brace_residues(structure: Any, protein_chain: str = DEFAULT_PROTEIN_CHAIN) -> list[Any]:
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


def _get_atom_by_exact_name(residue: Any, atom_name: str) -> Any | None:
    """Get atom by exact atom name from residue."""
    for atom in residue:
        if atom.name.strip() == atom_name:
            return atom
    return None


def _select_cu_his_gate_nitrogens(
    candidate_n_atoms: list[dict[str, Any]],
    cu_pos: np.ndarray,
    max_search_a: float,
) -> tuple[list[dict[str, Any]], list[str]]:
    """Select the 3 Cu-His gate nitrogens using locked Stage 4 rules.

    Rules:
    - His1:N and His1:ND1 are mandatory
    - His1:NE2 is excluded from selected gate atoms
    - Third atom must be histidine N from residue != 1
    """
    errors: list[str] = []

    def _find_one(resnum: int, atom_name: str) -> dict[str, Any] | None:
        for rec in candidate_n_atoms:
            if int(rec.get("his_resnum", -1)) == resnum and rec.get("atom_name") == atom_name:
                return rec
        return None

    his1_n = _find_one(1, "N")
    his1_nd1 = _find_one(1, "ND1")

    if his1_n is None:
        errors.append("Cu-His selection failed: mandatory His1:N not found")
    if his1_nd1 is None:
        errors.append("Cu-His selection failed: mandatory His1:ND1 not found")

    selected: list[dict[str, Any]] = []
    for rec in [his1_n, his1_nd1]:
        if rec is None:
            continue
        pos = np.asarray(rec["position"], dtype=float)
        dist = float(np.linalg.norm(cu_pos - pos))
        selected.append({**rec, "distance_angstrom": dist})

    third_candidates: list[dict[str, Any]] = []
    for rec in candidate_n_atoms:
        resnum = int(rec.get("his_resnum", -1))
        atom_name = str(rec.get("atom_name", ""))
        if resnum == 1:
            continue
        if atom_name not in {"ND1", "NE2"}:
            continue
        pos = np.asarray(rec["position"], dtype=float)
        dist = float(np.linalg.norm(cu_pos - pos))
        if dist <= max_search_a:
            third_candidates.append({**rec, "distance_angstrom": dist})

    if third_candidates:
        third = min(third_candidates, key=lambda r: float(r["distance_angstrom"]))
        selected.append(third)
    else:
        errors.append(
            "Cu-His selection failed: no non-His1 histidine ND1/NE2 within "
            f"{max_search_a:.2f} Å"
        )

    return selected, errors


def _compute_his_brace_angle(
    measurements: list[CuHisMeasurement],
    cu_pos: np.ndarray,
) -> float | None:
    """Compute angle between the two His-N–Cu vectors (His-brace angle).

    Groups measurements by His residue, picks the closest coordinating N
    atom from each, and returns the N1–Cu–N2 angle in degrees.

    Returns None if fewer than 2 distinct His residues are measured.
    """
    best_per_his: dict[int, CuHisMeasurement] = {}
    for m in measurements:
        prev = best_per_his.get(m.his_resnum)
        if prev is None or m.distance_angstrom < prev.distance_angstrom:
            best_per_his[m.his_resnum] = m

    if len(best_per_his) < 2:
        return None

    ordered = sorted(best_per_his.values(), key=lambda m: m.distance_angstrom)
    m1, m2 = ordered[0], ordered[1]

    v1 = np.asarray(m1.his_atom_position) - cu_pos
    v2 = np.asarray(m2.his_atom_position) - cu_pos

    norm1 = float(np.linalg.norm(v1))
    norm2 = float(np.linalg.norm(v2))
    if norm1 < 1e-9 or norm2 < 1e-9:
        return None

    cos_angle = float(np.dot(v1, v2) / (norm1 * norm2))
    cos_angle = max(-1.0, min(1.0, cos_angle))
    return float(np.degrees(np.arccos(cos_angle)))
