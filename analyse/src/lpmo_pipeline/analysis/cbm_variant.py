# src/lpmo_pipeline/analysis/cbm_variant.py
"""
Responsibility: CBM-specific analysis for DEL B (full-length LPMO + CBM).
Input:  Full-length predictions (LPMO + CBM + ligand)
Output: Dual IFP (IFP_LPMO + IFP_CBM), CBM proximity metrics

STEP 12 in masterplan.
RQ3: Does CBM modeling shift ligand poses / alter binding contacts?
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from lpmo_pipeline.config import load_defaults_config

logger = logging.getLogger(__name__)

_DEFAULTS_CONFIG = load_defaults_config()
_CHAIN_SCHEMA = _DEFAULTS_CONFIG.get("chain_schema") or {}
DEFAULT_PROTEIN_CHAIN = str(_CHAIN_SCHEMA.get("protein") or "A")
DEFAULT_GLYCAN_CHAINS = tuple(str(chain) for chain in (_CHAIN_SCHEMA.get("glycans") or ["B", "C", "D"]))
DEFAULT_CU_CHAIN = str(_CHAIN_SCHEMA.get("metal") or "E")


@dataclass
class CBMProximityMetrics:
    """CBM-specific proximity measurements for a single pose."""

    pose_id: str = ""
    cbm_ligand_min_dist: float = float("inf")  # Å (heavy atoms)
    cbm_cu_min_dist: float = float("inf")       # Å
    cbm_contacts_ligand: bool = False            # True if cbm_ligand_min_dist < threshold
    cbm_contacts_near_cu: bool = False           # True if close to Cu-bound ligand


@dataclass
class CBMAnalysisResult:
    """Full CBM analysis result for a pose."""

    pose_id: str = ""
    protein_id: str = ""
    ligand_id: str = ""

    # Dual IFP
    ifp_lpmo_bitvector: list[int] = field(default_factory=list)
    ifp_cbm_bitvector: list[int] = field(default_factory=list)

    # Proximity
    proximity: CBMProximityMetrics = field(default_factory=CBMProximityMetrics)

    # Residue-level: which LPMO residues contact ligand, which CBM residues
    lpmo_contact_residues: list[str] = field(default_factory=list)
    cbm_contact_residues: list[str] = field(default_factory=list)


CBM_LIGAND_CONTACT_THRESHOLD = 4.5  # Å for "CBM contacts ligand"


def compute_cbm_analysis(
    complex_pdb: Path,
    ligand_pdb: Path,
    pose_id: str = "",
    lpmo_residues: list[int] | None = None,
    cbm_residues: list[int] | None = None,
    protein_chain: str = DEFAULT_PROTEIN_CHAIN,
    glycan_chains: list[str] | None = None,
    cu_chain: str = DEFAULT_CU_CHAIN,
) -> CBMAnalysisResult:
    """Compute CBM-specific metrics for a full-length prediction.

    Key difference from standard analysis: separate LPMO and CBM contacts.

    Args:
        complex_pdb: Path to non-protonated full-length complex PDB.
        ligand_pdb: Path to ligand-only PDB.
        pose_id: Identifier.
        lpmo_residues: Residue numbers belonging to LPMO domain.
        cbm_residues: Residue numbers belonging to CBM domain.
        protein_chain: Chain for protein.
        glycan_chains: Chains for glycan.
        cu_chain: Chain for Cu.

    Returns:
        CBMAnalysisResult with dual IFP and proximity metrics.
    """
    if glycan_chains is None:
        glycan_chains = list(DEFAULT_GLYCAN_CHAINS)

    result = CBMAnalysisResult(pose_id=pose_id)
    logger.info("CBM analysis for pose %s", pose_id)

    logger.info(
        "Pose-level CBM dual-IFP/proximity backend is not active for %s; "
        "condition-level CBM paired analysis is implemented in cbm_comparison.py",
        pose_id,
    )

    return result


def write_cbm_analysis(result: CBMAnalysisResult, output_path: Path) -> None:
    """Write CBM analysis result to JSON."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "pose_id": result.pose_id,
        "protein_id": result.protein_id,
        "ligand_id": result.ligand_id,
        "ifp_lpmo_n_active": sum(result.ifp_lpmo_bitvector),
        "ifp_cbm_n_active": sum(result.ifp_cbm_bitvector),
        "cbm_ligand_min_dist": result.proximity.cbm_ligand_min_dist,
        "cbm_cu_min_dist": result.proximity.cbm_cu_min_dist,
        "cbm_contacts_ligand": result.proximity.cbm_contacts_ligand,
        "lpmo_contact_residues": result.lpmo_contact_residues,
        "cbm_contact_residues": result.cbm_contact_residues,
    }
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
