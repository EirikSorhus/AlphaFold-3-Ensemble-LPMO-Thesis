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

import numpy as np

logger = logging.getLogger(__name__)


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
    ligand_mol2: Path,
    pose_id: str = "",
    lpmo_residues: list[int] | None = None,
    cbm_residues: list[int] | None = None,
    protein_chain: str = "A",
    glycan_chains: list[str] | None = None,
    cu_chain: str = "E",
) -> CBMAnalysisResult:
    """Compute CBM-specific metrics for a full-length prediction.

    Key difference from standard analysis: separate LPMO and CBM contacts.

    Args:
        complex_pdb: Path to protonated full-length complex.
        ligand_mol2: Path to ligand MOL2.
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
        glycan_chains = ["B", "C", "D"]

    result = CBMAnalysisResult(pose_id=pose_id)
    logger.info("CBM analysis for pose %s", pose_id)

    # --- Step 1: Compute IFP for LPMO domain residues only ---
    # PSEUDOCODE: Use ProLIF with selection restricted to LPMO residues
    # import prolif as plf
    # u = mda.Universe(str(complex_pdb))
    # if lpmo_residues:
    #     lpmo_sel = "resid " + " ".join(str(r) for r in lpmo_residues)
    #     prot_lpmo = plf.Molecule.from_mda(u, selection=f"segid {protein_chain} and ({lpmo_sel})")
    # lig = plf.Molecule.from_file(str(ligand_mol2))
    # fp_lpmo = plf.Fingerprint()
    # fp_lpmo.run(lig, prot_lpmo)
    # result.ifp_lpmo_bitvector = fp_lpmo.to_bitvectors()[0].tolist()

    # --- Step 2: Compute IFP for CBM domain residues only ---
    # if cbm_residues:
    #     cbm_sel = "resid " + " ".join(str(r) for r in cbm_residues)
    #     prot_cbm = plf.Molecule.from_mda(u, selection=f"segid {protein_chain} and ({cbm_sel})")
    #     fp_cbm = plf.Fingerprint()
    #     fp_cbm.run(lig, prot_cbm)
    #     result.ifp_cbm_bitvector = fp_cbm.to_bitvectors()[0].tolist()

    # --- Step 3: CBM proximity metrics ---
    # PSEUDOCODE: MDAnalysis distance calculations
    # cbm_atoms = u.select_atoms(f"segid {protein_chain} and ({cbm_sel}) and not name H*")
    # lig_atoms = u.select_atoms(f"segid {' '.join(glycan_chains)} and not name H*")
    # cu_atoms = u.select_atoms(f"segid {cu_chain} and element Cu")
    #
    # dists_cbm_lig = mda.lib.distances.distance_array(cbm_atoms.positions, lig_atoms.positions)
    # result.proximity.cbm_ligand_min_dist = float(dists_cbm_lig.min())
    # result.proximity.cbm_contacts_ligand = result.proximity.cbm_ligand_min_dist < CBM_LIGAND_CONTACT_THRESHOLD
    #
    # if len(cu_atoms) > 0:
    #     dists_cbm_cu = mda.lib.distances.distance_array(cbm_atoms.positions, cu_atoms.positions)
    #     result.proximity.cbm_cu_min_dist = float(dists_cbm_cu.min())

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
