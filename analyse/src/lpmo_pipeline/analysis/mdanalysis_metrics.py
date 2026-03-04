# src/lpmo_pipeline/analysis/mdanalysis_metrics.py
"""
Responsibility: Compute 3D structural metrics using MDAnalysis.
Input:  PDB file (complex_H.pdb from protonation step)
Output: geometry_metrics.json per pose

Metrics computed:
  - Cu–C1, Cu–C4 distances (per glycan residue)
  - Cu–His-brace distances
  - His-brace angle / tilt
  - Core RMSD (align on protein core, exclude loops)
  - H-bonds (optional, if trajectory)
  - SASA (optional)

RQ relevance:
  RQ1: Cu–C1 vs Cu–C4 geometry → regioselectivity
  RQ2: contact patterns → substrate specificity
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Result model
# ---------------------------------------------------------------------------
@dataclass
class PoseGeometryMetrics:
    """All geometry metrics for a single pose."""

    pose_id: str = ""
    model: str = ""
    protein_id: str = ""
    ligand_id: str = ""

    # Cu–substrate distances
    cu_c1_distances: list[dict[str, Any]] = field(default_factory=list)
    cu_c4_distances: list[dict[str, Any]] = field(default_factory=list)
    min_cu_c1: float = float("inf")
    min_cu_c4: float = float("inf")

    # Cu–His-brace
    cu_his_distances: list[dict[str, Any]] = field(default_factory=list)
    his_brace_angle_deg: float | None = None

    # Alignment
    core_rmsd_vs_reference: float | None = None  # If reference provided
    pocket_rmsd_vs_crystal: float | None = None   # If crystal provided

    # Optional
    hbonds: list[dict[str, Any]] = field(default_factory=list)
    sasa_protein: float | None = None
    sasa_ligand: float | None = None


# ---------------------------------------------------------------------------
# Core metric computation
# ---------------------------------------------------------------------------
def compute_pose_metrics(
    pdb_path: Path,
    pose_id: str = "",
    protein_chain: str = "A",
    glycan_chains: list[str] | None = None,
    cu_chain: str = "E",
    reference_pdb: Path | None = None,
    crystal_pdb: Path | None = None,
    core_selection: str = "protein and not (resid 1:5 or name H*)",
    pocket_residues: list[int] | None = None,
) -> PoseGeometryMetrics:
    """Compute all geometry metrics for a single pose.

    PSEUDOCODE — relies on MDAnalysis Universe.

    Args:
        pdb_path: Path to protonated complex PDB.
        pose_id: Identifier for this pose.
        protein_chain: Chain ID for protein (default "A").
        glycan_chains: Chain IDs for glycans (default ["B","C","D"]).
        cu_chain: Chain ID for Cu (default "E").
        reference_pdb: Optional reference for core RMSD alignment.
        crystal_pdb: Optional crystal structure for pocket RMSD.
        core_selection: MDAnalysis selection string for core atoms.
        pocket_residues: Residue numbers defining the binding pocket.

    Returns:
        PoseGeometryMetrics.
    """
    if glycan_chains is None:
        glycan_chains = ["B", "C", "D"]

    metrics = PoseGeometryMetrics(pose_id=pose_id)
    logger.info("Computing geometry metrics for pose %s: %s", pose_id, pdb_path)

    # --- Step 1: Load structure ---
    # import MDAnalysis as mda
    # u = mda.Universe(str(pdb_path))

    # --- Step 2: Find Cu ---
    # cu_sel = u.select_atoms(f"segid {cu_chain} and element Cu")
    # if len(cu_sel) == 0:
    #     logger.error("Cu not found in chain %s", cu_chain)
    #     return metrics
    # cu_pos = cu_sel.positions[0]

    # --- Step 3: Cu–C1 and Cu–C4 distances ---
    # for gly_chain in glycan_chains:
    #     gly_residues = u.select_atoms(f"segid {gly_chain}").residues
    #     for res in gly_residues:
    #         c1_atoms = res.atoms.select_atoms("name C1")
    #         c4_atoms = res.atoms.select_atoms("name C4")
    #         if len(c1_atoms) > 0:
    #             d_c1 = np.linalg.norm(c1_atoms.positions[0] - cu_pos)
    #             metrics.cu_c1_distances.append({
    #                 "chain": gly_chain, "resnum": res.resid,
    #                 "resname": res.resname, "distance": float(d_c1),
    #             })
    #             metrics.min_cu_c1 = min(metrics.min_cu_c1, d_c1)
    #         if len(c4_atoms) > 0:
    #             d_c4 = np.linalg.norm(c4_atoms.positions[0] - cu_pos)
    #             metrics.cu_c4_distances.append({
    #                 "chain": gly_chain, "resnum": res.resid,
    #                 "resname": res.resname, "distance": float(d_c4),
    #             })
    #             metrics.min_cu_c4 = min(metrics.min_cu_c4, d_c4)

    # --- Step 4: Cu–His distances ---
    # his_residues = u.select_atoms(f"segid {protein_chain} and resname HIS").residues
    # for his in his_residues:
    #     for atom_name in ["NE2", "ND1"]:
    #         atoms = his.atoms.select_atoms(f"name {atom_name}")
    #         if len(atoms) > 0:
    #             d = np.linalg.norm(atoms.positions[0] - cu_pos)
    #             metrics.cu_his_distances.append({
    #                 "resnum": his.resid, "atom": atom_name, "distance": float(d),
    #             })

    # --- Step 5: His-brace angle ---
    # PSEUDOCODE: compute angle N-Cu-N between two Cu-coordinating His
    # metrics.his_brace_angle_deg = _compute_ncu_n_angle(cu_pos, his_n_positions)

    # --- Step 6: Core RMSD vs reference ---
    # if reference_pdb:
    #     ref = mda.Universe(str(reference_pdb))
    #     from MDAnalysis.analysis import align
    #     align.alignto(u, ref, select=core_selection)
    #     core_u = u.select_atoms(core_selection)
    #     core_ref = ref.select_atoms(core_selection)
    #     metrics.core_rmsd_vs_reference = float(
    #         np.sqrt(np.mean(np.sum((core_u.positions - core_ref.positions)**2, axis=1)))
    #     )

    # --- Step 7: Pocket RMSD vs crystal ---
    # if crystal_pdb and pocket_residues:
    #     xtal = mda.Universe(str(crystal_pdb))
    #     pocket_sel = "resid " + " ".join(str(r) for r in pocket_residues)
    #     pocket_u = u.select_atoms(f"segid {protein_chain} and ({pocket_sel}) and name CA")
    #     pocket_xtal = xtal.select_atoms(f"segid A and ({pocket_sel}) and name CA")
    #     if len(pocket_u) > 0 and len(pocket_xtal) > 0:
    #         metrics.pocket_rmsd_vs_crystal = float(
    #             np.sqrt(np.mean(np.sum(
    #                 (pocket_u.positions - pocket_xtal.positions)**2, axis=1
    #             )))
    #         )

    return metrics


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def write_geometry_metrics(metrics: PoseGeometryMetrics, output_path: Path) -> None:
    """Write geometry metrics to JSON."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "pose_id": metrics.pose_id,
        "model": metrics.model,
        "protein_id": metrics.protein_id,
        "ligand_id": metrics.ligand_id,
        "min_cu_c1": metrics.min_cu_c1,
        "min_cu_c4": metrics.min_cu_c4,
        "cu_c1_distances": metrics.cu_c1_distances,
        "cu_c4_distances": metrics.cu_c4_distances,
        "cu_his_distances": metrics.cu_his_distances,
        "his_brace_angle_deg": metrics.his_brace_angle_deg,
        "core_rmsd_vs_reference": metrics.core_rmsd_vs_reference,
        "pocket_rmsd_vs_crystal": metrics.pocket_rmsd_vs_crystal,
        "hbonds": metrics.hbonds,
        "sasa_protein": metrics.sasa_protein,
        "sasa_ligand": metrics.sasa_ligand,
    }
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    logger.info("Wrote geometry metrics to %s", output_path)
