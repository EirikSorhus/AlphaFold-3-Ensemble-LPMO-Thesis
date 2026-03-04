# src/lpmo_pipeline/analysis/crystal_anchoring.py
"""
Responsibility: Compare predicted clusters against crystal structures (if available).
Input:  Predicted cluster IFPs + geometry, crystal complex PDB
Output: crystal_metrics.json with per-cluster similarity scores

STEP 10 in masterplan.
RQ4: Model credibility via crystallographic comparison.

Metrics:
  - Tanimoto(cluster_medoid_IFP, crystal_IFP) — interaction fingerprint similarity
  - Pocket RMSD (local, pocket residues only) vs crystal
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from lpmo_pipeline.analysis.prolif_ifp import compute_ifp_single, compute_tanimoto_similarity

logger = logging.getLogger(__name__)


CRYSTAL_SIM_MODERATE_THRESHOLD = 0.5  # Empirical; flag below this
POCKET_RMSD_THRESHOLD = 2.5  # Å


@dataclass
class CrystalComparisonResult:
    """Comparison of a single cluster against crystal."""

    cluster_id: int
    medoid_pose_id: str = ""
    ifp_tanimoto: float = 0.0
    pocket_rmsd: float | None = None
    ifp_above_threshold: bool = False
    pocket_rmsd_below_threshold: bool = False


@dataclass
class CrystalAnchoringReport:
    """Full crystal anchoring report for a protein×ligand pair."""

    protein_id: str = ""
    ligand_id: str = ""
    crystal_pdb: str = ""
    has_crystal_ligand: bool = False
    comparisons: list[CrystalComparisonResult] = field(default_factory=list)
    best_cluster_id: int = -1
    best_tanimoto: float = 0.0


def run_crystal_anchoring(
    crystal_pdb: Path,
    crystal_ligand_mol2: Path | None,
    cluster_medoid_data: list[dict[str, Any]],
    protein_id: str = "",
    ligand_id: str = "",
    pocket_residues: list[int] | None = None,
    protein_chain: str = "A",
) -> CrystalAnchoringReport:
    """Compare predicted clusters to crystal structure.

    Args:
        crystal_pdb: Path to crystal complex PDB.
        crystal_ligand_mol2: Path to crystal ligand MOL2 (None if apo).
        cluster_medoid_data: List of {
            "cluster_id": int,
            "medoid_pose_id": str,
            "medoid_pdb": Path,
            "medoid_mol2": Path,
            "medoid_ifp": list[int],  # flat bitvector
        }
        protein_id: For metadata.
        ligand_id: For metadata.
        pocket_residues: Residue numbers for pocket RMSD.
        protein_chain: Protein chain.

    Returns:
        CrystalAnchoringReport.
    """
    report = CrystalAnchoringReport(
        protein_id=protein_id,
        ligand_id=ligand_id,
        crystal_pdb=str(crystal_pdb),
    )

    # --- Case 1: Crystal has ligand → compute IFP similarity ---
    if crystal_ligand_mol2 is not None and crystal_ligand_mol2.exists():
        report.has_crystal_ligand = True

        # Compute crystal IFP
        crystal_ifp_result = compute_ifp_single(
            complex_pdb=crystal_pdb,
            ligand_mol2=crystal_ligand_mol2,
            pose_id="crystal",
            protein_chain=protein_chain,
        )
        crystal_bitvec = crystal_ifp_result.flat_bitvector

        # Compare each cluster medoid to crystal
        for cdata in cluster_medoid_data:
            tanimoto = compute_tanimoto_similarity(
                cdata["medoid_ifp"], crystal_bitvec
            )
            comparison = CrystalComparisonResult(
                cluster_id=cdata["cluster_id"],
                medoid_pose_id=cdata["medoid_pose_id"],
                ifp_tanimoto=tanimoto,
                ifp_above_threshold=tanimoto >= CRYSTAL_SIM_MODERATE_THRESHOLD,
            )

            # Optional: pocket RMSD
            if pocket_residues:
                rmsd = _compute_pocket_rmsd(
                    pred_pdb=cdata["medoid_pdb"],
                    crystal_pdb=crystal_pdb,
                    pocket_residues=pocket_residues,
                    protein_chain=protein_chain,
                )
                comparison.pocket_rmsd = rmsd
                comparison.pocket_rmsd_below_threshold = (
                    rmsd is not None and rmsd < POCKET_RMSD_THRESHOLD
                )

            report.comparisons.append(comparison)
    else:
        # --- Case 2: Apo crystal → pocket RMSD only ---
        report.has_crystal_ligand = False
        if pocket_residues:
            for cdata in cluster_medoid_data:
                rmsd = _compute_pocket_rmsd(
                    pred_pdb=cdata["medoid_pdb"],
                    crystal_pdb=crystal_pdb,
                    pocket_residues=pocket_residues,
                    protein_chain=protein_chain,
                )
                comparison = CrystalComparisonResult(
                    cluster_id=cdata["cluster_id"],
                    medoid_pose_id=cdata["medoid_pose_id"],
                    pocket_rmsd=rmsd,
                    pocket_rmsd_below_threshold=(
                        rmsd is not None and rmsd < POCKET_RMSD_THRESHOLD
                    ),
                )
                report.comparisons.append(comparison)

    # Best cluster
    if report.comparisons:
        best = max(report.comparisons, key=lambda c: c.ifp_tanimoto)
        report.best_cluster_id = best.cluster_id
        report.best_tanimoto = best.ifp_tanimoto

    logger.info(
        "Crystal anchoring %s×%s: %d clusters compared, best Tanimoto=%.3f (cluster %d)",
        protein_id, ligand_id, len(report.comparisons),
        report.best_tanimoto, report.best_cluster_id,
    )
    return report


def _compute_pocket_rmsd(
    pred_pdb: Path,
    crystal_pdb: Path,
    pocket_residues: list[int],
    protein_chain: str = "A",
) -> float | None:
    """Compute pocket RMSD between predicted and crystal structures.

    PSEUDOCODE — uses MDAnalysis align.
    """
    # import MDAnalysis as mda
    # from MDAnalysis.analysis.rms import rmsd as compute_rmsd
    #
    # u_pred = mda.Universe(str(pred_pdb))
    # u_xtal = mda.Universe(str(crystal_pdb))
    #
    # pocket_sel = "resid " + " ".join(str(r) for r in pocket_residues)
    # sel = f"segid {protein_chain} and ({pocket_sel}) and name CA"
    #
    # pred_atoms = u_pred.select_atoms(sel)
    # xtal_atoms = u_xtal.select_atoms(sel)
    #
    # if len(pred_atoms) == 0 or len(xtal_atoms) == 0:
    #     return None
    # if len(pred_atoms) != len(xtal_atoms):
    #     logger.warning("Pocket atom count mismatch: pred=%d, xtal=%d",
    #                    len(pred_atoms), len(xtal_atoms))
    #     return None
    #
    # return float(compute_rmsd(pred_atoms.positions, xtal_atoms.positions))
    return None  # PSEUDOCODE placeholder


def write_crystal_anchoring_report(
    report: CrystalAnchoringReport,
    output_path: Path,
) -> None:
    """Write crystal anchoring report to JSON."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "protein_id": report.protein_id,
        "ligand_id": report.ligand_id,
        "crystal_pdb": report.crystal_pdb,
        "has_crystal_ligand": report.has_crystal_ligand,
        "best_cluster_id": report.best_cluster_id,
        "best_tanimoto": report.best_tanimoto,
        "comparisons": [
            {
                "cluster_id": c.cluster_id,
                "medoid_pose_id": c.medoid_pose_id,
                "ifp_tanimoto": c.ifp_tanimoto,
                "pocket_rmsd": c.pocket_rmsd,
                "ifp_above_threshold": c.ifp_above_threshold,
                "pocket_rmsd_below_threshold": c.pocket_rmsd_below_threshold,
            }
            for c in report.comparisons
        ],
    }
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    logger.info("Wrote crystal anchoring report to %s", output_path)
