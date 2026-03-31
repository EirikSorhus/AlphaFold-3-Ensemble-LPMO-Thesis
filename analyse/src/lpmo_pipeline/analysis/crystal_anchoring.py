# src/lpmo_pipeline/analysis/crystal_anchoring.py
"""
Responsibility: Compare predicted clusters against crystal structures (if available).
Input:  Predicted cluster IFPs + geometry, crystal complex PDB
Output: crystal_metrics.json with per-cluster similarity scores

STEP 12 in masterplan.
RQ4: Model credibility via crystallographic comparison.

Alignment method: PyMOL `pair_fit` for optimal local superposition.
  - Align on: histidine-brace (Cα/Nε2), Cu-coordinating residues,
    and substrate-recognition surface residues.
  - Substrate-recognition residues: literature-based (preferred) or
    proximity-based fallback (all protein residues within cutoff of ligand).
  - Always document that RMSD is measured after optimized local alignment.

Metrics:
  - Tanimoto(cluster_medoid_IFP, crystal_IFP) — interaction fingerprint similarity
  - Pocket RMSD after optimal local alignment via PyMOL pair_fit
  - Per-residue deviations for alignment atoms

Adapted patterns from:
  - PoseBench (MIT): general structure comparison workflow
  - benchmarking-af3 (MIT): pocket residue identification via proximity cutoff
  See ATTRIBUTION.md for full credits.
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

    Uses PyMOL pair_fit for optimal local superposition on:
      1. Histidine-brace atoms (Cα, Nε2)
      2. Cu-coordinating residues
      3. Substrate-recognition surface residues (from pocket_residues)

    The RMSD is measured AFTER optimized local alignment (not global).
    This gives better results for local pocket comparison but must be
    clearly documented as an optimized alignment metric.

    Inspired by benchmarking-af3 pocket residue approach (MIT license),
    but uses PyMOL pair_fit instead of APoc/Biopython NeighborSearch.
    """
    # TODO: implement with PyMOL pair_fit
    # from pymol import cmd
    #
    # cmd.load(str(crystal_pdb), "crystal")
    # cmd.load(str(pred_pdb), "predicted")
    #
    # # Build selection string for pair_fit atoms (CA of pocket residues)
    # pocket_sel = " or ".join(f"resi {r}" for r in pocket_residues)
    # crystal_sel = f"crystal and chain {protein_chain} and ({pocket_sel}) and name CA"
    # pred_sel = f"predicted and chain {protein_chain} and ({pocket_sel}) and name CA"
    #
    # # pair_fit returns RMSD after optimal superposition
    # rmsd = cmd.pair_fit(pred_sel, crystal_sel)
    #
    # cmd.delete("all")
    # return float(rmsd) if rmsd is not None else None
    return None  # placeholder — implementation pending


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


def identify_pocket_residues_by_proximity(
    complex_pdb: Path,
    ligand_chain: str = "B",
    protein_chain: str = "A",
    cutoff_a: float = 5.0,
) -> list[int]:
    """Identify protein residues near a ligand using spatial proximity.

    Fallback method when literature-based substrate-recognition residues
    are not available. Finds all protein residues with any heavy atom
    within `cutoff_a` Angstrom of any ligand heavy atom.

    Inspired by benchmarking-af3 `3_find_pocket_residues.py` pocket
    identification via NeighborSearch (MIT license). Our implementation
    uses gemmi instead of Biopython.

    Note: proximity-based selection may give different residue sets for
    different prediction models (AF3/RF3/Boltz-2) depending on how each
    model places the ligand. This must be documented in results.

    Args:
        complex_pdb: Path to PDB file with protein + ligand.
        ligand_chain: Chain ID for ligand.
        protein_chain: Chain ID for protein.
        cutoff_a: Distance cutoff in Angstrom.

    Returns:
        Sorted list of residue sequence numbers within cutoff.
    """
    # TODO: implement with gemmi
    # import gemmi
    #
    # st = gemmi.read_structure(str(complex_pdb))
    # model = st[0]
    # ns = gemmi.NeighborSearch(model, st.cell, cutoff_a).populate()
    #
    # ligand_atoms = []
    # protein_residues = set()
    #
    # for chain in model:
    #     if chain.name == ligand_chain:
    #         for res in chain:
    #             for atom in res:
    #                 ligand_atoms.append(atom)
    #
    # for latom in ligand_atoms:
    #     marks = ns.find_atoms(latom.pos, '\0', cutoff_a)
    #     for mark in marks:
    #         cra = mark.to_cra(model)
    #         if cra.chain.name == protein_chain:
    #             protein_residues.add(cra.residue.seqid.num)
    #
    # return sorted(protein_residues)
    return []  # placeholder — implementation pending
