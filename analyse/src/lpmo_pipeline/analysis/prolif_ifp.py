# src/lpmo_pipeline/analysis/prolif_ifp.py
"""
Responsibility: Generate interaction fingerprints (IFP) using ProLIF.
Input:  complex PDB + ligand MOL2 (with bond orders and charges)
Output: IFP matrix (binary: residues × interaction types)

RQ relevance:
  RQ1: IFP signatures distinguish C1 vs C4 binding modes
  RQ2: IFP per substrate type → specificity patterns
  RQ3: IFP_LPMO vs IFP_CBM in DEL B
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
class IFPResult:
    """Interaction fingerprint result for a single pose."""

    pose_id: str = ""
    n_residues: int = 0
    n_interaction_types: int = 0
    residue_names: list[str] = field(default_factory=list)
    interaction_types: list[str] = field(default_factory=list)
    fingerprint: list[list[int]] = field(default_factory=list)  # n_residues × n_types
    flat_bitvector: list[int] = field(default_factory=list)     # Flattened for clustering


@dataclass
class IFPBatch:
    """IFP results for multiple poses (same protein×ligand)."""

    protein_id: str = ""
    ligand_id: str = ""
    model: str = ""
    results: list[IFPResult] = field(default_factory=list)
    # Matrix: n_poses × n_features (for clustering input)
    matrix: list[list[int]] = field(default_factory=list)
    feature_names: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# IFP generation
# ---------------------------------------------------------------------------
def compute_ifp_single(
    complex_pdb: Path,
    ligand_mol2: Path,
    pose_id: str = "",
    protein_chain: str = "A",
    interaction_types: list[str] | None = None,
) -> IFPResult:
    """Compute ProLIF IFP for a single pose.

    PSEUDOCODE — wraps ProLIF Python API.

    Args:
        complex_pdb: Path to protonated complex PDB.
        ligand_mol2: Path to ligand MOL2 with bond orders & charges.
        pose_id: Identifier for this pose.
        protein_chain: Protein chain for IFP computation.
        interaction_types: Which interactions to compute.
            Default: HBDonor, HBAcceptor, Hydrophobic, PiStacking, VdWContact, etc.

    Returns:
        IFPResult with binary fingerprint matrix.
    """
    if interaction_types is None:
        interaction_types = [
            "HBDonor", "HBAcceptor", "Hydrophobic", "PiStacking",
            "PiCation", "Anionic", "Cationic", "VdWContact",
        ]

    logger.info("Computing IFP for pose %s: %s + %s", pose_id, complex_pdb, ligand_mol2)

    # --- Step 1: Load molecules ---
    # import prolif as plf
    # import MDAnalysis as mda
    #
    # u = mda.Universe(str(complex_pdb))
    # lig = plf.Molecule.from_file(str(ligand_mol2))
    # prot = plf.Molecule.from_mda(u, selection=f"segid {protein_chain}")

    # --- Step 2: Compute fingerprint ---
    # fp = plf.Fingerprint(interactions=interaction_types)
    # fp.run(lig, prot)
    # df = fp.to_dataframe()

    # --- Step 3: Extract binary matrix ---
    # PSEUDOCODE: df has columns like (residue_name, interaction_type) → 0/1
    # residue_names = [col[0] for col in df.columns]  # unique residues
    # interaction_cols = [col[1] for col in df.columns]
    # values = df.values[0]  # single frame
    residue_names: list[str] = []  # PSEUDOCODE: filled from df
    fingerprint: list[list[int]] = []  # PSEUDOCODE: filled from df
    flat_bitvector: list[int] = []  # PSEUDOCODE: flattened

    result = IFPResult(
        pose_id=pose_id,
        n_residues=len(residue_names),
        n_interaction_types=len(interaction_types),
        residue_names=residue_names,
        interaction_types=interaction_types,
        fingerprint=fingerprint,
        flat_bitvector=flat_bitvector,
    )

    logger.info(
        "IFP for %s: %d residues × %d types = %d features, %d active",
        pose_id, result.n_residues, result.n_interaction_types,
        len(flat_bitvector), sum(flat_bitvector),
    )
    return result


def compute_ifp_batch(
    pose_data: list[dict[str, Path]],
    protein_id: str = "",
    ligand_id: str = "",
    model: str = "",
    protein_chain: str = "A",
) -> IFPBatch:
    """Compute IFP for all poses in a protein×ligand×model combination.

    Args:
        pose_data: List of {"pose_id": str, "complex_pdb": Path, "ligand_mol2": Path}
        protein_id: For metadata.
        ligand_id: For metadata.
        model: For metadata.
        protein_chain: Protein chain for IFP.

    Returns:
        IFPBatch with aligned matrix ready for clustering.
    """
    results: list[IFPResult] = []
    for pd in pose_data:
        r = compute_ifp_single(
            complex_pdb=pd["complex_pdb"],
            ligand_mol2=pd["ligand_mol2"],
            pose_id=pd["pose_id"],
            protein_chain=protein_chain,
        )
        results.append(r)

    # Align all fingerprints to the same feature set
    all_feature_names: list[str] = []
    if results:
        # Use union of all residue×interaction features
        feature_set: set[str] = set()
        for r in results:
            for res in r.residue_names:
                for itype in r.interaction_types:
                    feature_set.add(f"{res}_{itype}")
        all_feature_names = sorted(feature_set)

    # Build aligned matrix
    matrix: list[list[int]] = []
    for r in results:
        row = [0] * len(all_feature_names)
        # PSEUDOCODE: map each r.flat_bitvector entry to the aligned feature index
        matrix.append(row)

    batch = IFPBatch(
        protein_id=protein_id,
        ligand_id=ligand_id,
        model=model,
        results=results,
        matrix=matrix,
        feature_names=all_feature_names,
    )

    logger.info(
        "IFP batch: %s × %s × %s → %d poses × %d features",
        protein_id, ligand_id, model, len(results), len(all_feature_names),
    )
    return batch


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def write_ifp_matrix(batch: IFPBatch, output_path: Path) -> None:
    """Write IFP matrix to CSV (columns = features, rows = poses)."""
    import csv

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pose_id"] + batch.feature_names)
        for r, row in zip(batch.results, batch.matrix):
            writer.writerow([r.pose_id] + row)
    logger.info("Wrote IFP matrix (%d×%d) to %s",
                len(batch.matrix), len(batch.feature_names), output_path)


def compute_tanimoto_similarity(bitvec_a: list[int], bitvec_b: list[int]) -> float:
    """Tanimoto similarity between two binary fingerprints.

    T(A,B) = |A∩B| / |A∪B|
    """
    a = np.array(bitvec_a, dtype=bool)
    b = np.array(bitvec_b, dtype=bool)
    intersection = np.sum(a & b)
    union = np.sum(a | b)
    if union == 0:
        return 0.0
    return float(intersection / union)
