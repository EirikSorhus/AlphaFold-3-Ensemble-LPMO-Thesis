# src/lpmo_pipeline/analysis/activity_mapping.py
"""
Responsibility: Map cluster signatures to enzyme activity annotations.
Input:  cluster_signatures.json + activity annotations (from dataset)
Output: activity_features.json

STEP 11 in masterplan.
RQ1: C1 vs C4 occupancy from geometry
RQ2: Substrate-specific binding patterns
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


@dataclass
class ActivityFeature:
    """Activity features extracted for a single protein."""

    protein_id: str = ""
    activity_class: str = ""  # "C1" | "C4" | "C1+C4" | "unknown"
    substrate_class: str = ""  # "chitin" | "cellulose" | "amylose" | "mixed"

    # Geometry-derived
    c1_occupancy: float = 0.0  # fraction of clusters where Cu-C1 < Cu-C4
    c4_occupancy: float = 0.0
    mean_cu_c1_best_cluster: float = 0.0
    mean_cu_c4_best_cluster: float = 0.0

    # Interaction-derived
    signature_interactions: list[dict[str, Any]] = field(default_factory=list)

    # Per-ligand breakdown
    per_ligand: dict[str, dict[str, Any]] = field(default_factory=dict)


def compute_activity_features(
    protein_id: str,
    cluster_signatures_per_ligand: dict[str, list[dict[str, Any]]],
    activity_annotation: dict[str, str] | None = None,
    cu_c1_c4_offset: float = 0.5,
) -> ActivityFeature:
    """Compute activity-related features from cluster signatures.

    C1 occupancy = sum(occupancy) for clusters where mean_Cu-C1 < mean_Cu-C4 - offset
    C4 occupancy = complement

    Args:
        protein_id: Protein identifier.
        cluster_signatures_per_ligand: {ligand_id: [cluster_signature_dict, ...]}
        activity_annotation: Optional {protein_id: "C1"|"C4"|"C1+C4"}
        cu_c1_c4_offset: Difference threshold for C1 vs C4 classification (Å).

    Returns:
        ActivityFeature.
    """
    feat = ActivityFeature(protein_id=protein_id)

    if activity_annotation and protein_id in activity_annotation:
        feat.activity_class = activity_annotation[protein_id]

    total_occupancy = 0.0
    c1_occ = 0.0
    c4_occ = 0.0

    for ligand_id, signatures in cluster_signatures_per_ligand.items():
        lig_c1 = 0.0
        lig_c4 = 0.0

        for sig in signatures:
            occ = sig.get("occupancy", 0.0)
            total_occupancy += occ
            geom = sig.get("geometry", {})
            mean_c1 = geom.get("mean_cu_c1", float("inf"))
            mean_c4 = geom.get("mean_cu_c4", float("inf"))

            if mean_c1 < mean_c4 - cu_c1_c4_offset:
                c1_occ += occ
                lig_c1 += occ
            elif mean_c4 < mean_c1 - cu_c1_c4_offset:
                c4_occ += occ
                lig_c4 += occ
            else:
                # Ambiguous: split evenly
                c1_occ += occ / 2
                c4_occ += occ / 2
                lig_c1 += occ / 2
                lig_c4 += occ / 2

        feat.per_ligand[ligand_id] = {
            "c1_occ": lig_c1,
            "c4_occ": lig_c4,
            "n_clusters": len(signatures),
        }

    if total_occupancy > 0:
        feat.c1_occupancy = c1_occ / total_occupancy
        feat.c4_occupancy = c4_occ / total_occupancy

    logger.info(
        "Activity features for %s: C1_occ=%.2f, C4_occ=%.2f, class=%s",
        protein_id, feat.c1_occupancy, feat.c4_occupancy, feat.activity_class,
    )
    return feat


def write_activity_features(
    features: list[ActivityFeature], output_path: Path,
) -> None:
    """Write activity features to JSON."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    data = [
        {
            "protein_id": f.protein_id,
            "activity_class": f.activity_class,
            "substrate_class": f.substrate_class,
            "c1_occupancy": f.c1_occupancy,
            "c4_occupancy": f.c4_occupancy,
            "per_ligand": f.per_ligand,
            "signature_interactions": f.signature_interactions,
        }
        for f in features
    ]
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    logger.info("Wrote %d activity features to %s", len(features), output_path)
