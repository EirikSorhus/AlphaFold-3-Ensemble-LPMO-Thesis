# src/lpmo_pipeline/analysis/activity_mapping.py
"""
Responsibility: Build activity-linked analysis tables from cluster annotations.
Input:  cluster_signatures/cluster_table + activity annotations (from metadata)
Output: predictive_cluster_table (primary) and optional enzyme summaries (secondary)

STEP 11 in masterplan.
RQ1: C1 vs C4 occupancy from geometry
RQ2: Substrate-specific binding patterns

Important design rule:
    - Main analyses are cluster-primary.
    - Enzyme-level aggregation is secondary sensitivity output only.
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


@dataclass
class ClusterActivityRow:
    """Primary predictive row: one row per cluster."""

    analysis_id: str = ""
    protein_id: str = ""
    family: str = ""
    cbm_status: str = "unknown"
    substrate_class: str = ""
    dp: int = 0
    cluster_id: str = ""
    occupancy: float = 0.0
    qc_factor: float = 1.0
    support_factor: float = 1.0
    cluster_weight: float = 0.0
    median_cu_c1: float = 0.0
    median_cu_c4: float = 0.0
    median_oxyl_h_c1: float = 0.0
    median_oxyl_h_c4: float = 0.0
    experimental_regio_label: str = "unknown"
    experimental_ligand_specificity_label: str = "unknown"
    extra_features: dict[str, Any] = field(default_factory=dict)


@dataclass
class ActivityFeature:
    """Secondary enzyme-level sensitivity summary."""

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


def build_predictive_cluster_rows(
    cluster_rows: list[dict[str, Any]],
    activity_annotation: dict[str, dict[str, str]],
) -> list[ClusterActivityRow]:
    """Build primary predictive rows from cluster-level annotations.

    Args:
        cluster_rows: Rows from cluster_table-like source.
        activity_annotation: {protein_id: {regio_label, ligand_specificity_label, ...}}

    Returns:
        List of ClusterActivityRow objects (one row per cluster).
    """
    out: list[ClusterActivityRow] = []
    for row in cluster_rows:
        protein_id = str(row.get("protein_id", ""))
        ann = activity_annotation.get(protein_id, {})
        occupancy = float(row.get("occupancy", 0.0))
        qc_factor = float(row.get("qc_factor", 1.0))
        support_factor = float(row.get("support_factor", 1.0))

        out.append(
            ClusterActivityRow(
                analysis_id=str(row.get("analysis_id", "")),
                protein_id=protein_id,
                family=str(row.get("family", "")),
                cbm_status=str(row.get("cbm_status", "unknown")),
                substrate_class=str(row.get("substrate_class", "")),
                dp=int(row.get("DP", row.get("dp", 0)) or 0),
                cluster_id=str(row.get("cluster_id", "")),
                occupancy=occupancy,
                qc_factor=qc_factor,
                support_factor=support_factor,
                cluster_weight=occupancy * qc_factor * support_factor,
                median_cu_c1=float(row.get("median_cu_c1", 0.0)),
                median_cu_c4=float(row.get("median_cu_c4", 0.0)),
                median_oxyl_h_c1=float(row.get("median_oxyl_h_c1", 0.0)),
                median_oxyl_h_c4=float(row.get("median_oxyl_h_c4", 0.0)),
                experimental_regio_label=str(ann.get("regio_label", "unknown")),
                experimental_ligand_specificity_label=str(
                    ann.get("ligand_specificity_label", "unknown")
                ),
                extra_features={
                    "cross_model_support": row.get("cross_model_support"),
                    "convergence_support": row.get("convergence_support"),
                    "plausible_geometry_fraction": row.get("plausible_geometry_fraction"),
                    "ifp_features_selected": row.get("ifp_features_selected"),
                },
            )
        )

    logger.info("Built %d predictive cluster rows", len(out))
    return out


def compute_activity_features(
    protein_id: str,
    cluster_signatures_per_ligand: dict[str, list[dict[str, Any]]],
    activity_annotation: dict[str, str] | None = None,
    cu_c1_c4_offset: float = 0.5,
) -> ActivityFeature:
    """Compute secondary enzyme-level summary features from cluster signatures.

    This function intentionally aggregates cluster-level signals and is only meant
    for sensitivity analyses or compact summaries. It is not the primary modeling
    representation.

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
