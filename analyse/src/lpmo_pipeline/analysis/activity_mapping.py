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

import csv
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

PREDICTIVE_CLUSTER_TABLE_COLUMNS = [
    "analysis_id",
    "protein_id",
    "family",
    "cbm_status",
    "condition_id",
    "construct_type",
    "ligand_id",
    "substrate_class",
    "DP",
    "cluster_id",
    "medoid_pose_id",
    "cluster_type",
    "occupancy",
    "cluster_size",
    "qc_factor",
    "support_factor",
    "cluster_weight",
    "convergence_support",
    "invalid_geometry_fraction",
    "plausible_geometry_fraction",
    "scorable_geometry_fraction",
    "median_oxyl_H_C1",
    "median_oxyl_H_C4",
    "median_Cu_C1",
    "median_Cu_C4",
    "median_attack_angle_C1",
    "median_attack_angle_C4",
    "face_orientation_summary",
    "ifp_features_selected",
    "mean_ranking_score",
    "median_ranking_score",
    "medoid_ranking_score",
    "mean_iptm",
    "median_iptm",
    "medoid_iptm",
    "mean_mean_plddt",
    "median_mean_plddt",
    "medoid_mean_plddt",
    "experimental_regio_label",
    "experimental_substrate_label",
    "experimental_activity_label",
    "activity_mapping_rule",
]


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


def _as_float(value: Any, default: float | None = None) -> float | None:
    if value in (None, ""):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _as_int(value: Any, default: int = 0) -> int:
    numeric = _as_float(value)
    return int(numeric) if numeric is not None else default


def _first_present(row: dict[str, Any], keys: tuple[str, ...], default: str = "") -> str:
    for key in keys:
        value = row.get(key)
        if value not in (None, ""):
            return str(value)
    return default


def _index_by_protein(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    indexed: dict[str, dict[str, Any]] = {}
    for row in rows:
        protein_id = str(row.get("protein_id", "") or row.get("uniprot_id", "")).strip()
        if protein_id and protein_id not in indexed:
            indexed[protein_id] = row
    return indexed


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_predictive_cluster_table(rows: list[dict[str, Any]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=PREDICTIVE_CLUSTER_TABLE_COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in PREDICTIVE_CLUSTER_TABLE_COLUMNS})


def build_activity_annotation_by_protein(
    activity_rows: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    """Index mapped EC/activity rows by protein_id."""

    return _index_by_protein(activity_rows)


def build_predictive_cluster_table_rows(
    cluster_rows: list[dict[str, Any]],
    *,
    protein_metadata_rows: list[dict[str, Any]] | None = None,
    activity_annotation_rows: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """Build the postprocess predictive cluster table from enriched cluster rows.

    The builder consumes the current Stage 7 column names directly, especially
    `Cu_C1_distance_median`/`Cu_C4_distance_median` and related geometry fields.
    """

    metadata_by_protein = _index_by_protein(protein_metadata_rows or [])
    activity_by_protein = build_activity_annotation_by_protein(activity_annotation_rows or [])

    out: list[dict[str, Any]] = []
    for row in cluster_rows:
        protein_id = str(row.get("protein_id", "")).strip()
        metadata = metadata_by_protein.get(protein_id, {})
        activity = activity_by_protein.get(protein_id, {})

        occupancy = _as_float(row.get("occupancy"), 0.0) or 0.0
        convergence_support = _as_float(row.get("convergent_fraction"))
        support_factor = convergence_support if convergence_support is not None else 1.0
        qc_factor = _as_float(row.get("qc_factor"), 1.0) or 1.0
        c1_computable = _as_float(row.get("c1_geometry_computable_fraction"), 0.0) or 0.0
        c4_computable = _as_float(row.get("c4_geometry_computable_fraction"), 0.0) or 0.0
        c1_plausible = _as_float(row.get("c1_geometry_plausible_fraction"))
        if c1_plausible is None:
            c1_plausible = _as_float(row.get("c1_plausible_fraction"), 0.0)
        c4_plausible = _as_float(row.get("c4_geometry_plausible_fraction"))
        if c4_plausible is None:
            c4_plausible = _as_float(row.get("c4_plausible_fraction"), 0.0)
        plausible_geometry_fraction = max(c1_plausible or 0.0, c4_plausible or 0.0)
        scorable_geometry_fraction = max(c1_computable, c4_computable)

        condition_id = str(row.get("condition_id", ""))
        cluster_id = str(row.get("cluster_id", ""))
        out.append(
            {
                "analysis_id": f"{condition_id}__cluster{cluster_id}",
                "protein_id": protein_id,
                "family": _first_present(metadata, ("family", "protein_family", "aa_family")),
                "cbm_status": _first_present(metadata, ("cbm_status", "has_cbm", "cbm_present"), "unknown"),
                "condition_id": condition_id,
                "construct_type": row.get("construct_type", ""),
                "ligand_id": row.get("ligand_id", ""),
                "substrate_class": row.get("substrate_class", ""),
                "DP": row.get("dp", ""),
                "cluster_id": cluster_id,
                "medoid_pose_id": row.get("medoid_pose_id", ""),
                "cluster_type": row.get("cluster_type", ""),
                "occupancy": occupancy,
                "cluster_size": row.get("cluster_size", row.get("n_poses", "")),
                "qc_factor": qc_factor,
                "support_factor": support_factor,
                "cluster_weight": occupancy * qc_factor * support_factor,
                "convergence_support": "" if convergence_support is None else convergence_support,
                "invalid_geometry_fraction": max(0.0, 1.0 - scorable_geometry_fraction),
                "plausible_geometry_fraction": plausible_geometry_fraction,
                "scorable_geometry_fraction": scorable_geometry_fraction,
                "median_oxyl_H_C1": row.get("oxyl_H_C1_distance_median", ""),
                "median_oxyl_H_C4": row.get("oxyl_H_C4_distance_median", ""),
                "median_Cu_C1": row.get("Cu_C1_distance_median", ""),
                "median_Cu_C4": row.get("Cu_C4_distance_median", ""),
                "median_attack_angle_C1": row.get("attack_angle_C1_median", ""),
                "median_attack_angle_C4": row.get("attack_angle_C4_median", ""),
                "face_orientation_summary": row.get("ring_normal_vs_brace_normal_median", ""),
                "ifp_features_selected": row.get("ifp_features_selected", ""),
                "mean_ranking_score": row.get("mean_ranking_score", ""),
                "median_ranking_score": row.get("median_ranking_score", ""),
                "medoid_ranking_score": row.get("medoid_ranking_score", ""),
                "mean_iptm": row.get("mean_iptm", ""),
                "median_iptm": row.get("median_iptm", ""),
                "medoid_iptm": row.get("medoid_iptm", ""),
                "mean_mean_plddt": row.get("mean_mean_plddt", ""),
                "median_mean_plddt": row.get("median_mean_plddt", ""),
                "medoid_mean_plddt": row.get("medoid_mean_plddt", ""),
                "experimental_regio_label": _first_present(
                    activity,
                    ("mapped_regio_class", "regio_class", "regio_label"),
                    "unknown",
                ),
                "experimental_substrate_label": _first_present(
                    activity,
                    ("mapped_substrate_class", "substrate_class", "substrate_label"),
                    "unknown",
                ),
                "experimental_activity_label": _first_present(
                    activity,
                    ("mapped_activity_label", "activity_label"),
                    "unknown",
                ),
                "activity_mapping_rule": activity.get("mapping_rule", ""),
            }
        )

    logger.info("Built %d predictive cluster table rows", len(out))
    return out


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
