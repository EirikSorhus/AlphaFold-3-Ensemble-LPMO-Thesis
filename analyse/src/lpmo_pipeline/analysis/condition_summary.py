"""Condition- and protein-level summary table builders."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np


_CLUSTER_WEIGHTED_FIELDS = [
    "c1_plausible_fraction",
    "c4_plausible_fraction",
    "Cu_C1_distance_median",
    "Cu_C4_distance_median",
    "oxyl_H_C1_distance_median",
    "oxyl_H_C4_distance_median",
    "Cu_oxyl_H_C1_angle_median",
    "Cu_oxyl_H_C4_angle_median",
    "ring_normal_vs_brace_normal_median",
    "oxyl_H_score_C1_median",
    "oxyl_H_score_C4_median",
    "his_brace_angle_deg_median",
]

_CONFIDENCE_AGG_FIELDS = [
    "ranking_score",
    "iptm",
    "ptm",
    "mean_plddt",
    "ligand_interface_confidence",
    "pae_summary",
]

CONDITION_TABLE_COLUMNS = [
    "condition_id",
    "protein_id",
    "construct_type",
    "ligand_id",
    "substrate_class",
    "dp",
    "n_generated",
    "n_prepared",
    "n_prep_error",
    "n_stage1_hard_fail",
    "n_stage1_soft_flag",
    "n_stage1_pass",
    "hard_fail_rate",
    "n_ifp_success",
    "n_contact_eligible",
    "contact_eligible_fraction",
    "minimum_clusterable_n",
    "clustering_status",
    "formal_clustering_allowed",
    "null_ifp_fraction",
    "vdw_only_fraction",
    "low_specific_contact_fraction",
    "median_n_non_vdw_interactions",
    "median_n_non_vdw_contact_residues",
    "n_ifp_clustered",
    "n_noise",
    "noise_fraction",
    "n_clusters",
    "top_cluster_occupancy",
    "cluster_entropy",
    "occupancy_gini",
    "convergence_reference_pose_id",
    "convergence_n_qc_pass_poses",
    "convergence_fraction",
    "median_ligand_rmsd",
    "iqr_ligand_rmsd",
    "low_convergence_flag",
    "any_valid_cluster",
    "n_clusters_considered",
    "total_nonnoise_cluster_occupancy",
    "weighted_contact_feature_mass",
    "aromatic_contact_fraction",
    "polar_contact_fraction",
    "charged_contact_fraction",
    "hbond_contact_fraction",
    "catalytic_surface_contact_fraction",
    "non_core_contact_fraction",
    "cbm_contact_fraction",
    "linker_contact_fraction",
    "n_confidence_rows",
    "n_confidence_available",
    *[
        column
        for field in _CONFIDENCE_AGG_FIELDS
        for column in (f"mean_{field}", f"median_{field}")
    ],
    "n_cluster_rows",
    "cluster_total_occupancy",
    "mean_cluster_size",
    *[f"occupancy_weighted_{field}" for field in _CLUSTER_WEIGHTED_FIELDS],
]

PROTEIN_SUMMARY_TABLE_COLUMNS = [
    "protein_id",
    "n_conditions",
    "n_construct_types",
    "construct_types",
    "n_ligands",
    "n_substrate_classes",
    "substrate_classes",
    "n_generated",
    "n_prepared",
    "n_stage1_hard_fail",
    "n_stage1_pass",
    "hard_fail_rate",
    "n_conditions_with_ifp_success",
    "n_conditions_with_contact_eligible",
    "n_conditions_with_clusters",
    "mean_contact_eligible_fraction",
    "mean_vdw_only_fraction",
    "mean_low_specific_contact_fraction",
    "mean_n_clusters",
    "mean_top_cluster_occupancy",
    "mean_cluster_entropy",
    "mean_convergence_fraction",
    "mean_median_ligand_rmsd",
    "mean_total_nonnoise_cluster_occupancy",
    "mean_aromatic_contact_fraction",
    "mean_polar_contact_fraction",
    "mean_charged_contact_fraction",
    "mean_hbond_contact_fraction",
    "mean_weighted_Cu_C1_distance_median",
    "mean_weighted_Cu_C4_distance_median",
    "mean_weighted_oxyl_H_C1_distance_median",
    "mean_weighted_oxyl_H_C4_distance_median",
    "mean_confidence_ranking_score",
    "mean_confidence_iptm",
    "mean_confidence_mean_plddt",
]


def read_tsv(path: Path | None) -> list[dict[str, str]]:
    if path is None or not Path(path).exists():
        return []
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_condition_table(rows: list[dict[str, Any]], output_path: Path) -> None:
    _write_tsv(output_path, CONDITION_TABLE_COLUMNS, rows)


def write_protein_summary_table(rows: list[dict[str, Any]], output_path: Path) -> None:
    _write_tsv(output_path, PROTEIN_SUMMARY_TABLE_COLUMNS, rows)


def _write_tsv(path: Path, columns: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def _as_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if np.isnan(numeric):
        return None
    return numeric


def _as_int(value: Any) -> int:
    numeric = _as_float(value)
    return int(numeric) if numeric is not None else 0


def _mean(values: list[float]) -> float | str:
    return float(np.mean(values)) if values else ""


def _median(values: list[float]) -> float | str:
    return float(np.median(values)) if values else ""


def _group_by(rows: list[dict[str, Any]], key: str) -> dict[str, list[dict[str, Any]]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row.get(key, ""))].append(row)
    return dict(grouped)


def _index_first(rows: list[dict[str, Any]], key: str) -> dict[str, dict[str, Any]]:
    indexed: dict[str, dict[str, Any]] = {}
    for row in rows:
        value = str(row.get(key, ""))
        if value and value not in indexed:
            indexed[value] = row
    return indexed


def _confidence_aggregates(rows: list[dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {
        "n_confidence_rows": len(rows),
        "n_confidence_available": sum(
            1 for row in rows if row.get("confidence_json_status") not in {"missing", ""}
        ),
    }
    for field in _CONFIDENCE_AGG_FIELDS:
        values = [
            numeric
            for row in rows
            if (numeric := _as_float(row.get(field))) is not None
        ]
        out[f"mean_{field}"] = _mean(values)
        out[f"median_{field}"] = _median(values)
    return out


def _cluster_aggregates(rows: list[dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {
        "n_cluster_rows": len(rows),
        "cluster_total_occupancy": "",
        "mean_cluster_size": "",
    }
    occupancies = [
        numeric
        for row in rows
        if (numeric := _as_float(row.get("occupancy"))) is not None
    ]
    sizes = [
        numeric
        for row in rows
        if (numeric := _as_float(row.get("n_poses"))) is not None
    ]
    total_occupancy = sum(occupancies)
    out["cluster_total_occupancy"] = total_occupancy if rows else ""
    out["mean_cluster_size"] = _mean(sizes)

    for field in _CLUSTER_WEIGHTED_FIELDS:
        weighted_sum = 0.0
        weight_sum = 0.0
        unweighted_values: list[float] = []
        for row in rows:
            value = _as_float(row.get(field))
            if value is None:
                continue
            unweighted_values.append(value)
            weight = _as_float(row.get("occupancy")) or 0.0
            weighted_sum += value * weight
            weight_sum += weight
        if weight_sum:
            out[f"occupancy_weighted_{field}"] = weighted_sum / weight_sum
        else:
            out[f"occupancy_weighted_{field}"] = _mean(unweighted_values)
    return out


def build_condition_table_rows(
    *,
    qc_attrition_rows: list[dict[str, Any]],
    condition_cluster_summary_rows: list[dict[str, Any]] | None = None,
    condition_convergence_summary_rows: list[dict[str, Any]] | None = None,
    condition_patch_summary_rows: list[dict[str, Any]] | None = None,
    pose_confidence_rows: list[dict[str, Any]] | None = None,
    cluster_table_rows: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """Build one condition summary row per QC attrition condition."""

    condition_cluster_summary_rows = condition_cluster_summary_rows or []
    condition_convergence_summary_rows = condition_convergence_summary_rows or []
    condition_patch_summary_rows = condition_patch_summary_rows or []
    pose_confidence_rows = pose_confidence_rows or []
    cluster_table_rows = cluster_table_rows or []

    cluster_summary_by_condition = _index_first(condition_cluster_summary_rows, "condition_id")
    convergence_by_condition = _index_first(condition_convergence_summary_rows, "condition_id")
    patch_by_condition = _index_first(condition_patch_summary_rows, "condition_id")
    confidence_by_condition = _group_by(pose_confidence_rows, "condition_id")
    clusters_by_condition = _group_by(cluster_table_rows, "condition_id")

    rows: list[dict[str, Any]] = []
    for qc_row in sorted(qc_attrition_rows, key=lambda row: str(row.get("condition_id", ""))):
        condition_id = str(qc_row.get("condition_id", ""))
        cluster_row = cluster_summary_by_condition.get(condition_id, {})
        convergence_row = convergence_by_condition.get(condition_id, {})
        patch_row = patch_by_condition.get(condition_id, {})

        row: dict[str, Any] = {
            "condition_id": condition_id,
            "protein_id": qc_row.get("protein_id", ""),
            "construct_type": qc_row.get("construct_type", ""),
            "ligand_id": qc_row.get("ligand_id", ""),
            "substrate_class": qc_row.get("substrate_class", ""),
            "dp": qc_row.get("dp", ""),
        }
        for field in [
            "n_generated",
            "n_prepared",
            "n_prep_error",
            "n_stage1_hard_fail",
            "n_stage1_soft_flag",
            "n_stage1_pass",
            "hard_fail_rate",
        ]:
            row[field] = qc_row.get(field, "")

        for field in [
            "n_ifp_success",
            "n_contact_eligible",
            "contact_eligible_fraction",
            "minimum_clusterable_n",
            "clustering_status",
            "formal_clustering_allowed",
            "null_ifp_fraction",
            "vdw_only_fraction",
            "low_specific_contact_fraction",
            "median_n_non_vdw_interactions",
            "median_n_non_vdw_contact_residues",
            "n_ifp_clustered",
            "n_noise",
            "noise_fraction",
            "n_clusters",
            "top_cluster_occupancy",
            "cluster_entropy",
            "occupancy_gini",
        ]:
            row[field] = cluster_row.get(field, 0 if field.startswith("n_") else "")

        row.update(
            {
                "convergence_reference_pose_id": convergence_row.get("reference_pose_id", ""),
                "convergence_n_qc_pass_poses": convergence_row.get("n_qc_pass_poses", ""),
                "convergence_fraction": convergence_row.get("convergence_fraction", ""),
                "median_ligand_rmsd": convergence_row.get("median_ligand_rmsd", ""),
                "iqr_ligand_rmsd": convergence_row.get("iqr_ligand_rmsd", ""),
                "low_convergence_flag": convergence_row.get("low_convergence_flag", ""),
            }
        )

        for field in [
            "any_valid_cluster",
            "n_clusters_considered",
            "total_nonnoise_cluster_occupancy",
            "weighted_contact_feature_mass",
            "aromatic_contact_fraction",
            "polar_contact_fraction",
            "charged_contact_fraction",
            "hbond_contact_fraction",
            "catalytic_surface_contact_fraction",
            "non_core_contact_fraction",
            "cbm_contact_fraction",
            "linker_contact_fraction",
        ]:
            row[field] = patch_row.get(field, "")

        row.update(_confidence_aggregates(confidence_by_condition.get(condition_id, [])))
        row.update(_cluster_aggregates(clusters_by_condition.get(condition_id, [])))
        rows.append(row)
    return rows


def build_protein_summary_rows(
    condition_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Build protein-level summaries as a pure groupby over condition rows."""

    grouped = _group_by(condition_rows, "protein_id")
    rows: list[dict[str, Any]] = []
    for protein_id in sorted(grouped):
        conditions = grouped[protein_id]
        constructs = sorted({str(row.get("construct_type", "")) for row in conditions if row.get("construct_type")})
        substrates = sorted({str(row.get("substrate_class", "")) for row in conditions if row.get("substrate_class")})
        ligands = sorted({str(row.get("ligand_id", "")) for row in conditions if row.get("ligand_id")})

        n_generated = sum(_as_int(row.get("n_generated")) for row in conditions)
        n_stage1_hard_fail = sum(_as_int(row.get("n_stage1_hard_fail")) for row in conditions)
        row: dict[str, Any] = {
            "protein_id": protein_id,
            "n_conditions": len(conditions),
            "n_construct_types": len(constructs),
            "construct_types": ",".join(constructs),
            "n_ligands": len(ligands),
            "n_substrate_classes": len(substrates),
            "substrate_classes": ",".join(substrates),
            "n_generated": n_generated,
            "n_prepared": sum(_as_int(row.get("n_prepared")) for row in conditions),
            "n_stage1_hard_fail": n_stage1_hard_fail,
            "n_stage1_pass": sum(_as_int(row.get("n_stage1_pass")) for row in conditions),
            "hard_fail_rate": n_stage1_hard_fail / n_generated if n_generated else 0.0,
            "n_conditions_with_ifp_success": sum(
                1 for row in conditions if _as_int(row.get("n_ifp_success")) > 0
            ),
            "n_conditions_with_contact_eligible": sum(
                1 for row in conditions if _as_int(row.get("n_contact_eligible")) > 0
            ),
            "n_conditions_with_clusters": sum(
                1 for row in conditions if _as_int(row.get("n_clusters")) > 0
            ),
        }

        mean_mappings = {
            "mean_contact_eligible_fraction": "contact_eligible_fraction",
            "mean_vdw_only_fraction": "vdw_only_fraction",
            "mean_low_specific_contact_fraction": "low_specific_contact_fraction",
            "mean_n_clusters": "n_clusters",
            "mean_top_cluster_occupancy": "top_cluster_occupancy",
            "mean_cluster_entropy": "cluster_entropy",
            "mean_convergence_fraction": "convergence_fraction",
            "mean_median_ligand_rmsd": "median_ligand_rmsd",
            "mean_total_nonnoise_cluster_occupancy": "total_nonnoise_cluster_occupancy",
            "mean_aromatic_contact_fraction": "aromatic_contact_fraction",
            "mean_polar_contact_fraction": "polar_contact_fraction",
            "mean_charged_contact_fraction": "charged_contact_fraction",
            "mean_hbond_contact_fraction": "hbond_contact_fraction",
            "mean_weighted_Cu_C1_distance_median": "occupancy_weighted_Cu_C1_distance_median",
            "mean_weighted_Cu_C4_distance_median": "occupancy_weighted_Cu_C4_distance_median",
            "mean_weighted_oxyl_H_C1_distance_median": "occupancy_weighted_oxyl_H_C1_distance_median",
            "mean_weighted_oxyl_H_C4_distance_median": "occupancy_weighted_oxyl_H_C4_distance_median",
            "mean_confidence_ranking_score": "mean_ranking_score",
            "mean_confidence_iptm": "mean_iptm",
            "mean_confidence_mean_plddt": "mean_mean_plddt",
        }
        for output_field, source_field in mean_mappings.items():
            values = [
                numeric
                for condition in conditions
                if (numeric := _as_float(condition.get(source_field))) is not None
            ]
            row[output_field] = _mean(values)

        rows.append(row)
    return rows
