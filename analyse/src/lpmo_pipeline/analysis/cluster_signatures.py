# src/lpmo_pipeline/analysis/cluster_signatures.py
"""
Responsibility: Extract geometry and interaction signatures per cluster.
Input:  cluster assignments (from HDBSCAN) + geometry_metrics + IFP data
Output: cluster_signatures.json conforming to schemas/cluster_signatures_schema.json

STEP 9 in masterplan: For each cluster, select medoid, aggregate geometry + interactions.
"""
from __future__ import annotations

import json
import logging
import csv
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from lpmo_pipeline.analysis.prolif_ifp import IFPResult, parse_ifp_feature_name
from lpmo_pipeline.config import load_runtime_paths_config

logger = logging.getLogger(__name__)

_RUNTIME_PATHS = load_runtime_paths_config()
_THRESHOLDS_PATH = _RUNTIME_PATHS.pipeline_assets.thresholds_config

CLUSTER_IFP_SIGNATURE_COLUMNS = [
    "condition_id",
    "cluster_id",
    "protein_id",
    "ligand_id",
    "cluster_type",
    "n_poses",
    "occupancy",
    "medoid_pose_id",
    "feature_name",
    "ligand_residue_label",
    "protein_residue_label",
    "interaction_type",
    "n_poses_with_contact",
    "contact_frequency",
]

CLUSTER_RESIDUE_SIGNATURE_COLUMNS = [
    "condition_id",
    "cluster_id",
    "protein_id",
    "ligand_id",
    "cluster_type",
    "n_poses",
    "occupancy",
    "medoid_pose_id",
    "residue_chain",
    "residue_number",
    "residue_name",
    "interaction_type",
    "ligand_residue_label",
    "n_poses_with_contact",
    "contact_frequency",
    "is_catalytic_surface_region",
    "is_cbm_region",
    "is_linker_region",
]

_GEOMETRY_STAT_FIELDS = [
    "Cu_C1_distance",
    "Cu_C4_distance",
    "oxyl_H_C1_distance",
    "oxyl_H_C4_distance",
    "attack_angle_C1",
    "attack_angle_C4",
    "ring_normal_vs_brace_normal",
    "oxyl_H_score_C1",
    "oxyl_H_score_C4",
    "his_brace_angle_deg",
    "core_rmsd_vs_reference",
    "pocket_rmsd_vs_crystal",
]

CLUSTER_TABLE_COLUMNS = [
    "condition_id",
    "cluster_id",
    "protein_id",
    "ligand_id",
    "cluster_type",
    "n_poses",
    "occupancy",
    "medoid_pose_id",
    "c1_plausible_fraction",
    "c4_plausible_fraction",
    *[
        column_name
        for field_name in _GEOMETRY_STAT_FIELDS
        for column_name in (f"{field_name}_median", f"{field_name}_iqr")
    ],
]

_PLAUSIBLE_GEOMETRY_STATUSES = {"geometry_plausible", "geometry_highly_plausible"}


def _geometry_value(metrics: dict[str, Any], *keys: str) -> float | None:
    """Read the first non-null value from legacy or pose_geometry row keys."""
    for key in keys:
        value = metrics.get(key)
        if value is not None:
            return value
    return None


@dataclass
class ClusterSignature:
    """Aggregated signature for a single cluster."""

    cluster_id: int
    n_poses: int = 0
    occupancy: float = 0.0  # Fraction of total poses in this cluster
    medoid_pose_id: str = ""

    # Geometry (mean ± std across cluster members)
    mean_cu_c1: float = 0.0
    std_cu_c1: float = 0.0
    mean_cu_c4: float = 0.0
    std_cu_c4: float = 0.0
    mean_his_brace_angle: float = 0.0
    std_his_brace_angle: float = 0.0

    # Dominant interactions (top-N by frequency)
    dominant_interactions: list[dict[str, Any]] = field(default_factory=list)
    # [{feature: str, frequency: float, residue: str, interaction_type: str}]

    # Raw member IDs for traceability
    member_pose_ids: list[str] = field(default_factory=list)


def compute_cluster_signatures(
    cluster_labels: dict[str, int],  # {pose_id: cluster_label}
    geometry_metrics: dict[str, dict[str, Any]],  # {pose_id: metrics}
    ifp_matrix: dict[str, list[int]],  # {pose_id: flat_bitvector}
    feature_names: list[str],
    total_poses: int,
    top_n_interactions: int = 5,
) -> list[ClusterSignature]:
    """Compute signatures for all clusters.

    Args:
        cluster_labels: {pose_id: cluster_id} where -1 = outlier.
        geometry_metrics: {pose_id: {min_cu_c1, min_cu_c4, his_brace_angle_deg, ...}}
        ifp_matrix: {pose_id: [0,1,0,...]} binary fingerprint.
        feature_names: Ordered list of feature names matching ifp_matrix columns.
        total_poses: Total poses for occupancy calculation.
        top_n_interactions: Number of dominant interactions to report per cluster.

    Returns:
        List of ClusterSignature (excludes cluster -1 = outliers).
    """
    # Group poses by cluster
    clusters: dict[int, list[str]] = {}
    for pose_id, label in cluster_labels.items():
        clusters.setdefault(label, []).append(pose_id)

    signatures: list[ClusterSignature] = []

    for cid, members in sorted(clusters.items()):
        if cid == -1:
            logger.info("Outlier cluster (-1): %d poses", len(members))
            continue

        sig = ClusterSignature(
            cluster_id=cid,
            n_poses=len(members),
            occupancy=len(members) / total_poses if total_poses > 0 else 0.0,
            member_pose_ids=members,
        )

        # --- Geometry aggregation ---
        cu_c1_vals = [
            _geometry_value(geometry_metrics[pid], "min_cu_c1", "Cu_C1_distance",)
            for pid in members
            if pid in geometry_metrics
        ]
        cu_c4_vals = [
            _geometry_value(geometry_metrics[pid], "min_cu_c4", "Cu_C4_distance",)
            for pid in members
            if pid in geometry_metrics
        ]
        angle_vals = [
            geometry_metrics[pid].get("his_brace_angle_deg", 0.0)
            for pid in members
            if pid in geometry_metrics
            and geometry_metrics[pid].get("his_brace_angle_deg") is not None
        ]

        if cu_c1_vals:
            arr = np.array([v for v in cu_c1_vals if v is not None])
            if len(arr) > 0:
                sig.mean_cu_c1 = float(np.mean(arr))
                sig.std_cu_c1 = float(np.std(arr))

        if cu_c4_vals:
            arr = np.array([v for v in cu_c4_vals if v is not None])
            if len(arr) > 0:
                sig.mean_cu_c4 = float(np.mean(arr))
                sig.std_cu_c4 = float(np.std(arr))

        if angle_vals:
            sig.mean_his_brace_angle = float(np.mean(angle_vals))
            sig.std_his_brace_angle = float(np.std(angle_vals))

        # --- Medoid: pose closest to cluster centroid in IFP space ---
        cluster_ifps = [
            np.array(ifp_matrix[pid], dtype=float)
            for pid in members
            if pid in ifp_matrix
        ]
        if cluster_ifps:
            centroid = np.mean(cluster_ifps, axis=0)
            dists_to_centroid = [np.linalg.norm(v - centroid) for v in cluster_ifps]
            medoid_idx = int(np.argmin(dists_to_centroid))
            filtered_members = [pid for pid in members if pid in ifp_matrix]
            sig.medoid_pose_id = filtered_members[medoid_idx]

        # --- Dominant interactions ---
        if cluster_ifps and feature_names:
            freq = np.mean(cluster_ifps, axis=0)  # Frequency of each feature across cluster
            top_indices = np.argsort(freq)[::-1][:top_n_interactions]
            for idx in top_indices:
                if freq[idx] > 0:
                    fname = feature_names[idx]
                    # Parse residue_interactiontype from feature name
                    parts = fname.rsplit("_", 1)
                    sig.dominant_interactions.append({
                        "feature": fname,
                        "frequency": float(freq[idx]),
                        "residue": parts[0] if len(parts) == 2 else fname,
                        "interaction_type": parts[1] if len(parts) == 2 else "",
                    })

        signatures.append(sig)

    logger.info("Computed %d cluster signatures", len(signatures))
    return signatures


@dataclass(frozen=True)
class ClusterSignatureTables:
    """Stage 7 cluster annotation tables built from production Stage 6 outputs."""

    cluster_summaries: list[dict[str, Any]]
    ifp_signature_rows: list[dict[str, Any]]
    residue_signature_rows: list[dict[str, Any]]


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


def _as_cluster_id(value: Any) -> int:
    if isinstance(value, (int, np.integer)):
        return int(value)
    return int(str(value).strip())


def _as_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(number):
        return None
    return number


def _median_iqr(values: list[Any]) -> tuple[float | None, float | None]:
    numeric = np.asarray(
        [value for value in (_as_float(value) for value in values) if value is not None],
        dtype=float,
    )
    if numeric.size == 0:
        return None, None
    q1, q3 = np.percentile(numeric, [25, 75])
    return float(np.median(numeric)), float(q3 - q1)


def _load_cluster_type_thresholds(config_path: Path | None = None) -> dict[str, float]:
    raw_config = yaml.safe_load((config_path or _THRESHOLDS_PATH).read_text()) or {}
    cluster_type = raw_config.get("cluster_type") or {}
    return {
        "c1_min": float(cluster_type.get("c1_compatible_min_plausible_fraction", 0.50)),
        "c4_min": float(cluster_type.get("c4_compatible_min_plausible_fraction", 0.50)),
        "mixed_min": float(cluster_type.get("mixed_compatible_min_plausible_fraction", 0.40)),
        "non_plausible_max": float(cluster_type.get("non_plausible_max_plausible_fraction", 0.10)),
        "uncertain_max_occupancy": float(
            cluster_type.get("uncertain_max_occupancy_for_classification", 0.05)
        ),
    }


def _plausible_fraction(rows: list[dict[str, Any]], key: str, denominator: int) -> float:
    if denominator <= 0:
        return 0.0
    return (
        sum(1 for row in rows if str(row.get(key, "")).strip() in _PLAUSIBLE_GEOMETRY_STATUSES)
        / denominator
    )


def _classify_cluster_type(
    *,
    occupancy: float,
    c1_plausible_fraction: float,
    c4_plausible_fraction: float,
    thresholds: dict[str, float],
) -> str:
    if occupancy <= thresholds["uncertain_max_occupancy"]:
        return "uncertain"
    if (
        c1_plausible_fraction <= thresholds["non_plausible_max"]
        and c4_plausible_fraction <= thresholds["non_plausible_max"]
    ):
        return "non_plausible"
    if (
        c1_plausible_fraction >= thresholds["c1_min"]
        and c4_plausible_fraction >= thresholds["c4_min"]
    ) or min(c1_plausible_fraction, c4_plausible_fraction) >= thresholds["mixed_min"]:
        return "mixed_compatible"
    if c1_plausible_fraction >= thresholds["c1_min"] and c1_plausible_fraction > c4_plausible_fraction:
        return "C1_compatible"
    if c4_plausible_fraction >= thresholds["c4_min"] and c4_plausible_fraction > c1_plausible_fraction:
        return "C4_compatible"
    return "uncertain"


def _cluster_members_from_assignments(
    cluster_assignment_rows: list[dict[str, Any]],
) -> tuple[dict[tuple[str, int], list[str]], dict[str, int]]:
    members_by_cluster: dict[tuple[str, int], list[str]] = {}
    total_clustered_by_condition: dict[str, int] = {}

    for row in cluster_assignment_rows:
        condition_id = str(row["condition_id"])
        cluster_id = _as_cluster_id(row["cluster_id"])
        total_clustered_by_condition[condition_id] = total_clustered_by_condition.get(condition_id, 0) + 1
        if cluster_id == -1 or not _as_bool(row.get("cluster_member_flag", cluster_id != -1)):
            continue
        members_by_cluster.setdefault((condition_id, cluster_id), []).append(str(row["pose_id"]))

    return members_by_cluster, total_clustered_by_condition


def _medoids_by_cluster(medoid_rows: list[dict[str, Any]]) -> dict[tuple[str, int], str]:
    return {
        (str(row["condition_id"]), _as_cluster_id(row["cluster_id"])): str(row["medoid_pose_id"])
        for row in medoid_rows
    }


def _ifp_frequency_rows(
    *,
    cluster_summary: dict[str, Any],
    members: list[str],
    ifp_by_pose_id: dict[str, IFPResult],
) -> list[dict[str, Any]]:
    feature_counts: dict[str, int] = {}
    n_ifp_poses = 0

    for pose_id in members:
        result = ifp_by_pose_id.get(pose_id)
        if result is None or result.status != "ok" or not result.feature_names:
            continue
        if len(result.feature_names) != len(result.flat_bitvector):
            raise ValueError(f"IFP feature/vector length mismatch for pose {pose_id!r}")
        n_ifp_poses += 1
        for feature_name, bit in zip(result.feature_names, result.flat_bitvector, strict=True):
            feature_counts.setdefault(feature_name, 0)
            feature_counts[feature_name] += int(bit)

    rows: list[dict[str, Any]] = []
    if n_ifp_poses == 0:
        return rows

    for feature_name, active_count in sorted(feature_counts.items()):
        ligand_residue, protein_residue, interaction_type = parse_ifp_feature_name(feature_name)
        rows.append(
            {
                **cluster_summary,
                "feature_name": feature_name,
                "ligand_residue_label": ligand_residue,
                "protein_residue_label": protein_residue,
                "interaction_type": interaction_type,
                "n_poses_with_contact": active_count,
                "contact_frequency": active_count / n_ifp_poses,
            }
        )
    return rows


def _residue_frequency_rows(
    *,
    cluster_summary: dict[str, Any],
    members: list[str],
    residue_contact_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    member_set = set(members)
    active_pose_ids_by_key: dict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    labels_by_key: dict[tuple[str, str, str, str], dict[str, set[str] | bool]] = {}

    for row in residue_contact_rows:
        pose_id = str(row.get("pose_id", ""))
        if pose_id not in member_set:
            continue
        key = (
            str(row.get("residue_chain", "")),
            str(row.get("residue_number", "")),
            str(row.get("residue_name", "")),
            str(row.get("condition_id", "")),
        )
        entry = labels_by_key.setdefault(
            key,
            {
                "interaction_types": set(),
                "ligand_residue_labels": set(),
                "is_catalytic_surface_region": False,
                "is_cbm_region": False,
                "is_linker_region": False,
            },
        )
        if int(row.get("contact_present") or 0):
            active_pose_ids_by_key[key].add(pose_id)
            interaction_type = str(row.get("interaction_type", ""))
            ligand_residue_label = str(row.get("ligand_residue_label", ""))
            if interaction_type:
                entry["interaction_types"].add(interaction_type)
            if ligand_residue_label:
                entry["ligand_residue_labels"].add(ligand_residue_label)
        entry["is_catalytic_surface_region"] = entry["is_catalytic_surface_region"] or _as_bool(
            row.get("is_catalytic_surface_region")
        )
        entry["is_cbm_region"] = entry["is_cbm_region"] or _as_bool(row.get("is_cbm_region"))
        entry["is_linker_region"] = entry["is_linker_region"] or _as_bool(row.get("is_linker_region"))

    rows: list[dict[str, Any]] = []
    for key, entry in sorted(labels_by_key.items()):
        residue_chain, residue_number, residue_name, _condition_id = key
        active_pose_ids = active_pose_ids_by_key.get(key, set())
        rows.append(
            {
                **cluster_summary,
                "residue_chain": residue_chain,
                "residue_number": residue_number,
                "residue_name": residue_name,
                "interaction_type": ",".join(sorted(entry["interaction_types"])) or "none",
                "ligand_residue_label": ",".join(sorted(entry["ligand_residue_labels"])) or "none",
                "n_poses_with_contact": len(active_pose_ids),
                "contact_frequency": len(active_pose_ids) / len(member_set) if member_set else 0.0,
                "is_catalytic_surface_region": entry["is_catalytic_surface_region"],
                "is_cbm_region": entry["is_cbm_region"],
                "is_linker_region": entry["is_linker_region"],
            }
        )
    return rows


def build_cluster_signature_tables(
    *,
    cluster_assignment_rows: list[dict[str, Any]],
    medoid_rows: list[dict[str, Any]],
    geometry_rows_by_pose_id: dict[str, dict[str, Any]],
    ifp_results: list[IFPResult],
    residue_contact_rows: list[dict[str, Any]],
    pose_metadata_by_id: dict[str, dict[str, str]],
    thresholds_path: Path | None = None,
) -> ClusterSignatureTables:
    """Build Stage 7 cluster-level IFP and residue-contact signature tables."""

    thresholds = _load_cluster_type_thresholds(thresholds_path)
    members_by_cluster, total_clustered_by_condition = _cluster_members_from_assignments(
        cluster_assignment_rows
    )
    medoid_by_cluster = _medoids_by_cluster(medoid_rows)
    ifp_by_pose_id = {result.pose_id: result for result in ifp_results}

    cluster_summaries: list[dict[str, Any]] = []
    ifp_signature_rows: list[dict[str, Any]] = []
    residue_signature_rows: list[dict[str, Any]] = []

    for (condition_id, cluster_id), members in sorted(members_by_cluster.items()):
        first_metadata = pose_metadata_by_id.get(members[0], {})
        n_poses = len(members)
        total_clustered = total_clustered_by_condition.get(condition_id, n_poses)
        occupancy = n_poses / total_clustered if total_clustered else 0.0
        geometry_rows = [
            geometry_rows_by_pose_id[pose_id]
            for pose_id in members
            if pose_id in geometry_rows_by_pose_id
        ]
        c1_plausible_fraction = _plausible_fraction(geometry_rows, "geometry_status_C1", n_poses)
        c4_plausible_fraction = _plausible_fraction(geometry_rows, "geometry_status_C4", n_poses)
        cluster_type = _classify_cluster_type(
            occupancy=occupancy,
            c1_plausible_fraction=c1_plausible_fraction,
            c4_plausible_fraction=c4_plausible_fraction,
            thresholds=thresholds,
        )

        summary = {
            "condition_id": condition_id,
            "cluster_id": cluster_id,
            "protein_id": first_metadata.get("protein_id", ""),
            "ligand_id": first_metadata.get("ligand_id", ""),
            "cluster_type": cluster_type,
            "n_poses": n_poses,
            "occupancy": occupancy,
            "medoid_pose_id": medoid_by_cluster.get((condition_id, cluster_id), ""),
            "member_pose_ids": members,
            "c1_plausible_fraction": c1_plausible_fraction,
            "c4_plausible_fraction": c4_plausible_fraction,
        }
        for field_name in _GEOMETRY_STAT_FIELDS:
            median, iqr = _median_iqr([row.get(field_name) for row in geometry_rows])
            summary[f"{field_name}_median"] = median
            summary[f"{field_name}_iqr"] = iqr

        cluster_summaries.append(summary)
        common_row = {
            "condition_id": condition_id,
            "cluster_id": cluster_id,
            "protein_id": summary["protein_id"],
            "ligand_id": summary["ligand_id"],
            "cluster_type": cluster_type,
            "n_poses": n_poses,
            "occupancy": occupancy,
            "medoid_pose_id": summary["medoid_pose_id"],
        }
        ifp_signature_rows.extend(
            _ifp_frequency_rows(
                cluster_summary=common_row,
                members=members,
                ifp_by_pose_id=ifp_by_pose_id,
            )
        )
        residue_signature_rows.extend(
            _residue_frequency_rows(
                cluster_summary=common_row,
                members=members,
                residue_contact_rows=residue_contact_rows,
            )
        )

    return ClusterSignatureTables(
        cluster_summaries=cluster_summaries,
        ifp_signature_rows=ifp_signature_rows,
        residue_signature_rows=residue_signature_rows,
    )


def _write_tsv(output_path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def write_cluster_ifp_signature_table(rows: list[dict[str, Any]], output_path: Path) -> None:
    """Write cluster_ifp_signature.tsv."""

    _write_tsv(output_path, CLUSTER_IFP_SIGNATURE_COLUMNS, rows)


def write_cluster_residue_signature_table(rows: list[dict[str, Any]], output_path: Path) -> None:
    """Write cluster_residue_signature.tsv."""

    _write_tsv(output_path, CLUSTER_RESIDUE_SIGNATURE_COLUMNS, rows)


def write_cluster_table_tsv(rows: list[dict[str, Any]], output_path: Path) -> None:
    """Write cluster_table.tsv."""

    _write_tsv(output_path, CLUSTER_TABLE_COLUMNS, rows)


def write_cluster_signature_summary_json(rows: list[dict[str, Any]], output_path: Path) -> None:
    """Write the per-cluster Stage 7 summary JSON."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps({"clusters": rows}, indent=2))


def write_cluster_signatures(
    signatures: list[ClusterSignature],
    output_path: Path,
) -> None:
    """Write cluster signatures to JSON."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    data = [
        {
            "cluster_id": s.cluster_id,
            "n_poses": s.n_poses,
            "occupancy": s.occupancy,
            "medoid_pose_id": s.medoid_pose_id,
            "geometry": {
                "mean_cu_c1": s.mean_cu_c1,
                "std_cu_c1": s.std_cu_c1,
                "mean_cu_c4": s.mean_cu_c4,
                "std_cu_c4": s.std_cu_c4,
                "mean_his_brace_angle": s.mean_his_brace_angle,
                "std_his_brace_angle": s.std_his_brace_angle,
            },
            "dominant_interactions": s.dominant_interactions,
            "member_pose_ids": s.member_pose_ids,
        }
        for s in signatures
    ]
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    logger.info("Wrote %d cluster signatures to %s", len(signatures), output_path)
