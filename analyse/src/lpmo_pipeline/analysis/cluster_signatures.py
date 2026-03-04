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
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)


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
            geometry_metrics[pid].get("min_cu_c1", float("inf"))
            for pid in members
            if pid in geometry_metrics
        ]
        cu_c4_vals = [
            geometry_metrics[pid].get("min_cu_c4", float("inf"))
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
            arr = np.array([v for v in cu_c1_vals if v < float("inf")])
            if len(arr) > 0:
                sig.mean_cu_c1 = float(np.mean(arr))
                sig.std_cu_c1 = float(np.std(arr))

        if cu_c4_vals:
            arr = np.array([v for v in cu_c4_vals if v < float("inf")])
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
