"""LPMO Pipeline: AF3-only clustering with HDBSCAN.

Responsibility: cluster one protein-ligand condition at a time from a binary
ProLIF IFP matrix, then build the raw Stage 6 clustering outputs.

Stage 6 contract:
  Input: binary IFP matrix for one protein-ligand condition
  Output: cluster assignments, medoids, and condition summary

Implementation notes:
  - The governing AF3-only plan removed the old within-model/cross-model split.
  - Clustering uses precomputed Jaccard distance on binary IFP vectors.
  - Medoids are chosen by minimum summed Jaccard distance within each cluster.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from pathlib import Path
from typing import Any

import numpy as np

from lpmo_pipeline.config import load_runtime_paths_config
from lpmo_pipeline.analysis.prolif_ifp import ContactEligibility

try:
    import hdbscan
except ImportError:
    hdbscan = None

from lpmo_pipeline.utils.logging import StructuredLogger


_RUNTIME_PATHS = load_runtime_paths_config()
_THRESHOLDS_PATH = _RUNTIME_PATHS.pipeline_assets.thresholds_config


@dataclass(frozen=True)
class HDBSCANConfig:
    """Condition-wise HDBSCAN settings loaded from thresholds.yaml."""

    min_cluster_size: int
    min_samples: int | None = None
    metric: str = "jaccard"
    cluster_selection_epsilon: float = 0.0
    cluster_selection_method: str = "eom"
    locked: bool = False


@dataclass
class ClusteringResult:
    """Outcome of one condition-wise clustering procedure."""

    n_clusters: int
    n_outliers: int
    outlier_rate: float
    cluster_sizes: dict[int, int]
    cluster_occupancy: dict[int, float]
    outlier_indices: list[int]
    cluster_labels: np.ndarray
    medoids: dict[int, int] = field(default_factory=dict)
    medoid_distance_sums: dict[int, float] = field(default_factory=dict)


@dataclass(frozen=True)
class ConditionClusterSummary:
    """Condition-level summary row for Stage 6 outputs."""

    condition_id: str
    n_qc_pass_poses: int
    n_ifp_success: int
    n_contact_eligible: int
    contact_eligible_fraction: float
    null_ifp_fraction: float
    vdw_only_fraction: float
    low_specific_contact_fraction: float
    median_n_non_vdw_interactions: float
    median_n_non_vdw_contact_residues: float
    n_ifp_clustered: int
    n_noise: int
    noise_fraction: float
    n_clusters: int
    top_cluster_occupancy: float
    cluster_entropy: float
    occupancy_gini: float

    def to_row(self) -> dict[str, Any]:
        return {
            "condition_id": self.condition_id,
            "n_qc_pass_poses": self.n_qc_pass_poses,
            "n_ifp_success": self.n_ifp_success,
            "n_contact_eligible": self.n_contact_eligible,
            "contact_eligible_fraction": self.contact_eligible_fraction,
            "null_ifp_fraction": self.null_ifp_fraction,
            "vdw_only_fraction": self.vdw_only_fraction,
            "low_specific_contact_fraction": self.low_specific_contact_fraction,
            "median_n_non_vdw_interactions": self.median_n_non_vdw_interactions,
            "median_n_non_vdw_contact_residues": self.median_n_non_vdw_contact_residues,
            "n_ifp_clustered": self.n_ifp_clustered,
            "n_noise": self.n_noise,
            "noise_fraction": self.noise_fraction,
            "n_clusters": self.n_clusters,
            "top_cluster_occupancy": self.top_cluster_occupancy,
            "cluster_entropy": self.cluster_entropy,
            "occupancy_gini": self.occupancy_gini,
        }


def load_hdbscan_config(config_path: Path | None = None) -> HDBSCANConfig:
    """Load the locked HDBSCAN settings from thresholds.yaml."""
    import yaml

    path = config_path or _THRESHOLDS_PATH
    raw_config = yaml.safe_load(path.read_text()) or {}
    hdbscan_config = raw_config.get("hdbscan") or {}
    clustering_config = raw_config.get("clustering") or {}

    min_samples_raw = hdbscan_config.get("min_samples", clustering_config.get("min_samples"))
    min_samples = None if min_samples_raw in {None, "", "null"} else int(min_samples_raw)

    return HDBSCANConfig(
        min_cluster_size=int(
            hdbscan_config.get("min_cluster_size", clustering_config.get("min_cluster_size", 3))
        ),
        min_samples=min_samples,
        metric=str(hdbscan_config.get("metric") or clustering_config.get("distance_metric") or "jaccard"),
        cluster_selection_epsilon=float(hdbscan_config.get("cluster_selection_epsilon", 0.0)),
        cluster_selection_method=str(hdbscan_config.get("cluster_selection_method", "eom")),
        locked=bool(hdbscan_config.get("locked", False)),
    )


def build_cluster_assignment_rows(
    condition_id: str,
    pose_ids: list[str],
    labels: np.ndarray,
    *,
    distance_matrix: np.ndarray | None = None,
    medoids: dict[int, int] | None = None,
) -> list[dict[str, Any]]:
    """Build the raw cluster_assignments.tsv rows for one condition."""
    if len(pose_ids) != len(labels):
        raise ValueError("pose_ids and labels must have the same length")

    rows: list[dict[str, Any]] = []
    for index, (pose_id, label) in enumerate(zip(pose_ids, labels, strict=True)):
        cluster_id = int(label)
        distance_to_representative: float | None = None
        if (
            cluster_id != -1
            and distance_matrix is not None
            and medoids is not None
            and cluster_id in medoids
        ):
            distance_to_representative = float(distance_matrix[index, medoids[cluster_id]])

        rows.append(
            {
                "pose_id": pose_id,
                "condition_id": condition_id,
                "cluster_id": cluster_id,
                "cluster_member_flag": cluster_id != -1,
                "noise_flag": cluster_id == -1,
                "distance_to_cluster_representative": distance_to_representative,
            }
        )
    return rows


def build_medoid_rows(
    condition_id: str,
    pose_ids: list[str],
    medoids: dict[int, int],
    medoid_distance_sums: dict[int, float],
    *,
    structure_paths: dict[str, str] | None = None,
) -> list[dict[str, Any]]:
    """Build the raw medoid_manifest.tsv rows for one condition."""
    rows: list[dict[str, Any]] = []
    structure_paths = structure_paths or {}
    for cluster_id, medoid_index in sorted(medoids.items()):
        medoid_pose_id = pose_ids[medoid_index]
        rows.append(
            {
                "condition_id": condition_id,
                "cluster_id": cluster_id,
                "medoid_pose_id": medoid_pose_id,
                "medoid_structure_path": structure_paths.get(medoid_pose_id, ""),
                "medoid_ifp_distance_sum": float(medoid_distance_sums.get(cluster_id, 0.0)),
            }
        )
    return rows


def build_condition_cluster_summary(
    condition_id: str,
    clustering_result: ClusteringResult,
    *,
    n_qc_pass_poses: int,
    n_ifp_success: int | None = None,
    contact_eligibilities: list[ContactEligibility] | None = None,
) -> ConditionClusterSummary:
    """Build the Stage 6 condition summary row.

    The cluster occupancy distribution for entropy and Gini is normalized over
    non-noise cluster members. Noise is reported separately via `n_noise` and
    `noise_fraction`.
    """
    contact_eligibilities = contact_eligibilities or []
    if n_ifp_success is None:
        n_ifp_success = len(contact_eligibilities) if contact_eligibilities else len(clustering_result.cluster_labels)

    if contact_eligibilities:
        n_contact_eligible = sum(1 for eligibility in contact_eligibilities if eligibility.eligible)
        null_ifp_count = sum(1 for eligibility in contact_eligibilities if eligibility.exclusion_class == "null_ifp")
        vdw_only_count = sum(1 for eligibility in contact_eligibilities if eligibility.exclusion_class == "vdw_only")
        low_specific_contact_count = sum(
            1 for eligibility in contact_eligibilities if eligibility.exclusion_class == "low_specific_contact"
        )
        non_vdw_interactions = np.asarray(
            [eligibility.n_non_vdw_interactions for eligibility in contact_eligibilities],
            dtype=float,
        )
        non_vdw_contact_residues = np.asarray(
            [eligibility.n_non_vdw_contact_residues for eligibility in contact_eligibilities],
            dtype=float,
        )
        median_n_non_vdw_interactions = float(np.median(non_vdw_interactions))
        median_n_non_vdw_contact_residues = float(np.median(non_vdw_contact_residues))
    else:
        n_contact_eligible = 0
        null_ifp_count = 0
        vdw_only_count = 0
        low_specific_contact_count = 0
        median_n_non_vdw_interactions = 0.0
        median_n_non_vdw_contact_residues = 0.0

    if n_ifp_success > 0:
        contact_eligible_fraction = n_contact_eligible / n_ifp_success
        null_ifp_fraction = null_ifp_count / n_ifp_success
        vdw_only_fraction = vdw_only_count / n_ifp_success
        low_specific_contact_fraction = low_specific_contact_count / n_ifp_success
    else:
        contact_eligible_fraction = 0.0
        null_ifp_fraction = 0.0
        vdw_only_fraction = 0.0
        low_specific_contact_fraction = 0.0

    non_noise_sizes = list(clustering_result.cluster_sizes.values())
    if non_noise_sizes:
        size_array = np.asarray(non_noise_sizes, dtype=float)
        probabilities = size_array / np.sum(size_array)
        top_cluster_occupancy = float(np.max(probabilities))
        cluster_entropy = float(
            -sum(probability * math.log2(probability) for probability in probabilities if probability > 0)
        )
        occupancy_gini = _gini(size_array)
    else:
        top_cluster_occupancy = 0.0
        cluster_entropy = 0.0
        occupancy_gini = 0.0

    return ConditionClusterSummary(
        condition_id=condition_id,
        n_qc_pass_poses=n_qc_pass_poses,
        n_ifp_success=n_ifp_success,
        n_contact_eligible=n_contact_eligible,
        contact_eligible_fraction=contact_eligible_fraction,
        null_ifp_fraction=null_ifp_fraction,
        vdw_only_fraction=vdw_only_fraction,
        low_specific_contact_fraction=low_specific_contact_fraction,
        median_n_non_vdw_interactions=median_n_non_vdw_interactions,
        median_n_non_vdw_contact_residues=median_n_non_vdw_contact_residues,
        n_ifp_clustered=len(clustering_result.cluster_labels),
        n_noise=clustering_result.n_outliers,
        noise_fraction=clustering_result.outlier_rate,
        n_clusters=clustering_result.n_clusters,
        top_cluster_occupancy=top_cluster_occupancy,
        cluster_entropy=cluster_entropy,
        occupancy_gini=occupancy_gini,
    )


def _gini(values: np.ndarray) -> float:
    """Compute the Gini coefficient for a positive-valued vector."""
    if values.size == 0:
        return 0.0
    sorted_values = np.sort(values.astype(float))
    total = float(np.sum(sorted_values))
    if total == 0.0:
        return 0.0
    n_values = sorted_values.size
    ranks = np.arange(1, n_values + 1, dtype=float)
    return float((2.0 * np.sum(ranks * sorted_values) / (n_values * total)) - ((n_values + 1) / n_values))


class HDBSCANClusterer:
    """Cluster one condition's IFP matrix using binary Jaccard distance."""

    def __init__(
        self,
        min_cluster_size: int | None = None,
        *,
        min_samples: int | None = None,
        metric: str | None = None,
        cluster_selection_epsilon: float | None = None,
        cluster_selection_method: str | None = None,
        output_dir: Path | None = None,
        config: HDBSCANConfig | None = None,
    ):
        resolved = config or load_hdbscan_config()
        self.config = HDBSCANConfig(
            min_cluster_size=min_cluster_size or resolved.min_cluster_size,
            min_samples=min_samples if min_samples is not None else resolved.min_samples,
            metric=metric or resolved.metric,
            cluster_selection_epsilon=(
                cluster_selection_epsilon
                if cluster_selection_epsilon is not None
                else resolved.cluster_selection_epsilon
            ),
            cluster_selection_method=cluster_selection_method or resolved.cluster_selection_method,
            locked=resolved.locked,
        )
        self.output_dir = Path(output_dir) if output_dir else Path(".")
        self.logger = StructuredLogger("clustering", self.output_dir)

        if hdbscan is None:
            raise ImportError("HDBSCAN not installed")
        if self.config.metric != "jaccard":
            raise ValueError(f"Unsupported clustering metric: {self.config.metric}")

    def cluster(self, ifp_matrix: np.ndarray, pose_ids: list[str]) -> ClusteringResult:
        """Cluster one condition's binary IFP matrix."""
        matrix = np.asarray(ifp_matrix, dtype=np.uint8)
        if matrix.ndim != 2:
            raise ValueError("ifp_matrix must be a 2D array")
        if matrix.shape[0] != len(pose_ids):
            raise ValueError("pose_ids must match the number of IFP rows")

        self.logger.log_step_start(
            "hdbscan_clustering",
            {
                "n_poses": matrix.shape[0],
                "ifp_dim": matrix.shape[1] if matrix.ndim == 2 else 0,
                "min_cluster_size": self.config.min_cluster_size,
                "min_samples": self.config.min_samples,
                "metric": self.config.metric,
            },
        )

        if matrix.shape[0] == 0:
            result = ClusteringResult(
                n_clusters=0,
                n_outliers=0,
                outlier_rate=0.0,
                cluster_sizes={},
                cluster_occupancy={},
                outlier_indices=[],
                cluster_labels=np.array([], dtype=int),
            )
            self.logger.log_step_end("hdbscan_clustering", "success", {"n_clusters": 0}, 0.0)
            return result

        distance_matrix = self.compute_jaccard_distances(matrix)
        labels = self._fit_labels(distance_matrix)
        medoids, medoid_distance_sums = self.select_medoids(distance_matrix, labels)
        result = self._aggregate_clusters(
            labels,
            medoids=medoids,
            medoid_distance_sums=medoid_distance_sums,
        )

        self.logger.log_step_end(
            "hdbscan_clustering",
            "success",
            {
                "n_clusters": result.n_clusters,
                "n_outliers": result.n_outliers,
                "outlier_rate": result.outlier_rate,
            },
            0.0,
        )
        return result

    def compute_jaccard_distances(self, ifp_matrix: np.ndarray) -> np.ndarray:
        """Compute pairwise Jaccard distances for binary IFP vectors."""
        matrix = np.asarray(ifp_matrix, dtype=np.uint8)
        n_rows = matrix.shape[0]
        distances = np.zeros((n_rows, n_rows), dtype=float)

        for row_index in range(n_rows):
            for col_index in range(row_index + 1, n_rows):
                row = matrix[row_index]
                col = matrix[col_index]
                intersection = int(np.dot(row, col))
                union = int(np.sum((row + col) > 0))
                jaccard_distance = 0.0 if union == 0 else 1.0 - (intersection / union)
                distances[row_index, col_index] = jaccard_distance
                distances[col_index, row_index] = jaccard_distance

        return distances

    def select_medoids(
        self,
        distance_matrix: np.ndarray,
        labels: np.ndarray,
    ) -> tuple[dict[int, int], dict[int, float]]:
        """Select exact medoids by minimum summed within-cluster distance."""
        medoids: dict[int, int] = {}
        medoid_distance_sums: dict[int, float] = {}

        for label in sorted({int(value) for value in labels.tolist()}):
            if label == -1:
                continue

            cluster_indices = np.where(labels == label)[0]
            cluster_distances = distance_matrix[np.ix_(cluster_indices, cluster_indices)]
            distance_sums = np.sum(cluster_distances, axis=1)
            medoid_position = int(np.argmin(distance_sums))
            medoid_index = int(cluster_indices[medoid_position])

            medoids[label] = medoid_index
            medoid_distance_sums[label] = float(distance_sums[medoid_position])

        return medoids, medoid_distance_sums

    def _fit_labels(self, distance_matrix: np.ndarray) -> np.ndarray:
        n_rows = distance_matrix.shape[0]

        if n_rows < self.config.min_cluster_size:
            return np.full(n_rows, -1, dtype=int)

        if np.all(distance_matrix == 0.0):
            return np.zeros(n_rows, dtype=int)

        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=self.config.min_cluster_size,
            min_samples=self.config.min_samples,
            metric="precomputed",
            cluster_selection_epsilon=self.config.cluster_selection_epsilon,
            cluster_selection_method=self.config.cluster_selection_method,
        )
        return clusterer.fit_predict(distance_matrix)

    def _aggregate_clusters(
        self,
        labels: np.ndarray,
        *,
        medoids: dict[int, int],
        medoid_distance_sums: dict[int, float],
    ) -> ClusteringResult:
        unique_labels = {int(label) for label in labels.tolist()}
        n_outliers = int(np.sum(labels == -1))
        total_rows = len(labels)
        cluster_sizes: dict[int, int] = {}
        cluster_occupancy: dict[int, float] = {}

        for label in sorted(unique_labels):
            if label == -1:
                continue
            cluster_size = int(np.sum(labels == label))
            cluster_sizes[label] = cluster_size
            cluster_occupancy[label] = cluster_size / total_rows if total_rows else 0.0

        return ClusteringResult(
            n_clusters=len(cluster_sizes),
            n_outliers=n_outliers,
            outlier_rate=(n_outliers / total_rows) if total_rows else 0.0,
            cluster_sizes=cluster_sizes,
            cluster_occupancy=cluster_occupancy,
            outlier_indices=np.where(labels == -1)[0].tolist(),
            cluster_labels=labels,
            medoids=medoids,
            medoid_distance_sums=medoid_distance_sums,
        )


class CrossRunClusterer:
    """Legacy compatibility shim for the retired cross-run clustering design."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        del args, kwargs

    def aggregate_medoid_ifps(self, *args: Any, **kwargs: Any) -> np.ndarray:
        del args, kwargs
        raise NotImplementedError(
            "Cross-run clustering is retired in the AF3-only pipeline; cluster one condition at a time."
        )

    def cluster_cross_run(self, *args: Any, **kwargs: Any) -> ClusteringResult:
        del args, kwargs
        raise NotImplementedError(
            "Cross-run clustering is retired in the AF3-only pipeline; cluster one condition at a time."
        )


HDBANSCANClusterer = HDBSCANClusterer
