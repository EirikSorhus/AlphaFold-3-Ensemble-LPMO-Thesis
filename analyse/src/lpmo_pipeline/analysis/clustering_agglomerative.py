"""Agglomerative Jaccard clustering for pilot method comparison.

This module keeps the implementation additive: it returns the same
``ClusteringResult`` contract as the production HDBSCAN path so pilot outputs
can compare methods without changing the standard Stage 6 surface.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from lpmo_pipeline.analysis.clustering_hdbscan import ClusteringResult
from lpmo_pipeline.utils.logging import StructuredLogger

try:
    from sklearn.cluster import AgglomerativeClustering
except ImportError:
    AgglomerativeClustering = None


@dataclass(frozen=True)
class AgglomerativeJaccardConfig:
    """Settings for pilot agglomerative clustering on binary IFP matrices."""

    linkage: str = "average"
    distance_threshold: float = 0.5
    min_cluster_size: int = 10


class AgglomerativeJaccardClusterer:
    """Cluster one condition's IFP matrix using Jaccard distance thresholding."""

    def __init__(
        self,
        *,
        linkage: str | None = None,
        distance_threshold: float | None = None,
        min_cluster_size: int | None = None,
        output_dir: Path | None = None,
        config: AgglomerativeJaccardConfig | None = None,
    ) -> None:
        resolved = config or AgglomerativeJaccardConfig()
        self.config = AgglomerativeJaccardConfig(
            linkage=linkage or resolved.linkage,
            distance_threshold=(
                distance_threshold if distance_threshold is not None else resolved.distance_threshold
            ),
            min_cluster_size=min_cluster_size or resolved.min_cluster_size,
        )
        self.output_dir = Path(output_dir) if output_dir else Path(".")
        self.logger = StructuredLogger("agglomerative_clustering", self.output_dir)

        if AgglomerativeClustering is None:
            raise ImportError("scikit-learn not installed")
        if self.config.linkage not in {"average", "complete", "single"}:
            raise ValueError(f"Unsupported agglomerative linkage: {self.config.linkage}")
        if not 0.0 <= self.config.distance_threshold <= 1.0:
            raise ValueError("distance_threshold must be between 0.0 and 1.0")
        if self.config.min_cluster_size < 1:
            raise ValueError("min_cluster_size must be at least 1")

    def cluster(self, ifp_matrix: np.ndarray, pose_ids: list[str]) -> ClusteringResult:
        """Cluster one condition's binary IFP matrix."""
        matrix = np.asarray(ifp_matrix, dtype=np.uint8)
        if matrix.ndim != 2:
            raise ValueError("ifp_matrix must be a 2D array")
        if matrix.shape[0] != len(pose_ids):
            raise ValueError("pose_ids must match the number of IFP rows")

        self.logger.log_step_start(
            "agglomerative_clustering",
            {
                "n_poses": matrix.shape[0],
                "ifp_dim": matrix.shape[1] if matrix.ndim == 2 else 0,
                "linkage": self.config.linkage,
                "distance_threshold": self.config.distance_threshold,
                "min_cluster_size": self.config.min_cluster_size,
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
                medoids={},
                medoid_distance_sums={},
            )
            self.logger.log_step_end("agglomerative_clustering", "success", {"n_clusters": 0}, 0.0)
            return result

        distance_matrix = self.compute_jaccard_distances(matrix)
        labels = self._fit_labels(distance_matrix)
        labels = self._relabel_small_clusters_as_noise(labels)
        labels = self._renumber_labels(labels)
        medoids, medoid_distance_sums = self.select_medoids(distance_matrix, labels)
        result = self._aggregate_clusters(
            labels,
            medoids=medoids,
            medoid_distance_sums=medoid_distance_sums,
        )

        self.logger.log_step_end(
            "agglomerative_clustering",
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
        if n_rows == 1 or np.all(distance_matrix == 0.0):
            return np.zeros(n_rows, dtype=int)

        clusterer = AgglomerativeClustering(
            n_clusters=None,
            metric="precomputed",
            linkage=self.config.linkage,
            distance_threshold=self.config.distance_threshold,
            compute_full_tree=True,
        )
        return np.asarray(clusterer.fit_predict(distance_matrix), dtype=int)

    def _relabel_small_clusters_as_noise(self, labels: np.ndarray) -> np.ndarray:
        relabeled = np.asarray(labels, dtype=int).copy()
        for label in sorted({int(value) for value in relabeled.tolist()}):
            if label == -1:
                continue
            cluster_size = int(np.sum(relabeled == label))
            if cluster_size < self.config.min_cluster_size:
                relabeled[relabeled == label] = -1
        return relabeled

    def _renumber_labels(self, labels: np.ndarray) -> np.ndarray:
        original = np.asarray(labels, dtype=int)
        label_map: dict[int, int] = {}
        next_label = 0
        for label in original.tolist():
            label = int(label)
            if label == -1:
                continue
            if label not in label_map:
                label_map[label] = next_label
                next_label += 1

        renumbered = np.full(original.shape, -1, dtype=int)
        for index, label in enumerate(original.tolist()):
            label = int(label)
            if label != -1:
                renumbered[index] = label_map[label]
        return renumbered

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