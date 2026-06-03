from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from lpmo_pipeline.analysis import clustering_hdbscan as clustering_hdbscan_module
from lpmo_pipeline.analysis.clustering_agglomerative import AgglomerativeJaccardClusterer
from lpmo_pipeline.analysis.clustering_hdbscan import (
    HDBANSCANClusterer,
    HDBSCANClusterer,
    build_cluster_assignment_rows,
    build_condition_cluster_summary,
    load_hdbscan_config,
)


def test_config_from_yaml() -> None:
    config = load_hdbscan_config()

    assert config.min_cluster_size == 5
    assert config.min_samples is None
    assert config.metric == "jaccard"
    assert config.cluster_selection_method == "eom"
    assert config.allow_single_cluster is True
    assert config.locked is True


def test_production_clusterer_uses_loaded_config_without_mutation() -> None:
    config = load_hdbscan_config()
    clusterer = HDBSCANClusterer(config=config, output_dir=Path("."))

    result = clusterer.cluster(np.zeros((2, 4), dtype=np.uint8), ["p0", "p1"])

    assert clusterer.config == config
    assert result.n_clusters == 0
    assert result.n_outliers == 2
    assert HDBANSCANClusterer is HDBSCANClusterer


def test_hdbscan_constructor_receives_allow_single_cluster(monkeypatch: pytest.MonkeyPatch) -> None:
    captured: dict[str, object] = {}

    class FakeHDBSCAN:
        def __init__(self, **kwargs: object) -> None:
            captured.update(kwargs)

        def fit_predict(self, distance_matrix: np.ndarray) -> np.ndarray:
            return np.full(distance_matrix.shape[0], -1, dtype=int)

    monkeypatch.setattr(clustering_hdbscan_module.hdbscan, "HDBSCAN", FakeHDBSCAN)

    clusterer = HDBSCANClusterer(min_cluster_size=2, min_samples=1, allow_single_cluster=True)
    distance_matrix = np.array(
        [
            [0.0, 0.5],
            [0.5, 0.0],
        ],
        dtype=float,
    )

    labels = clusterer._fit_labels(distance_matrix)

    assert labels.tolist() == [-1, -1]
    assert captured["metric"] == "precomputed"
    assert captured["allow_single_cluster"] is True


def test_single_cluster_from_identical_vectors() -> None:
    matrix = np.ones((20, 50), dtype=np.uint8)
    pose_ids = [f"p{i}" for i in range(20)]

    clusterer = HDBSCANClusterer(min_cluster_size=5, min_samples=1)
    result = clusterer.cluster(matrix, pose_ids)

    assert result.n_clusters == 1
    assert result.n_outliers == 0
    assert result.cluster_labels.tolist() == [0] * 20
    assert result.cluster_sizes == {0: 20}
    assert result.medoids == {0: 0}


def test_outliers_labeled_minus_one() -> None:
    matrix = np.array(
        [
            [1, 1, 0, 0],
            [1, 1, 0, 0],
            [1, 1, 0, 0],
            [0, 0, 1, 1],
            [0, 0, 1, 1],
            [0, 0, 1, 1],
            [0, 0, 0, 0],
        ],
        dtype=np.uint8,
    )
    pose_ids = [f"p{i}" for i in range(len(matrix))]

    clusterer = HDBSCANClusterer(min_cluster_size=3, min_samples=1)
    result = clusterer.cluster(matrix, pose_ids)
    distance_matrix = clusterer.compute_jaccard_distances(matrix)
    rows = build_cluster_assignment_rows(
        "P1__domain_only__chitin_DP4",
        pose_ids,
        result.cluster_labels,
        distance_matrix=distance_matrix,
        medoids=result.medoids,
    )
    summary = build_condition_cluster_summary(
        "P1__domain_only__chitin_DP4",
        result,
        n_qc_pass_poses=9,
    )

    assert result.n_clusters == 2
    assert result.n_outliers == 1
    assert sorted(result.cluster_sizes.values()) == [3, 3]
    assert result.cluster_labels[-1] == -1
    assert rows[-1]["noise_flag"] is True
    assert rows[-1]["cluster_member_flag"] is False
    assert rows[-1]["distance_to_cluster_representative"] is None
    assert summary.n_ifp_success == 7
    assert summary.n_contact_eligible == 0
    assert summary.n_ifp_clustered == 7
    assert summary.n_noise == 1
    assert summary.top_cluster_occupancy == pytest.approx(0.5)
    assert summary.cluster_entropy == pytest.approx(1.0)
    assert summary.occupancy_gini == pytest.approx(0.0)


def test_too_few_poses_returns_all_noise() -> None:
    matrix = np.array(
        [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
        ],
        dtype=np.uint8,
    )

    clusterer = HDBSCANClusterer(min_cluster_size=3, min_samples=1)
    result = clusterer.cluster(matrix, ["p0", "p1"])

    assert result.n_clusters == 0
    assert result.n_outliers == 2
    assert result.cluster_labels.tolist() == [-1, -1]
    assert result.medoids == {}


def test_select_medoids_uses_minimum_summed_jaccard_distance() -> None:
    matrix = np.array(
        [
            [1, 1, 0, 0],
            [1, 1, 0, 1],
            [1, 1, 1, 1],
        ],
        dtype=np.uint8,
    )
    labels = np.array([0, 0, 0], dtype=int)
    clusterer = HDBSCANClusterer(min_cluster_size=2, min_samples=1)

    medoids, medoid_distance_sums = clusterer.select_medoids(
        clusterer.compute_jaccard_distances(matrix),
        labels,
    )

    assert medoids == {0: 1}
    assert medoid_distance_sums[0] == pytest.approx((1.0 / 3.0) + (1.0 / 4.0))


def test_agglomerative_finds_two_separable_clusters() -> None:
    matrix = np.array(
        [
            [1, 1, 0, 0],
            [1, 1, 0, 0],
            [1, 1, 0, 0],
            [0, 0, 1, 1],
            [0, 0, 1, 1],
            [0, 0, 1, 1],
        ],
        dtype=np.uint8,
    )
    pose_ids = [f"p{i}" for i in range(len(matrix))]

    clusterer = AgglomerativeJaccardClusterer(distance_threshold=0.2, min_cluster_size=2)
    result = clusterer.cluster(matrix, pose_ids)

    assert result.n_clusters == 2
    assert result.n_outliers == 0
    assert sorted(result.cluster_sizes.values()) == [3, 3]
    assert result.cluster_labels[:3].tolist() == [0, 0, 0]
    assert result.cluster_labels[3:].tolist() == [1, 1, 1]
    assert result.medoids == {0: 0, 1: 3}


def test_agglomerative_distance_threshold_changes_cluster_membership() -> None:
    matrix = np.array(
        [
            [1, 1, 0, 0],
            [1, 1, 0, 1],
            [0, 0, 1, 1],
        ],
        dtype=np.uint8,
    )
    pose_ids = ["p0", "p1", "p2"]

    permissive = AgglomerativeJaccardClusterer(distance_threshold=0.4, min_cluster_size=2)
    strict = AgglomerativeJaccardClusterer(distance_threshold=0.2, min_cluster_size=2)

    permissive_result = permissive.cluster(matrix, pose_ids)
    strict_result = strict.cluster(matrix, pose_ids)

    assert permissive_result.cluster_labels.tolist() == [0, 0, -1]
    assert permissive_result.n_clusters == 1
    assert strict_result.cluster_labels.tolist() == [-1, -1, -1]
    assert strict_result.n_clusters == 0


def test_agglomerative_too_few_poses_becomes_noise_after_min_cluster_filter() -> None:
    matrix = np.array(
        [
            [1, 0, 0, 0],
            [1, 0, 0, 0],
        ],
        dtype=np.uint8,
    )

    clusterer = AgglomerativeJaccardClusterer(distance_threshold=0.1, min_cluster_size=3)
    result = clusterer.cluster(matrix, ["p0", "p1"])

    assert result.n_clusters == 0
    assert result.n_outliers == 2
    assert result.cluster_labels.tolist() == [-1, -1]


def test_agglomerative_selects_medoid_by_minimum_summed_distance() -> None:
    matrix = np.array(
        [
            [1, 1, 0, 0],
            [1, 1, 0, 1],
            [1, 1, 1, 1],
        ],
        dtype=np.uint8,
    )

    clusterer = AgglomerativeJaccardClusterer(distance_threshold=0.8, min_cluster_size=2)
    result = clusterer.cluster(matrix, ["p0", "p1", "p2"])

    assert result.n_clusters == 1
    assert result.medoids == {0: 1}
    assert result.medoid_distance_sums[0] == pytest.approx((1.0 / 3.0) + (1.0 / 4.0))
