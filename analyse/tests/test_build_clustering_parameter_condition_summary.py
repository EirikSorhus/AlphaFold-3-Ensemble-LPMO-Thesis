from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path


SCRIPT_PATH = (
    Path(__file__).parent
    / "run_tests_scripts"
    / "build_clustering_parameter_condition_summary.py"
)


def _load_module():
    spec = importlib.util.spec_from_file_location(
        "build_clustering_parameter_condition_summary",
        SCRIPT_PATH,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_build_condition_summary_rows_computes_profile_stats() -> None:
    module = _load_module()

    grid_rows = [
        {
            "condition_id": "P1__domain_only__chitin_DP4",
            "construct_type": "domain_only",
            "protein_id": "P1",
            "target": "chitin_DP4",
            "n_selected_poses": "21",
            "main_feature_count": "22",
            "original_n_qc_pass_poses": "75",
            "original_n_ifp_success": "75",
            "original_n_contact_eligible": "21",
            "original_clustering_status": "ok",
            "parameter_label": "hdbscan_jaccard__min_cluster_size_5__min_samples_none",
            "n_clusters": "2",
            "n_noise": "5",
            "noise_fraction": "0.2380952381",
            "top_cluster_occupancy": "0.5",
            "cluster_entropy": "1.0",
            "occupancy_gini": "0.1",
            "cluster_sizes_json": '{"0": 10, "1": 6}',
            "top_medoid_pose_id": "pose_a",
            "max_single_seed_fraction": "0.75",
            "seed_artifact_flag": "True",
            "cluster_assignments_tsv": "/tmp/assignments_primary.tsv",
            "medoid_manifest_tsv": "/tmp/medoids_primary.tsv",
        },
        {
            "condition_id": "P1__domain_only__chitin_DP4",
            "construct_type": "domain_only",
            "protein_id": "P1",
            "target": "chitin_DP4",
            "n_selected_poses": "21",
            "main_feature_count": "22",
            "original_n_qc_pass_poses": "75",
            "original_n_ifp_success": "75",
            "original_n_contact_eligible": "21",
            "original_clustering_status": "ok",
            "parameter_label": "hdbscan_jaccard__min_cluster_size_3__min_samples_none",
            "n_clusters": "0",
            "n_noise": "21",
            "noise_fraction": "1.0",
            "top_cluster_occupancy": "0.0",
            "cluster_entropy": "0.0",
            "occupancy_gini": "0.0",
            "cluster_sizes_json": "{}",
            "top_medoid_pose_id": "",
            "max_single_seed_fraction": "0.0",
            "seed_artifact_flag": "False",
            "cluster_assignments_tsv": "/tmp/assignments_hdbscan.tsv",
            "medoid_manifest_tsv": "/tmp/medoids_hdbscan.tsv",
        },
    ]

    fieldnames, summary_rows = module.build_condition_summary_rows(
        grid_rows,
        [
            (
                "primary_hdbscan_min_5",
                "hdbscan_jaccard__min_cluster_size_5__min_samples_none",
            ),
            (
                "sensitivity_hdbscan_min_3",
                "hdbscan_jaccard__min_cluster_size_3__min_samples_none",
            ),
        ],
    )

    assert len(summary_rows) == 1
    row = summary_rows[0]
    assert row["condition_id"] == "P1__domain_only__chitin_DP4"
    assert row["primary_hdbscan_min_5_cluster_size_mean"] == 8.0
    assert row["primary_hdbscan_min_5_cluster_size_sd"] == 2.0
    assert row["primary_hdbscan_min_5_cluster_size_min"] == 6
    assert row["primary_hdbscan_min_5_cluster_size_max"] == 10
    assert row["primary_hdbscan_min_5_cluster_size_median"] == 8.0
    assert row["primary_hdbscan_min_5_seed_artifact_flag"] == "True"
    assert row["sensitivity_hdbscan_min_3_cluster_sizes_json"] == "{}"
    assert row["sensitivity_hdbscan_min_3_cluster_size_mean"] == ""
    assert row["sensitivity_hdbscan_min_3_seed_artifact_flag"] == "False"
    assert "primary_hdbscan_min_5_cluster_assignments_tsv" in fieldnames
