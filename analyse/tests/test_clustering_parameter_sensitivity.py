from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import numpy as np

from lpmo_pipeline.analysis.clustering_hdbscan import ClusteringResult


SCRIPT_PATH = Path(__file__).parent / "run_tests_scripts" / "run_clustering_parameter_sensitivity.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("run_clustering_parameter_sensitivity", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_condition_and_parameter_outputs_match_main_analysis_shape(tmp_path: Path) -> None:
    module = _load_module()

    matrix_path = tmp_path / "input" / "main_contact_eligible_ifp_matrix.csv"
    matrix_path.parent.mkdir(parents=True, exist_ok=True)
    matrix_path.write_text(
        "\n".join(
            [
                "pose_id,feat_a,feat_b",
                "P1_NAG4_seed-1_sample-0_model,1,0",
                "P1_NAG4_seed-1_sample-1_model,1,1",
                "P1_NAG4_seed-2_sample-0_model,0,0",
            ]
        )
        + "\n"
    )

    matrix_input = module.MatrixInput(
        condition_id="P1__domain_only__chitin_DP4",
        construct_type="domain_only",
        protein_id="P1",
        target="chitin_DP4",
        matrix_path=matrix_path,
        pose_ids=[
            "P1_NAG4_seed-1_sample-0_model",
            "P1_NAG4_seed-1_sample-1_model",
            "P1_NAG4_seed-2_sample-0_model",
        ],
        feature_names=["feat_a", "feat_b"],
        matrix=np.asarray([[1, 0], [1, 1], [0, 0]], dtype=np.uint8),
        original_summary={
            "minimum_clusterable_n": "3",
            "formal_clustering_allowed": "true",
            "n_qc_pass_poses": "3",
            "n_ifp_success": "3",
            "n_contact_eligible": "3",
            "contact_eligible_fraction": "1.0",
            "null_ifp_fraction": "0.0",
            "vdw_only_fraction": "0.0",
            "low_specific_contact_fraction": "0.0",
            "median_n_non_vdw_interactions": "2",
            "median_n_non_vdw_contact_residues": "1",
            "raw_feature_count": "2",
        },
    )
    result = ClusteringResult(
        n_clusters=1,
        n_outliers=1,
        outlier_rate=1.0 / 3.0,
        cluster_sizes={0: 2},
        cluster_occupancy={0: 2.0 / 3.0},
        outlier_indices=[2],
        cluster_labels=np.asarray([0, 0, -1], dtype=int),
        medoids={0: 0},
        medoid_distance_sums={0: 0.25},
    )
    distance_matrix = np.asarray(
        [
            [0.0, 0.25, 1.0],
            [0.25, 0.0, 1.0],
            [1.0, 1.0, 0.0],
        ]
    )
    parameter_label = "agglomerative_jaccard__threshold_0p55__min_cluster_size_3"
    bundle = module._write_condition_outputs(
        matrix_input=matrix_input,
        method_output_dir=tmp_path / "condition_output",
        result=result,
        distance_matrix=distance_matrix,
        metadata=module._metadata_for(
            matrix_input,
            method="agglomerative_jaccard",
            parameter_label=parameter_label,
            distance_threshold=0.55,
            min_cluster_size=3,
        ),
    )

    assert Path(bundle.grid_row["cluster_assignments_tsv"]).exists()
    assert Path(bundle.grid_row["medoid_manifest_tsv"]).exists()
    assert Path(bundle.grid_row["condition_cluster_summary_tsv"]).exists()
    assert Path(bundle.grid_row["cluster_table_tsv"]).exists()
    assert Path(bundle.grid_row["condition_table_tsv"]).exists()

    condition_table_rows = _read_tsv(Path(bundle.grid_row["condition_table_tsv"]))
    assert len(condition_table_rows) == 1
    assert condition_table_rows[0]["condition_id"] == "P1__domain_only__chitin_DP4"
    assert condition_table_rows[0]["n_clusters"] == "1"
    assert condition_table_rows[0]["cluster_total_occupancy"] == "0.6666666666666666"

    index_rows = module._write_parameter_run_outputs(
        tmp_path / "aggregated_output",
        {parameter_label: [bundle]},
    )

    assert len(index_rows) == 1
    aggregate_condition_rows = _read_tsv(Path(index_rows[0]["condition_table_tsv"]))
    aggregate_cluster_rows = _read_tsv(Path(index_rows[0]["cluster_table_tsv"]))
    assert len(aggregate_condition_rows) == 1
    assert aggregate_condition_rows[0]["condition_id"] == "P1__domain_only__chitin_DP4"
    assert len(aggregate_cluster_rows) == 1
    assert aggregate_cluster_rows[0]["cluster_id"] == "0"
    assert aggregate_cluster_rows[0]["occupancy"] == "0.6666666666666666"
