from __future__ import annotations

import csv
from pathlib import Path

import pytest

from lpmo_pipeline.analysis.activity_mapping import (
    build_predictive_cluster_table_rows,
    write_predictive_cluster_table,
)


def test_build_predictive_cluster_table_uses_current_stage7_fields(tmp_path: Path) -> None:
    rows = build_predictive_cluster_table_rows(
        [
            {
                "condition_id": "P1__domain_only__chitin_DP4",
                "cluster_id": "0",
                "protein_id": "P1",
                "construct_type": "domain_only",
                "ligand_id": "NAG4",
                "substrate_class": "chitin",
                "dp": "4",
                "cluster_type": "C1_compatible",
                "n_poses": "3",
                "cluster_size": "3",
                "occupancy": "0.75",
                "medoid_pose_id": "pose_a",
                "c1_geometry_computable_fraction": "1.0",
                "c4_geometry_computable_fraction": "0.5",
                "c1_geometry_plausible_fraction": "0.67",
                "c4_geometry_plausible_fraction": "0.0",
                "Cu_C1_distance_median": "3.1",
                "Cu_C4_distance_median": "6.2",
                "oxyl_H_C1_distance_median": "2.2",
                "oxyl_H_C4_distance_median": "4.8",
                "Cu_oxyl_H_C1_angle_median": "142.0",
                "Cu_oxyl_H_C4_angle_median": "91.0",
                "ring_normal_vs_brace_normal_median": "36.0",
                "convergent_fraction": "0.8",
                "mean_ranking_score": "0.9",
                "median_iptm": "0.7",
                "medoid_mean_plddt": "91.0",
            }
        ],
        protein_metadata_rows=[
            {"protein_id": "P1", "family": "AA9", "cbm_status": "present"}
        ],
        activity_annotation_rows=[
            {
                "protein_id": "P1",
                "mapped_regio_class": "C1",
                "mapped_substrate_class": "chitin",
                "mapped_activity_label": "chitin_C1_hydroxylating",
                "mapping_rule": "exact_1.14.99.53",
            }
        ],
    )

    assert len(rows) == 1
    row = rows[0]
    assert row["analysis_id"] == "P1__domain_only__chitin_DP4__cluster0"
    assert row["family"] == "AA9"
    assert row["cbm_status"] == "present"
    assert row["DP"] == "4"
    assert row["median_Cu_C1"] == "3.1"
    assert row["median_oxyl_H_C4"] == "4.8"
    assert row["convergence_support"] == pytest.approx(0.8)
    assert row["support_factor"] == pytest.approx(0.8)
    assert row["cluster_weight"] == pytest.approx(0.6)
    assert row["invalid_geometry_fraction"] == pytest.approx(0.0)
    assert row["plausible_geometry_fraction"] == pytest.approx(0.67)
    assert row["experimental_regio_label"] == "C1"
    assert row["experimental_substrate_label"] == "chitin"

    output_path = tmp_path / "predictive_cluster_table.tsv"
    write_predictive_cluster_table(rows, output_path)
    with output_path.open(newline="") as handle:
        written = list(csv.DictReader(handle, delimiter="\t"))
    assert written[0]["analysis_id"] == "P1__domain_only__chitin_DP4__cluster0"
    assert written[0]["experimental_activity_label"] == "chitin_C1_hydroxylating"
