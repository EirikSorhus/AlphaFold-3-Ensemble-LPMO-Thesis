from __future__ import annotations

import csv
from pathlib import Path

import pytest

from lpmo_pipeline.analysis.condition_summary import (
    build_condition_table_rows,
    build_protein_summary_rows,
    write_condition_table,
    write_protein_summary_table,
)


def test_condition_table_uses_qc_attrition_as_master_and_keeps_null_cluster_conditions(
    tmp_path: Path,
) -> None:
    condition_rows = build_condition_table_rows(
        qc_attrition_rows=[
            {
                "condition_id": "P1__domain_only__chitin_DP4",
                "protein_id": "P1",
                "construct_type": "domain_only",
                "ligand_id": "NAG4",
                "substrate_class": "chitin",
                "dp": "4",
                "n_generated": "75",
                "n_prepared": "74",
                "n_prep_error": "1",
                "n_stage1_hard_fail": "2",
                "n_stage1_soft_flag": "3",
                "n_stage1_pass": "72",
                "hard_fail_rate": "0.0266666667",
            },
            {
                "condition_id": "P1__domain_only__cellulose_DP6",
                "protein_id": "P1",
                "construct_type": "domain_only",
                "ligand_id": "CEL6",
                "substrate_class": "cellulose",
                "dp": "6",
                "n_generated": "75",
                "n_prepared": "75",
                "n_stage1_hard_fail": "0",
                "n_stage1_pass": "75",
                "hard_fail_rate": "0.0",
            },
        ],
        condition_cluster_summary_rows=[
            {
                "condition_id": "P1__domain_only__chitin_DP4",
                "n_ifp_success": "70",
                "n_contact_eligible": "20",
                "contact_eligible_fraction": "0.2857142857",
                "vdw_only_fraction": "0.2",
                "low_specific_contact_fraction": "0.1",
                "n_clusters": "2",
                "top_cluster_occupancy": "0.4",
                "cluster_entropy": "0.9",
            }
        ],
        condition_convergence_summary_rows=[
            {
                "condition_id": "P1__domain_only__chitin_DP4",
                "reference_pose_id": "pose_a",
                "n_qc_pass_poses": "72",
                "convergence_fraction": "0.75",
                "median_ligand_rmsd": "1.2",
                "iqr_ligand_rmsd": "0.3",
                "low_convergence_flag": "False",
            }
        ],
        condition_patch_summary_rows=[
            {
                "condition_id": "P1__domain_only__chitin_DP4",
                "any_valid_cluster": "True",
                "total_nonnoise_cluster_occupancy": "0.6",
                "aromatic_contact_fraction": "0.25",
                "polar_contact_fraction": "0.5",
            }
        ],
        pose_confidence_rows=[
            {
                "condition_id": "P1__domain_only__chitin_DP4",
                "confidence_json_status": "ok",
                "ranking_score": "0.8",
                "iptm": "0.7",
                "mean_plddt": "80",
            },
            {
                "condition_id": "P1__domain_only__chitin_DP4",
                "confidence_json_status": "ok",
                "ranking_score": "0.6",
                "iptm": "0.5",
                "mean_plddt": "70",
            },
        ],
        cluster_table_rows=[
            {
                "condition_id": "P1__domain_only__chitin_DP4",
                "n_poses": "3",
                "occupancy": "0.25",
                "c1_plausible_fraction": "1.0",
                "Cu_C1_distance_median": "3.0",
            },
            {
                "condition_id": "P1__domain_only__chitin_DP4",
                "n_poses": "6",
                "occupancy": "0.75",
                "c1_plausible_fraction": "0.0",
                "Cu_C1_distance_median": "5.0",
            },
        ],
    )

    assert [row["condition_id"] for row in condition_rows] == [
        "P1__domain_only__cellulose_DP6",
        "P1__domain_only__chitin_DP4",
    ]
    null_cluster_row = condition_rows[0]
    assert null_cluster_row["n_generated"] == "75"
    assert null_cluster_row["n_clusters"] == 0
    assert null_cluster_row["n_cluster_rows"] == 0

    clustered_row = condition_rows[1]
    assert clustered_row["n_clusters"] == "2"
    assert clustered_row["convergence_reference_pose_id"] == "pose_a"
    assert clustered_row["mean_ranking_score"] == pytest.approx(0.7)
    assert clustered_row["median_iptm"] == pytest.approx(0.6)
    assert clustered_row["cluster_total_occupancy"] == pytest.approx(1.0)
    assert clustered_row["mean_cluster_size"] == pytest.approx(4.5)
    assert clustered_row["occupancy_weighted_c1_plausible_fraction"] == pytest.approx(0.25)
    assert clustered_row["occupancy_weighted_Cu_C1_distance_median"] == pytest.approx(4.5)

    condition_table = tmp_path / "condition_table.tsv"
    write_condition_table(condition_rows, condition_table)
    with condition_table.open(newline="") as handle:
        written = list(csv.DictReader(handle, delimiter="\t"))
    assert written[0]["condition_id"] == "P1__domain_only__cellulose_DP6"


def test_protein_summary_is_groupby_over_condition_table(tmp_path: Path) -> None:
    protein_rows = build_protein_summary_rows(
        [
            {
                "protein_id": "P1",
                "construct_type": "domain_only",
                "ligand_id": "NAG4",
                "substrate_class": "chitin",
                "n_generated": "75",
                "n_prepared": "74",
                "n_stage1_hard_fail": "1",
                "n_stage1_pass": "74",
                "n_ifp_success": "70",
                "n_contact_eligible": "10",
                "n_clusters": "2",
                "contact_eligible_fraction": "0.2",
                "top_cluster_occupancy": "0.4",
                "mean_ranking_score": "0.8",
            },
            {
                "protein_id": "P1",
                "construct_type": "full_length",
                "ligand_id": "CEL4",
                "substrate_class": "cellulose",
                "n_generated": "75",
                "n_prepared": "75",
                "n_stage1_hard_fail": "0",
                "n_stage1_pass": "75",
                "n_ifp_success": "75",
                "n_contact_eligible": "0",
                "n_clusters": "0",
                "contact_eligible_fraction": "0.0",
                "top_cluster_occupancy": "0.0",
                "mean_ranking_score": "0.6",
            },
        ]
    )

    assert len(protein_rows) == 1
    row = protein_rows[0]
    assert row["protein_id"] == "P1"
    assert row["n_conditions"] == 2
    assert row["n_construct_types"] == 2
    assert row["construct_types"] == "domain_only,full_length"
    assert row["n_conditions_with_clusters"] == 1
    assert row["hard_fail_rate"] == pytest.approx(1 / 150)
    assert row["mean_contact_eligible_fraction"] == pytest.approx(0.1)
    assert row["mean_confidence_ranking_score"] == pytest.approx(0.7)

    protein_table = tmp_path / "protein_summary_table.tsv"
    write_protein_summary_table(protein_rows, protein_table)
    with protein_table.open(newline="") as handle:
        written = list(csv.DictReader(handle, delimiter="\t"))
    assert written[0]["protein_id"] == "P1"
