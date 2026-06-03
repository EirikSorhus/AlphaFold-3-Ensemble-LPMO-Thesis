from __future__ import annotations

import pytest

from lpmo_pipeline.analysis.predictive_postprocess import (
    build_c1_c4_modeling_rows,
    build_substrate_modeling_rows,
    run_predictive_postprocess,
)


def test_build_c1_c4_modeling_rows_prefers_domain_only_and_aggregates_dp_rows() -> None:
    condition_rows = [
        {
            "protein_id": "P1",
            "construct_type": "domain_only",
            "substrate_class": "chitin",
            "dp": "4",
            "noise_fraction": "0.20",
            "top_cluster_occupancy": "0.80",
            "median_n_non_vdw_interactions": "6.0",
            "occupancy_weighted_c1_plausible_fraction": "0.8",
            "occupancy_weighted_c4_plausible_fraction": "0.2",
            "occupancy_weighted_Cu_oxyl_H_C1_angle_median": "140.0",
            "occupancy_weighted_Cu_oxyl_H_C4_angle_median": "90.0",
            "occupancy_weighted_Cu_C1_distance_median": "2.8",
            "occupancy_weighted_Cu_C4_distance_median": "5.8",
            "occupancy_weighted_ring_normal_vs_brace_normal_median": "30.0",
            "hbond_contact_fraction": "0.40",
            "aromatic_contact_fraction": "0.30",
            "contact_eligible_fraction": "0.70",
        },
        {
            "protein_id": "P1",
            "construct_type": "domain_only",
            "substrate_class": "chitin",
            "dp": "6",
            "noise_fraction": "0.30",
            "top_cluster_occupancy": "0.70",
            "median_n_non_vdw_interactions": "8.0",
            "occupancy_weighted_c1_plausible_fraction": "0.6",
            "occupancy_weighted_c4_plausible_fraction": "0.1",
            "occupancy_weighted_Cu_oxyl_H_C1_angle_median": "150.0",
            "occupancy_weighted_Cu_oxyl_H_C4_angle_median": "100.0",
            "occupancy_weighted_Cu_C1_distance_median": "3.0",
            "occupancy_weighted_Cu_C4_distance_median": "6.0",
            "occupancy_weighted_ring_normal_vs_brace_normal_median": "40.0",
            "hbond_contact_fraction": "0.50",
            "aromatic_contact_fraction": "0.20",
            "contact_eligible_fraction": "0.60",
        },
        {
            "protein_id": "P1",
            "construct_type": "catalytic_domain",
            "substrate_class": "chitin",
            "dp": "8",
            "noise_fraction": "0.95",
            "top_cluster_occupancy": "0.05",
            "median_n_non_vdw_interactions": "1.0",
            "occupancy_weighted_c1_plausible_fraction": "0.1",
            "occupancy_weighted_c4_plausible_fraction": "0.9",
            "occupancy_weighted_Cu_C1_distance_median": "7.0",
            "occupancy_weighted_Cu_C4_distance_median": "2.0",
            "occupancy_weighted_ring_normal_vs_brace_normal_median": "80.0",
            "hbond_contact_fraction": "0.90",
            "aromatic_contact_fraction": "0.90",
            "contact_eligible_fraction": "0.10",
        },
    ]
    protein_metadata_rows = [
        {"protein_id": "P1", "experimental_regio_label": "C1/C4"}
    ]

    rows = build_c1_c4_modeling_rows(condition_rows, protein_metadata_rows)

    assert len(rows) == 1
    row = rows[0]
    assert row["model_row_id"] == "P1__chitin"
    assert row["construct_type"] == "domain_only"
    assert row["n_source_rows"] == 2
    assert row["has_C1_activity"] == 1
    assert row["has_C4_activity"] == 1
    assert row["occupancy_weighted_any_plausibility"] == pytest.approx(0.7)
    assert row["occupancy_weighted_C1_minus_C4_plausibility"] == pytest.approx(0.55)
    assert row["contact_eligible_fraction"] == pytest.approx(0.65)
    assert row["noise_fraction"] == pytest.approx(0.25)
    assert row["top_cluster_occupancy"] == pytest.approx(0.75)
    assert row["median_n_non_vdw_interactions"] == pytest.approx(7.0)
    assert row["hbond_contact_fraction"] == pytest.approx(0.45)
    assert row["aromatic_contact_fraction"] == pytest.approx(0.25)


def test_build_substrate_modeling_rows_aggregates_by_protein_and_encodes_target() -> None:
    condition_rows = [
        {
            "protein_id": "P1",
            "construct_type": "domain_only",
            "substrate_class": "cellulose",
            "dp": "4",
            "contact_eligible_fraction": "0.80",
            "noise_fraction": "0.10",
            "top_cluster_occupancy": "0.60",
            "hbond_contact_fraction": "0.40",
            "aromatic_contact_fraction": "0.30",
            "median_n_non_vdw_interactions": "6.0",
            "occupancy_weighted_c1_plausible_fraction": "0.60",
            "occupancy_weighted_c4_plausible_fraction": "0.20",
            "occupancy_weighted_Cu_C1_distance_median": "3.0",
            "occupancy_weighted_Cu_C4_distance_median": "5.0",
            "occupancy_weighted_ring_normal_vs_brace_normal_median": "20.0",
        },
        {
            "protein_id": "P1",
            "construct_type": "domain_only",
            "substrate_class": "cellulose",
            "dp": "6",
            "contact_eligible_fraction": "0.60",
            "noise_fraction": "0.20",
            "top_cluster_occupancy": "0.50",
            "hbond_contact_fraction": "0.50",
            "aromatic_contact_fraction": "0.20",
            "median_n_non_vdw_interactions": "8.0",
            "occupancy_weighted_c1_plausible_fraction": "0.40",
            "occupancy_weighted_c4_plausible_fraction": "0.10",
            "occupancy_weighted_Cu_C1_distance_median": "2.0",
            "occupancy_weighted_Cu_C4_distance_median": "6.0",
            "occupancy_weighted_ring_normal_vs_brace_normal_median": "30.0",
        },
        {
            "protein_id": "P2",
            "construct_type": "domain_only",
            "substrate_class": "cellulose",
            "dp": "4",
            "contact_eligible_fraction": "0.30",
            "noise_fraction": "0.40",
            "top_cluster_occupancy": "0.25",
            "hbond_contact_fraction": "0.10",
            "aromatic_contact_fraction": "0.05",
            "median_n_non_vdw_interactions": "2.0",
            "occupancy_weighted_c1_plausible_fraction": "0.10",
            "occupancy_weighted_c4_plausible_fraction": "0.05",
            "occupancy_weighted_Cu_C1_distance_median": "6.0",
            "occupancy_weighted_Cu_C4_distance_median": "4.0",
            "occupancy_weighted_ring_normal_vs_brace_normal_median": "70.0",
        },
    ]
    protein_metadata_rows = [
        {"protein_id": "P1", "experimental_substrate_label": "cellulose+chitin"},
        {"protein_id": "P2", "experimental_substrate_label": "chitin"},
    ]

    rows = build_substrate_modeling_rows(
        condition_rows,
        protein_metadata_rows,
        substrate_class="cellulose",
    )

    assert len(rows) == 2
    p1_row = next(row for row in rows if row["protein_id"] == "P1")
    p2_row = next(row for row in rows if row["protein_id"] == "P2")
    assert p1_row["is_active_on_substrate"] == 1
    assert p2_row["is_active_on_substrate"] == 0
    assert p1_row["occupancy_weighted_any_plausibility"] == pytest.approx(0.5)
    assert p1_row["occupancy_weighted_C1_minus_C4_plausibility"] == pytest.approx(0.35)
    assert p1_row["median_n_non_vdw_interactions"] == pytest.approx(7.0)
    assert p1_row["hbond_contact_fraction"] == pytest.approx(0.45)
    assert p1_row["aromatic_contact_fraction"] == pytest.approx(0.25)


def test_build_modeling_rows_use_actual_metadata_columns_and_ec_mapping() -> None:
    condition_rows = [
        {
            "protein_id": "PX1",
            "construct_type": "domain_only",
            "substrate_class": "cellulose",
            "dp": "4",
            "contact_eligible_fraction": "0.80",
            "noise_fraction": "0.10",
            "top_cluster_occupancy": "0.70",
            "hbond_contact_fraction": "0.30",
            "aromatic_contact_fraction": "0.20",
            "median_n_non_vdw_interactions": "5.0",
            "occupancy_weighted_c1_plausible_fraction": "0.60",
            "occupancy_weighted_c4_plausible_fraction": "0.20",
        },
        {
            "protein_id": "PX1",
            "construct_type": "domain_only",
            "substrate_class": "cellulose",
            "dp": "6",
            "contact_eligible_fraction": "0.60",
            "noise_fraction": "0.30",
            "top_cluster_occupancy": "0.50",
            "hbond_contact_fraction": "0.20",
            "aromatic_contact_fraction": "0.10",
            "median_n_non_vdw_interactions": "4.0",
            "occupancy_weighted_c1_plausible_fraction": "0.40",
            "occupancy_weighted_c4_plausible_fraction": "0.30",
        },
    ]
    protein_metadata_rows = [
        {
            "UniProt_ID": "PX1",
            "CAZy_family": "AA9",
            "EC_Number": "1.14.99.54 1.14.99.56",
        }
    ]

    c1_c4_rows = build_c1_c4_modeling_rows(condition_rows, protein_metadata_rows)
    substrate_rows = build_substrate_modeling_rows(
        condition_rows,
        protein_metadata_rows,
        substrate_class="cellulose",
    )

    assert len(c1_c4_rows) == 1
    assert c1_c4_rows[0]["has_C1_activity"] == 1
    assert c1_c4_rows[0]["has_C4_activity"] == 1
    assert c1_c4_rows[0]["occupancy_weighted_any_plausibility"] == pytest.approx(0.5)
    assert len(substrate_rows) == 1
    assert substrate_rows[0]["is_active_on_substrate"] == 1


def test_build_substrate_rows_use_chitin_mapping_from_actual_metadata_columns() -> None:
    condition_rows = [
        {
            "protein_id": "PX2",
            "construct_type": "domain_only",
            "substrate_class": "chitin",
            "dp": "4",
            "contact_eligible_fraction": "0.70",
            "noise_fraction": "0.20",
            "top_cluster_occupancy": "0.60",
            "hbond_contact_fraction": "0.40",
            "aromatic_contact_fraction": "0.30",
            "median_n_non_vdw_interactions": "7.0",
            "occupancy_weighted_c1_plausible_fraction": "0.70",
            "occupancy_weighted_c4_plausible_fraction": "0.50",
        }
    ]
    protein_metadata_rows = [
        {
            "UniProt_ID": "PX2",
            "CAZy_family": "AA10",
            "EC_Number": "1.14.99.53",
        }
    ]

    rows = build_substrate_modeling_rows(
        condition_rows,
        protein_metadata_rows,
        substrate_class="chitin",
    )

    assert len(rows) == 1
    assert rows[0]["is_active_on_substrate"] == 1


def test_build_substrate_rows_treat_amylose_predictions_as_starch_activity() -> None:
    condition_rows = [
        {
            "protein_id": "PX3",
            "construct_type": "domain_only",
            "substrate_class": "amylose",
            "dp": "4",
            "contact_eligible_fraction": "0.80",
            "noise_fraction": "0.10",
            "top_cluster_occupancy": "0.70",
            "hbond_contact_fraction": "0.30",
            "aromatic_contact_fraction": "0.20",
            "median_n_non_vdw_interactions": "5.0",
            "occupancy_weighted_c1_plausible_fraction": "0.60",
            "occupancy_weighted_c4_plausible_fraction": "0.20",
        },
        {
            "protein_id": "PX4",
            "construct_type": "domain_only",
            "substrate_class": "amylose",
            "dp": "6",
            "contact_eligible_fraction": "0.50",
            "noise_fraction": "0.20",
            "top_cluster_occupancy": "0.40",
            "hbond_contact_fraction": "0.10",
            "aromatic_contact_fraction": "0.10",
            "median_n_non_vdw_interactions": "3.0",
            "occupancy_weighted_c1_plausible_fraction": "0.10",
            "occupancy_weighted_c4_plausible_fraction": "0.20",
        },
    ]
    protein_metadata_rows = [
        {
            "UniProt_ID": "PX3",
            "CAZy_family": "AA13",
            "EC_Number": "1.14.99.55",
        },
        {
            "UniProt_ID": "PX4",
            "CAZy_family": "AA10",
            "EC_Number": "1.14.99.53",
        },
    ]

    rows = build_substrate_modeling_rows(
        condition_rows,
        protein_metadata_rows,
        substrate_class="starch",
    )

    assert len(rows) == 2
    assert {row["prediction_substrate"] for row in rows} == {"starch"}
    assert {row["protein_id"]: row["is_active_on_substrate"] for row in rows} == {
        "PX3": 1,
        "PX4": 0,
    }


def test_build_rows_prefer_mapped_labels_over_experimental_labels() -> None:
    condition_rows = [
        {
            "protein_id": "PX3",
            "construct_type": "domain_only",
            "substrate_class": "cellulose",
            "dp": "4",
            "contact_eligible_fraction": "0.70",
            "noise_fraction": "0.20",
            "top_cluster_occupancy": "0.60",
            "hbond_contact_fraction": "0.40",
            "aromatic_contact_fraction": "0.30",
            "median_n_non_vdw_interactions": "7.0",
            "occupancy_weighted_c1_plausible_fraction": "0.70",
            "occupancy_weighted_c4_plausible_fraction": "0.10",
        }
    ]
    protein_metadata_rows = [
        {
            "UniProt_ID": "PX3",
            "mapped_regio_class": "C4",
            "mapped_substrate_class": "cellulose",
            "mapped_activity_label": "cellulose_C4_dehydrogenating",
            "experimental_regio_label": "C1",
            "experimental_substrate_label": "chitin",
        }
    ]

    c1_c4_rows = build_c1_c4_modeling_rows(condition_rows, protein_metadata_rows)
    substrate_rows = build_substrate_modeling_rows(
        condition_rows,
        protein_metadata_rows,
        substrate_class="cellulose",
    )

    assert len(c1_c4_rows) == 1
    assert c1_c4_rows[0]["has_C1_activity"] == 0
    assert c1_c4_rows[0]["has_C4_activity"] == 1
    assert len(substrate_rows) == 1
    assert substrate_rows[0]["is_active_on_substrate"] == 1


def test_build_rows_merge_partial_preannotations_with_ec_fallback() -> None:
    condition_rows = [
        {
            "protein_id": "PX4",
            "construct_type": "domain_only",
            "substrate_class": "cellulose",
            "dp": "4",
            "contact_eligible_fraction": "0.70",
            "noise_fraction": "0.20",
            "top_cluster_occupancy": "0.60",
            "hbond_contact_fraction": "0.40",
            "aromatic_contact_fraction": "0.30",
            "median_n_non_vdw_interactions": "7.0",
            "occupancy_weighted_c1_plausible_fraction": "0.70",
            "occupancy_weighted_c4_plausible_fraction": "0.10",
        }
    ]
    protein_metadata_rows = [
        {
            "UniProt_ID": "PX4",
            "mapped_regio_class": "C1",
            "EC_Number": "1.14.99.54",
        }
    ]

    substrate_rows = build_substrate_modeling_rows(
        condition_rows,
        protein_metadata_rows,
        substrate_class="cellulose",
    )

    assert len(substrate_rows) == 1
    assert substrate_rows[0]["is_active_on_substrate"] == 1


def test_predictive_postprocess_summary_records_locked_backend_and_validation(tmp_path) -> None:
    condition_table = tmp_path / "condition_table.tsv"
    metadata_table = tmp_path / "metadata.tsv"
    output_dir = tmp_path / "output"
    condition_table.write_text(
        "\t".join(
            [
                "condition_id",
                "protein_id",
                "construct_type",
                "substrate_class",
                "dp",
                "contact_eligible_fraction",
                "noise_fraction",
                "top_cluster_occupancy",
                "hbond_contact_fraction",
                "aromatic_contact_fraction",
                "median_n_non_vdw_interactions",
                "occupancy_weighted_c1_plausible_fraction",
                "occupancy_weighted_c4_plausible_fraction",
            ]
        )
        + "\n"
        + "P1__domain__cellulose_DP4\tP1\tdomain_only\tcellulose\t4\t0.8\t0.1\t0.7\t0.3\t0.2\t5\t0.7\t0.1\n"
        + "P2__domain__cellulose_DP4\tP2\tdomain_only\tcellulose\t4\t0.4\t0.2\t0.5\t0.2\t0.1\t3\t0.2\t0.6\n"
    )
    metadata_table.write_text(
        "UniProt_ID\tCAZy_family\tEC_Number\n"
        "P1\tAA9\t1.14.99.54\n"
        "P2\tAA9\t1.14.99.56\n"
    )

    result = run_predictive_postprocess(
        condition_table_path=condition_table,
        protein_metadata_path=metadata_table,
        output_dir=output_dir,
        task="all",
        n_folds=2,
        random_state=2,
    )

    import json

    summary = json.loads(result.summary_path.read_text())
    assert summary["implementation"]["primary_backend"] == "sklearn_logistic_regression_l2"
    assert summary["implementation"]["status"] == "final_compact_exploratory_backend"
    assert summary["validation"]["passed"] is True
    assert summary["validation"]["n_condition_rows"] == 2
    assert summary["validation"]["modeling_tables"]["c1_c4"]["n_rows"] == 2
    assert summary["validation"]["models"]["c1_activity"]["status"] in {
        "validated_cv_outputs",
        "single_class_no_cv",
        "no_evaluable_cv_folds",
    }
