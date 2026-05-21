from __future__ import annotations

import csv
from pathlib import Path

import pytest

from lpmo_pipeline.analysis.family_residue_enrichment import (
    compute_family_residue_enrichment_outputs,
    write_family_aligned_residue_table,
    write_family_residue_enrichment,
)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_compute_family_residue_enrichment_outputs_joins_alignment_and_stage16b_rows() -> None:
    outputs = compute_family_residue_enrichment_outputs(
        protein_condition_residue_score_rows=[
            {
                "protein_id": "P1",
                "condition_id": "P1__domain_only__chitin_DP4",
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASN",
                "residue_label": "ASN10.A",
                "residue_contact_score": 0.8,
            },
            {
                "protein_id": "P1",
                "condition_id": "P1__domain_only__chitin_DP6",
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASN",
                "residue_label": "ASN10.A",
                "residue_contact_score": 0.4,
            },
            {
                "protein_id": "P2",
                "condition_id": "P2__domain_only__cellulose_DP4",
                "residue_chain": "A",
                "residue_number": 11,
                "residue_name": "ASP",
                "residue_label": "ASP11.A",
                "residue_contact_score": 0.5,
            },
            {
                "protein_id": "P1",
                "condition_id": "P1__domain_only__chitin_DP4",
                "residue_chain": "A",
                "residue_number": 20,
                "residue_name": "TYR",
                "residue_label": "TYR20.A",
                "residue_contact_score": 0.25,
            },
        ],
        protein_residue_regio_delta_rows=[
            {
                "protein_id": "P1",
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASN",
                "residue_label": "ASN10.A",
                "n_conditions_with_valid_clusters": 2,
                "n_conditions_with_contact": 2,
                "mean_c1_weighted_residue_score": 0.4,
                "mean_c4_weighted_residue_score": 0.05,
                "c1_minus_c4_weighted_delta": 0.35,
                "is_catalytic_surface_region": True,
                "is_cbm_region": False,
                "is_linker_region": False,
            },
            {
                "protein_id": "P2",
                "residue_chain": "A",
                "residue_number": 11,
                "residue_name": "ASP",
                "residue_label": "ASP11.A",
                "n_conditions_with_valid_clusters": 1,
                "n_conditions_with_contact": 1,
                "mean_c1_weighted_residue_score": 0.05,
                "mean_c4_weighted_residue_score": 0.2,
                "c1_minus_c4_weighted_delta": -0.15,
                "is_catalytic_surface_region": True,
                "is_cbm_region": False,
                "is_linker_region": False,
            },
            {
                "protein_id": "P1",
                "residue_chain": "A",
                "residue_number": 20,
                "residue_name": "TYR",
                "residue_label": "TYR20.A",
                "n_conditions_with_valid_clusters": 1,
                "n_conditions_with_contact": 1,
                "mean_c1_weighted_residue_score": 0.1,
                "mean_c4_weighted_residue_score": 0.0,
                "c1_minus_c4_weighted_delta": 0.1,
                "is_catalytic_surface_region": False,
                "is_cbm_region": False,
                "is_linker_region": False,
            },
        ],
        residue_alignment_rows=[
            {
                "protein_id": "P1",
                "alignment_column": 5,
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASN",
                "residue_label": "ASN10.A",
            },
            {
                "protein_id": "P2",
                "alignment_column": 5,
                "residue_chain": "A",
                "residue_number": 11,
                "residue_name": "ASP",
                "residue_label": "ASP11.A",
            },
            {
                "protein_id": "P1",
                "alignment_column": 9,
                "residue_chain": "A",
                "residue_number": 20,
                "residue_name": "TYR",
                "residue_label": "TYR20.A",
            },
        ],
        protein_metadata_rows=[
            {"protein_id": "P1", "CAZy_family": "AA10"},
            {"protein_id": "P2", "CAZy_family": "AA10"},
        ],
    )

    assert outputs.family_aligned_residue_table == [
        {
            "family_label": "AA10",
            "family_aggregation_id": "P1",
            "alignment_column": "5",
            "protein_id": "P1",
            "construct_type": "",
            "catalytic_core_residue_index": 10,
            "alignment_residue": "",
            "residue_chain": "A",
            "residue_number": 10,
            "residue_name": "ASN",
            "residue_label": "ASN10.A",
            "n_conditions_with_valid_clusters": 2,
            "n_conditions_with_contact": 2,
            "mean_condition_residue_contact_score": pytest.approx(0.6),
            "max_condition_residue_contact_score": pytest.approx(0.8),
            "mean_c1_weighted_residue_score": pytest.approx(0.4),
            "mean_c4_weighted_residue_score": pytest.approx(0.05),
            "c1_minus_c4_weighted_delta": pytest.approx(0.35),
            "is_catalytic_surface_region": True,
            "is_cbm_region": False,
            "is_linker_region": False,
        },
        {
            "family_label": "AA10",
            "family_aggregation_id": "P2",
            "alignment_column": "5",
            "protein_id": "P2",
            "construct_type": "",
            "catalytic_core_residue_index": 11,
            "alignment_residue": "",
            "residue_chain": "A",
            "residue_number": 11,
            "residue_name": "ASP",
            "residue_label": "ASP11.A",
            "n_conditions_with_valid_clusters": 1,
            "n_conditions_with_contact": 1,
            "mean_condition_residue_contact_score": pytest.approx(0.5),
            "max_condition_residue_contact_score": pytest.approx(0.5),
            "mean_c1_weighted_residue_score": pytest.approx(0.05),
            "mean_c4_weighted_residue_score": pytest.approx(0.2),
            "c1_minus_c4_weighted_delta": pytest.approx(-0.15),
            "is_catalytic_surface_region": True,
            "is_cbm_region": False,
            "is_linker_region": False,
        },
        {
            "family_label": "AA10",
            "family_aggregation_id": "P1",
            "alignment_column": "9",
            "protein_id": "P1",
            "construct_type": "",
            "catalytic_core_residue_index": 20,
            "alignment_residue": "",
            "residue_chain": "A",
            "residue_number": 20,
            "residue_name": "TYR",
            "residue_label": "TYR20.A",
            "n_conditions_with_valid_clusters": 1,
            "n_conditions_with_contact": 1,
            "mean_condition_residue_contact_score": pytest.approx(0.25),
            "max_condition_residue_contact_score": pytest.approx(0.25),
            "mean_c1_weighted_residue_score": pytest.approx(0.1),
            "mean_c4_weighted_residue_score": pytest.approx(0.0),
            "c1_minus_c4_weighted_delta": pytest.approx(0.1),
            "is_catalytic_surface_region": False,
            "is_cbm_region": False,
            "is_linker_region": False,
        },
    ]

    assert outputs.family_residue_enrichment == [
        {
            "family_label": "AA10",
            "alignment_column": "5",
            "n_proteins_observed": 2,
            "n_proteins_with_contact": 2,
            "n_family_aggregation_units_observed": 2,
            "n_family_aggregation_units_with_contact": 2,
            "mean_condition_residue_contact_score": pytest.approx(0.55),
            "max_condition_residue_contact_score": pytest.approx(0.8),
            "mean_c1_weighted_residue_score": pytest.approx(0.225),
            "mean_c4_weighted_residue_score": pytest.approx(0.125),
            "mean_c1_minus_c4_weighted_delta": pytest.approx(0.1),
            "n_positive_c1_minus_c4_delta": 1,
            "n_negative_c1_minus_c4_delta": 1,
        },
        {
            "family_label": "AA10",
            "alignment_column": "9",
            "n_proteins_observed": 1,
            "n_proteins_with_contact": 1,
            "n_family_aggregation_units_observed": 1,
            "n_family_aggregation_units_with_contact": 1,
            "mean_condition_residue_contact_score": pytest.approx(0.25),
            "max_condition_residue_contact_score": pytest.approx(0.25),
            "mean_c1_weighted_residue_score": pytest.approx(0.1),
            "mean_c4_weighted_residue_score": pytest.approx(0.0),
            "mean_c1_minus_c4_weighted_delta": pytest.approx(0.1),
            "n_positive_c1_minus_c4_delta": 1,
            "n_negative_c1_minus_c4_delta": 0,
        },
    ]


def test_write_family_residue_outputs(tmp_path: Path) -> None:
    outputs = compute_family_residue_enrichment_outputs(
        protein_condition_residue_score_rows=[],
        protein_residue_regio_delta_rows=[
            {
                "protein_id": "P1",
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASN",
                "residue_label": "ASN10.A",
                "n_conditions_with_valid_clusters": 1,
                "n_conditions_with_contact": 0,
                "mean_c1_weighted_residue_score": 0.0,
                "mean_c4_weighted_residue_score": 0.0,
                "c1_minus_c4_weighted_delta": 0.0,
                "is_catalytic_surface_region": False,
                "is_cbm_region": False,
                "is_linker_region": False,
            }
        ],
        residue_alignment_rows=[
            {
                "protein_id": "P1",
                "alignment_column": 5,
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASN",
                "residue_label": "ASN10.A",
                "family_label": "AA10",
            }
        ],
    )

    aligned_path = tmp_path / "family_aligned_residue_table.tsv"
    enrichment_path = tmp_path / "family_residue_enrichment.tsv"
    write_family_aligned_residue_table(outputs.family_aligned_residue_table, aligned_path)
    write_family_residue_enrichment(outputs.family_residue_enrichment, enrichment_path)

    assert _read_tsv(aligned_path) == [
        {
            "family_label": "AA10",
            "family_aggregation_id": "P1",
            "alignment_column": "5",
            "protein_id": "P1",
            "construct_type": "",
            "catalytic_core_residue_index": "10",
            "alignment_residue": "",
            "residue_chain": "A",
            "residue_number": "10",
            "residue_name": "ASN",
            "residue_label": "ASN10.A",
            "n_conditions_with_valid_clusters": "1",
            "n_conditions_with_contact": "0",
            "mean_condition_residue_contact_score": "",
            "max_condition_residue_contact_score": "",
            "mean_c1_weighted_residue_score": "0.0",
            "mean_c4_weighted_residue_score": "0.0",
            "c1_minus_c4_weighted_delta": "0.0",
            "is_catalytic_surface_region": "False",
            "is_cbm_region": "False",
            "is_linker_region": "False",
        }
    ]
    assert _read_tsv(enrichment_path) == [
        {
            "family_label": "AA10",
            "alignment_column": "5",
            "n_proteins_observed": "1",
            "n_proteins_with_contact": "0",
            "n_family_aggregation_units_observed": "1",
            "n_family_aggregation_units_with_contact": "0",
            "mean_condition_residue_contact_score": "",
            "max_condition_residue_contact_score": "",
            "mean_c1_weighted_residue_score": "0.0",
            "mean_c4_weighted_residue_score": "0.0",
            "mean_c1_minus_c4_weighted_delta": "0.0",
            "n_positive_c1_minus_c4_delta": "0",
            "n_negative_c1_minus_c4_delta": "0",
        }
    ]


def test_compute_family_residue_enrichment_outputs_aggregates_family_enrichment_by_sequence_group() -> None:
    outputs = compute_family_residue_enrichment_outputs(
        protein_condition_residue_score_rows=[
            {
                "protein_id": "P1",
                "condition_id": "P1__domain_only__chitin_DP4",
                "residue_label": "ASN10.A",
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASN",
                "residue_contact_score": 0.8,
            },
            {
                "protein_id": "P2",
                "condition_id": "P2__domain_only__chitin_DP4",
                "residue_label": "ASN10.A",
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASN",
                "residue_contact_score": 0.4,
            },
            {
                "protein_id": "P3",
                "condition_id": "P3__domain_only__chitin_DP4",
                "residue_label": "ASP10.A",
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASP",
                "residue_contact_score": 0.2,
            },
        ],
        protein_residue_regio_delta_rows=[
            {
                "protein_id": "P1",
                "residue_label": "ASN10.A",
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASN",
                "mean_c1_weighted_residue_score": 0.6,
                "mean_c4_weighted_residue_score": 0.1,
                "c1_minus_c4_weighted_delta": 0.5,
            },
            {
                "protein_id": "P2",
                "residue_label": "ASN10.A",
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASN",
                "mean_c1_weighted_residue_score": 0.2,
                "mean_c4_weighted_residue_score": 0.0,
                "c1_minus_c4_weighted_delta": 0.2,
            },
            {
                "protein_id": "P3",
                "residue_label": "ASP10.A",
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASP",
                "mean_c1_weighted_residue_score": 0.1,
                "mean_c4_weighted_residue_score": 0.4,
                "c1_minus_c4_weighted_delta": -0.3,
            },
        ],
        residue_alignment_rows=[
            {
                "protein_id": "P1",
                "family_label": "AA10",
                "family_aggregation_id": "SEQ_A",
                "alignment_column": 10,
                "residue_label": "ASN10.A",
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASN",
            },
            {
                "protein_id": "P2",
                "family_label": "AA10",
                "family_aggregation_id": "SEQ_A",
                "alignment_column": 10,
                "residue_label": "ASN10.A",
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASN",
            },
            {
                "protein_id": "P3",
                "family_label": "AA10",
                "family_aggregation_id": "SEQ_B",
                "alignment_column": 10,
                "residue_label": "ASP10.A",
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASP",
            },
        ],
    )

    assert len(outputs.family_residue_enrichment) == 1
    row = outputs.family_residue_enrichment[0]
    assert row["family_label"] == "AA10"
    assert row["alignment_column"] == "10"
    assert row["n_proteins_observed"] == 3
    assert row["n_proteins_with_contact"] == 3
    assert row["n_family_aggregation_units_observed"] == 2
    assert row["n_family_aggregation_units_with_contact"] == 2
    assert row["mean_condition_residue_contact_score"] == pytest.approx(0.4)
    assert row["max_condition_residue_contact_score"] == pytest.approx(0.8)
    assert row["mean_c1_weighted_residue_score"] == pytest.approx(0.25)
    assert row["mean_c4_weighted_residue_score"] == pytest.approx(0.225)
    assert row["mean_c1_minus_c4_weighted_delta"] == pytest.approx(0.025)
    assert row["n_positive_c1_minus_c4_delta"] == 1
    assert row["n_negative_c1_minus_c4_delta"] == 1