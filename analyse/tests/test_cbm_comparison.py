from __future__ import annotations

import csv
from pathlib import Path

import pytest

from lpmo_pipeline.analysis.cbm_comparison import (
    build_cbm_construct_condition_summary_rows,
    build_cbm_paired_comparison_rows,
    write_cbm_construct_condition_summary,
    write_cbm_paired_comparison_table,
)


def test_cbm_condition_summary_and_pairs_are_condition_level(tmp_path: Path) -> None:
    condition_rows = [
        {
            "condition_id": "P1__domain_only__chitin_DP4",
            "protein_id": "P1",
            "construct_type": "domain_only",
            "substrate_class": "chitin",
            "dp": "4",
            "n_generated": "10",
            "n_stage1_pass": "8",
            "n_ifp_success": "6",
            "n_clusters": "1",
            "noise_fraction": "0.1",
            "top_cluster_occupancy": "0.7",
            "cluster_entropy": "0.2",
            "occupancy_weighted_c1_plausible_fraction": "0.8",
            "occupancy_weighted_c4_plausible_fraction": "0.1",
            "catalytic_surface_contact_fraction": "0.5",
            "aromatic_contact_fraction": "0.25",
            "polar_contact_fraction": "0.3",
        },
        {
            "condition_id": "P1__full_length__chitin_DP4",
            "protein_id": "P1",
            "construct_type": "full_length",
            "substrate_class": "chitin",
            "dp": "4",
            "n_generated": "10",
            "n_stage1_pass": "7",
            "n_ifp_success": "7",
            "n_clusters": "2",
            "noise_fraction": "0.2",
            "top_cluster_occupancy": "0.6",
            "cluster_entropy": "0.4",
            "occupancy_weighted_c1_plausible_fraction": "0.2",
            "occupancy_weighted_c4_plausible_fraction": "0.6",
            "cbm_contact_fraction": "0.45",
            "linker_contact_fraction": "0.15",
            "bridge_fraction": "0.25",
            "catalytic_surface_contact_fraction": "0.4",
            "aromatic_contact_fraction": "0.35",
            "polar_contact_fraction": "0.2",
        },
        {
            "condition_id": "P2__domain_only__chitin_DP4",
            "protein_id": "P2",
            "construct_type": "domain_only",
            "substrate_class": "chitin",
            "dp": "4",
            "n_generated": "4",
            "n_stage1_pass": "4",
            "n_ifp_success": "0",
            "n_clusters": "0",
        },
    ]
    cluster_rows = [
        {
            "condition_id": "P1__domain_only__chitin_DP4",
            "cluster_type": "C1_compatible",
            "occupancy": "0.7",
        },
        {
            "condition_id": "P1__full_length__chitin_DP4",
            "cluster_type": "C4_compatible",
            "occupancy": "0.6",
        },
    ]

    construct_rows = build_cbm_construct_condition_summary_rows(
        condition_rows,
        cluster_rows=cluster_rows,
        protein_metadata_rows=[{"protein_id": "P1", "family": "AA9", "cbm_type": "CBM1"}],
    )
    assert len(construct_rows) == 3
    domain = next(row for row in construct_rows if row["construct_type"] == "domain_only" and row["protein_id"] == "P1")
    full = next(row for row in construct_rows if row["construct_type"] == "full_length")
    assert domain["qc_pass_fraction"] == pytest.approx(0.8)
    assert domain["ifp_success_fraction"] == pytest.approx(0.75)
    assert domain["C1_compatible_fraction"] == pytest.approx(1.0)
    assert full["C4_compatible_fraction"] == pytest.approx(1.0)
    assert full["cbm_recruitment_score"] == pytest.approx(0.45)

    pair_rows = build_cbm_paired_comparison_rows(construct_rows)
    assert len(pair_rows) == 1
    pair = pair_rows[0]
    assert pair["protein_id"] == "P1"
    assert pair["domain_only_condition_id"] == "P1__domain_only__chitin_DP4"
    assert pair["full_length_condition_id"] == "P1__full_length__chitin_DP4"
    assert pair["delta_qc_pass_fraction"] == pytest.approx(-0.1)
    assert pair["delta_C4_minus_C1_geometry_bias"] == pytest.approx(2.0)
    assert pair["bridge_fraction_full_length"] == pytest.approx(0.25)

    construct_path = tmp_path / "cbm_construct_condition_summary.tsv"
    pair_path = tmp_path / "cbm_paired_comparison_table.tsv"
    write_cbm_construct_condition_summary(construct_rows, construct_path)
    write_cbm_paired_comparison_table(pair_rows, pair_path)
    with pair_path.open(newline="") as handle:
        written_pairs = list(csv.DictReader(handle, delimiter="\t"))
    assert written_pairs[0]["protein_id"] == "P1"
    assert construct_path.exists()
