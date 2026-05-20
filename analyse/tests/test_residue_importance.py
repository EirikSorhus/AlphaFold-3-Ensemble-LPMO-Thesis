from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from lpmo_pipeline.analysis.residue_importance import (
    compute_residue_importance_outputs,
    write_condition_patch_summary,
    write_protein_condition_residue_scores,
    write_protein_patch_summary,
    write_protein_residue_regio_delta,
)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_compute_residue_importance_outputs_consumes_stage16_signatures() -> None:
    outputs = compute_residue_importance_outputs(
        [
            {
                "condition_id": "Q7SCE9__domain_only__chitin_DP4",
                "cluster_id": 0,
                "protein_id": "Q7SCE9",
                "ligand_id": "NAG4",
                "cluster_type": "C1_compatible",
                "n_poses": 2,
                "occupancy": 2 / 3,
                "medoid_pose_id": "pose-1",
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASN",
                "interaction_type": "HBDonor",
                "ligand_residue_label": "NAG1.B",
                "n_poses_with_contact": 2,
                "contact_frequency": 1.0,
                "is_catalytic_surface_region": False,
                "is_cbm_region": False,
                "is_linker_region": False,
            },
            {
                "condition_id": "Q7SCE9__domain_only__chitin_DP4",
                "cluster_id": 1,
                "protein_id": "Q7SCE9",
                "ligand_id": "NAG4",
                "cluster_type": "C4_compatible",
                "n_poses": 1,
                "occupancy": 1 / 3,
                "medoid_pose_id": "pose-3",
                "residue_chain": "A",
                "residue_number": 20,
                "residue_name": "TYR",
                "interaction_type": "PiStacking",
                "ligand_residue_label": "NAG2.B",
                "n_poses_with_contact": 1,
                "contact_frequency": 1.0,
                "is_catalytic_surface_region": False,
                "is_cbm_region": False,
                "is_linker_region": False,
            },
        ],
        [
            {
                "condition_id": "Q7SCE9__domain_only__chitin_DP4",
                "cluster_id": 0,
                "protein_id": "Q7SCE9",
                "ligand_id": "NAG4",
                "cluster_type": "C1_compatible",
                "n_poses": 2,
                "occupancy": 2 / 3,
                "medoid_pose_id": "pose-1",
                "member_pose_ids": ["pose-1", "pose-2"],
                "c1_plausible_fraction": 0.5,
                "c4_plausible_fraction": 0.0,
            },
            {
                "condition_id": "Q7SCE9__domain_only__chitin_DP4",
                "cluster_id": 1,
                "protein_id": "Q7SCE9",
                "ligand_id": "NAG4",
                "cluster_type": "C4_compatible",
                "n_poses": 1,
                "occupancy": 1 / 3,
                "medoid_pose_id": "pose-3",
                "member_pose_ids": ["pose-3"],
                "c1_plausible_fraction": 0.0,
                "c4_plausible_fraction": 1.0,
            },
        ],
        cluster_ifp_signature_rows=[
            {
                "condition_id": "Q7SCE9__domain_only__chitin_DP4",
                "cluster_id": 0,
                "protein_id": "Q7SCE9",
                "ligand_id": "NAG4",
                "cluster_type": "C1_compatible",
                "n_poses": 2,
                "occupancy": 2 / 3,
                "medoid_pose_id": "pose-1",
                "feature_name": "NAG1.B|ASN10.A|HBDonor",
                "ligand_residue_label": "NAG1.B",
                "protein_residue_label": "ASN10.A",
                "interaction_type": "HBDonor",
                "n_poses_with_contact": 2,
                "contact_frequency": 1.0,
            },
            {
                "condition_id": "Q7SCE9__domain_only__chitin_DP4",
                "cluster_id": 1,
                "protein_id": "Q7SCE9",
                "ligand_id": "NAG4",
                "cluster_type": "C4_compatible",
                "n_poses": 1,
                "occupancy": 1 / 3,
                "medoid_pose_id": "pose-3",
                "feature_name": "NAG2.B|TYR20.A|PiStacking",
                "ligand_residue_label": "NAG2.B",
                "protein_residue_label": "TYR20.A",
                "interaction_type": "PiStacking",
                "n_poses_with_contact": 1,
                "contact_frequency": 1.0,
            },
        ],
        condition_metadata_by_id={
            "Q7SCE9__domain_only__chitin_DP6": {
                "protein_id": "Q7SCE9",
                "construct_type": "domain_only",
                "substrate_class": "chitin",
                "dp": 6,
            },
        },
    )

    assert [row["residue_label"] for row in outputs.protein_condition_residue_scores] == [
        "ASN10.A",
        "TYR20.A",
    ]

    asn_row = outputs.protein_condition_residue_scores[0]
    tyr_row = outputs.protein_condition_residue_scores[1]
    assert asn_row["residue_contact_score"] == pytest.approx(2 / 3)
    assert asn_row["c1_weighted_residue_score"] == pytest.approx(1 / 3)
    assert asn_row["c4_weighted_residue_score"] == pytest.approx(0.0)
    assert asn_row["c1_minus_c4_weighted_delta"] == pytest.approx(1 / 3)
    assert asn_row["cluster_support_count"] == 1

    assert tyr_row["residue_contact_score"] == pytest.approx(1 / 3)
    assert tyr_row["c1_weighted_residue_score"] == pytest.approx(0.0)
    assert tyr_row["c4_weighted_residue_score"] == pytest.approx(1 / 3)
    assert tyr_row["c1_minus_c4_weighted_delta"] == pytest.approx(-1 / 3)

    assert outputs.protein_residue_regio_delta == [
        {
            "protein_id": "Q7SCE9",
            "residue_chain": "A",
            "residue_number": 10,
            "residue_name": "ASN",
            "residue_label": "ASN10.A",
            "n_conditions_with_valid_clusters": 1,
            "n_conditions_with_contact": 1,
            "mean_c1_weighted_residue_score": pytest.approx(1 / 3),
            "mean_c4_weighted_residue_score": pytest.approx(0.0),
            "c1_minus_c4_weighted_delta": pytest.approx(1 / 3),
            "is_catalytic_surface_region": False,
            "is_cbm_region": False,
            "is_linker_region": False,
        },
        {
            "protein_id": "Q7SCE9",
            "residue_chain": "A",
            "residue_number": 20,
            "residue_name": "TYR",
            "residue_label": "TYR20.A",
            "n_conditions_with_valid_clusters": 1,
            "n_conditions_with_contact": 1,
            "mean_c1_weighted_residue_score": pytest.approx(0.0),
            "mean_c4_weighted_residue_score": pytest.approx(1 / 3),
            "c1_minus_c4_weighted_delta": pytest.approx(-1 / 3),
            "is_catalytic_surface_region": False,
            "is_cbm_region": False,
            "is_linker_region": False,
        },
    ]

    assert outputs.condition_patch_summary == [
        {
            "protein_id": "Q7SCE9",
            "condition_id": "Q7SCE9__domain_only__chitin_DP4",
            "construct_type": "domain_only",
            "substrate_class": "chitin",
            "dp": 4,
            "any_valid_cluster": True,
            "n_clusters_considered": 2,
            "total_nonnoise_cluster_occupancy": 1.0,
            "weighted_contact_feature_mass": 1.0,
            "aromatic_contact_fraction": pytest.approx(1 / 3),
            "polar_contact_fraction": pytest.approx(2 / 3),
            "charged_contact_fraction": pytest.approx(0.0),
            "hydrophobic_contact_fraction": pytest.approx(0.0),
            "hbond_contact_fraction": pytest.approx(2 / 3),
            "catalytic_surface_contact_fraction": pytest.approx(0.0),
            "cbm_contact_fraction": pytest.approx(0.0),
            "linker_contact_fraction": pytest.approx(0.0),
        },
        {
            "protein_id": "Q7SCE9",
            "condition_id": "Q7SCE9__domain_only__chitin_DP6",
            "construct_type": "domain_only",
            "substrate_class": "chitin",
            "dp": 6,
            "any_valid_cluster": False,
            "n_clusters_considered": 0,
            "total_nonnoise_cluster_occupancy": 0,
            "weighted_contact_feature_mass": 0.0,
            "aromatic_contact_fraction": 0.0,
            "polar_contact_fraction": 0.0,
            "charged_contact_fraction": 0.0,
            "hydrophobic_contact_fraction": 0.0,
            "hbond_contact_fraction": 0.0,
            "catalytic_surface_contact_fraction": 0.0,
            "cbm_contact_fraction": 0.0,
            "linker_contact_fraction": 0.0,
        },
    ]

    assert outputs.protein_patch_summary == [
        {
            "protein_id": "Q7SCE9",
            "n_conditions_total": 2,
            "n_conditions_with_valid_clusters": 1,
            "mean_aromatic_contact_fraction": pytest.approx(1 / 6),
            "mean_polar_contact_fraction": pytest.approx(1 / 3),
            "mean_charged_contact_fraction": pytest.approx(0.0),
            "mean_hydrophobic_contact_fraction": pytest.approx(0.0),
            "mean_hbond_contact_fraction": pytest.approx(1 / 3),
            "mean_catalytic_surface_contact_fraction": pytest.approx(0.0),
            "mean_cbm_contact_fraction": pytest.approx(0.0),
            "mean_linker_contact_fraction": pytest.approx(0.0),
        }
    ]


def test_compute_residue_importance_outputs_backfills_zero_rows_without_valid_clusters() -> None:
    outputs = compute_residue_importance_outputs(
        cluster_residue_signature_rows=[],
        cluster_summary_rows=[],
        observed_residue_contact_rows=[
            {
                "pose_id": "pose-1",
                "protein_id": "Q7SCE9",
                "condition_id": "Q7SCE9__domain_only__chitin_DP4",
                "residue_chain": "A",
                "residue_number": 29,
                "residue_name": "GLU",
                "interaction_type": "VdWContact",
                "contact_present": 1,
                "ligand_residue_label": "NAG1.B",
                "distance_if_available": "",
                "is_catalytic_surface_region": True,
                "is_cbm_region": False,
                "is_linker_region": False,
            },
            {
                "pose_id": "pose-1",
                "protein_id": "Q7SCE9",
                "condition_id": "Q7SCE9__domain_only__chitin_DP4",
                "residue_chain": "A",
                "residue_number": 88,
                "residue_name": "ASN",
                "interaction_type": "HBDonor",
                "contact_present": 1,
                "ligand_residue_label": "NAG3.B",
                "distance_if_available": "",
                "is_catalytic_surface_region": False,
                "is_cbm_region": False,
                "is_linker_region": False,
            },
        ],
    )

    assert [row["residue_label"] for row in outputs.protein_condition_residue_scores] == [
        "GLU29.A",
        "ASN88.A",
    ]
    assert all(row["residue_contact_score"] == 0.0 for row in outputs.protein_condition_residue_scores)
    assert all(row["n_conditions_with_valid_clusters"] == 0 for row in outputs.protein_residue_regio_delta)
    assert outputs.condition_patch_summary[0]["any_valid_cluster"] is False
    assert outputs.protein_patch_summary[0]["n_conditions_with_valid_clusters"] == 0


def test_write_residue_importance_outputs_write_headers_and_rows(tmp_path: Path) -> None:
    outputs = compute_residue_importance_outputs(
        cluster_residue_signature_rows=[],
        cluster_summary_rows=[],
        condition_metadata_by_id={
            "Q7SCE9__domain_only__chitin_DP4": {
                "protein_id": "Q7SCE9",
                "construct_type": "domain_only",
                "substrate_class": "chitin",
                "dp": 4,
            }
        },
    )

    protein_condition_path = tmp_path / "protein_condition_residue_scores.tsv"
    protein_regio_path = tmp_path / "protein_residue_regio_delta.tsv"
    condition_patch_path = tmp_path / "condition_patch_summary.tsv"
    protein_patch_path = tmp_path / "protein_patch_summary.tsv"

    write_protein_condition_residue_scores(outputs.protein_condition_residue_scores, protein_condition_path)
    write_protein_residue_regio_delta(outputs.protein_residue_regio_delta, protein_regio_path)
    write_condition_patch_summary(outputs.condition_patch_summary, condition_patch_path)
    write_protein_patch_summary(outputs.protein_patch_summary, protein_patch_path)

    with protein_condition_path.open(newline="") as handle:
        protein_condition_rows = list(csv.DictReader(handle, delimiter="\t"))
    with protein_regio_path.open(newline="") as handle:
        protein_regio_rows = list(csv.DictReader(handle, delimiter="\t"))
    with condition_patch_path.open(newline="") as handle:
        condition_patch_rows = list(csv.DictReader(handle, delimiter="\t"))
    with protein_patch_path.open(newline="") as handle:
        protein_patch_rows = list(csv.DictReader(handle, delimiter="\t"))

    assert protein_condition_rows == []
    assert protein_regio_rows == []
    assert condition_patch_rows == [
        {
            "protein_id": "Q7SCE9",
            "condition_id": "Q7SCE9__domain_only__chitin_DP4",
            "construct_type": "domain_only",
            "substrate_class": "chitin",
            "dp": "4",
            "any_valid_cluster": "False",
            "n_clusters_considered": "0",
            "total_nonnoise_cluster_occupancy": "0",
            "weighted_contact_feature_mass": "0.0",
            "aromatic_contact_fraction": "0.0",
            "polar_contact_fraction": "0.0",
            "charged_contact_fraction": "0.0",
            "hydrophobic_contact_fraction": "0.0",
            "hbond_contact_fraction": "0.0",
            "catalytic_surface_contact_fraction": "0.0",
            "cbm_contact_fraction": "0.0",
            "linker_contact_fraction": "0.0",
        }
    ]
    assert protein_patch_rows == [
        {
            "protein_id": "Q7SCE9",
            "n_conditions_total": "1",
            "n_conditions_with_valid_clusters": "0",
            "mean_aromatic_contact_fraction": "0.0",
            "mean_polar_contact_fraction": "0.0",
            "mean_charged_contact_fraction": "0.0",
            "mean_hydrophobic_contact_fraction": "0.0",
            "mean_hbond_contact_fraction": "0.0",
            "mean_catalytic_surface_contact_fraction": "0.0",
            "mean_cbm_contact_fraction": "0.0",
            "mean_linker_contact_fraction": "0.0",
        }
    ]
