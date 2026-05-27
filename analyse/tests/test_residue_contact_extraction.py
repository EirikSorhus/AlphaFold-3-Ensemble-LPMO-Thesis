from __future__ import annotations

import csv
from pathlib import Path

import pytest

from lpmo_pipeline.analysis.prolif_ifp import IFPResult
from lpmo_pipeline.analysis.residue_contact_extraction import (
    build_pose_residue_contact_rows,
    write_pose_residue_contact_table,
)


def test_build_pose_residue_contact_rows_tracks_ifp_features_exactly() -> None:
    result = IFPResult(
        pose_id="pose-1",
        status="ok",
        feature_names=[
            "NAG1.B|ASN10.A|HBDonor",
            "NAG1.B|ASN10.A|HBAcceptor",
            "NAG2.B|TYR25.A|Hydrophobic",
        ],
        flat_bitvector=[1, 0, 1],
    )

    rows = build_pose_residue_contact_rows(
        [result],
        {"pose-1": {"protein_id": "Q7SCE9", "condition_id": "Q7SCE9__domain_only__chitin_DP4"}},
    )

    assert rows == [
        {
            "pose_id": "pose-1",
            "protein_id": "Q7SCE9",
            "condition_id": "Q7SCE9__domain_only__chitin_DP4",
            "residue_chain": "A",
            "residue_number": 10,
            "residue_name": "ASN",
            "interaction_type": "HBDonor",
            "contact_present": 1,
            "ligand_residue_label": "NAG1.B",
            "distance_if_available": "",
            "is_core_region": False,
            "is_catalytic_surface_region": False,
            "is_non_core_region": False,
            "is_cbm_region": False,
            "is_linker_region": False,
        },
        {
            "pose_id": "pose-1",
            "protein_id": "Q7SCE9",
            "condition_id": "Q7SCE9__domain_only__chitin_DP4",
            "residue_chain": "A",
            "residue_number": 10,
            "residue_name": "ASN",
            "interaction_type": "HBAcceptor",
            "contact_present": 0,
            "ligand_residue_label": "NAG1.B",
            "distance_if_available": "",
            "is_core_region": False,
            "is_catalytic_surface_region": False,
            "is_non_core_region": False,
            "is_cbm_region": False,
            "is_linker_region": False,
        },
        {
            "pose_id": "pose-1",
            "protein_id": "Q7SCE9",
            "condition_id": "Q7SCE9__domain_only__chitin_DP4",
            "residue_chain": "A",
            "residue_number": 25,
            "residue_name": "TYR",
            "interaction_type": "Hydrophobic",
            "contact_present": 1,
            "ligand_residue_label": "NAG2.B",
            "distance_if_available": "",
            "is_core_region": False,
            "is_catalytic_surface_region": False,
            "is_non_core_region": False,
            "is_cbm_region": False,
            "is_linker_region": False,
        },
    ]


def test_build_pose_residue_contact_rows_rejects_unexpected_feature_shape() -> None:
    result = IFPResult(
        pose_id="pose-1",
        status="ok",
        feature_names=["bad-feature-name"],
        flat_bitvector=[1],
    )

    with pytest.raises(ValueError, match="Unsupported IFP feature format"):
        build_pose_residue_contact_rows(
            [result],
            {"pose-1": {"protein_id": "Q7SCE9", "condition_id": "cond"}},
        )


def test_write_pose_residue_contact_table_writes_header_and_rows(tmp_path: Path) -> None:
    output_path = tmp_path / "pose_residue_contact_table.tsv"
    rows = [
        {
            "pose_id": "pose-1",
            "protein_id": "Q7SCE9",
            "condition_id": "cond",
            "residue_chain": "A",
            "residue_number": 10,
            "residue_name": "ASN",
            "interaction_type": "HBDonor",
            "contact_present": 1,
            "ligand_residue_label": "NAG1.B",
            "distance_if_available": "",
            "is_core_region": False,
            "is_catalytic_surface_region": False,
            "is_non_core_region": False,
            "is_cbm_region": False,
            "is_linker_region": False,
        }
    ]

    write_pose_residue_contact_table(rows, output_path)

    with output_path.open(newline="") as handle:
        written_rows = list(csv.DictReader(handle, delimiter="\t"))

    assert written_rows == [
        {
            "pose_id": "pose-1",
            "protein_id": "Q7SCE9",
            "condition_id": "cond",
            "residue_chain": "A",
            "residue_number": "10",
            "residue_name": "ASN",
            "interaction_type": "HBDonor",
            "contact_present": "1",
            "ligand_residue_label": "NAG1.B",
            "distance_if_available": "",
            "is_core_region": "False",
            "is_catalytic_surface_region": "False",
            "is_non_core_region": "False",
            "is_cbm_region": "False",
            "is_linker_region": "False",
        }
    ]
