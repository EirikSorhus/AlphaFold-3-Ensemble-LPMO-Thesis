from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from lpmo_pipeline.analysis.cluster_signatures import (
    build_cluster_signature_tables,
    write_cluster_table_tsv,
    write_cluster_ifp_signature_table,
    write_cluster_residue_signature_table,
    write_cluster_signature_summary_json,
)
from lpmo_pipeline.analysis.prolif_ifp import IFPResult


def _thresholds_path(tmp_path: Path) -> Path:
    path = tmp_path / "thresholds.yaml"
    path.write_text(
        "\n".join(
            [
                "cluster_type:",
                "  c1_compatible_min_plausible_fraction: 0.50",
                "  c4_compatible_min_plausible_fraction: 0.50",
                "  mixed_compatible_min_plausible_fraction: 0.40",
                "  non_plausible_max_plausible_fraction: 0.10",
                "  uncertain_max_occupancy_for_classification: 0.05",
            ]
        )
        + "\n"
    )
    return path


def test_build_cluster_signature_tables_aggregates_cluster_features_and_geometry(
    tmp_path: Path,
) -> None:
    tables = build_cluster_signature_tables(
        cluster_assignment_rows=[
            {
                "pose_id": "pose_a",
                "condition_id": "P1__domain_only__chitin_DP4",
                "cluster_id": 0,
                "cluster_member_flag": True,
            },
            {
                "pose_id": "pose_b",
                "condition_id": "P1__domain_only__chitin_DP4",
                "cluster_id": 0,
                "cluster_member_flag": True,
            },
            {
                "pose_id": "pose_noise",
                "condition_id": "P1__domain_only__chitin_DP4",
                "cluster_id": -1,
                "cluster_member_flag": False,
            },
        ],
        medoid_rows=[
            {
                "condition_id": "P1__domain_only__chitin_DP4",
                "cluster_id": 0,
                "medoid_pose_id": "pose_a",
            }
        ],
        geometry_rows_by_pose_id={
            "pose_a": {
                "Cu_C1_distance": 2.0,
                "Cu_C4_distance": 6.0,
                "geometry_status_C1": "geometry_plausible",
                "geometry_status_C4": "geometry_computable_implausible",
            },
            "pose_b": {
                "Cu_C1_distance": 4.0,
                "Cu_C4_distance": 7.0,
                "geometry_status_C1": "geometry_plausible",
                "geometry_status_C4": "geometry_computable_implausible",
            },
            "pose_noise": {
                "Cu_C1_distance": 1.0,
                "geometry_status_C1": "geometry_highly_plausible",
            },
        },
        ifp_results=[
            IFPResult(
                pose_id="pose_a",
                status="ok",
                feature_names=[
                    "NAG1.B|ASN10.A|HBDonor",
                    "NAG1.B|TYR20.A|Hydrophobic",
                ],
                flat_bitvector=[1, 0],
            ),
            IFPResult(
                pose_id="pose_b",
                status="ok",
                feature_names=[
                    "NAG1.B|ASN10.A|HBDonor",
                    "NAG1.B|TYR20.A|Hydrophobic",
                ],
                flat_bitvector=[1, 1],
            ),
            IFPResult(
                pose_id="pose_noise",
                status="ok",
                feature_names=["NAG1.B|HIS1.A|HBDonor"],
                flat_bitvector=[1],
            ),
        ],
        residue_contact_rows=[
            {
                "pose_id": "pose_a",
                "condition_id": "P1__domain_only__chitin_DP4",
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASN",
                "interaction_type": "HBDonor",
                "ligand_residue_label": "NAG1.B",
                "contact_present": 1,
                "is_catalytic_surface_region": True,
                "is_cbm_region": False,
                "is_linker_region": False,
            },
            {
                "pose_id": "pose_b",
                "condition_id": "P1__domain_only__chitin_DP4",
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASN",
                "interaction_type": "HBDonor",
                "ligand_residue_label": "NAG1.B",
                "contact_present": 1,
                "is_catalytic_surface_region": True,
                "is_cbm_region": False,
                "is_linker_region": False,
            },
            {
                "pose_id": "pose_noise",
                "condition_id": "P1__domain_only__chitin_DP4",
                "residue_chain": "A",
                "residue_number": 1,
                "residue_name": "HIS",
                "interaction_type": "HBDonor",
                "ligand_residue_label": "NAG1.B",
                "contact_present": 1,
                "is_catalytic_surface_region": False,
                "is_cbm_region": False,
                "is_linker_region": False,
            },
        ],
        pose_metadata_by_id={
            "pose_a": {
                "protein_id": "P1",
                "ligand_id": "NAG4",
                "construct_type": "domain_only",
            },
            "pose_b": {
                "protein_id": "P1",
                "ligand_id": "NAG4",
                "construct_type": "domain_only",
            },
            "pose_noise": {
                "protein_id": "P1",
                "ligand_id": "NAG4",
                "construct_type": "domain_only",
            },
        },
        condition_metadata_by_id={
            "P1__domain_only__chitin_DP4": {
                "protein_id": "P1",
                "construct_type": "domain_only",
                "substrate_class": "chitin",
                "dp": 4,
            }
        },
        pose_confidence_rows_by_id={
            "pose_a": {"ranking_score": 0.9, "iptm": 0.8, "mean_plddt": 90.0},
            "pose_b": {"ranking_score": 0.7, "iptm": 0.6, "mean_plddt": 80.0},
        },
        pose_convergence_rows_by_id={
            "pose_a": {
                "ligand_rmsd_to_reference": 0.0,
                "convergent_flag": True,
            },
            "pose_b": {
                "ligand_rmsd_to_reference": 1.0,
                "convergent_flag": True,
            },
        },
        thresholds_path=_thresholds_path(tmp_path),
    )

    assert len(tables.cluster_summaries) == 1
    summary = tables.cluster_summaries[0]
    assert summary["cluster_id"] == 0
    assert summary["n_poses"] == 2
    assert summary["occupancy"] == pytest.approx(2 / 3)
    assert summary["medoid_pose_id"] == "pose_a"
    assert summary["cluster_type"] == "C1_compatible"
    assert summary["Cu_C1_distance_median"] == pytest.approx(3.0)
    assert summary["Cu_C1_distance_iqr"] == pytest.approx(1.0)
    assert summary["construct_type"] == "domain_only"
    assert summary["substrate_class"] == "chitin"
    assert summary["dp"] == 4
    assert summary["cluster_size"] == 2
    assert summary["c1_geometry_computable_fraction"] == pytest.approx(1.0)
    assert summary["c4_geometry_computable_fraction"] == pytest.approx(1.0)
    assert summary["c1_geometry_plausible_fraction"] == pytest.approx(1.0)
    assert summary["c4_geometry_plausible_fraction"] == pytest.approx(0.0)
    assert summary["mean_ranking_score"] == pytest.approx(0.8)
    assert summary["median_iptm"] == pytest.approx(0.7)
    assert summary["medoid_mean_plddt"] == 90.0
    assert summary["median_ligand_rmsd_to_reference"] == pytest.approx(0.5)
    assert summary["convergent_fraction"] == pytest.approx(1.0)
    assert summary["medoid_ligand_rmsd_to_reference"] == 0.0

    hbond_rows = [
        row
        for row in tables.ifp_signature_rows
        if row["feature_name"] == "NAG1.B|ASN10.A|HBDonor"
    ]
    assert len(hbond_rows) == 1
    assert hbond_rows[0]["ligand_residue_label"] == "NAG1.B"
    assert hbond_rows[0]["protein_residue_label"] == "ASN10.A"
    assert hbond_rows[0]["interaction_type"] == "HBDonor"
    assert hbond_rows[0]["contact_frequency"] == pytest.approx(1.0)

    assert all(row["cluster_id"] != -1 for row in tables.ifp_signature_rows)
    assert all(row["cluster_id"] != -1 for row in tables.residue_signature_rows)
    assert tables.residue_signature_rows[0]["contact_frequency"] == pytest.approx(1.0)


def test_cluster_signature_writers_emit_contract_files(tmp_path: Path) -> None:
    cluster_table_path = tmp_path / "cluster_table.tsv"
    ifp_path = tmp_path / "cluster_ifp_signature.tsv"
    residue_path = tmp_path / "cluster_residue_signature.tsv"
    json_path = tmp_path / "cluster_signatures.json"

    write_cluster_table_tsv(
        [
            {
                "condition_id": "condition",
                "cluster_id": 0,
                "protein_id": "P1",
                "ligand_id": "NAG4",
                "cluster_type": "C1_compatible",
                "n_poses": 2,
                "occupancy": 1.0,
                "medoid_pose_id": "pose_a",
                "c1_plausible_fraction": 1.0,
                "c4_plausible_fraction": 0.0,
                "Cu_C1_distance_median": 3.0,
                "Cu_C1_distance_iqr": 1.0,
            }
        ],
        cluster_table_path,
    )
    write_cluster_ifp_signature_table(
        [
            {
                "condition_id": "condition",
                "cluster_id": 0,
                "protein_id": "P1",
                "ligand_id": "NAG4",
                "cluster_type": "C1_compatible",
                "n_poses": 2,
                "occupancy": 1.0,
                "medoid_pose_id": "pose_a",
                "feature_name": "NAG1.B|ASN10.A|HBDonor",
                "ligand_residue_label": "NAG1.B",
                "protein_residue_label": "ASN10.A",
                "interaction_type": "HBDonor",
                "n_poses_with_contact": 2,
                "contact_frequency": 1.0,
            }
        ],
        ifp_path,
    )
    write_cluster_residue_signature_table([], residue_path)
    write_cluster_signature_summary_json([{"cluster_id": 0}], json_path)

    with open(cluster_table_path, newline="") as handle:
        cluster_rows = list(csv.DictReader(handle, delimiter="\t"))
    with open(ifp_path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert cluster_rows[0]["cluster_type"] == "C1_compatible"
    assert cluster_rows[0]["Cu_C1_distance_median"] == "3.0"
    assert rows[0]["feature_name"] == "NAG1.B|ASN10.A|HBDonor"
    assert residue_path.read_text().splitlines()[0].startswith("condition_id\tcluster_id")
    assert json.loads(json_path.read_text()) == {"clusters": [{"cluster_id": 0}]}
