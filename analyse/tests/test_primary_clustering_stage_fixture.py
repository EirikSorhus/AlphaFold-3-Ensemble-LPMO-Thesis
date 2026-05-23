from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).parent / "run_tests_scripts" / "build_primary_clustering_stage_fixture.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("build_primary_clustering_stage_fixture", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_ifp_result_from_pose_ifp_row_parses_json_payloads() -> None:
    module = _load_module()

    result = module.ifp_result_from_pose_ifp_row(
        {
            "pose_id": "P1_NAG4_seed-1_sample-0_model",
            "ifp_generation_status": "ok",
            "ifp_vector": "[1, 0]",
            "ifp_feature_names": '["NAG1.B|ASN10.A|ImplicitHBDonor", "NAG1.B|ASN10.A|VdWContact"]',
            "n_total_contacts": "1",
            "ifp_interaction_counts": '{"ImplicitHBDonor": 1, "VdWContact": 0}',
            "ifp_error": "",
        }
    )

    assert result.pose_id == "P1_NAG4_seed-1_sample-0_model"
    assert result.status == "ok"
    assert result.flat_bitvector == [1, 0]
    assert result.feature_names == ["NAG1.B|ASN10.A|ImplicitHBDonor", "NAG1.B|ASN10.A|VdWContact"]
    assert result.interaction_counts == {"ImplicitHBDonor": 1, "VdWContact": 0}


def test_build_fixture_from_synthetic_condition_is_deterministic(tmp_path: Path) -> None:
    module = _load_module()
    condition_id = "P1__domain_only__chitin_DP4"
    output_root = tmp_path / "clustering_stage_outputs"
    production_root = tmp_path / "pilot" / "domain_only_shards" / "shard_0001" / "production_output"
    matrix_path = (
        production_root
        / "clustering_pilot"
        / "pilot"
        / "conditions"
        / condition_id
        / "main_contact_eligible_ifp_matrix.csv"
    )
    matrix_path.parent.mkdir(parents=True, exist_ok=True)
    matrix_path.write_text(
        "\n".join(
            [
                "pose_id,NAG1.B|ASN10.A|ImplicitHBDonor",
                "P1_NAG4_seed-1_sample-0_model,1",
                "P1_NAG4_seed-1_sample-1_model,1",
                "P1_NAG4_seed-1_sample-2_model,1",
                "P1_NAG4_seed-2_sample-0_model,1",
                "P1_NAG4_seed-2_sample-1_model,1",
            ]
        )
        + "\n"
    )

    pose_ids = [
        "P1_NAG4_seed-1_sample-0_model",
        "P1_NAG4_seed-1_sample-1_model",
        "P1_NAG4_seed-1_sample-2_model",
        "P1_NAG4_seed-2_sample-0_model",
        "P1_NAG4_seed-2_sample-1_model",
    ]
    _write_tsv(
        production_root / "pose_ifp_table.tsv",
        [
            "pose_id",
            "ifp_generation_status",
            "ifp_vector",
            "ifp_feature_names",
            "n_total_contacts",
            "n_implicit_hbond_acceptor",
            "n_implicit_hbond_donor",
            "n_vdw_contact",
            "ifp_interaction_counts",
            "ifp_interaction_occurrence_counts",
            "ifp_error",
        ],
        [
            {
                "pose_id": pose_id,
                "ifp_generation_status": "ok",
                "ifp_vector": "[1]",
                "ifp_feature_names": '["NAG1.B|ASN10.A|ImplicitHBDonor"]',
                "n_total_contacts": 1,
                "n_implicit_hbond_acceptor": 0,
                "n_implicit_hbond_donor": 1,
                "n_vdw_contact": 0,
                "ifp_interaction_counts": '{"ImplicitHBDonor": 1}',
                "ifp_interaction_occurrence_counts": '{"ImplicitHBDonor": 1}',
                "ifp_error": "",
            }
            for pose_id in pose_ids
        ],
    )
    _write_tsv(
        production_root / "pose_residue_contact_table.tsv",
        [
            "pose_id",
            "protein_id",
            "condition_id",
            "residue_chain",
            "residue_number",
            "residue_name",
            "interaction_type",
            "contact_present",
            "ligand_residue_label",
            "distance_if_available",
            "is_catalytic_surface_region",
            "is_cbm_region",
            "is_linker_region",
        ],
        [
            {
                "pose_id": pose_id,
                "protein_id": "P1",
                "condition_id": condition_id,
                "residue_chain": "A",
                "residue_number": 10,
                "residue_name": "ASN",
                "interaction_type": "ImplicitHBDonor",
                "contact_present": 1,
                "ligand_residue_label": "NAG1.B",
                "distance_if_available": "",
                "is_catalytic_surface_region": False,
                "is_cbm_region": False,
                "is_linker_region": False,
            }
            for pose_id in pose_ids
        ],
    )
    _write_tsv(
        production_root / "pose_geometry.tsv",
        [
            "pose_id",
            "model",
            "protein_id",
            "ligand_id",
            "geometry_status_C1",
            "geometry_status_C4",
            "Cu_C1_distance",
            "Cu_C4_distance",
        ],
        [
            {
                "pose_id": pose_id,
                "model": "af3",
                "protein_id": "P1",
                "ligand_id": "NAG4",
                "geometry_status_C1": "geometry_plausible",
                "geometry_status_C4": "geometry_computable_implausible",
                "Cu_C1_distance": 4.0,
                "Cu_C4_distance": 8.0,
            }
            for pose_id in pose_ids
        ],
    )
    _write_tsv(
        production_root / "condition_cluster_summary.tsv",
        [
            "condition_id",
            "n_qc_pass_poses",
            "n_ifp_success",
            "n_contact_eligible",
            "contact_eligible_fraction",
            "null_ifp_fraction",
            "vdw_only_fraction",
            "low_specific_contact_fraction",
            "median_n_non_vdw_interactions",
            "median_n_non_vdw_contact_residues",
            "n_ifp_clustered",
            "n_noise",
            "noise_fraction",
            "n_clusters",
            "top_cluster_occupancy",
            "cluster_entropy",
            "occupancy_gini",
        ],
        [
            {
                "condition_id": condition_id,
                "n_qc_pass_poses": 3,
                "n_ifp_success": 3,
                "n_contact_eligible": 3,
                "contact_eligible_fraction": 1.0,
                "null_ifp_fraction": 0.0,
                "vdw_only_fraction": 0.0,
                "low_specific_contact_fraction": 0.0,
                "median_n_non_vdw_interactions": 1.0,
                "median_n_non_vdw_contact_residues": 1.0,
                "n_ifp_clustered": 3,
                "n_noise": 0,
                "noise_fraction": 0.0,
                "n_clusters": 1,
                "top_cluster_occupancy": 1.0,
                "cluster_entropy": 0.0,
                "occupancy_gini": 0.0,
            }
        ],
    )
    (production_root / "analysis_core_summary.json").write_text(
        json.dumps(
            {
                "cases": [
                    {
                        "pose_id": pose_id,
                        "condition_id": condition_id,
                        "normalized_cif": str(tmp_path / f"{pose_id}.cif"),
                    }
                    for pose_id in pose_ids
                ]
            }
        )
    )

    provenance_dir = output_root / "conditions" / condition_id / module.SELECTED_PARAMETER_LABEL
    provenance_dir.mkdir(parents=True)
    (provenance_dir / "cluster_assignments.tsv").write_text("placeholder\n")
    (provenance_dir / "medoid_manifest.tsv").write_text("placeholder\n")
    selected_grid = output_root / "selected_primary_clustering_grid_results.tsv"
    _write_tsv(
        selected_grid,
        [
            "condition_id",
            "construct_type",
            "protein_id",
            "target",
            "method",
            "parameter_label",
            "distance_threshold",
            "min_cluster_size",
            "original_formal_clustering_allowed",
            "original_n_qc_pass_poses",
            "matrix_path",
            "cluster_assignments_tsv",
            "medoid_manifest_tsv",
        ],
        [
            {
                "condition_id": condition_id,
                "construct_type": "domain_only",
                "protein_id": "P1",
                "target": "chitin_DP4",
                "method": module.SELECTED_METHOD,
                "parameter_label": module.SELECTED_PARAMETER_LABEL,
                "distance_threshold": module.SELECTED_DISTANCE_THRESHOLD,
                "min_cluster_size": module.SELECTED_MIN_CLUSTER_SIZE,
                "original_formal_clustering_allowed": "True",
                "original_n_qc_pass_poses": 5,
                "matrix_path": str(matrix_path),
                "cluster_assignments_tsv": str(provenance_dir / "cluster_assignments.tsv"),
                "medoid_manifest_tsv": str(provenance_dir / "medoid_manifest.tsv"),
            }
        ],
    )

    module.build_fixture(
        selected_grid_results=selected_grid,
        output_root=output_root,
        allow_condition_count_mismatch=True,
    )
    first_checksums = {
        name: _sha256(output_root / name)
        for name in [
            "cluster_assignments.tsv",
            "medoid_manifest.tsv",
            "condition_cluster_summary.tsv",
            "cluster_signatures.json",
        ]
    }
    module.build_fixture(
        selected_grid_results=output_root / "selected_primary_clustering_grid_results.tsv",
        output_root=output_root,
        allow_condition_count_mismatch=True,
    )
    second_checksums = {name: _sha256(output_root / name) for name in first_checksums}

    assert first_checksums == second_checksums
    assert all((output_root / name).exists() for name in module.ROOT_OUTPUT_FILES)
    assert json.loads((output_root / "cluster_signatures.json").read_text())["clusters"]
    summary = json.loads((output_root / "fixture_build_summary.json").read_text())
    assert summary["n_conditions_total"] == 1
    assert summary["n_clusters"] == 1
