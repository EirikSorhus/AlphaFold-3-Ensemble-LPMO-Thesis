from __future__ import annotations

import csv
import importlib.util
from pathlib import Path


SCRIPT_PATH = Path(__file__).parent / "run_tests_scripts" / "run_analysis_core_real_cifs.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("run_analysis_core_real_cifs", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys()) if rows else ["condition_id"]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _write_matrix(path: Path, pose_ids: list[str], features: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["pose_id", *features])
        for pose_id in pose_ids:
            writer.writerow([pose_id, *([1] * len(features))])


def test_validate_stage6_outputs_checks_matrix_rows_and_status(tmp_path: Path) -> None:
    module = _load_module()
    production_output = tmp_path / "production_output"
    condition_id = "P1__domain_only__chitin_DP4"
    condition_slug = "P1__domain_only__chitin_DP4"

    condition_row = {
        "condition_id": condition_id,
        "n_qc_pass_poses": 3,
        "n_ifp_success": 3,
        "n_contact_eligible": 3,
        "contact_eligible_fraction": 1.0,
        "minimum_clusterable_n": 10,
        "clustering_status": "insufficient_clusterable_signal",
        "formal_clustering_allowed": False,
        "n_ifp_clustered": 0,
        "n_noise": 0,
        "noise_fraction": 0.0,
        "n_clusters": 0,
    }
    _write_tsv(production_output / "condition_cluster_summary.tsv", [condition_row])
    _write_tsv(production_output / "condition_table.tsv", [condition_row])
    _write_tsv(production_output / "cluster_assignments.tsv", [])
    _write_tsv(production_output / "medoid_manifest.tsv", [])
    _write_matrix(
        production_output / "ifp_matrices" / condition_slug / "ifp_matrix.csv",
        ["pose_1", "pose_2", "pose_3"],
        ["NAG1.B|ASN10.A|HBDonor"],
    )
    _write_matrix(
        production_output / "ifp_matrices" / condition_slug / "main_clustering_ifp_matrix.csv",
        ["pose_1", "pose_2", "pose_3"],
        ["NAG1.B|ASN10.A|HBDonor"],
    )

    validation = module._validate_stage6_outputs(production_output)

    assert validation["passed"] is True
    assert validation["status_counts"] == {"insufficient_clusterable_signal": 1}
    assert validation["conditions"][0]["main_matrix"]["n_rows"] == 3
    assert validation["conditions"][0]["formal_clustering_allowed"] is False


def test_validate_stage6_outputs_fails_on_status_mismatch(tmp_path: Path) -> None:
    module = _load_module()
    production_output = tmp_path / "production_output"
    condition_id = "P1__domain_only__chitin_DP4"

    condition_row = {
        "condition_id": condition_id,
        "n_contact_eligible": 3,
        "minimum_clusterable_n": 10,
        "clustering_status": "ok",
        "formal_clustering_allowed": True,
        "n_ifp_clustered": 3,
        "n_clusters": 1,
    }
    _write_tsv(production_output / "condition_cluster_summary.tsv", [condition_row])
    _write_tsv(production_output / "condition_table.tsv", [condition_row])
    _write_tsv(production_output / "cluster_assignments.tsv", [])
    _write_tsv(production_output / "medoid_manifest.tsv", [])
    _write_matrix(
        production_output / "ifp_matrices" / condition_id / "ifp_matrix.csv",
        ["pose_1", "pose_2", "pose_3"],
        ["NAG1.B|ASN10.A|HBDonor"],
    )
    _write_matrix(
        production_output / "ifp_matrices" / condition_id / "main_clustering_ifp_matrix.csv",
        ["pose_1", "pose_2", "pose_3"],
        ["NAG1.B|ASN10.A|HBDonor"],
    )

    validation = module._validate_stage6_outputs(production_output)

    assert validation["passed"] is False
    assert any("clustering_status" in error for error in validation["errors"])
