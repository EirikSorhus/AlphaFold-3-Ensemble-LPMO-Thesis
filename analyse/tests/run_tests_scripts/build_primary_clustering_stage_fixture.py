#!/usr/bin/env python3
"""Build a production-like clustering fixture from selected pilot IFP matrices.

The builder reuses existing pilot artifacts and does not rerun discovery, QC,
ProLIF, geometry, or crystal anchoring. It reclusters the selected pilot main
IFP matrices with the locked primary HDBSCAN Jaccard settings, then
regenerates the flat downstream Stage 6/7/16b surfaces used by analysis tests.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from lpmo_pipeline.analysis.cluster_signatures import (
    build_cluster_signature_tables,
    write_cluster_ifp_signature_table,
    write_cluster_residue_signature_table,
    write_cluster_signature_summary_json,
)
from lpmo_pipeline.analysis.clustering_hdbscan import (
    HDBSCANClusterer,
    HDBSCANConfig,
    build_cluster_assignment_rows,
    build_condition_cluster_summary,
    build_medoid_rows,
)
from lpmo_pipeline.analysis.prolif_ifp import IFPResult, write_pose_ifp_table
from lpmo_pipeline.analysis.residue_importance import (
    compute_residue_importance_outputs,
    write_condition_patch_summary,
    write_protein_condition_residue_scores,
    write_protein_patch_summary,
    write_protein_residue_regio_delta,
)


SELECTED_PARAMETER_LABEL = "hdbscan_jaccard__min_cluster_size_5__min_samples_none"
SELECTED_METHOD = "hdbscan_jaccard"
SELECTED_DISTANCE_THRESHOLD = ""
SELECTED_MIN_CLUSTER_SIZE = 5
SELECTED_MIN_SAMPLES = None
EXPECTED_CONDITION_COUNT = 145

ROOT_OUTPUT_FILES = (
    "cluster_assignments.tsv",
    "medoid_manifest.tsv",
    "condition_cluster_summary.tsv",
    "cluster_ifp_signature.tsv",
    "cluster_residue_signature.tsv",
    "cluster_signatures.json",
    "protein_condition_residue_scores.tsv",
    "protein_residue_regio_delta.tsv",
    "condition_patch_summary.tsv",
    "protein_patch_summary.tsv",
    "pose_ifp_table.tsv",
    "pose_residue_contact_table.tsv",
    "pose_geometry.tsv",
    "fixture_build_summary.json",
)

CLUSTER_ASSIGNMENT_COLUMNS = [
    "pose_id",
    "condition_id",
    "cluster_id",
    "cluster_member_flag",
    "noise_flag",
    "distance_to_cluster_representative",
]
MEDOID_COLUMNS = [
    "condition_id",
    "cluster_id",
    "medoid_pose_id",
    "medoid_structure_path",
    "medoid_ifp_distance_sum",
]
CONDITION_CLUSTER_SUMMARY_COLUMNS = [
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
]


@dataclass(frozen=True)
class MatrixInput:
    condition_id: str
    construct_type: str
    protein_id: str
    target: str
    matrix_path: Path
    pose_ids: list[str]
    feature_names: list[str]
    matrix: np.ndarray
    selected_row: dict[str, str]
    source_output_root: Path


@dataclass
class SourceOutputData:
    pose_ifp_results_by_id: dict[str, IFPResult]
    pose_ifp_rows_by_id: dict[str, dict[str, str]]
    pose_residue_contact_rows: list[dict[str, str]]
    pose_geometry_rows_by_id: dict[str, dict[str, str]]
    condition_summary_by_id: dict[str, dict[str, str]]
    pose_condition_by_id: dict[str, str]
    normalized_path_by_pose_id: dict[str, str]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build production-like clustering stage fixture from selected pilot matrices.",
    )
    parser.add_argument(
        "--selected-grid-results",
        default="tests/fixtures/clustering_stage_outputs/selected_primary_clustering_grid_results.tsv",
        help="Selected primary parameter result TSV.",
    )
    parser.add_argument(
        "--output-root",
        default="tests/fixtures/clustering_stage_outputs",
        help="Fixture output root to replace atomically after staging.",
    )
    parser.add_argument(
        "--allow-condition-count-mismatch",
        action="store_true",
        help="Allow selected condition count to differ from 145.",
    )
    return parser.parse_args()


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def _read_matrix(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if not header or header[0] != "pose_id":
            raise ValueError(f"IFP matrix must start with pose_id column: {path}")
        pose_ids: list[str] = []
        rows: list[list[int]] = []
        for row in reader:
            if not row:
                continue
            pose_ids.append(row[0])
            rows.append([int(value) for value in row[1:]])
    matrix = np.asarray(rows, dtype=np.uint8)
    if not rows:
        matrix = np.zeros((0, len(header) - 1), dtype=np.uint8)
    return pose_ids, header[1:], matrix


def _json_list(value: str) -> list[Any]:
    if not value:
        return []
    parsed = json.loads(value)
    if not isinstance(parsed, list):
        raise ValueError(f"Expected JSON list, got: {value[:80]}")
    return parsed


def _json_dict(value: str) -> dict[str, Any]:
    if not value:
        return {}
    parsed = json.loads(value)
    if not isinstance(parsed, dict):
        raise ValueError(f"Expected JSON object, got: {value[:80]}")
    return parsed


def ifp_result_from_pose_ifp_row(row: dict[str, str]) -> IFPResult:
    """Reconstruct the subset of IFPResult needed by downstream signatures."""
    feature_names = [str(value) for value in _json_list(row.get("ifp_feature_names", ""))]
    flat_bitvector = [int(value) for value in _json_list(row.get("ifp_vector", ""))]
    interaction_counts = {
        str(key): int(value)
        for key, value in _json_dict(row.get("ifp_interaction_counts", "")).items()
    }
    return IFPResult(
        pose_id=str(row.get("pose_id", "")),
        status=str(row.get("ifp_generation_status", "")) or "ok",
        error=str(row.get("ifp_error", "")) or None,
        feature_names=feature_names,
        flat_bitvector=flat_bitvector,
        n_total_contacts=int(float(row.get("n_total_contacts") or 0)),
        interaction_counts=interaction_counts,
    )


def _production_output_root(path: Path) -> Path:
    for parent in path.resolve().parents:
        if parent.name == "production_output":
            return parent
    raise ValueError(f"Could not find production_output ancestor for {path}")


def _condition_metadata(condition_id: str, selected_row: dict[str, str]) -> dict[str, Any]:
    target = selected_row.get("target", "")
    substrate_class = target
    dp = ""
    if "_DP" in target:
        substrate_class, dp = target.rsplit("_DP", 1)
    return {
        "protein_id": selected_row.get("protein_id", condition_id.split("__", 1)[0]),
        "construct_type": selected_row.get("construct_type", ""),
        "substrate_class": substrate_class,
        "dp": int(dp) if str(dp).isdigit() else dp,
    }


def _ligand_id_from_pose_id(pose_id: str, protein_id: str, target: str) -> str:
    prefix = f"{protein_id}_"
    if pose_id.startswith(prefix) and "_seed-" in pose_id:
        return pose_id[len(prefix) :].split("_seed-", 1)[0]
    if "_DP" in target:
        substrate, dp = target.rsplit("_DP", 1)
        prefix_by_substrate = {"chitin": "NAG", "cellulose": "CEL", "amylose": "STA"}
        return f"{prefix_by_substrate.get(substrate, substrate.upper())}{dp}"
    return target


def _load_source_output_data(root: Path) -> SourceOutputData:
    pose_ifp_rows = {
        str(row.get("pose_id", "")): row
        for row in _read_tsv(root / "pose_ifp_table.tsv")
        if row.get("pose_id")
    }
    pose_condition_by_id: dict[str, str] = {}
    normalized_path_by_pose_id: dict[str, str] = {}
    summary_path = root / "analysis_core_summary.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text())
        for case in summary.get("cases", []):
            pose_id = str(case.get("pose_id", ""))
            if not pose_id:
                continue
            if case.get("condition_id"):
                pose_condition_by_id[pose_id] = str(case["condition_id"])
            normalized_path = str(case.get("normalized_cif") or case.get("normalized_cif_path") or "")
            if normalized_path:
                normalized_path_by_pose_id[pose_id] = normalized_path

    return SourceOutputData(
        pose_ifp_results_by_id={
            pose_id: ifp_result_from_pose_ifp_row(row)
            for pose_id, row in pose_ifp_rows.items()
        },
        pose_ifp_rows_by_id=pose_ifp_rows,
        pose_residue_contact_rows=_read_tsv(root / "pose_residue_contact_table.tsv"),
        pose_geometry_rows_by_id={
            str(row.get("pose_id", "")): row
            for row in _read_tsv(root / "pose_geometry.tsv")
            if row.get("pose_id")
        },
        condition_summary_by_id={
            str(row.get("condition_id", "")): row
            for row in _read_tsv(root / "condition_cluster_summary.tsv")
            if row.get("condition_id")
        },
        pose_condition_by_id=pose_condition_by_id,
        normalized_path_by_pose_id=normalized_path_by_pose_id,
    )


def _validate_selected_rows(rows: list[dict[str, str]], allow_condition_count_mismatch: bool) -> None:
    if not rows:
        raise ValueError("No selected rows found")
    if not allow_condition_count_mismatch and len(rows) != EXPECTED_CONDITION_COUNT:
        raise ValueError(f"Expected {EXPECTED_CONDITION_COUNT} conditions, found {len(rows)}")
    for row in rows:
        if row.get("method") != SELECTED_METHOD:
            raise ValueError(f"Unexpected method for {row.get('condition_id')}: {row.get('method')}")
        if row.get("parameter_label") != SELECTED_PARAMETER_LABEL:
            raise ValueError(
                f"Unexpected parameter label for {row.get('condition_id')}: {row.get('parameter_label')}"
            )
        for key in ("matrix_path", "cluster_assignments_tsv", "medoid_manifest_tsv"):
            path = Path(row.get(key, ""))
            if not path.exists():
                raise FileNotFoundError(f"Missing selected {key} for {row.get('condition_id')}: {path}")


def _load_matrix_inputs(rows: list[dict[str, str]]) -> list[MatrixInput]:
    matrices: list[MatrixInput] = []
    for row in rows:
        matrix_path = Path(row["matrix_path"])
        pose_ids, feature_names, matrix = _read_matrix(matrix_path)
        matrices.append(
            MatrixInput(
                condition_id=str(row["condition_id"]),
                construct_type=str(row.get("construct_type", "")),
                protein_id=str(row.get("protein_id", "")),
                target=str(row.get("target", "")),
                matrix_path=matrix_path,
                pose_ids=pose_ids,
                feature_names=feature_names,
                matrix=matrix,
                selected_row=row,
                source_output_root=_production_output_root(matrix_path),
            )
        )
    return matrices


def _cluster_field_overrides(summary_row: dict[str, Any]) -> dict[str, Any]:
    return {
        key: summary_row.get(key, "")
        for key in (
            "n_ifp_clustered",
            "n_noise",
            "noise_fraction",
            "n_clusters",
            "top_cluster_occupancy",
            "cluster_entropy",
            "occupancy_gini",
        )
    }


def _copy_selected_provenance(rows: list[dict[str, str]], staging_root: Path) -> None:
    for row in rows:
        condition_id = str(row["condition_id"])
        destination = staging_root / "conditions" / condition_id / SELECTED_PARAMETER_LABEL
        destination.mkdir(parents=True, exist_ok=True)
        for key in ("cluster_assignments_tsv", "medoid_manifest_tsv"):
            source = Path(row[key]).parent
            if source.exists():
                for path in source.iterdir():
                    if path.is_file():
                        shutil.copy2(path, destination / path.name)
                break


def _write_readme(staging_root: Path) -> None:
    (staging_root / "README.md").write_text(
        "\n".join(
            [
                "# Clustering Stage Outputs Fixture",
                "",
                "This directory is a production-like clustering output fixture built from the",
                "selected pilot primary clustering method.",
                "",
                "Primary method and parameters:",
                "",
                "- Method: `hdbscan_jaccard`",
                "- Minimum cluster size: `5`",
                "- Minimum samples: `null`",
                "- Cluster selection method: `eom`",
                "",
                "The root-level TSV/JSON files mirror the main analysis output surface for",
                "clustering and downstream cluster-dependent stages. The `conditions/`",
                "subtree keeps the selected per-condition pilot clustering artifacts as",
                "provenance.",
                "",
                "This fixture is regenerated by:",
                "",
                "```bash",
                "sbatch tests/run_tests_scripts/build_primary_clustering_stage_fixture.slurm",
                "```",
                "",
                "See `fixture_build_summary.json` for source roots, counts, and build metadata.",
                "",
            ]
        )
    )


def _write_selected_rows_for_output(rows: list[dict[str, str]], output_root: Path, staging_root: Path) -> None:
    rewritten: list[dict[str, str]] = []
    for row in rows:
        condition_id = str(row["condition_id"])
        method_dir = output_root / "conditions" / condition_id / SELECTED_PARAMETER_LABEL
        rewritten.append(
            {
                **row,
                "cluster_assignments_tsv": str(method_dir / "cluster_assignments.tsv"),
                "medoid_manifest_tsv": str(method_dir / "medoid_manifest.tsv"),
            }
        )
    _write_tsv(staging_root / "selected_primary_clustering_grid_results.tsv", list(rows[0].keys()), rewritten)


def _write_root_outputs(
    *,
    staging_root: Path,
    output_root: Path,
    selected_grid_results: Path,
    matrices: list[MatrixInput],
    rows: list[dict[str, str]],
    source_by_root: dict[Path, SourceOutputData],
) -> dict[str, Any]:
    clusterer = HDBSCANClusterer(
        output_dir=staging_root,
        config=HDBSCANConfig(
            min_cluster_size=SELECTED_MIN_CLUSTER_SIZE,
            min_samples=SELECTED_MIN_SAMPLES,
            metric="jaccard",
            cluster_selection_method="eom",
            locked=True,
        ),
    )

    selected_condition_ids = {matrix.condition_id for matrix in matrices}
    cluster_assignment_rows: list[dict[str, Any]] = []
    medoid_rows: list[dict[str, Any]] = []
    condition_summary_rows: list[dict[str, Any]] = []
    pose_ids_to_include: set[str] = set()
    missing_optional_metadata: set[str] = set()
    condition_metadata_by_id: dict[str, dict[str, Any]] = {}

    for matrix in matrices:
        source = source_by_root[matrix.source_output_root]
        condition_metadata_by_id[matrix.condition_id] = _condition_metadata(matrix.condition_id, matrix.selected_row)
        pose_ids_to_include.update(
            pose_id
            for pose_id, condition_id in source.pose_condition_by_id.items()
            if condition_id == matrix.condition_id
        )
        pose_ids_to_include.update(matrix.pose_ids)

        result = clusterer.cluster(matrix.matrix, matrix.pose_ids)
        distance_matrix = clusterer.compute_jaccard_distances(matrix.matrix)
        cluster_assignment_rows.extend(
            build_cluster_assignment_rows(
                matrix.condition_id,
                matrix.pose_ids,
                result.cluster_labels,
                distance_matrix=distance_matrix,
                medoids=result.medoids,
            )
        )
        structure_paths = {
            pose_id: source.normalized_path_by_pose_id.get(pose_id, "")
            for pose_id in matrix.pose_ids
        }
        if any(not value for value in structure_paths.values()):
            missing_optional_metadata.add("medoid_structure_path")
        medoid_rows.extend(
            build_medoid_rows(
                matrix.condition_id,
                matrix.pose_ids,
                result.medoids,
                result.medoid_distance_sums,
                structure_paths=structure_paths,
            )
        )

        base_summary = dict(source.condition_summary_by_id.get(matrix.condition_id, {}))
        if not base_summary:
            missing_optional_metadata.add("condition_cluster_summary_source_row")
        n_qc_pass = int(float(base_summary.get("n_qc_pass_poses") or matrix.selected_row.get("original_n_qc_pass_poses") or 0))
        cluster_summary = build_condition_cluster_summary(
            matrix.condition_id,
            result,
            n_qc_pass_poses=n_qc_pass,
        ).to_row()
        condition_summary_rows.append(
            {
                **cluster_summary,
                **{
                    key: base_summary.get(key, cluster_summary.get(key, ""))
                    for key in (
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
                    )
                },
                **_cluster_field_overrides(cluster_summary),
            }
        )

    pose_ifp_results: list[IFPResult] = []
    pose_ifp_rows: list[dict[str, str]] = []
    pose_geometry_rows: list[dict[str, str]] = []
    pose_residue_contact_rows: list[dict[str, str]] = []
    pose_metadata_by_id: dict[str, dict[str, str]] = {}

    selected_rows_by_condition = {str(row["condition_id"]): row for row in rows}
    for source in source_by_root.values():
        for pose_id in sorted(pose_ids_to_include):
            condition_id = source.pose_condition_by_id.get(pose_id)
            if condition_id not in selected_condition_ids:
                continue
            if pose_id in source.pose_ifp_results_by_id:
                pose_ifp_results.append(source.pose_ifp_results_by_id[pose_id])
                pose_ifp_rows.append(source.pose_ifp_rows_by_id[pose_id])
            if pose_id in source.pose_geometry_rows_by_id:
                pose_geometry_rows.append(source.pose_geometry_rows_by_id[pose_id])
            selected_row = selected_rows_by_condition[condition_id]
            pose_metadata_by_id[pose_id] = {
                "protein_id": selected_row.get("protein_id", ""),
                "ligand_id": _ligand_id_from_pose_id(pose_id, selected_row.get("protein_id", ""), selected_row.get("target", "")),
                "condition_id": condition_id,
            }
        pose_residue_contact_rows.extend(
            row
            for row in source.pose_residue_contact_rows
            if row.get("condition_id") in selected_condition_ids
        )

    geometry_rows_by_pose_id = {
        str(row.get("pose_id", "")): row
        for row in pose_geometry_rows
        if row.get("pose_id")
    }
    cluster_signature_tables = build_cluster_signature_tables(
        cluster_assignment_rows=cluster_assignment_rows,
        medoid_rows=medoid_rows,
        geometry_rows_by_pose_id=geometry_rows_by_pose_id,
        ifp_results=pose_ifp_results,
        residue_contact_rows=pose_residue_contact_rows,
        pose_metadata_by_id=pose_metadata_by_id,
    )
    residue_importance = compute_residue_importance_outputs(
        cluster_signature_tables.residue_signature_rows,
        cluster_signature_tables.cluster_summaries,
        cluster_ifp_signature_rows=cluster_signature_tables.ifp_signature_rows,
        observed_residue_contact_rows=pose_residue_contact_rows,
        condition_metadata_by_id=condition_metadata_by_id,
    )

    _write_tsv(staging_root / "cluster_assignments.tsv", CLUSTER_ASSIGNMENT_COLUMNS, cluster_assignment_rows)
    _write_tsv(staging_root / "medoid_manifest.tsv", MEDOID_COLUMNS, medoid_rows)
    _write_tsv(staging_root / "condition_cluster_summary.tsv", CONDITION_CLUSTER_SUMMARY_COLUMNS, condition_summary_rows)
    write_cluster_ifp_signature_table(cluster_signature_tables.ifp_signature_rows, staging_root / "cluster_ifp_signature.tsv")
    write_cluster_residue_signature_table(
        cluster_signature_tables.residue_signature_rows,
        staging_root / "cluster_residue_signature.tsv",
    )
    write_cluster_signature_summary_json(cluster_signature_tables.cluster_summaries, staging_root / "cluster_signatures.json")
    write_protein_condition_residue_scores(
        residue_importance.protein_condition_residue_scores,
        staging_root / "protein_condition_residue_scores.tsv",
    )
    write_protein_residue_regio_delta(
        residue_importance.protein_residue_regio_delta,
        staging_root / "protein_residue_regio_delta.tsv",
    )
    write_condition_patch_summary(residue_importance.condition_patch_summary, staging_root / "condition_patch_summary.tsv")
    write_protein_patch_summary(residue_importance.protein_patch_summary, staging_root / "protein_patch_summary.tsv")
    write_pose_ifp_table(pose_ifp_results, staging_root / "pose_ifp_table.tsv")
    _write_tsv(
        staging_root / "pose_residue_contact_table.tsv",
        list(pose_residue_contact_rows[0].keys()) if pose_residue_contact_rows else [],
        pose_residue_contact_rows,
    )
    _write_tsv(
        staging_root / "pose_geometry.tsv",
        list(pose_geometry_rows[0].keys()) if pose_geometry_rows else [],
        pose_geometry_rows,
    )

    n_conditions_with_nonnoise_clusters = sum(1 for row in condition_summary_rows if int(float(row.get("n_clusters") or 0)) > 0)
    n_clusters = len(cluster_signature_tables.cluster_summaries)
    summary = {
        "input_selected_grid_results": str(selected_grid_results.resolve()),
        "input_pilot_source_roots": [str(root) for root in sorted(source_by_root)],
        "output_root": str(output_root),
        "parameter_label": SELECTED_PARAMETER_LABEL,
        "method": SELECTED_METHOD,
        "linkage": "average",
        "distance_threshold": SELECTED_DISTANCE_THRESHOLD,
        "min_cluster_size": SELECTED_MIN_CLUSTER_SIZE,
        "n_conditions_total": len(matrices),
        "n_conditions_formal_clusterable_original": sum(
            1 for row in rows if str(row.get("original_formal_clustering_allowed", "")).lower() == "true"
        ),
        "n_conditions_with_nonnoise_clusters": n_conditions_with_nonnoise_clusters,
        "n_clusters": n_clusters,
        "n_conditions_all_noise_or_empty": len(matrices) - n_conditions_with_nonnoise_clusters,
        "n_pose_ifp_rows": len(pose_ifp_results),
        "n_pose_residue_contact_rows": len(pose_residue_contact_rows),
        "n_pose_geometry_rows": len(pose_geometry_rows),
        "missing_optional_metadata": sorted(missing_optional_metadata),
        "git_commit": _git_commit(),
    }
    (staging_root / "fixture_build_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    _write_readme(staging_root)
    _write_selected_rows_for_output(rows, output_root, staging_root)
    if (output_root / "pilot_clustering_parameter_grid_summary.tsv").exists():
        shutil.copy2(
            output_root / "pilot_clustering_parameter_grid_summary.tsv",
            staging_root / "pilot_clustering_parameter_grid_summary.tsv",
        )
    _copy_selected_provenance(rows, staging_root)
    return summary


def _git_commit() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            cwd=Path(__file__).resolve().parents[2],
        )
    except Exception:
        return ""
    return result.stdout.strip()


def _validate_output(root: Path) -> None:
    missing = [name for name in ROOT_OUTPUT_FILES if not (root / name).exists()]
    if missing:
        raise FileNotFoundError(f"Fixture output missing required files: {missing}")
    clusters = json.loads((root / "cluster_signatures.json").read_text()).get("clusters", [])
    if not clusters:
        raise ValueError("cluster_signatures.json contains no non-noise clusters")


def build_fixture(
    *,
    selected_grid_results: Path,
    output_root: Path,
    allow_condition_count_mismatch: bool = False,
) -> dict[str, Any]:
    rows = _read_tsv(selected_grid_results)
    _validate_selected_rows(rows, allow_condition_count_mismatch)
    matrices = _load_matrix_inputs(rows)
    source_by_root = {
        root: _load_source_output_data(root)
        for root in sorted({matrix.source_output_root for matrix in matrices})
    }

    staging_root = output_root.with_name(f"{output_root.name}.__tmp__")
    backup_root = output_root.with_name(f"{output_root.name}.__previous__")
    if staging_root.exists():
        shutil.rmtree(staging_root)
    staging_root.mkdir(parents=True)
    summary = _write_root_outputs(
        staging_root=staging_root,
        output_root=output_root,
        selected_grid_results=selected_grid_results,
        matrices=matrices,
        rows=rows,
        source_by_root=source_by_root,
    )
    _validate_output(staging_root)

    if backup_root.exists():
        shutil.rmtree(backup_root)
    if output_root.exists():
        output_root.rename(backup_root)
    staging_root.rename(output_root)
    if backup_root.exists():
        shutil.rmtree(backup_root)
    return summary


def main() -> int:
    args = _parse_args()
    summary = build_fixture(
        selected_grid_results=Path(args.selected_grid_results),
        output_root=Path(args.output_root),
        allow_condition_count_mismatch=args.allow_condition_count_mismatch,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
