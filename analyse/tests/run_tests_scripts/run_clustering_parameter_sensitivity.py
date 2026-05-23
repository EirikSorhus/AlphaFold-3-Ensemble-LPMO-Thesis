#!/usr/bin/env python3
"""Run clustering parameter sensitivity from existing pilot IFP matrices.

This script intentionally does not rerun discovery, staging, ProLIF, geometry,
or contact eligibility. It reads the already-written
``main_contact_eligible_ifp_matrix.csv`` files from a staged clustering pilot
and writes a separate parameter-sensitivity result tree.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from lpmo_pipeline.analysis.cluster_signatures import write_cluster_table_tsv
from lpmo_pipeline.analysis.clustering_agglomerative import (
    AgglomerativeJaccardClusterer,
    AgglomerativeJaccardConfig,
)
from lpmo_pipeline.analysis.clustering_hdbscan import (
    HDBSCANClusterer,
    HDBSCANConfig,
    build_cluster_assignment_rows,
    build_medoid_rows,
)
from lpmo_pipeline.analysis.condition_summary import build_condition_table_rows, write_condition_table


DEFAULT_INPUT_RUN_ROOT = Path("tests/tests_results/clustering_pilot_staged")
DEFAULT_HDBSCAN_MIN_CLUSTER_SIZES = (3, 5, 10)
DEFAULT_AGGLOMERATIVE_DISTANCE_THRESHOLDS = (0.35, 0.45, 0.55)
DEFAULT_AGGLOMERATIVE_MIN_CLUSTER_SIZES = (3, 5, 10)

CONDITION_CLUSTER_SUMMARY_COLUMNS = [
    "method",
    "parameter_label",
    "condition_id",
    "n_selected_poses",
    "raw_feature_count",
    "main_feature_count",
    "minimum_clusterable_n",
    "clustering_status",
    "formal_clustering_allowed",
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
    "cluster_assignments_tsv",
    "medoid_manifest_tsv",
]

PARAMETER_RUN_INDEX_COLUMNS = [
    "method",
    "parameter_label",
    "cluster_assignments_tsv",
    "medoid_manifest_tsv",
    "condition_cluster_summary_tsv",
    "cluster_table_tsv",
    "condition_table_tsv",
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
    original_summary: dict[str, str]


@dataclass(frozen=True)
class ConditionOutputBundle:
    grid_row: dict[str, Any]
    assignment_rows: list[dict[str, Any]]
    medoid_rows: list[dict[str, Any]]
    condition_cluster_summary_row: dict[str, Any]
    cluster_table_rows: list[dict[str, Any]]
    qc_attrition_row: dict[str, Any]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run HDBSCAN/agglomerative parameter sensitivity from existing pilot matrices.",
    )
    parser.add_argument(
        "--input-run-root",
        default=str(DEFAULT_INPUT_RUN_ROOT),
        help="Existing staged clustering pilot root.",
    )
    parser.add_argument(
        "--output-root",
        required=True,
        help="New output directory for parameter sensitivity artifacts.",
    )
    parser.add_argument(
        "--hdbscan-min-cluster-sizes",
        nargs="+",
        type=int,
        default=list(DEFAULT_HDBSCAN_MIN_CLUSTER_SIZES),
        help="HDBSCAN min_cluster_size values.",
    )
    parser.add_argument(
        "--hdbscan-min-samples",
        default="none",
        help="Fixed HDBSCAN min_samples value, or 'none'.",
    )
    parser.add_argument(
        "--agglomerative-distance-thresholds",
        nargs="+",
        type=float,
        default=list(DEFAULT_AGGLOMERATIVE_DISTANCE_THRESHOLDS),
        help="Agglomerative Jaccard distance thresholds.",
    )
    parser.add_argument(
        "--agglomerative-min-cluster-sizes",
        nargs="+",
        type=int,
        default=list(DEFAULT_AGGLOMERATIVE_MIN_CLUSTER_SIZES),
        help="Agglomerative post-clustering minimum cluster sizes.",
    )
    parser.add_argument(
        "--allow-existing-empty-output-dir",
        action="store_true",
        help="Allow output-root if it exists and is empty. Non-empty directories always fail.",
    )
    return parser.parse_args()


def _slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def _float_slug(value: float) -> str:
    return str(value).replace(".", "p")


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
            writer.writerow(row)


def _read_matrix(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if not header or header[0] != "pose_id":
            raise ValueError(f"IFP matrix must start with pose_id column: {path}")
        feature_names = header[1:]
        pose_ids: list[str] = []
        rows: list[list[int]] = []
        for row in reader:
            if not row:
                continue
            pose_ids.append(row[0])
            values = [int(value) for value in row[1:]]
            if len(values) != len(feature_names):
                raise ValueError(f"Matrix row has wrong feature count in {path}: {row[0]}")
            rows.append(values)
    matrix = np.asarray(rows, dtype=np.uint8)
    if not rows:
        matrix = np.zeros((0, len(feature_names)), dtype=np.uint8)
    return pose_ids, feature_names, matrix


def _load_original_method_rows(input_run_root: Path) -> dict[str, dict[str, str]]:
    """Return one metadata row per condition from original pilot summaries."""
    rows_by_condition: dict[str, dict[str, str]] = {}
    for path in sorted(input_run_root.glob("*_shards/shard_*/production_output/clustering_pilot/*/pilot_clustering_method_summary.tsv")):
        for row in _read_tsv(path):
            condition_id = row.get("condition_id", "")
            if not condition_id:
                continue
            # Prefer the HDBSCAN row just to avoid duplicate method-specific rows.
            if condition_id not in rows_by_condition or row.get("method") == "hdbscan":
                rows_by_condition[condition_id] = row
    return rows_by_condition


def _construct_from_path(path: Path) -> str:
    for part in path.parts:
        if part in {"domain_only_shards", "full_length_shards"}:
            return part.removesuffix("_shards")
    return ""


def _condition_parts(condition_id: str, construct_type: str) -> tuple[str, str]:
    parts = condition_id.split("__")
    protein_id = parts[0] if parts else ""
    target = parts[2] if len(parts) >= 3 else ""
    if len(parts) >= 2 and not construct_type:
        construct_type = parts[1]
    return protein_id, target


def _discover_matrices(input_run_root: Path) -> list[MatrixInput]:
    original_rows = _load_original_method_rows(input_run_root)
    matrices: list[MatrixInput] = []
    for matrix_path in sorted(
        input_run_root.glob("*_shards/shard_*/production_output/clustering_pilot/*/conditions/*/main_contact_eligible_ifp_matrix.csv")
    ):
        condition_id = matrix_path.parent.name
        construct_type = _construct_from_path(matrix_path)
        protein_id, target = _condition_parts(condition_id, construct_type)
        pose_ids, feature_names, matrix = _read_matrix(matrix_path)
        matrices.append(
            MatrixInput(
                condition_id=condition_id,
                construct_type=construct_type,
                protein_id=protein_id,
                target=target,
                matrix_path=matrix_path,
                pose_ids=pose_ids,
                feature_names=feature_names,
                matrix=matrix,
                original_summary=original_rows.get(condition_id, {}),
            )
        )
    return matrices


def _parse_min_samples(value: str) -> int | None:
    normalized = value.strip().lower()
    if normalized in {"", "none", "null"}:
        return None
    return int(normalized)


def _gini(values: list[int]) -> float:
    if not values:
        return 0.0
    sorted_values = np.sort(np.asarray(values, dtype=float))
    total = float(np.sum(sorted_values))
    if total == 0.0:
        return 0.0
    n_values = sorted_values.size
    ranks = np.arange(1, n_values + 1, dtype=float)
    return float((2.0 * np.sum(ranks * sorted_values) / (n_values * total)) - ((n_values + 1) / n_values))


def _as_int(value: Any, default: int = 0) -> int:
    if value in {None, ""}:
        return default
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _as_float(value: Any, default: float = 0.0) -> float:
    if value in {None, ""}:
        return default
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return default
    if math.isnan(numeric):
        return default
    return numeric


def _condition_identity(matrix_input: MatrixInput) -> dict[str, Any]:
    match = re.match(r"(?P<substrate>.+)_DP(?P<dp>[0-9]+)$", matrix_input.target)
    substrate_class = match.group("substrate") if match else ""
    dp = int(match.group("dp")) if match else ""
    return {
        "condition_id": matrix_input.condition_id,
        "protein_id": matrix_input.protein_id,
        "construct_type": matrix_input.construct_type,
        "ligand_id": matrix_input.target,
        "substrate_class": substrate_class,
        "dp": dp,
    }


def _clustering_status(metrics: dict[str, Any]) -> str:
    n_ifp_clustered = int(metrics["n_ifp_clustered"])
    n_clusters = int(metrics["n_clusters"])
    n_noise = int(metrics["n_noise"])
    if n_ifp_clustered == 0:
        return "empty_matrix"
    if n_clusters > 0:
        return "clustered"
    if n_noise == n_ifp_clustered:
        return "all_noise"
    return "no_clusters"


def _build_condition_cluster_summary_row(
    matrix_input: MatrixInput,
    metadata: dict[str, Any],
    metrics: dict[str, Any],
    *,
    assignment_path: Path,
    medoid_path: Path,
) -> dict[str, Any]:
    original = matrix_input.original_summary
    n_selected_poses = len(matrix_input.pose_ids)
    n_qc_pass_poses = _as_int(original.get("n_qc_pass_poses"), n_selected_poses)
    n_ifp_success = _as_int(original.get("n_ifp_success"), n_selected_poses)
    n_contact_eligible = _as_int(original.get("n_contact_eligible"), n_selected_poses)
    contact_eligible_fraction = _as_float(
        original.get("contact_eligible_fraction"),
        (n_contact_eligible / n_ifp_success) if n_ifp_success else 0.0,
    )
    return {
        "method": metadata["method"],
        "parameter_label": metadata["parameter_label"],
        "condition_id": matrix_input.condition_id,
        "n_selected_poses": n_selected_poses,
        "raw_feature_count": _as_int(original.get("raw_feature_count"), len(matrix_input.feature_names)),
        "main_feature_count": len(matrix_input.feature_names),
        "minimum_clusterable_n": original.get("minimum_clusterable_n", ""),
        "clustering_status": _clustering_status(metrics),
        "formal_clustering_allowed": original.get("formal_clustering_allowed", ""),
        "n_qc_pass_poses": n_qc_pass_poses,
        "n_ifp_success": n_ifp_success,
        "n_contact_eligible": n_contact_eligible,
        "contact_eligible_fraction": contact_eligible_fraction,
        "null_ifp_fraction": _as_float(original.get("null_ifp_fraction"), 0.0),
        "vdw_only_fraction": _as_float(original.get("vdw_only_fraction"), 0.0),
        "low_specific_contact_fraction": _as_float(original.get("low_specific_contact_fraction"), 0.0),
        "median_n_non_vdw_interactions": _as_float(original.get("median_n_non_vdw_interactions"), 0.0),
        "median_n_non_vdw_contact_residues": _as_float(original.get("median_n_non_vdw_contact_residues"), 0.0),
        "n_ifp_clustered": metrics["n_ifp_clustered"],
        "n_noise": metrics["n_noise"],
        "noise_fraction": metrics["noise_fraction"],
        "n_clusters": metrics["n_clusters"],
        "top_cluster_occupancy": metrics["top_cluster_occupancy"],
        "cluster_entropy": metrics["cluster_entropy"],
        "occupancy_gini": metrics["occupancy_gini"],
        "cluster_assignments_tsv": str(assignment_path),
        "medoid_manifest_tsv": str(medoid_path),
    }


def _build_cluster_table_rows(matrix_input: MatrixInput, result: Any) -> list[dict[str, Any]]:
    identity = _condition_identity(matrix_input)
    total_poses = len(matrix_input.pose_ids)
    rows: list[dict[str, Any]] = []
    for cluster_id in sorted(result.cluster_sizes):
        medoid_pose_id = ""
        if cluster_id in result.medoids:
            medoid_pose_id = matrix_input.pose_ids[result.medoids[cluster_id]]
        cluster_size = int(result.cluster_sizes[cluster_id])
        occupancy = result.cluster_occupancy.get(cluster_id)
        if occupancy is None:
            occupancy = (cluster_size / total_poses) if total_poses else 0.0
        rows.append(
            {
                **identity,
                "cluster_id": cluster_id,
                "cluster_type": "",
                "n_poses": cluster_size,
                "occupancy": occupancy,
                "medoid_pose_id": medoid_pose_id,
                "cluster_size": cluster_size,
            }
        )
    return rows


def _build_qc_attrition_row(matrix_input: MatrixInput) -> dict[str, Any]:
    identity = _condition_identity(matrix_input)
    n_generated = len(matrix_input.pose_ids)
    n_stage1_pass = _as_int(matrix_input.original_summary.get("n_qc_pass_poses"), n_generated)
    n_stage1_hard_fail = max(n_generated - n_stage1_pass, 0)
    return {
        **identity,
        "n_generated": n_generated,
        "n_prepared": n_generated,
        "n_prep_error": 0,
        "n_stage1_hard_fail": n_stage1_hard_fail,
        "n_stage1_soft_flag": 0,
        "n_stage1_pass": n_stage1_pass,
        "hard_fail_rate": (n_stage1_hard_fail / n_generated) if n_generated else 0.0,
    }


def _cluster_metrics(labels: np.ndarray) -> dict[str, Any]:
    total = int(len(labels))
    n_noise = int(np.sum(labels == -1))
    cluster_sizes = {
        int(label): int(np.sum(labels == label))
        for label in sorted({int(value) for value in labels.tolist()})
        if int(label) != -1
    }
    sizes = list(cluster_sizes.values())
    if sizes:
        size_array = np.asarray(sizes, dtype=float)
        probabilities = size_array / np.sum(size_array)
        top_cluster_occupancy = float(np.max(probabilities))
        cluster_entropy = float(
            -sum(float(probability) * math.log2(float(probability)) for probability in probabilities if probability > 0)
        )
        occupancy_gini = _gini(sizes)
        top_cluster_id = max(sorted(cluster_sizes), key=lambda cluster_id: cluster_sizes[cluster_id])
    else:
        top_cluster_occupancy = 0.0
        cluster_entropy = 0.0
        occupancy_gini = 0.0
        top_cluster_id = ""
    return {
        "n_ifp_clustered": total,
        "n_noise": n_noise,
        "noise_fraction": (n_noise / total) if total else 0.0,
        "n_clusters": len(cluster_sizes),
        "cluster_sizes_json": json.dumps(cluster_sizes, sort_keys=True),
        "top_cluster_id": top_cluster_id,
        "top_cluster_occupancy": top_cluster_occupancy,
        "cluster_entropy": cluster_entropy,
        "occupancy_gini": occupancy_gini,
    }


def _seed_from_pose_id(pose_id: str) -> str:
    match = re.search(r"seed-([0-9]+)", pose_id)
    return match.group(1) if match else ""


def _seed_mixing(labels: np.ndarray, pose_ids: list[str]) -> dict[str, Any]:
    max_single_seed_fraction = 0.0
    worst_cluster_id: int | str = ""
    for label in sorted({int(value) for value in labels.tolist()}):
        if label == -1:
            continue
        indices = np.where(labels == label)[0].tolist()
        seeds = [_seed_from_pose_id(pose_ids[index]) for index in indices]
        counts = Counter(seeds)
        if not counts:
            continue
        fraction = max(counts.values()) / len(indices)
        if fraction > max_single_seed_fraction:
            max_single_seed_fraction = fraction
            worst_cluster_id = label
    return {
        "max_single_seed_fraction": max_single_seed_fraction,
        "seed_artifact_flag": max_single_seed_fraction > 0.70,
        "seed_artifact_cluster_id": worst_cluster_id,
    }


def _all_noise_result(n_rows: int):
    from lpmo_pipeline.analysis.clustering_hdbscan import ClusteringResult

    labels = np.full(n_rows, -1, dtype=int)
    return ClusteringResult(
        n_clusters=0,
        n_outliers=n_rows,
        outlier_rate=1.0 if n_rows else 0.0,
        cluster_sizes={},
        cluster_occupancy={},
        outlier_indices=list(range(n_rows)),
        cluster_labels=labels,
        medoids={},
        medoid_distance_sums={},
    )


def _run_hdbscan(matrix_input: MatrixInput, min_cluster_size: int, min_samples: int | None, output_dir: Path):
    if len(matrix_input.feature_names) == 0 and len(matrix_input.pose_ids) > 0:
        return _all_noise_result(len(matrix_input.pose_ids)), np.zeros((len(matrix_input.pose_ids), len(matrix_input.pose_ids)))
    clusterer = HDBSCANClusterer(
        output_dir=output_dir,
        config=HDBSCANConfig(
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            metric="jaccard",
            cluster_selection_epsilon=0.0,
            cluster_selection_method="eom",
            locked=False,
        ),
    )
    result = clusterer.cluster(matrix_input.matrix, matrix_input.pose_ids)
    return result, clusterer.compute_jaccard_distances(matrix_input.matrix)


def _run_agglomerative(
    matrix_input: MatrixInput,
    distance_threshold: float,
    min_cluster_size: int,
    output_dir: Path,
):
    if len(matrix_input.feature_names) == 0 and len(matrix_input.pose_ids) > 0:
        return _all_noise_result(len(matrix_input.pose_ids)), np.zeros((len(matrix_input.pose_ids), len(matrix_input.pose_ids)))
    clusterer = AgglomerativeJaccardClusterer(
        output_dir=output_dir,
        config=AgglomerativeJaccardConfig(
            linkage="average",
            distance_threshold=distance_threshold,
            min_cluster_size=min_cluster_size,
        ),
    )
    result = clusterer.cluster(matrix_input.matrix, matrix_input.pose_ids)
    return result, clusterer.compute_jaccard_distances(matrix_input.matrix)


def _write_condition_outputs(
    *,
    matrix_input: MatrixInput,
    method_output_dir: Path,
    result: Any,
    distance_matrix: np.ndarray,
    metadata: dict[str, Any],
) -> ConditionOutputBundle:
    labels = result.cluster_labels
    assignment_rows = build_cluster_assignment_rows(
        matrix_input.condition_id,
        matrix_input.pose_ids,
        labels,
        distance_matrix=distance_matrix,
        medoids=result.medoids,
    )
    medoid_rows = build_medoid_rows(
        matrix_input.condition_id,
        matrix_input.pose_ids,
        result.medoids,
        result.medoid_distance_sums,
    )

    assignment_path = method_output_dir / "cluster_assignments.tsv"
    medoid_path = method_output_dir / "medoid_manifest.tsv"
    _write_tsv(
        assignment_path,
        [
            "pose_id",
            "condition_id",
            "cluster_id",
            "cluster_member_flag",
            "noise_flag",
            "distance_to_cluster_representative",
        ],
        assignment_rows,
    )
    _write_tsv(
        medoid_path,
        [
            "condition_id",
            "cluster_id",
            "medoid_pose_id",
            "medoid_structure_path",
            "medoid_ifp_distance_sum",
        ],
        medoid_rows,
    )

    metrics = _cluster_metrics(labels)
    seed_metrics = _seed_mixing(labels, matrix_input.pose_ids)
    condition_cluster_summary_row = _build_condition_cluster_summary_row(
        matrix_input,
        metadata,
        metrics,
        assignment_path=assignment_path,
        medoid_path=medoid_path,
    )
    cluster_table_rows = _build_cluster_table_rows(matrix_input, result)
    qc_attrition_row = _build_qc_attrition_row(matrix_input)
    condition_table_rows = build_condition_table_rows(
        qc_attrition_rows=[qc_attrition_row],
        condition_cluster_summary_rows=[condition_cluster_summary_row],
        cluster_table_rows=cluster_table_rows,
    )
    condition_cluster_summary_path = method_output_dir / "condition_cluster_summary.tsv"
    cluster_table_path = method_output_dir / "cluster_table.tsv"
    condition_table_path = method_output_dir / "condition_table.tsv"
    _write_tsv(
        condition_cluster_summary_path,
        CONDITION_CLUSTER_SUMMARY_COLUMNS,
        [condition_cluster_summary_row],
    )
    write_cluster_table_tsv(cluster_table_rows, cluster_table_path)
    write_condition_table(condition_table_rows, condition_table_path)

    top_cluster_id = metrics["top_cluster_id"]
    top_medoid_pose_id = ""
    if top_cluster_id != "" and int(top_cluster_id) in result.medoids:
        top_medoid_pose_id = matrix_input.pose_ids[result.medoids[int(top_cluster_id)]]
    row = {
        **metadata,
        **metrics,
        **seed_metrics,
        "top_medoid_pose_id": top_medoid_pose_id,
        "cluster_assignments_tsv": str(assignment_path),
        "medoid_manifest_tsv": str(medoid_path),
        "condition_cluster_summary_tsv": str(condition_cluster_summary_path),
        "cluster_table_tsv": str(cluster_table_path),
        "condition_table_tsv": str(condition_table_path),
    }
    summary_path = method_output_dir / "clustering_summary.json"
    summary_path.write_text(json.dumps(row, indent=2))
    return ConditionOutputBundle(
        grid_row=row,
        assignment_rows=assignment_rows,
        medoid_rows=medoid_rows,
        condition_cluster_summary_row=condition_cluster_summary_row,
        cluster_table_rows=cluster_table_rows,
        qc_attrition_row=qc_attrition_row,
    )


def _metadata_for(matrix_input: MatrixInput, *, method: str, parameter_label: str, **params: Any) -> dict[str, Any]:
    original = matrix_input.original_summary
    return {
        "condition_id": matrix_input.condition_id,
        "construct_type": matrix_input.construct_type,
        "protein_id": matrix_input.protein_id,
        "target": matrix_input.target,
        "method": method,
        "parameter_label": parameter_label,
        "distance_threshold": params.get("distance_threshold", ""),
        "min_cluster_size": params.get("min_cluster_size", ""),
        "min_samples": params.get("min_samples", ""),
        "n_selected_poses": len(matrix_input.pose_ids),
        "main_feature_count": len(matrix_input.feature_names),
        "matrix_path": str(matrix_input.matrix_path),
        "original_minimum_clusterable_n": original.get("minimum_clusterable_n", ""),
        "original_clustering_status": original.get("clustering_status", ""),
        "original_formal_clustering_allowed": original.get("formal_clustering_allowed", ""),
        "original_n_qc_pass_poses": original.get("n_qc_pass_poses", ""),
        "original_n_ifp_success": original.get("n_ifp_success", ""),
        "original_n_contact_eligible": original.get("n_contact_eligible", ""),
    }


def _write_parameter_run_outputs(
    output_root: Path,
    bundles_by_parameter_label: dict[str, list[ConditionOutputBundle]],
) -> list[dict[str, Any]]:
    index_rows: list[dict[str, Any]] = []
    for parameter_label, bundles in sorted(bundles_by_parameter_label.items()):
        parameter_dir = output_root / "parameter_runs" / _slug(parameter_label)
        assignment_path = parameter_dir / "cluster_assignments.tsv"
        medoid_path = parameter_dir / "medoid_manifest.tsv"
        condition_cluster_summary_path = parameter_dir / "condition_cluster_summary.tsv"
        cluster_table_path = parameter_dir / "cluster_table.tsv"
        condition_table_path = parameter_dir / "condition_table.tsv"

        assignment_rows = [row for bundle in bundles for row in bundle.assignment_rows]
        medoid_rows = [row for bundle in bundles for row in bundle.medoid_rows]
        condition_cluster_summary_rows = [bundle.condition_cluster_summary_row for bundle in bundles]
        cluster_table_rows = [row for bundle in bundles for row in bundle.cluster_table_rows]
        qc_attrition_rows = [bundle.qc_attrition_row for bundle in bundles]
        condition_table_rows = build_condition_table_rows(
            qc_attrition_rows=qc_attrition_rows,
            condition_cluster_summary_rows=condition_cluster_summary_rows,
            cluster_table_rows=cluster_table_rows,
        )

        _write_tsv(
            assignment_path,
            [
                "pose_id",
                "condition_id",
                "cluster_id",
                "cluster_member_flag",
                "noise_flag",
                "distance_to_cluster_representative",
            ],
            assignment_rows,
        )
        _write_tsv(
            medoid_path,
            [
                "condition_id",
                "cluster_id",
                "medoid_pose_id",
                "medoid_structure_path",
                "medoid_ifp_distance_sum",
            ],
            medoid_rows,
        )
        _write_tsv(
            condition_cluster_summary_path,
            CONDITION_CLUSTER_SUMMARY_COLUMNS,
            condition_cluster_summary_rows,
        )
        write_cluster_table_tsv(cluster_table_rows, cluster_table_path)
        write_condition_table(condition_table_rows, condition_table_path)

        index_rows.append(
            {
                "method": bundles[0].grid_row["method"],
                "parameter_label": parameter_label,
                "cluster_assignments_tsv": str(assignment_path),
                "medoid_manifest_tsv": str(medoid_path),
                "condition_cluster_summary_tsv": str(condition_cluster_summary_path),
                "cluster_table_tsv": str(cluster_table_path),
                "condition_table_tsv": str(condition_table_path),
            }
        )
    return index_rows


def _summary_rows(grid_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in grid_rows:
        grouped[str(row["parameter_label"])].append(row)

    summaries: list[dict[str, Any]] = []
    for parameter_label, rows in sorted(grouped.items()):
        n_conditions = len(rows)
        n_clusters = [int(row["n_clusters"]) for row in rows]
        noise_fractions = [float(row["noise_fraction"]) for row in rows]
        top_occupancies = [float(row["top_cluster_occupancy"]) for row in rows]
        summaries.append(
            {
                "method": rows[0]["method"],
                "parameter_label": parameter_label,
                "n_conditions": n_conditions,
                "n_conditions_with_clusters": sum(value > 0 for value in n_clusters),
                "total_clusters": sum(n_clusters),
                "median_n_clusters": float(np.median(n_clusters)) if n_clusters else 0.0,
                "median_noise_fraction": float(np.median(noise_fractions)) if noise_fractions else 0.0,
                "median_top_cluster_occupancy": float(np.median(top_occupancies)) if top_occupancies else 0.0,
                "n_conditions_high_noise_gt_0p75": sum(value > 0.75 for value in noise_fractions),
                "n_conditions_seed_artifact_flag": sum(bool(row["seed_artifact_flag"]) for row in rows),
            }
        )
    return summaries


def _write_selection_criteria(output_root: Path) -> None:
    text = """# Clustering Parameter Selection Criteria

These criteria were written before the parameter-sensitivity run and should be
used to interpret the grid without picking parameters condition-by-condition.

1. Choose one global primary clustering method and one global parameter set.
2. Do not choose the parameter set that simply maximizes the number of clusters.
3. Reject parameter sets with very high noise in many conditions, many tiny/seed-specific clusters, or strong instability across adjacent parameter values.
4. Prefer HDBSCAN min_cluster_size=5 if it balances cluster recovery, noise handling, and cluster granularity.
5. Use HDBSCAN min_cluster_size=3 as lenient sensitivity for low-support modes.
6. Use agglomerative Jaccard 0.55/min5 as the orthogonal method check.
7. Treat conditions with empty main IFP matrices or too few contact-eligible poses as weak evidence even if a small min_cluster_size creates clusters.
"""
    (output_root / "parameter_selection_criteria.md").write_text(text)


def main() -> int:
    args = _parse_args()
    input_run_root = Path(args.input_run_root).resolve()
    output_root = Path(args.output_root).resolve()
    min_samples = _parse_min_samples(str(args.hdbscan_min_samples))

    if not input_run_root.exists():
        raise FileNotFoundError(f"Input run root does not exist: {input_run_root}")
    if output_root.exists():
        if any(output_root.iterdir()):
            raise FileExistsError(f"Refusing to write into non-empty output directory: {output_root}")
        if not args.allow_existing_empty_output_dir:
            raise FileExistsError(
                f"Output directory already exists: {output_root}. "
                "Use --allow-existing-empty-output-dir only for an empty directory."
            )
    output_root.mkdir(parents=True, exist_ok=True)
    _write_selection_criteria(output_root)

    matrices = _discover_matrices(input_run_root)
    if not matrices:
        raise RuntimeError(f"No main_contact_eligible_ifp_matrix.csv files found under {input_run_root}")

    grid_rows: list[dict[str, Any]] = []
    bundles_by_parameter_label: dict[str, list[ConditionOutputBundle]] = defaultdict(list)
    for matrix_input in matrices:
        condition_dir = output_root / "conditions" / _slug(matrix_input.condition_id)

        for min_cluster_size in args.hdbscan_min_cluster_sizes:
            parameter_label = (
                f"hdbscan_jaccard__min_cluster_size_{min_cluster_size}__"
                f"min_samples_{'none' if min_samples is None else min_samples}"
            )
            method_output_dir = condition_dir / parameter_label
            result, distance_matrix = _run_hdbscan(
                matrix_input,
                min_cluster_size=min_cluster_size,
                min_samples=min_samples,
                output_dir=method_output_dir,
            )
            bundle = _write_condition_outputs(
                matrix_input=matrix_input,
                method_output_dir=method_output_dir,
                result=result,
                distance_matrix=distance_matrix,
                metadata=_metadata_for(
                    matrix_input,
                    method="hdbscan_jaccard",
                    parameter_label=parameter_label,
                    min_cluster_size=min_cluster_size,
                    min_samples=("none" if min_samples is None else min_samples),
                ),
            )
            grid_rows.append(bundle.grid_row)
            bundles_by_parameter_label[parameter_label].append(bundle)

        for distance_threshold in args.agglomerative_distance_thresholds:
            for min_cluster_size in args.agglomerative_min_cluster_sizes:
                parameter_label = (
                    "agglomerative_jaccard__"
                    f"threshold_{_float_slug(distance_threshold)}__"
                    f"min_cluster_size_{min_cluster_size}"
                )
                method_output_dir = condition_dir / parameter_label
                result, distance_matrix = _run_agglomerative(
                    matrix_input,
                    distance_threshold=distance_threshold,
                    min_cluster_size=min_cluster_size,
                    output_dir=method_output_dir,
                )
                bundle = _write_condition_outputs(
                    matrix_input=matrix_input,
                    method_output_dir=method_output_dir,
                    result=result,
                    distance_matrix=distance_matrix,
                    metadata=_metadata_for(
                        matrix_input,
                        method="agglomerative_jaccard",
                        parameter_label=parameter_label,
                        distance_threshold=distance_threshold,
                        min_cluster_size=min_cluster_size,
                    ),
                )
                grid_rows.append(bundle.grid_row)
                bundles_by_parameter_label[parameter_label].append(bundle)

    parameter_run_index_rows = _write_parameter_run_outputs(output_root, bundles_by_parameter_label)
    _write_tsv(output_root / "parameter_run_index.tsv", PARAMETER_RUN_INDEX_COLUMNS, parameter_run_index_rows)

    grid_fieldnames = [
        "condition_id",
        "construct_type",
        "protein_id",
        "target",
        "method",
        "parameter_label",
        "distance_threshold",
        "min_cluster_size",
        "min_samples",
        "n_selected_poses",
        "main_feature_count",
        "n_ifp_clustered",
        "n_noise",
        "noise_fraction",
        "n_clusters",
        "cluster_sizes_json",
        "top_cluster_id",
        "top_cluster_occupancy",
        "cluster_entropy",
        "occupancy_gini",
        "top_medoid_pose_id",
        "max_single_seed_fraction",
        "seed_artifact_flag",
        "seed_artifact_cluster_id",
        "original_minimum_clusterable_n",
        "original_clustering_status",
        "original_formal_clustering_allowed",
        "original_n_qc_pass_poses",
        "original_n_ifp_success",
        "original_n_contact_eligible",
        "matrix_path",
        "cluster_assignments_tsv",
        "medoid_manifest_tsv",
        "condition_cluster_summary_tsv",
        "cluster_table_tsv",
        "condition_table_tsv",
    ]
    _write_tsv(output_root / "pilot_clustering_parameter_grid_results.tsv", grid_fieldnames, grid_rows)

    summary_fieldnames = [
        "method",
        "parameter_label",
        "n_conditions",
        "n_conditions_with_clusters",
        "total_clusters",
        "median_n_clusters",
        "median_noise_fraction",
        "median_top_cluster_occupancy",
        "n_conditions_high_noise_gt_0p75",
        "n_conditions_seed_artifact_flag",
    ]
    summary_rows = _summary_rows(grid_rows)
    _write_tsv(output_root / "pilot_clustering_parameter_grid_summary.tsv", summary_fieldnames, summary_rows)

    run_summary = {
        "input_run_root": str(input_run_root),
        "output_root": str(output_root),
        "n_conditions": len(matrices),
        "n_grid_rows": len(grid_rows),
        "hdbscan_min_cluster_sizes": list(args.hdbscan_min_cluster_sizes),
        "hdbscan_min_samples": "none" if min_samples is None else min_samples,
        "agglomerative_distance_thresholds": list(args.agglomerative_distance_thresholds),
        "agglomerative_min_cluster_sizes": list(args.agglomerative_min_cluster_sizes),
        "grid_results_tsv": str(output_root / "pilot_clustering_parameter_grid_results.tsv"),
        "grid_summary_tsv": str(output_root / "pilot_clustering_parameter_grid_summary.tsv"),
        "parameter_run_index_tsv": str(output_root / "parameter_run_index.tsv"),
        "parameter_selection_criteria": str(output_root / "parameter_selection_criteria.md"),
        "derived_ligand_id_source": "condition target token from condition_id",
    }
    (output_root / "run_summary.json").write_text(json.dumps(run_summary, indent=2))
    print(json.dumps(run_summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
