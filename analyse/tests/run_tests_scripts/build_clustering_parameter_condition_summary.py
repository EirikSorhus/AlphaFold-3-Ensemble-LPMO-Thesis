#!/usr/bin/env python3
"""Build a wide one-row-per-condition summary from parameter sensitivity outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


DEFAULT_PROFILE_SPECS = (
    (
        "primary_hdbscan_min_5",
        "hdbscan_jaccard__min_cluster_size_5__min_samples_none",
    ),
    (
        "lenient_hdbscan_min_3",
        "hdbscan_jaccard__min_cluster_size_3__min_samples_none",
    ),
    (
        "orthogonal_agglomerative_thr_0p55_min_5",
        "agglomerative_jaccard__threshold_0p55__min_cluster_size_5",
    ),
)

BASE_COLUMNS = [
    "condition_id",
    "construct_type",
    "protein_id",
    "target",
    "n_selected_poses",
    "main_feature_count",
    "original_n_qc_pass_poses",
    "original_n_ifp_success",
    "original_n_contact_eligible",
    "original_clustering_status",
]

PROFILE_VALUE_FIELDS = [
    "n_clusters",
    "n_noise",
    "noise_fraction",
    "top_cluster_occupancy",
    "cluster_entropy",
    "occupancy_gini",
    "cluster_sizes_json",
    "cluster_size_mean",
    "cluster_size_sd",
    "cluster_size_min",
    "cluster_size_max",
    "cluster_size_median",
    "top_medoid_pose_id",
    "max_single_seed_fraction",
    "seed_artifact_flag",
    "cluster_assignments_tsv",
    "medoid_manifest_tsv",
]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a wide condition-level summary from clustering parameter sensitivity outputs.",
    )
    parser.add_argument(
        "--sensitivity-root",
        required=True,
        help="Root directory containing pilot_clustering_parameter_grid_results.tsv.",
    )
    parser.add_argument(
        "--output-path",
        required=True,
        help="TSV path for the wide condition summary.",
    )
    parser.add_argument(
        "--profile",
        action="append",
        default=[],
        help=(
            "Profile alias and parameter label as alias=parameter_label. "
            "May be repeated. Defaults to primary HDBSCAN min5, lenient HDBSCAN min3, "
            "and orthogonal agglomerative 0.55/min5 rows."
        ),
    )
    return parser.parse_args()


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def _parse_profile_specs(raw_specs: list[str]) -> list[tuple[str, str]]:
    if not raw_specs:
        return list(DEFAULT_PROFILE_SPECS)
    parsed: list[tuple[str, str]] = []
    for raw_spec in raw_specs:
        alias, separator, parameter_label = raw_spec.partition("=")
        if not separator or not alias or not parameter_label:
            raise ValueError(f"Invalid --profile value: {raw_spec!r}. Expected alias=parameter_label")
        parsed.append((alias, parameter_label))
    return parsed


def _cluster_sizes(cluster_sizes_json: str) -> list[int]:
    if not cluster_sizes_json:
        return []
    decoded = json.loads(cluster_sizes_json)
    if isinstance(decoded, dict):
        return [int(value) for value in decoded.values()]
    raise ValueError(f"Expected cluster_sizes_json to decode to an object, got: {type(decoded)!r}")


def _cluster_size_stats(cluster_sizes_json: str) -> dict[str, Any]:
    sizes = _cluster_sizes(cluster_sizes_json)
    if not sizes:
        return {
            "cluster_size_mean": "",
            "cluster_size_sd": "",
            "cluster_size_min": "",
            "cluster_size_max": "",
            "cluster_size_median": "",
        }
    values = np.asarray(sizes, dtype=float)
    return {
        "cluster_size_mean": float(np.mean(values)),
        "cluster_size_sd": float(np.std(values, ddof=0)),
        "cluster_size_min": int(np.min(values)),
        "cluster_size_max": int(np.max(values)),
        "cluster_size_median": float(np.median(values)),
    }


def _normalize_bool(value: Any) -> str:
    if isinstance(value, bool):
        return "True" if value else "False"
    text = str(value).strip().lower()
    if text in {"true", "1"}:
        return "True"
    if text in {"false", "0", ""}:
        return "False" if text else ""
    return str(value)


def build_condition_summary_rows(
    grid_rows: list[dict[str, str]],
    profile_specs: list[tuple[str, str]],
) -> tuple[list[str], list[dict[str, Any]]]:
    rows_by_condition: dict[str, dict[str, str]] = {}
    rows_by_condition_and_parameter: dict[tuple[str, str], dict[str, str]] = {}

    for row in grid_rows:
        condition_id = row["condition_id"]
        parameter_label = row["parameter_label"]
        rows_by_condition.setdefault(condition_id, row)
        rows_by_condition_and_parameter[(condition_id, parameter_label)] = row

    fieldnames = list(BASE_COLUMNS)
    for alias, _parameter_label in profile_specs:
        fieldnames.extend([f"{alias}_{field}" for field in PROFILE_VALUE_FIELDS])

    summary_rows: list[dict[str, Any]] = []
    for condition_id in sorted(rows_by_condition):
        exemplar = rows_by_condition[condition_id]
        summary_row: dict[str, Any] = {column: exemplar.get(column, "") for column in BASE_COLUMNS}

        for alias, parameter_label in profile_specs:
            parameter_row = rows_by_condition_and_parameter.get((condition_id, parameter_label), {})
            cluster_sizes_json = parameter_row.get("cluster_sizes_json", "")
            stats = _cluster_size_stats(cluster_sizes_json) if parameter_row else _cluster_size_stats("")

            summary_row[f"{alias}_n_clusters"] = parameter_row.get("n_clusters", "")
            summary_row[f"{alias}_n_noise"] = parameter_row.get("n_noise", "")
            summary_row[f"{alias}_noise_fraction"] = parameter_row.get("noise_fraction", "")
            summary_row[f"{alias}_top_cluster_occupancy"] = parameter_row.get("top_cluster_occupancy", "")
            summary_row[f"{alias}_cluster_entropy"] = parameter_row.get("cluster_entropy", "")
            summary_row[f"{alias}_occupancy_gini"] = parameter_row.get("occupancy_gini", "")
            summary_row[f"{alias}_cluster_sizes_json"] = cluster_sizes_json
            summary_row[f"{alias}_cluster_size_mean"] = stats["cluster_size_mean"]
            summary_row[f"{alias}_cluster_size_sd"] = stats["cluster_size_sd"]
            summary_row[f"{alias}_cluster_size_min"] = stats["cluster_size_min"]
            summary_row[f"{alias}_cluster_size_max"] = stats["cluster_size_max"]
            summary_row[f"{alias}_cluster_size_median"] = stats["cluster_size_median"]
            summary_row[f"{alias}_top_medoid_pose_id"] = parameter_row.get("top_medoid_pose_id", "")
            summary_row[f"{alias}_max_single_seed_fraction"] = parameter_row.get(
                "max_single_seed_fraction", ""
            )
            summary_row[f"{alias}_seed_artifact_flag"] = _normalize_bool(
                parameter_row.get("seed_artifact_flag", "")
            )
            summary_row[f"{alias}_cluster_assignments_tsv"] = parameter_row.get(
                "cluster_assignments_tsv", ""
            )
            summary_row[f"{alias}_medoid_manifest_tsv"] = parameter_row.get("medoid_manifest_tsv", "")

        summary_rows.append(summary_row)

    return fieldnames, summary_rows


def main() -> int:
    args = _parse_args()
    sensitivity_root = Path(args.sensitivity_root).resolve()
    output_path = Path(args.output_path).resolve()
    grid_results_path = sensitivity_root / "pilot_clustering_parameter_grid_results.tsv"

    if not grid_results_path.exists():
        raise FileNotFoundError(f"Missing grid results TSV: {grid_results_path}")

    profile_specs = _parse_profile_specs(args.profile)
    grid_rows = _read_tsv(grid_results_path)
    fieldnames, summary_rows = build_condition_summary_rows(grid_rows, profile_specs)
    _write_tsv(output_path, fieldnames, summary_rows)

    print(
        json.dumps(
            {
                "sensitivity_root": str(sensitivity_root),
                "output_path": str(output_path),
                "n_conditions": len(summary_rows),
                "profiles": [
                    {"alias": alias, "parameter_label": parameter_label}
                    for alias, parameter_label in profile_specs
                ],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
