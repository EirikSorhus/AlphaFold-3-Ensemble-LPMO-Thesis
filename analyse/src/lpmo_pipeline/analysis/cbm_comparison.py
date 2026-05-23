# src/lpmo_pipeline/analysis/cbm_comparison.py
"""
Responsibility: Statistical comparison of DEL A (catalytic/core domain only)
                vs DEL B (full-length with non-core module(s)).
Input:  Metrics from DEL A and DEL B for the same protein×ligand pairs
Output: Paired statistical tests, effect sizes, comparison summary

RQ3: Does adding non-core sequence/module context shift ligand poses?
"""
from __future__ import annotations

import csv
import json
import logging
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import binomtest, wilcoxon

logger = logging.getLogger(__name__)

CBM_CONSTRUCT_CONDITION_COLUMNS = [
    "condition_id",
    "protein_id",
    "family_label",
    "cbm_type",
    "construct_type",
    "substrate_class",
    "dp",
    "n_generated",
    "n_qc_pass",
    "qc_pass_fraction",
    "n_ifp_success",
    "ifp_success_fraction",
    "n_clusters",
    "any_valid_cluster",
    "no_valid_cluster_flag",
    "noise_fraction",
    "top_cluster_occupancy",
    "cluster_entropy",
    "geometry_plausible_fraction",
    "C1_compatible_fraction",
    "C4_compatible_fraction",
    "C4_minus_C1_geometry_bias",
    "non_core_ligand_contact_fraction",
    "bridge_fraction",
    "non_core_recruitment_score",
    "active_site_ifp_weighted_vector",
    "catalytic_surface_contact_fraction",
    "aromatic_contact_fraction",
    "polar_contact_fraction",
]

CBM_PAIRED_COMPARISON_COLUMNS = [
    "protein_id",
    "family_label",
    "cbm_type",
    "substrate_class",
    "dp",
    "domain_only_condition_id",
    "full_length_condition_id",
    "n_generated_domain_only",
    "n_generated_full_length",
    "n_qc_pass_domain_only",
    "n_qc_pass_full_length",
    "any_valid_cluster_domain_only",
    "any_valid_cluster_full_length",
    "non_core_ligand_contact_fraction_full_length",
    "bridge_fraction_full_length",
    "non_core_recruitment_score_full_length",
    "catalytic_domain_ifp_jaccard_distance",
    "delta_qc_pass_fraction",
    "delta_ifp_success_fraction",
    "delta_n_clusters",
    "delta_noise_fraction",
    "delta_top_cluster_occupancy",
    "delta_cluster_entropy",
    "delta_C1_compatible_fraction",
    "delta_C4_compatible_fraction",
    "delta_C4_minus_C1_geometry_bias",
    "delta_geometry_plausible_fraction",
    "delta_catalytic_surface_contact_fraction",
    "delta_aromatic_contact_fraction",
    "delta_polar_contact_fraction",
]

CBM_PRIMARY_METRIC_COLUMNS = [
    "metric_name",
    "source_column",
    "endpoint_type",
    "n_pairs",
    "n_nonzero_pairs",
    "sample_size_label",
    "statistical_test",
    "median",
    "mean",
    "std",
    "q1",
    "q3",
    "iqr",
    "minimum",
    "maximum",
    "n_positive_delta",
    "n_negative_delta",
    "n_zero_delta",
    "bootstrap_ci_low",
    "bootstrap_ci_high",
    "sign_test_p_value",
    "wilcoxon_p_value",
    "effect_size_rank_biserial",
    "bh_fdr_p_value",
    "significant_fdr_0_05",
    "interpretation_scope",
]

CBM_SECONDARY_METRIC_COLUMNS = [
    "metric_name",
    "source_column",
    "endpoint_type",
    "n_pairs",
    "median",
    "q1",
    "q3",
    "iqr",
    "minimum",
    "maximum",
    "n_positive_delta",
    "n_negative_delta",
    "n_zero_delta",
]

CBM_STRATIFIED_SUMMARY_COLUMNS = [
    "stratum_name",
    "stratum_value",
    "metric_name",
    "source_column",
    "n_pairs",
    "median",
    "q1",
    "q3",
    "iqr",
    "n_positive_delta",
    "n_negative_delta",
    "n_zero_delta",
]

CBM_REPRESENTATIVE_EXAMPLE_COLUMNS = [
    "example_type",
    "protein_id",
    "family_label",
    "cbm_type",
    "substrate_class",
    "dp",
    "domain_only_condition_id",
    "full_length_condition_id",
    "metric_name",
    "metric_value",
    "selection_reason",
]

CBM_FIGURE_MANIFEST_COLUMNS = [
    "figure_name",
    "description",
    "source_table",
    "required_columns",
    "plot_unit",
    "status",
]

PRIMARY_ENDPOINTS = [
    {
        "metric_name": "bridge_fraction",
        "source_column": "bridge_fraction_full_length",
        "endpoint_type": "full_length_only",
    },
    {
        "metric_name": "catalytic_domain_ifp_jaccard_distance",
        "source_column": "catalytic_domain_ifp_jaccard_distance",
        "endpoint_type": "paired_distance",
    },
    {
        "metric_name": "delta_C4_minus_C1_geometry_bias",
        "source_column": "delta_C4_minus_C1_geometry_bias",
        "endpoint_type": "paired_delta",
    },
    {
        "metric_name": "delta_qc_pass_fraction",
        "source_column": "delta_qc_pass_fraction",
        "endpoint_type": "paired_delta",
    },
    {
        "metric_name": "delta_cluster_entropy",
        "source_column": "delta_cluster_entropy",
        "endpoint_type": "paired_delta",
    },
]

SECONDARY_ENDPOINTS = [
    {"metric_name": "non_core_ligand_contact_fraction", "source_column": "non_core_ligand_contact_fraction_full_length", "endpoint_type": "full_length_only"},
    {"metric_name": "non_core_recruitment_score", "source_column": "non_core_recruitment_score_full_length", "endpoint_type": "full_length_only"},
    {"metric_name": "delta_C1_compatible_fraction", "source_column": "delta_C1_compatible_fraction", "endpoint_type": "paired_delta"},
    {"metric_name": "delta_C4_compatible_fraction", "source_column": "delta_C4_compatible_fraction", "endpoint_type": "paired_delta"},
    {"metric_name": "delta_geometry_plausible_fraction", "source_column": "delta_geometry_plausible_fraction", "endpoint_type": "paired_delta"},
    {"metric_name": "delta_ifp_success_fraction", "source_column": "delta_ifp_success_fraction", "endpoint_type": "paired_delta"},
    {"metric_name": "delta_noise_fraction", "source_column": "delta_noise_fraction", "endpoint_type": "paired_delta"},
    {"metric_name": "delta_top_cluster_occupancy", "source_column": "delta_top_cluster_occupancy", "endpoint_type": "paired_delta"},
    {"metric_name": "delta_n_clusters", "source_column": "delta_n_clusters", "endpoint_type": "paired_delta"},
    {"metric_name": "delta_catalytic_surface_contact_fraction", "source_column": "delta_catalytic_surface_contact_fraction", "endpoint_type": "paired_delta"},
    {"metric_name": "delta_aromatic_contact_fraction", "source_column": "delta_aromatic_contact_fraction", "endpoint_type": "paired_delta"},
    {"metric_name": "delta_polar_contact_fraction", "source_column": "delta_polar_contact_fraction", "endpoint_type": "paired_delta"},
]

CBM_FIGURE_MANIFEST_ROWS = [
    {
        "figure_name": "paired_delta_plot_primary_metrics.png",
        "description": "One point per matched pair for primary full_length - domain_only delta metrics.",
        "source_table": "cbm_paired_comparison_table.tsv; cbm_primary_metric_summary.tsv",
        "required_columns": "delta_C4_minus_C1_geometry_bias,delta_qc_pass_fraction,delta_cluster_entropy",
        "plot_unit": "matched protein_id x substrate_class x dp pair",
        "status": "planned_downstream_plot",
    },
    {
        "figure_name": "domain_only_vs_full_length_paired_lines.png",
        "description": "Paired line plots for selected construct-level metrics before and after adding non-core sequence/module context.",
        "source_table": "cbm_construct_condition_summary.tsv",
        "required_columns": "protein_id,substrate_class,dp,construct_type,qc_pass_fraction,cluster_entropy,geometry_plausible_fraction",
        "plot_unit": "matched construct condition",
        "status": "planned_downstream_plot",
    },
    {
        "figure_name": "cbm_bridge_fraction_by_substrate_dp.png",
        "description": "Full-length catalytic/core-to-non-core bridge fraction by substrate_class and dp with individual points and medians.",
        "source_table": "cbm_paired_comparison_table.tsv",
        "required_columns": "substrate_class,dp,bridge_fraction_full_length",
        "plot_unit": "matched pair full-length condition",
        "status": "planned_downstream_plot",
    },
    {
        "figure_name": "active_site_ifp_change_heatmap.png",
        "description": "Heatmap of catalytic-domain IFP Jaccard distance per protein and ligand condition when vectors are available.",
        "source_table": "cbm_paired_comparison_table.tsv",
        "required_columns": "protein_id,substrate_class,dp,catalytic_domain_ifp_jaccard_distance",
        "plot_unit": "matched pair",
        "status": "planned_downstream_plot",
    },
    {
        "figure_name": "geometry_bias_delta_heatmap.png",
        "description": "Heatmap of delta_C4_minus_C1_geometry_bias per protein and ligand condition.",
        "source_table": "cbm_paired_comparison_table.tsv",
        "required_columns": "protein_id,substrate_class,dp,delta_C4_minus_C1_geometry_bias",
        "plot_unit": "matched pair",
        "status": "planned_downstream_plot",
    },
]


@dataclass
class PairedComparison:
    """Result of comparing DEL A vs DEL B for one metric."""

    metric_name: str
    n_pairs: int = 0
    del_a_mean: float = 0.0
    del_b_mean: float = 0.0
    difference_mean: float = 0.0
    difference_std: float = 0.0
    p_value: float = 1.0  # Wilcoxon signed-rank
    effect_size: float = 0.0  # Rank-biserial correlation or similar
    significant: bool = False


@dataclass
class CBMComparisonReport:
    """Full DEL A vs DEL B comparison report."""

    protein_id: str = ""
    ligand_id: str = ""
    n_proteins_compared: int = 0
    comparisons: list[PairedComparison] = field(default_factory=list)


@dataclass
class CBMPairedAnalysisResult:
    """Paths written by the CBM paired postprocess."""

    output_dir: Path
    summary_path: Path
    table_paths: dict[str, Path] = field(default_factory=dict)


def _as_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(number):
        return None
    return number


def _as_int(value: Any) -> int:
    number = _as_float(value)
    return int(number) if number is not None else 0


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


def _first_present(row: dict[str, Any], *keys: str) -> str:
    for key in keys:
        value = row.get(key)
        if value not in (None, ""):
            text = str(value).strip()
            if text:
                return text
    return ""


def _safe_divide(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator else 0.0


def _numeric_values(rows: list[dict[str, Any]], column_name: str) -> list[float]:
    values = []
    for row in rows:
        value = _as_float(row.get(column_name))
        if value is not None:
            values.append(value)
    return values


def _sample_size_label(n_pairs: int) -> str:
    if n_pairs < 5:
        return "extremely_small_n"
    if n_pairs < 8:
        return "very_small_n"
    if n_pairs < 12:
        return "small_n"
    return "adequate_for_simple_nonparametric_paired_summary"


def _summary_stats(values: list[float]) -> dict[str, float | int | str]:
    if not values:
        return {
            "n_pairs": 0,
            "median": "",
            "mean": "",
            "std": "",
            "q1": "",
            "q3": "",
            "iqr": "",
            "minimum": "",
            "maximum": "",
            "n_positive_delta": 0,
            "n_negative_delta": 0,
            "n_zero_delta": 0,
        }
    array = np.asarray(values, dtype=float)
    q1 = float(np.percentile(array, 25))
    q3 = float(np.percentile(array, 75))
    return {
        "n_pairs": len(values),
        "median": float(np.median(array)),
        "mean": float(np.mean(array)),
        "std": float(np.std(array, ddof=1)) if len(values) > 1 else 0.0,
        "q1": q1,
        "q3": q3,
        "iqr": q3 - q1,
        "minimum": float(np.min(array)),
        "maximum": float(np.max(array)),
        "n_positive_delta": int(np.sum(array > 0)),
        "n_negative_delta": int(np.sum(array < 0)),
        "n_zero_delta": int(np.sum(array == 0)),
    }


def _bootstrap_median_ci(
    values: list[float],
    *,
    n_resamples: int = 10000,
    random_state: int = 42,
) -> tuple[float | str, float | str]:
    if len(values) < 8:
        return "", ""
    array = np.asarray(values, dtype=float)
    rng = np.random.default_rng(random_state)
    sample_indices = rng.integers(0, len(array), size=(n_resamples, len(array)))
    medians = np.median(array[sample_indices], axis=1)
    return float(np.percentile(medians, 2.5)), float(np.percentile(medians, 97.5))


def _rank_biserial_from_wilcoxon(values: list[float]) -> float | str:
    nonzero_abs = [abs(value) for value in values if value != 0.0]
    if not nonzero_abs:
        return ""
    ranks = _average_ranks(nonzero_abs)
    positive_rank_sum = 0.0
    negative_rank_sum = 0.0
    rank_index = 0
    for value in values:
        if value == 0.0:
            continue
        if value > 0:
            positive_rank_sum += ranks[rank_index]
        else:
            negative_rank_sum += ranks[rank_index]
        rank_index += 1
    denominator = positive_rank_sum + negative_rank_sum
    if denominator == 0.0:
        return ""
    return float((positive_rank_sum - negative_rank_sum) / denominator)


def _average_ranks(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0 for _ in values]
    cursor = 0
    while cursor < len(indexed):
        end = cursor + 1
        while end < len(indexed) and indexed[end][1] == indexed[cursor][1]:
            end += 1
        average_rank = (cursor + 1 + end) / 2.0
        for original_index, _ in indexed[cursor:end]:
            ranks[original_index] = average_rank
        cursor = end
    return ranks


def _benjamini_hochberg(p_values: list[float | str]) -> list[float | str]:
    indexed = [
        (index, float(value))
        for index, value in enumerate(p_values)
        if value not in ("", None) and math.isfinite(float(value))
    ]
    adjusted: list[float | str] = ["" for _ in p_values]
    if not indexed:
        return adjusted
    sorted_values = sorted(indexed, key=lambda item: item[1])
    n_tests = len(sorted_values)
    running_min = 1.0
    for rank_from_end, (original_index, p_value) in enumerate(reversed(sorted_values), start=1):
        rank = n_tests - rank_from_end + 1
        running_min = min(running_min, p_value * n_tests / rank)
        adjusted[original_index] = float(min(running_min, 1.0))
    return adjusted


def _parse_weighted_vector(value: Any) -> list[float] | None:
    if value in (None, ""):
        return None
    if isinstance(value, (list, tuple)):
        raw_values = value
    else:
        text = str(value).strip()
        try:
            decoded = json.loads(text)
        except (TypeError, ValueError, json.JSONDecodeError):
            decoded = None
        if isinstance(decoded, dict):
            raw_values = [decoded[key] for key in sorted(decoded)]
        elif isinstance(decoded, list):
            raw_values = decoded
        else:
            raw_values = [part for part in text.replace(",", " ").split() if part]
    parsed: list[float] = []
    for raw_value in raw_values:
        numeric = _as_float(raw_value)
        if numeric is None:
            return None
        parsed.append(numeric)
    return parsed or None


def _jaccard_distance_from_vectors(
    domain_only: dict[str, Any],
    full_length: dict[str, Any],
    *,
    threshold: float = 0.10,
) -> float | str:
    domain_vector = _parse_weighted_vector(domain_only.get("active_site_ifp_weighted_vector"))
    full_vector = _parse_weighted_vector(full_length.get("active_site_ifp_weighted_vector"))
    if not domain_vector or not full_vector or len(domain_vector) != len(full_vector):
        return ""
    domain_bits = [value >= threshold for value in domain_vector]
    full_bits = [value >= threshold for value in full_vector]
    union = sum(domain or full for domain, full in zip(domain_bits, full_bits, strict=True))
    if union == 0:
        return 0.0
    intersection = sum(domain and full for domain, full in zip(domain_bits, full_bits, strict=True))
    return float(1.0 - (intersection / union))


def _index_protein_metadata(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    indexed: dict[str, dict[str, Any]] = {}
    for row in rows:
        protein_id = _first_present(row, "protein_id", "uniprot_id", "UniProt_ID")
        if protein_id and protein_id not in indexed:
            indexed[protein_id] = row
    return indexed


def _write_tsv(path: Path, columns: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def read_tsv(path: Path | None) -> list[dict[str, str]]:
    if path is None or not Path(path).exists():
        return []
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _cluster_type_fractions(cluster_rows: list[dict[str, Any]]) -> dict[str, float]:
    total_occupancy = sum(_as_float(row.get("occupancy")) or 0.0 for row in cluster_rows)
    if not total_occupancy:
        return {"C1_compatible_fraction": 0.0, "C4_compatible_fraction": 0.0}
    c1 = sum(
        (_as_float(row.get("occupancy")) or 0.0)
        for row in cluster_rows
        if str(row.get("cluster_type", "")) == "C1_compatible"
    )
    c4 = sum(
        (_as_float(row.get("occupancy")) or 0.0)
        for row in cluster_rows
        if str(row.get("cluster_type", "")) == "C4_compatible"
    )
    return {
        "C1_compatible_fraction": c1 / total_occupancy,
        "C4_compatible_fraction": c4 / total_occupancy,
    }


def build_cbm_construct_condition_summary_rows(
    condition_rows: list[dict[str, Any]],
    *,
    cluster_rows: list[dict[str, Any]] | None = None,
    protein_metadata_rows: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """Build one CBM side-analysis row per construct condition."""

    clusters_by_condition: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in cluster_rows or []:
        clusters_by_condition[str(row.get("condition_id", ""))].append(row)
    metadata_by_protein = _index_protein_metadata(protein_metadata_rows or [])

    out: list[dict[str, Any]] = []
    for condition in sorted(condition_rows, key=lambda row: str(row.get("condition_id", ""))):
        condition_id = str(condition.get("condition_id", ""))
        protein_id = str(condition.get("protein_id", ""))
        metadata = metadata_by_protein.get(protein_id, {})
        cluster_type_fractions = _cluster_type_fractions(clusters_by_condition.get(condition_id, []))

        n_generated = _as_int(condition.get("n_generated"))
        n_qc_pass = _as_int(condition.get("n_stage1_pass"))
        n_ifp_success = _as_int(condition.get("n_ifp_success"))
        n_clusters = _as_int(condition.get("n_clusters"))
        c1_fraction = cluster_type_fractions["C1_compatible_fraction"]
        c4_fraction = cluster_type_fractions["C4_compatible_fraction"]
        geometry_plausible = max(
            _as_float(condition.get("occupancy_weighted_c1_plausible_fraction")) or 0.0,
            _as_float(condition.get("occupancy_weighted_c4_plausible_fraction")) or 0.0,
        )
        non_core_fraction = _as_float(condition.get("non_core_contact_fraction"))
        if non_core_fraction is None:
            cbm_fraction = _as_float(condition.get("cbm_contact_fraction")) or 0.0
            linker_fraction = _as_float(condition.get("linker_contact_fraction")) or 0.0
            non_core_fraction = cbm_fraction + linker_fraction
        bridge_fraction = _as_float(condition.get("bridge_fraction"))
        if bridge_fraction is None:
            bridge_fraction = 0.0

        out.append(
            {
                "condition_id": condition_id,
                "protein_id": protein_id,
                "family_label": _first_present(
                    metadata,
                    "family",
                    "family_label",
                    "protein_family",
                    "aa_family",
                    "CAZy_family",
                ),
                "cbm_type": _first_present(
                    metadata,
                    "cbm_type",
                    "cbm_status",
                    "Binding_Modules",
                ),
                "construct_type": condition.get("construct_type", ""),
                "substrate_class": condition.get("substrate_class", ""),
                "dp": condition.get("dp", ""),
                "n_generated": n_generated,
                "n_qc_pass": n_qc_pass,
                "qc_pass_fraction": _safe_divide(n_qc_pass, n_generated),
                "n_ifp_success": n_ifp_success,
                "ifp_success_fraction": _safe_divide(n_ifp_success, n_qc_pass),
                "n_clusters": n_clusters,
                "any_valid_cluster": _as_bool(condition.get("any_valid_cluster")) or n_clusters > 0,
                "no_valid_cluster_flag": n_clusters == 0,
                "noise_fraction": condition.get("noise_fraction", ""),
                "top_cluster_occupancy": condition.get("top_cluster_occupancy", ""),
                "cluster_entropy": condition.get("cluster_entropy", ""),
                "geometry_plausible_fraction": geometry_plausible,
                "C1_compatible_fraction": c1_fraction,
                "C4_compatible_fraction": c4_fraction,
                "C4_minus_C1_geometry_bias": c4_fraction - c1_fraction,
                "non_core_ligand_contact_fraction": non_core_fraction,
                "bridge_fraction": bridge_fraction,
                "non_core_recruitment_score": max(non_core_fraction, bridge_fraction),
                "active_site_ifp_weighted_vector": condition.get(
                    "active_site_ifp_weighted_vector",
                    condition.get("catalytic_domain_ifp_weighted_vector", ""),
                ),
                "catalytic_surface_contact_fraction": condition.get(
                    "catalytic_surface_contact_fraction", ""
                ),
                "aromatic_contact_fraction": condition.get("aromatic_contact_fraction", ""),
                "polar_contact_fraction": condition.get("polar_contact_fraction", ""),
            }
        )
    return out


def _pair_key(row: dict[str, Any]) -> tuple[str, str, str]:
    return (
        str(row.get("protein_id", "")),
        str(row.get("substrate_class", "")),
        str(row.get("dp", "")),
    )


def _delta(full_length: dict[str, Any], domain_only: dict[str, Any], field: str) -> float | str:
    full_value = _as_float(full_length.get(field))
    domain_value = _as_float(domain_only.get(field))
    if full_value is None or domain_value is None:
        return ""
    return full_value - domain_value


def build_cbm_paired_comparison_rows(
    construct_condition_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Build full-length versus domain-only condition pairs."""

    by_key: dict[tuple[str, str, str], dict[str, dict[str, Any]]] = defaultdict(dict)
    for row in construct_condition_rows:
        construct_type = str(row.get("construct_type", ""))
        if construct_type in {"domain_only", "catalytic_domain"}:
            by_key[_pair_key(row)]["domain_only"] = row
        elif construct_type == "full_length":
            by_key[_pair_key(row)]["full_length"] = row

    out: list[dict[str, Any]] = []
    for key in sorted(by_key):
        pair = by_key[key]
        if "domain_only" not in pair or "full_length" not in pair:
            continue
        domain_only = pair["domain_only"]
        full_length = pair["full_length"]
        protein_id, substrate_class, dp = key
        out.append(
            {
                "protein_id": protein_id,
                "family_label": full_length.get("family_label") or domain_only.get("family_label", ""),
                "cbm_type": full_length.get("cbm_type") or domain_only.get("cbm_type", ""),
                "substrate_class": substrate_class,
                "dp": dp,
                "domain_only_condition_id": domain_only.get("condition_id", ""),
                "full_length_condition_id": full_length.get("condition_id", ""),
                "n_generated_domain_only": domain_only.get("n_generated", ""),
                "n_generated_full_length": full_length.get("n_generated", ""),
                "n_qc_pass_domain_only": domain_only.get("n_qc_pass", ""),
                "n_qc_pass_full_length": full_length.get("n_qc_pass", ""),
                "any_valid_cluster_domain_only": domain_only.get("any_valid_cluster", ""),
                "any_valid_cluster_full_length": full_length.get("any_valid_cluster", ""),
                "non_core_ligand_contact_fraction_full_length": full_length.get(
                    "non_core_ligand_contact_fraction", ""
                ),
                "bridge_fraction_full_length": full_length.get("bridge_fraction", ""),
                "non_core_recruitment_score_full_length": full_length.get(
                    "non_core_recruitment_score", ""
                ),
                "catalytic_domain_ifp_jaccard_distance": _jaccard_distance_from_vectors(
                    domain_only,
                    full_length,
                ),
                "delta_qc_pass_fraction": _delta(full_length, domain_only, "qc_pass_fraction"),
                "delta_ifp_success_fraction": _delta(full_length, domain_only, "ifp_success_fraction"),
                "delta_n_clusters": _delta(full_length, domain_only, "n_clusters"),
                "delta_noise_fraction": _delta(full_length, domain_only, "noise_fraction"),
                "delta_top_cluster_occupancy": _delta(
                    full_length, domain_only, "top_cluster_occupancy"
                ),
                "delta_cluster_entropy": _delta(full_length, domain_only, "cluster_entropy"),
                "delta_C1_compatible_fraction": _delta(
                    full_length, domain_only, "C1_compatible_fraction"
                ),
                "delta_C4_compatible_fraction": _delta(
                    full_length, domain_only, "C4_compatible_fraction"
                ),
                "delta_C4_minus_C1_geometry_bias": _delta(
                    full_length, domain_only, "C4_minus_C1_geometry_bias"
                ),
                "delta_geometry_plausible_fraction": _delta(
                    full_length, domain_only, "geometry_plausible_fraction"
                ),
                "delta_catalytic_surface_contact_fraction": _delta(
                    full_length, domain_only, "catalytic_surface_contact_fraction"
                ),
                "delta_aromatic_contact_fraction": _delta(
                    full_length, domain_only, "aromatic_contact_fraction"
                ),
                "delta_polar_contact_fraction": _delta(full_length, domain_only, "polar_contact_fraction"),
            }
        )
    return out


def write_cbm_construct_condition_summary(rows: list[dict[str, Any]], output_path: Path) -> None:
    _write_tsv(output_path, CBM_CONSTRUCT_CONDITION_COLUMNS, rows)


def write_cbm_paired_comparison_table(rows: list[dict[str, Any]], output_path: Path) -> None:
    _write_tsv(output_path, CBM_PAIRED_COMPARISON_COLUMNS, rows)


def build_cbm_primary_metric_summary_rows(
    paired_rows: list[dict[str, Any]],
    *,
    alpha: float = 0.05,
    random_state: int = 42,
) -> list[dict[str, Any]]:
    """Summarize primary CBM endpoints using the sample-size rules in the plan."""

    summary_rows: list[dict[str, Any]] = []
    wilcoxon_p_values: list[float | str] = []
    for endpoint in PRIMARY_ENDPOINTS:
        values = _numeric_values(paired_rows, endpoint["source_column"])
        stats = _summary_stats(values)
        n_pairs = int(stats["n_pairs"])
        n_nonzero = int(stats["n_positive_delta"]) + int(stats["n_negative_delta"])
        ci_low, ci_high = _bootstrap_median_ci(values, random_state=random_state)
        sign_p_value: float | str = ""
        wilcoxon_p_value: float | str = ""
        effect_size: float | str = ""
        statistical_test = "none"

        if endpoint["endpoint_type"] == "full_length_only":
            statistical_test = "descriptive_full_length_only"
        elif n_pairs >= 8:
            sign_p_value = float(
                binomtest(
                    int(stats["n_positive_delta"]),
                    n=int(stats["n_positive_delta"]) + int(stats["n_negative_delta"]),
                    p=0.5,
                    alternative="two-sided",
                ).pvalue
            ) if n_nonzero else ""
            statistical_test = "sign_test"
        if endpoint["endpoint_type"] != "full_length_only" and n_nonzero >= 8:
            try:
                wilcoxon_p_value = float(
                    wilcoxon(values, zero_method="wilcox", alternative="two-sided").pvalue
                )
                effect_size = _rank_biserial_from_wilcoxon(values)
                statistical_test = "paired_wilcoxon_signed_rank"
            except ValueError:
                wilcoxon_p_value = ""
                effect_size = ""

        row = {
            "metric_name": endpoint["metric_name"],
            "source_column": endpoint["source_column"],
            "endpoint_type": endpoint["endpoint_type"],
            "n_pairs": n_pairs,
            "n_nonzero_pairs": n_nonzero,
            "sample_size_label": _sample_size_label(n_pairs),
            "statistical_test": statistical_test,
            "median": stats["median"],
            "mean": stats["mean"],
            "std": stats["std"],
            "q1": stats["q1"],
            "q3": stats["q3"],
            "iqr": stats["iqr"],
            "minimum": stats["minimum"],
            "maximum": stats["maximum"],
            "n_positive_delta": stats["n_positive_delta"],
            "n_negative_delta": stats["n_negative_delta"],
            "n_zero_delta": stats["n_zero_delta"],
            "bootstrap_ci_low": ci_low,
            "bootstrap_ci_high": ci_high,
            "sign_test_p_value": sign_p_value,
            "wilcoxon_p_value": wilcoxon_p_value,
            "effect_size_rank_biserial": effect_size,
            "bh_fdr_p_value": "",
            "significant_fdr_0_05": "",
            "interpretation_scope": "exploratory_af3_construct_comparison",
        }
        summary_rows.append(row)
        wilcoxon_p_values.append(wilcoxon_p_value)

    adjusted_p_values = _benjamini_hochberg(wilcoxon_p_values)
    for row, adjusted_p_value in zip(summary_rows, adjusted_p_values, strict=True):
        row["bh_fdr_p_value"] = adjusted_p_value
        row["significant_fdr_0_05"] = (
            bool(adjusted_p_value <= alpha) if adjusted_p_value != "" else ""
        )
    return summary_rows


def build_cbm_secondary_descriptive_summary_rows(
    paired_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Build descriptive summaries for secondary CBM endpoints."""

    out: list[dict[str, Any]] = []
    for endpoint in SECONDARY_ENDPOINTS:
        stats = _summary_stats(_numeric_values(paired_rows, endpoint["source_column"]))
        out.append(
            {
                "metric_name": endpoint["metric_name"],
                "source_column": endpoint["source_column"],
                "endpoint_type": endpoint["endpoint_type"],
                "n_pairs": stats["n_pairs"],
                "median": stats["median"],
                "q1": stats["q1"],
                "q3": stats["q3"],
                "iqr": stats["iqr"],
                "minimum": stats["minimum"],
                "maximum": stats["maximum"],
                "n_positive_delta": stats["n_positive_delta"],
                "n_negative_delta": stats["n_negative_delta"],
                "n_zero_delta": stats["n_zero_delta"],
            }
        )
    return out


def build_cbm_stratified_summary_rows(
    paired_rows: list[dict[str, Any]],
    *,
    stratum_name: str,
    endpoints: list[dict[str, str]] | None = None,
) -> list[dict[str, Any]]:
    """Build descriptive CBM summaries within one stratification variable."""

    if endpoints is None:
        endpoints = PRIMARY_ENDPOINTS
    rows_by_value: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in paired_rows:
        rows_by_value[str(row.get(stratum_name, ""))].append(row)

    out: list[dict[str, Any]] = []
    for stratum_value in sorted(rows_by_value):
        stratum_rows = rows_by_value[stratum_value]
        for endpoint in endpoints:
            stats = _summary_stats(_numeric_values(stratum_rows, endpoint["source_column"]))
            out.append(
                {
                    "stratum_name": stratum_name,
                    "stratum_value": stratum_value,
                    "metric_name": endpoint["metric_name"],
                    "source_column": endpoint["source_column"],
                    "n_pairs": stats["n_pairs"],
                    "median": stats["median"],
                    "q1": stats["q1"],
                    "q3": stats["q3"],
                    "iqr": stats["iqr"],
                    "n_positive_delta": stats["n_positive_delta"],
                    "n_negative_delta": stats["n_negative_delta"],
                    "n_zero_delta": stats["n_zero_delta"],
                }
            )
    return out


def build_cbm_representative_example_rows(
    paired_rows: list[dict[str, Any]],
    *,
    top_n: int = 3,
) -> list[dict[str, Any]]:
    """Select deterministic rows worth inspecting structurally after the table analysis."""

    selections = [
        (
            "full_length_medoid_with_ligand_bridging_catalytic_domain_and_cbm",
            "bridge_fraction_full_length",
            "largest full-length catalytic/core plus non-core bridge fraction",
            False,
        ),
        (
            "example_with_large_active_site_ifp_change",
            "catalytic_domain_ifp_jaccard_distance",
            "largest catalytic-domain IFP Jaccard distance between constructs",
            False,
        ),
        (
            "example_with_large_C1_to_C4_geometry_bias_shift",
            "delta_C4_minus_C1_geometry_bias",
            "largest absolute shift in C4-minus-C1 geometry bias",
            True,
        ),
    ]
    out: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str, str]] = set()
    for example_type, column_name, reason, use_abs in selections:
        candidate_rows = []
        for row in paired_rows:
            value = _as_float(row.get(column_name))
            if value is None:
                continue
            candidate_rows.append((abs(value) if use_abs else value, value, row))
        for _, metric_value, row in sorted(
            candidate_rows,
            key=lambda item: (
                -item[0],
                str(item[2].get("protein_id", "")),
                str(item[2].get("substrate_class", "")),
                str(item[2].get("dp", "")),
            ),
        )[:top_n]:
            key = (
                example_type,
                str(row.get("protein_id", "")),
                str(row.get("substrate_class", "")),
                str(row.get("dp", "")),
            )
            if key in seen:
                continue
            seen.add(key)
            out.append(
                {
                    "example_type": example_type,
                    "protein_id": row.get("protein_id", ""),
                    "family_label": row.get("family_label", ""),
                    "cbm_type": row.get("cbm_type", ""),
                    "substrate_class": row.get("substrate_class", ""),
                    "dp": row.get("dp", ""),
                    "domain_only_condition_id": row.get("domain_only_condition_id", ""),
                    "full_length_condition_id": row.get("full_length_condition_id", ""),
                    "metric_name": column_name,
                    "metric_value": metric_value,
                    "selection_reason": reason,
                }
            )
    return out


def build_cbm_figure_manifest_rows() -> list[dict[str, Any]]:
    """Return the documented downstream CBM figure contract."""

    return [dict(row) for row in CBM_FIGURE_MANIFEST_ROWS]


def write_cbm_primary_metric_summary(rows: list[dict[str, Any]], output_path: Path) -> None:
    _write_tsv(output_path, CBM_PRIMARY_METRIC_COLUMNS, rows)


def write_cbm_secondary_descriptive_summary(rows: list[dict[str, Any]], output_path: Path) -> None:
    _write_tsv(output_path, CBM_SECONDARY_METRIC_COLUMNS, rows)


def write_cbm_stratified_summary(rows: list[dict[str, Any]], output_path: Path) -> None:
    _write_tsv(output_path, CBM_STRATIFIED_SUMMARY_COLUMNS, rows)


def write_cbm_representative_examples(rows: list[dict[str, Any]], output_path: Path) -> None:
    _write_tsv(output_path, CBM_REPRESENTATIVE_EXAMPLE_COLUMNS, rows)


def write_cbm_figure_manifest(rows: list[dict[str, Any]], output_path: Path) -> None:
    _write_tsv(output_path, CBM_FIGURE_MANIFEST_COLUMNS, rows)


def run_cbm_paired_analysis(
    *,
    condition_table_path: Path,
    output_dir: Path,
    cluster_table_path: Path | None = None,
    protein_metadata_path: Path | None = None,
    random_state: int = 42,
) -> CBMPairedAnalysisResult:
    """Run the condition-level CBM paired side analysis from production TSVs."""

    cbm_root = output_dir / "15_cbm_paired_analysis"
    summary_path = cbm_root / "cbm_paired_analysis_summary.json"
    result = CBMPairedAnalysisResult(output_dir=cbm_root, summary_path=summary_path)

    condition_rows = read_tsv(condition_table_path)
    cluster_rows = read_tsv(cluster_table_path)
    protein_metadata_rows = read_tsv(protein_metadata_path)

    construct_rows = build_cbm_construct_condition_summary_rows(
        condition_rows,
        cluster_rows=cluster_rows,
        protein_metadata_rows=protein_metadata_rows,
    )
    paired_rows = build_cbm_paired_comparison_rows(construct_rows)
    primary_rows = build_cbm_primary_metric_summary_rows(
        paired_rows,
        random_state=random_state,
    )
    secondary_rows = build_cbm_secondary_descriptive_summary_rows(paired_rows)
    stratified_by_substrate_rows = build_cbm_stratified_summary_rows(
        paired_rows,
        stratum_name="substrate_class",
    )
    stratified_by_dp_rows = build_cbm_stratified_summary_rows(paired_rows, stratum_name="dp")
    stratified_by_cbm_type_rows = build_cbm_stratified_summary_rows(
        paired_rows,
        stratum_name="cbm_type",
    )
    representative_rows = build_cbm_representative_example_rows(paired_rows)
    figure_rows = build_cbm_figure_manifest_rows()

    tables = {
        "construct_condition_summary": (
            cbm_root / "cbm_construct_condition_summary.tsv",
            write_cbm_construct_condition_summary,
            construct_rows,
        ),
        "paired_comparison": (
            cbm_root / "cbm_paired_comparison_table.tsv",
            write_cbm_paired_comparison_table,
            paired_rows,
        ),
        "primary_metric_summary": (
            cbm_root / "cbm_primary_metric_summary.tsv",
            write_cbm_primary_metric_summary,
            primary_rows,
        ),
        "secondary_descriptive_summary": (
            cbm_root / "cbm_secondary_descriptive_summary.tsv",
            write_cbm_secondary_descriptive_summary,
            secondary_rows,
        ),
        "stratified_summary_by_substrate": (
            cbm_root / "cbm_stratified_summary_by_substrate.tsv",
            write_cbm_stratified_summary,
            stratified_by_substrate_rows,
        ),
        "stratified_summary_by_dp": (
            cbm_root / "cbm_stratified_summary_by_dp.tsv",
            write_cbm_stratified_summary,
            stratified_by_dp_rows,
        ),
        "stratified_summary_by_cbm_type": (
            cbm_root / "cbm_stratified_summary_by_cbm_type.tsv",
            write_cbm_stratified_summary,
            stratified_by_cbm_type_rows,
        ),
        "representative_examples": (
            cbm_root / "cbm_representative_examples.tsv",
            write_cbm_representative_examples,
            representative_rows,
        ),
        "figure_manifest": (
            cbm_root / "cbm_figure_manifest.tsv",
            write_cbm_figure_manifest,
            figure_rows,
        ),
    }
    for table_name, (path, writer, rows) in tables.items():
        writer(rows, path)
        result.table_paths[table_name] = path

    summary_data = {
        "condition_table": str(condition_table_path),
        "cluster_table": str(cluster_table_path) if cluster_table_path else "",
        "protein_metadata": str(protein_metadata_path) if protein_metadata_path else "",
        "n_construct_condition_rows": len(construct_rows),
        "n_paired_rows": len(paired_rows),
        "n_primary_metrics": len(primary_rows),
        "n_secondary_metrics": len(secondary_rows),
        "tables": {name: str(path) for name, path in sorted(result.table_paths.items())},
        "figures_to_generate": [row["figure_name"] for row in figure_rows],
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w") as handle:
        json.dump(summary_data, handle, indent=2)

    return result


def compare_del_a_vs_del_b(
    del_a_metrics: list[dict[str, Any]],
    del_b_metrics: list[dict[str, Any]],
    protein_id: str = "",
    ligand_id: str = "",
    metrics_to_compare: list[str] | None = None,
    alpha: float = 0.05,
) -> CBMComparisonReport:
    """Run paired comparison between DEL A and DEL B.

    Args:
        del_a_metrics: List of metric dicts from DEL A.
        del_b_metrics: List of metric dicts from DEL B (paired by run).
        protein_id: For metadata.
        ligand_id: For metadata.
        metrics_to_compare: Which metrics to test. Defaults to:
            [min_cu_c1, min_cu_c4, n_ifp_contacts_lpmo, cbm_ligand_min_dist,
             cluster_occupancy_best]
        alpha: Significance level.

    Returns:
        CBMComparisonReport with paired tests.
    """
    if metrics_to_compare is None:
        metrics_to_compare = [
            "min_cu_c1",
            "min_cu_c4",
            "n_ifp_contacts_lpmo",
            "cluster_occupancy_best",
        ]

    report = CBMComparisonReport(
        protein_id=protein_id,
        ligand_id=ligand_id,
        n_proteins_compared=min(len(del_a_metrics), len(del_b_metrics)),
    )

    for metric_name in metrics_to_compare:
        # Extract paired values
        pairs_a = [m.get(metric_name) for m in del_a_metrics]
        pairs_b = [m.get(metric_name) for m in del_b_metrics]

        # Filter None/missing
        valid_pairs = [
            (a, b) for a, b in zip(pairs_a, pairs_b)
            if a is not None and b is not None
        ]

        if len(valid_pairs) < 3:
            logger.warning(
                "Not enough valid pairs for %s: %d (need ≥3)",
                metric_name, len(valid_pairs),
            )
            continue

        a_vals = np.array([p[0] for p in valid_pairs])
        b_vals = np.array([p[1] for p in valid_pairs])
        diffs = b_vals - a_vals

        nonzero_diffs = [float(diff) for diff in diffs if diff != 0.0]
        if len(nonzero_diffs) >= 3:
            p_value = float(
                wilcoxon(diffs, zero_method="wilcox", alternative="two-sided").pvalue
            )
            rank_biserial = _rank_biserial_from_wilcoxon([float(diff) for diff in diffs])
            effect_size = float(rank_biserial) if rank_biserial != "" else 0.0
        else:
            p_value = 1.0
            effect_size = 0.0

        comparison = PairedComparison(
            metric_name=metric_name,
            n_pairs=len(valid_pairs),
            del_a_mean=float(np.mean(a_vals)),
            del_b_mean=float(np.mean(b_vals)),
            difference_mean=float(np.mean(diffs)),
            difference_std=float(np.std(diffs)),
            p_value=p_value,
            effect_size=effect_size,
            significant=p_value < alpha,
        )
        report.comparisons.append(comparison)

    logger.info(
        "CBM comparison %s×%s: %d metrics tested, %d significant",
        protein_id, ligand_id,
        len(report.comparisons),
        sum(1 for c in report.comparisons if c.significant),
    )
    return report


def write_cbm_comparison_report(
    report: CBMComparisonReport, output_path: Path,
) -> None:
    """Write comparison report to JSON."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "protein_id": report.protein_id,
        "ligand_id": report.ligand_id,
        "n_proteins_compared": report.n_proteins_compared,
        "comparisons": [
            {
                "metric": c.metric_name,
                "n_pairs": c.n_pairs,
                "del_a_mean": c.del_a_mean,
                "del_b_mean": c.del_b_mean,
                "diff_mean": c.difference_mean,
                "diff_std": c.difference_std,
                "p_value": c.p_value,
                "effect_size": c.effect_size,
                "significant": c.significant,
            }
            for c in report.comparisons
        ],
    }
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
