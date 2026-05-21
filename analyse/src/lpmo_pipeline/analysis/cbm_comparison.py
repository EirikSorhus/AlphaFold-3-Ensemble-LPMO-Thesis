# src/lpmo_pipeline/analysis/cbm_comparison.py
"""
Responsibility: Statistical comparison of DEL A (catalytic domain only)
                vs DEL B (full-length with CBM).
Input:  Metrics from DEL A and DEL B for the same protein×ligand pairs
Output: Paired statistical tests, effect sizes, comparison summary

RQ3: Does CBM shift ligand poses toward CBM-assisted positioning?
"""
from __future__ import annotations

import csv
import json
import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

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
    "cbm_ligand_contact_fraction",
    "linker_ligand_contact_fraction",
    "bridge_fraction",
    "cbm_recruitment_score",
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
    "cbm_ligand_contact_fraction_full_length",
    "linker_ligand_contact_fraction_full_length",
    "bridge_fraction_full_length",
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


def _safe_divide(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator else 0.0


def _index_protein_metadata(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    indexed: dict[str, dict[str, Any]] = {}
    for row in rows:
        protein_id = str(row.get("protein_id", "") or row.get("uniprot_id", "")).strip()
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
        cbm_fraction = _as_float(condition.get("cbm_contact_fraction")) or 0.0
        linker_fraction = _as_float(condition.get("linker_contact_fraction")) or 0.0
        bridge_fraction = _as_float(condition.get("bridge_fraction"))
        if bridge_fraction is None:
            bridge_fraction = 0.0

        out.append(
            {
                "condition_id": condition_id,
                "protein_id": protein_id,
                "family_label": metadata.get("family", metadata.get("family_label", "")),
                "cbm_type": metadata.get("cbm_type", metadata.get("cbm_status", "")),
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
                "cbm_ligand_contact_fraction": cbm_fraction,
                "linker_ligand_contact_fraction": linker_fraction,
                "bridge_fraction": bridge_fraction,
                "cbm_recruitment_score": max(cbm_fraction, bridge_fraction),
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
                "cbm_ligand_contact_fraction_full_length": full_length.get(
                    "cbm_ligand_contact_fraction", ""
                ),
                "linker_ligand_contact_fraction_full_length": full_length.get(
                    "linker_ligand_contact_fraction", ""
                ),
                "bridge_fraction_full_length": full_length.get("bridge_fraction", ""),
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

        # --- Wilcoxon signed-rank test ---
        # from scipy.stats import wilcoxon
        # stat, p_value = wilcoxon(diffs, alternative='two-sided')
        p_value = 1.0  # PSEUDOCODE placeholder

        # --- Effect size (rank-biserial) ---
        # r = 1 - (2 * stat) / (n * (n + 1))
        effect_size = float(np.mean(diffs) / np.std(diffs)) if np.std(diffs) > 0 else 0.0

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
