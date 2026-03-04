# src/lpmo_pipeline/analysis/cbm_comparison.py
"""
Responsibility: Statistical comparison of DEL A (catalytic domain only)
                vs DEL B (full-length with CBM).
Input:  Metrics from DEL A and DEL B for the same protein×ligand pairs
Output: Paired statistical tests, effect sizes, comparison summary

RQ3: Does CBM shift ligand poses toward CBM-assisted positioning?
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)


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
