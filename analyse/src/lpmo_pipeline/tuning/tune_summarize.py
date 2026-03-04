# src/lpmo_pipeline/tuning/tune_summarize.py
"""
Responsibility: Post-hoc tuning summary and visualization data.
Input:  tuning_summary.json from tune_orchestrator
Output: Formatted tuning report data for report/build_report_html.py
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def format_tuning_summary_for_report(
    tuning_summary_path: Path,
) -> dict[str, Any]:
    """Load and format tuning summary for inclusion in the HTML report.

    Returns dict with:
      - table_data: list of rows for a parameter-comparison table
      - best_params: the locked parameters
      - decision_rationale: text explanation
    """
    with open(tuning_summary_path) as f:
        summary = json.load(f)

    model = summary["model"]
    best = summary["best_params"]

    # Build table rows for each phase
    table_data: list[dict[str, Any]] = []

    for phase_key in ["refinement", "diversity"]:
        phase = summary.get(phase_key, {})
        best_hash = phase.get("best_hash", "")
        for score in phase.get("scores", []):
            table_data.append({
                "phase": phase_key,
                "param_hash": score["hash"],
                "qc_pass_rate": score["qc_pass"],
                "outlier_rate": score["outlier_rate"],
                "mean_n_clusters": score["n_clusters"],
                "crystal_sim": score["crystal_sim"],
                "is_best": score["hash"] == best_hash,
            })

    rationale = (
        f"Model {model}: selected by decision rule "
        f"(max QC pass → min outlier → stable clusters → crystal sanity)."
    )

    return {
        "model": model,
        "best_params": best,
        "table_data": table_data,
        "decision_rationale": rationale,
    }
