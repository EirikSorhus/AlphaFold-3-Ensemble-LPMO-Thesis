# src/lpmo_pipeline/tuning/tune_orchestrator.py
"""
Responsibility: Top-level tuning orchestration.
  1. Load param grid
  2. Run sweep for each point (optionally in parallel)
  3. Apply decision rule to pick best params
  4. Lock best params → configs/best_{model}.yaml

Input:  configs/tuning_{model}.yaml + tuning dataset
Output: best_params.yaml, tuning_summary.json, tuning_report artifacts
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import yaml

from lpmo_pipeline.tuning.parameter_grid import ParamPoint, load_param_grid, resolve_dependencies
from lpmo_pipeline.tuning.sweep_runner import SweepPointResult, TuningTestCase, run_sweep_point
from lpmo_pipeline.utils.hashing import hash_dict

logger = logging.getLogger(__name__)


def run_tuning(
    model: str,
    tuning_config_path: Path,
    test_cases: list[TuningTestCase],
    output_dir: Path,
    pipeline_config: dict[str, Any] | None = None,
    max_parallel: int = 1,
) -> dict[str, Any]:
    """Run full tuning pipeline for one model.

    Orchestration:
      Phase 1 (refinement): sweep + pick best
      Phase 2 (diversity): resolve dependencies from phase 1, sweep + pick best

    Args:
        model: "AF3" | "RF3" | "Boltz2"
        tuning_config_path: Path to tuning YAML.
        test_cases: List of test (protein, ligand) pairs.
        output_dir: Root output directory.
        pipeline_config: Additional pipeline settings.
        max_parallel: Max parallel sweep points (1 = sequential).

    Returns:
        Dict with best params and tuning summary.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== TUNING %s ===", model)

    # --- Step 1: Load param grid ---
    all_points = load_param_grid(tuning_config_path)

    # --- Step 2: Hash each point for identification ---
    for point in all_points:
        merged = {**point.held_fixed, **point.params, "phase": point.phase}
        point.param_hash = hash_dict(merged)[:12]

    # --- Step 3: Run refinement phase ---
    refinement_points = [p for p in all_points if p.phase == "refinement"]
    refinement_results: list[SweepPointResult] = []

    for point in refinement_points:
        point_dir = output_dir / f"refinement_{point.param_hash}"
        result = run_sweep_point(point, test_cases, point_dir, pipeline_config)
        refinement_results.append(result)

    # Pick best refinement params
    best_refinement = _apply_decision_rule(refinement_results)
    best_refine_params = {
        **(best_refinement.param_point.held_fixed if best_refinement.param_point else {}),
        **(best_refinement.param_point.params if best_refinement.param_point else {}),
    }
    logger.info("Best refinement params: %s", best_refine_params)

    # --- Step 4: Resolve diversity phase dependencies ---
    diversity_points = [p for p in all_points if p.phase == "diversity"]
    diversity_points = resolve_dependencies(diversity_points, best_refine_params)

    # --- Step 5: Run diversity phase ---
    diversity_results: list[SweepPointResult] = []
    for point in diversity_points:
        point_dir = output_dir / f"diversity_{point.param_hash}"
        result = run_sweep_point(point, test_cases, point_dir, pipeline_config)
        diversity_results.append(result)

    best_diversity = _apply_decision_rule(diversity_results)
    best_div_params = {
        **(best_diversity.param_point.held_fixed if best_diversity.param_point else {}),
        **(best_diversity.param_point.params if best_diversity.param_point else {}),
    }

    # --- Step 6: Merge best params ---
    final_best = {**best_refine_params, **best_div_params}

    # --- Step 7: Write locked best params ---
    best_yaml_path = output_dir / f"best_{model.lower()}.yaml"
    with open(best_yaml_path, "w") as f:
        yaml.dump({"model": model, "locked_params": final_best}, f, default_flow_style=False)
    logger.info("Locked best params to %s", best_yaml_path)

    # --- Step 8: Write tuning summary ---
    summary = {
        "model": model,
        "best_params": final_best,
        "refinement": {
            "n_points": len(refinement_results),
            "best_hash": best_refinement.param_hash,
            "scores": [
                {
                    "hash": r.param_hash,
                    "qc_pass": r.posebusters_pass_rate,
                    "outlier_rate": r.mean_outlier_rate,
                    "n_clusters": r.mean_n_clusters,
                    "crystal_sim": r.mean_crystal_ifp_sim,
                }
                for r in refinement_results
            ],
        },
        "diversity": {
            "n_points": len(diversity_results),
            "best_hash": best_diversity.param_hash,
            "scores": [
                {
                    "hash": r.param_hash,
                    "qc_pass": r.posebusters_pass_rate,
                    "outlier_rate": r.mean_outlier_rate,
                    "n_clusters": r.mean_n_clusters,
                    "crystal_sim": r.mean_crystal_ifp_sim,
                }
                for r in diversity_results
            ],
        },
    }
    summary_path = output_dir / "tuning_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    logger.info("=== TUNING %s COMPLETE ===", model)
    return summary


def _apply_decision_rule(results: list[SweepPointResult]) -> SweepPointResult:
    """Apply the decision rule to select best parameter point.

    Priority (descending):
      1. Highest QC pass-rate (PoseBusters + Privateer)
      2. Lowest outlier-rate (HDBSCAN)
      3. Most stable cluster landscape (fewer, clearer clusters)
      4. Crystal sanity-check (IFP sim ≥ baseline; do not degrade)

    Returns:
        Best SweepPointResult.
    """
    if not results:
        logger.warning("No results to apply decision rule on")
        return SweepPointResult()

    def sort_key(r: SweepPointResult) -> tuple:
        # Maximize QC pass, minimize outlier, maximize crystal sim
        return (
            r.posebusters_pass_rate + r.privateer_pass_rate,  # Higher = better
            -r.mean_outlier_rate,                              # Lower = better (negate)
            -r.mean_n_clusters,                                # Fewer = better (negate)
            r.mean_crystal_ifp_sim,                            # Higher = better
        )

    results_sorted = sorted(results, key=sort_key, reverse=True)
    best = results_sorted[0]

    logger.info(
        "Decision rule: best=%s (QC=%.2f, outlier=%.3f, clusters=%.1f, xtal=%.3f)",
        best.param_hash, best.posebusters_pass_rate,
        best.mean_outlier_rate, best.mean_n_clusters, best.mean_crystal_ifp_sim,
    )
    return best
