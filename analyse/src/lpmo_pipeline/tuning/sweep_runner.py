# src/lpmo_pipeline/tuning/sweep_runner.py
"""
Responsibility: Execute a single parameter sweep point through the pipeline.
Input:  ParamPoint + tuning dataset (protein×ligand list)
Output: Per-point metrics: {qc_pass_rate, outlier_rate, n_clusters, crystal_sim, ...}

This module runs STEPS 1–10 for a single parameter configuration and aggregates
the tuning metrics needed for the decision rule.
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from lpmo_pipeline.tuning.parameter_grid import ParamPoint

logger = logging.getLogger(__name__)


@dataclass
class SweepPointResult:
    """Aggregated result for a single param-grid point across all test cases."""

    param_point: ParamPoint | None = None
    param_hash: str = ""

    # QC aggregates
    posebusters_pass_rate: float = 0.0
    privateer_pass_rate: float = 0.0
    cu_his_pass_rate: float = 0.0

    # Clustering aggregates
    mean_n_clusters: float = 0.0
    mean_outlier_rate: float = 0.0
    cluster_stability_score: float = 0.0  # Higher = more stable

    # Crystal comparison (sanity check)
    mean_crystal_ifp_sim: float = 0.0
    mean_pocket_rmsd: float = 0.0

    # Per-test-case breakdown
    per_case: list[dict[str, Any]] = field(default_factory=list)


@dataclass
class TuningTestCase:
    """A single (protein, ligand) test case for tuning."""

    protein_id: str
    ligand_id: str
    protein_input: Path = Path(".")
    ligand_input: Path = Path(".")
    crystal_pdb: Path | None = None
    crystal_ligand_mol2: Path | None = None
    activity_class: str = ""  # For reference


def run_sweep_point(
    point: ParamPoint,
    test_cases: list[TuningTestCase],
    output_dir: Path,
    pipeline_config: dict[str, Any] | None = None,
) -> SweepPointResult:
    """Run the full pipeline for one parameter configuration.

    For each (protein, ligand) in test_cases:
      1. Generate predictions (with point.params)
      2. Ingest → Normalize → Privateer prep → Protonate → PLACER
      3. QC (PoseBusters + Privateer + Cu-geom)
      4. Analysis (ProLIF IFP + MDAnalysis metrics)
      5. Clustering (HDBSCAN)
      6. Crystal anchoring (if crystal available)
      7. Aggregate metrics for this test case

    Args:
        point: Parameter configuration to evaluate.
        test_cases: Tuning dataset.
        output_dir: Directory for this sweep point's outputs.
        pipeline_config: Additional pipeline settings.

    Returns:
        SweepPointResult with aggregated metrics.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    result = SweepPointResult(param_point=point, param_hash=point.param_hash)

    all_qc_pass: list[bool] = []
    all_priv_pass: list[bool] = []
    all_cuhi_pass: list[bool] = []
    all_n_clusters: list[int] = []
    all_outlier_rates: list[float] = []
    all_crystal_sims: list[float] = []

    for tc in test_cases:
        logger.info(
            "Sweep point %s: running %s × %s",
            point.param_hash, tc.protein_id, tc.ligand_id,
        )
        case_dir = output_dir / f"{tc.protein_id}_{tc.ligand_id}"
        case_dir.mkdir(parents=True, exist_ok=True)

        # --- PSEUDOCODE: Run pipeline steps 1–10 ---
        # merged_params = {**point.held_fixed, **point.params}
        #
        # # Step 1: Ingest
        # from lpmo_pipeline.io.mmcif_ingest import ingest_mmcif
        # ingest_result = ingest_mmcif(raw_cif_path, case_dir)
        #
        # # Step 2: Normalize
        # from lpmo_pipeline.io.normalize_mmcif import normalize
        # norm_result = normalize(ingest_result.structure, case_dir)
        #
        # # Step 3: Glycan expansion
        # (validate CCD monosaccharides)
        #
        # # Step 4: Protonation
        # from lpmo_pipeline.io.protonate_export import protonate_and_export
        # prot_result = protonate_and_export(norm_result.output_cif, case_dir)
        #
        # # Step 5: PLACER
        # from lpmo_pipeline.placer.run_placer import run_placer
        # placer_result = run_placer(norm_result.output_cif, case_dir / "placer_ensemble")
        #
        # # Step 6: QC
        # from lpmo_pipeline.qc.posebusters_runner import run_posebusters_batch
        # from lpmo_pipeline.qc.privateer_runner import run_privateer
        # from lpmo_pipeline.qc.custom_geometry_checks import check_geometry
        # pb_result = run_posebusters_batch(...)
        # priv_result = run_privateer(...)
        # geom_results = [check_geometry(...) for pose in placer_result.poses]
        #
        # # Step 7: Analysis
        # from lpmo_pipeline.analysis.prolif_ifp import compute_ifp_batch
        # from lpmo_pipeline.analysis.mdanalysis_metrics import compute_pose_metrics
        # ifp_batch = compute_ifp_batch(...)
        # geo_metrics = [compute_pose_metrics(...) for pose in passed_poses]
        #
        # # Step 8: Clustering
        # from lpmo_pipeline.analysis.clustering_hdbscan import run_hdbscan_clustering
        # cluster_result = run_hdbscan_clustering(ifp_batch.matrix, ...)
        #
        # # Step 10: Crystal anchoring
        # if tc.crystal_pdb:
        #     from lpmo_pipeline.analysis.crystal_anchoring import run_crystal_anchoring
        #     crystal_result = run_crystal_anchoring(tc.crystal_pdb, ...)

        # --- Collect per-case metrics ---
        case_metrics: dict[str, Any] = {
            "protein_id": tc.protein_id,
            "ligand_id": tc.ligand_id,
            "params": {**point.held_fixed, **point.params},
            # "posebusters_pass_rate": pb_result.pass_rate,
            # "privateer_pass_rate": priv_result.recognition_rate,
            # "n_clusters": cluster_result.n_clusters,
            # "outlier_rate": cluster_result.outlier_rate,
        }
        result.per_case.append(case_metrics)

        # Accumulate for aggregation
        # all_qc_pass.append(pb_result.pass_rate >= 1.0)
        # all_priv_pass.append(priv_result.all_pass)
        # all_n_clusters.append(cluster_result.n_clusters)
        # all_outlier_rates.append(cluster_result.outlier_rate)

    # --- Aggregate ---
    import numpy as np

    if all_qc_pass:
        result.posebusters_pass_rate = sum(all_qc_pass) / len(all_qc_pass)
    if all_priv_pass:
        result.privateer_pass_rate = sum(all_priv_pass) / len(all_priv_pass)
    if all_n_clusters:
        result.mean_n_clusters = float(np.mean(all_n_clusters))
    if all_outlier_rates:
        result.mean_outlier_rate = float(np.mean(all_outlier_rates))
    if all_crystal_sims:
        result.mean_crystal_ifp_sim = float(np.mean(all_crystal_sims))

    # Write per-point result
    result_path = output_dir / "sweep_result.json"
    with open(result_path, "w") as f:
        json.dump({
            "param_hash": result.param_hash,
            "params": {**point.held_fixed, **point.params},
            "posebusters_pass_rate": result.posebusters_pass_rate,
            "privateer_pass_rate": result.privateer_pass_rate,
            "mean_n_clusters": result.mean_n_clusters,
            "mean_outlier_rate": result.mean_outlier_rate,
            "mean_crystal_ifp_sim": result.mean_crystal_ifp_sim,
            "per_case": result.per_case,
        }, f, indent=2)

    logger.info(
        "Sweep point %s: QC=%.2f, outlier=%.3f, clusters=%.1f",
        result.param_hash, result.posebusters_pass_rate,
        result.mean_outlier_rate, result.mean_n_clusters,
    )
    return result
