# src/lpmo_pipeline/qc/hard_qc_orchestrator.py
"""
Responsibility: Orchestrate the hard QC sequence for a batch of poses.

Sequence per pose:
    1. Run active-site proximity pre-gate (Stage 8)
    2. Run Cu-His + substrate geometry checks
    3. Run PoseBusters
    4. Run Privateer when input is available
    5. Aggregate into unified QC verdict via compute_verdict()

INVARIANT — "Ikke-slett regel":
  All numeric metrics are preserved in the verdict even when the pose is
    dropped by a hard gate. Nothing computed is discarded, so thresholds can
    be revisited without rerunning downstream QC.
"""
from __future__ import annotations

import logging
from time import perf_counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from lpmo_pipeline.config import load_defaults_config, load_runtime_paths_config
from lpmo_pipeline.qc.active_site_proximity import check_active_site_proximity
from lpmo_pipeline.qc.custom_geometry_checks import GeometryResult, check_geometry
from lpmo_pipeline.qc.gates import GateConfig, load_gate_config_from_yaml
from lpmo_pipeline.qc.posebusters_runner import PoseBustersSingleResult, run_posebusters_single
from lpmo_pipeline.qc.privateer_runner import (
    PrivateerBatchInput,
    PrivateerResult,
    run_privateer_batch,
)
from lpmo_pipeline.qc.qc_report import (
    PoseQCVerdict,
    QCReport,
    build_qc_report,
    compute_verdict,
)

logger = logging.getLogger(__name__)

_RUNTIME_PATHS = load_runtime_paths_config()
_DEFAULTS_CONFIG = load_defaults_config()
_CHAIN_SCHEMA = _DEFAULTS_CONFIG.get("chain_schema") or {}
DEFAULT_CU_CHAIN = str(_CHAIN_SCHEMA.get("metal") or "E")
DEFAULT_GLYCAN_CHAINS = tuple(str(chain) for chain in (_CHAIN_SCHEMA.get("glycans") or ["B", "C", "D"]))


def _default_qc_config_path() -> Path:
    return _RUNTIME_PATHS.pipeline_assets.thresholds_config


@dataclass
class _PoseEvaluation:
    pose: HardQCInput
    proximity_result: Any
    pb_result: PoseBustersSingleResult | None
    geom_result: GeometryResult | None
    timing_events: list[dict[str, Any]] = field(default_factory=list)


def _posebusters_runner_error_result(pose_id: str, exc: Exception) -> PoseBustersSingleResult:
    return PoseBustersSingleResult(
        pose_id=pose_id,
        passed=False,
        critical_errors=["posebusters_runner_error"],
        details={"runner_error": str(exc)},
    )


def _privateer_runner_error_result(pose_id: str, error: str) -> PrivateerResult:
    return PrivateerResult(
        pose_id=pose_id,
        all_pass=False,
        runner_error=error,
    )


# ---------------------------------------------------------------------------
# Input model
# ---------------------------------------------------------------------------
@dataclass
class HardQCInput:
    """Input for a single pose entering hard QC.

    Attributes:
        pose_id: Unique identifier for this pose.
        mol_pred_path: Path to the PDB file (for PoseBusters).
        structure: Parsed structure object (gemmi) for geometry checks.
        protein_path: Optional protein PDB for PoseBusters dock/redock mode.
        reference_path: Optional reference ligand for PoseBusters redock RMSD.
        privateer_cif_path: Optional path to privateer_input.cif.
        cu_chain: Expected Cu chain (default ``"E"`` per normalisation).
        glycan_chains: Expected glycan chains (default ``["B", "C", "D"]``).
    """

    pose_id: str
    mol_pred_path: Path
    structure: Any  # gemmi.Structure
    protein_path: Path | None = None
    reference_path: Path | None = None
    privateer_cif_path: Path | None = None
    cu_chain: str = DEFAULT_CU_CHAIN
    glycan_chains: list[str] | None = field(default_factory=lambda: list(DEFAULT_GLYCAN_CHAINS))


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------
def run_hard_qc(
    poses: list[HardQCInput],
    run_id: str = "",
    config: GateConfig | None = None,
    config_path: Path | None = None,
    *,
    run_posebusters: bool = True,
    run_privateer: bool = True,
    max_workers: int | None = None,
    collect_timing_events: bool = False,
) -> QCReport:
    """Run the full hard-QC sequence on a list of poses.

    For each pose the following checks are executed:
        1. **Active-site proximity** — pre-QC gate (Stage 8).
            If this fails, the pose is dropped from further QC and downstream
            analysis.
        2. **Cu-His + substrate geometry** — LPMO-specific active-site check.
            If this fails, the pose is dropped from further QC and downstream
            analysis.
        3. **PoseBusters** — chemical/stereochemical validation.
            If the PoseBusters backend raises an unexpected exception, the
            pose gets an explicit hard-fail PoseBusters result.
        4. **Privateer** — if ``privateer_cif_path`` is provided.
            If the Privateer batch backend raises an unexpected exception, all
            Privateer-eligible poses get explicit hard-fail Privateer results.

    All numeric metrics are retained regardless of pass/fail status
    ("Ikke-slett regel").

    Args:
        poses: List of :class:`HardQCInput` objects.
        run_id: Identifier for this QC run (used in the report).
        run_posebusters: Whether to run PoseBusters for QC-eligible poses.
        run_privateer: Whether to run Privateer for QC-eligible poses.

    Returns:
        :class:`QCReport` with per-pose verdicts and aggregate counts.
    """
    gate_config = config or load_gate_config_from_yaml(config_path or _default_qc_config_path())
    verdicts: list[PoseQCVerdict] = []
    worker_count = max(1, int(max_workers or 1))
    total_start = perf_counter()

    def _timing_event(
        *,
        step: str,
        start: float,
        pose_id: str = "",
        status: str = "ok",
        detail: str = "",
    ) -> dict[str, Any]:
        return {
            "step": step,
            "pose_id": pose_id,
            "condition_id": "",
            "worker_count": worker_count,
            "wall_time_s": round(perf_counter() - start, 6),
            "status": status,
            "detail": detail,
        }

    def _evaluate_pose(pose: HardQCInput) -> _PoseEvaluation:
        logger.info("Hard QC: processing pose %s", pose.pose_id)
        pose_timing_events: list[dict[str, Any]] = []

        def _append_pose_timing(event: dict[str, Any]) -> None:
            if collect_timing_events:
                pose_timing_events.append(event)

        # --- Step 1: Active-site proximity pre-QC ---
        proximity_start = perf_counter()
        proximity_result = check_active_site_proximity(
            structure=pose.structure,
            pose_id=pose.pose_id,
            cu_chain=pose.cu_chain,
            glycan_chains=pose.glycan_chains,
            hard_cutoff_a=gate_config.active_site_proximity_max_a,
            cu_c_soft_flag_a=gate_config.cu_c_proximity_threshold_a,
        )
        _append_pose_timing(
            _timing_event(
                step="hard_qc.active_site_proximity",
                start=proximity_start,
                pose_id=pose.pose_id,
                status="passed" if proximity_result.passed else "failed",
                detail=";".join(proximity_result.failure_reasons or []),
            )
        )

        # --- Step 2: Geometry (Cu-His + Cu-substrate) ---
        geom_result: GeometryResult | None = None
        if proximity_result.passed:
            geometry_start = perf_counter()
            try:
                geom_result = check_geometry(
                    structure=pose.structure,
                    pose_id=pose.pose_id,
                    cu_chain=pose.cu_chain,
                    glycan_chains=pose.glycan_chains,
                    hard_cu_his_min_a=gate_config.cu_his_dist_min,
                    hard_cu_his_max_a=gate_config.cu_his_dist_max,
                    soft_cu_his_min_a=gate_config.cu_his_soft_min_a,
                    soft_cu_his_max_a=gate_config.cu_his_soft_max_a,
                    cu_c_soft_flag_a=gate_config.cu_c_proximity_threshold_a,
                    his_brace_max_search_a=gate_config.his_brace_max_search_a,
                )
                _append_pose_timing(
                    _timing_event(
                        step="hard_qc.geometry",
                        start=geometry_start,
                        pose_id=pose.pose_id,
                        status="passed" if geom_result.passed else "failed",
                        detail=";".join(geom_result.failure_reasons or []),
                    )
                )
            except Exception:
                logger.exception(
                    "Geometry check raised an unexpected error for pose %s",
                    pose.pose_id,
                )
                geom_result = GeometryResult(
                    pose_id=pose.pose_id,
                    passed=False,
                    failure_reasons=["geometry_runner_error"],
                )
                _append_pose_timing(
                    _timing_event(
                        step="hard_qc.geometry",
                        start=geometry_start,
                        pose_id=pose.pose_id,
                        status="error",
                        detail="geometry_runner_error",
                    )
                )
        else:
            logger.info(
                "Pre-QC hard gate failed for %s; skipping geometry, PoseBusters, and Privateer",
                pose.pose_id,
            )
            _append_pose_timing(
                {
                    "step": "hard_qc.geometry",
                    "pose_id": pose.pose_id,
                    "condition_id": "",
                    "worker_count": worker_count,
                    "wall_time_s": 0.0,
                    "status": "skipped",
                    "detail": "active_site_proximity_failed",
                }
            )

        qc_eligible = proximity_result.passed and geom_result is not None and geom_result.passed

        # --- Step 3: PoseBusters ---
        pb_result: PoseBustersSingleResult | None = None
        if qc_eligible and run_posebusters:
            posebusters_start = perf_counter()
            try:
                pb_result = run_posebusters_single(
                    pdb_path=pose.mol_pred_path,
                    pose_id=pose.pose_id,
                    protein_path=pose.protein_path,
                    reference_path=pose.reference_path,
                )
                _append_pose_timing(
                    _timing_event(
                        step="hard_qc.posebusters",
                        start=posebusters_start,
                        pose_id=pose.pose_id,
                        status="passed" if pb_result.passed else "failed",
                        detail=";".join(pb_result.critical_errors or []),
                    )
                )
            except Exception as exc:
                logger.exception(
                    "PoseBusters raised an unexpected error for pose %s; failing pose closed",
                    pose.pose_id,
                )
                pb_result = _posebusters_runner_error_result(pose.pose_id, exc)
                _append_pose_timing(
                    _timing_event(
                        step="hard_qc.posebusters",
                        start=posebusters_start,
                        pose_id=pose.pose_id,
                        status="error",
                        detail="posebusters_runner_error",
                    )
                )
        elif not qc_eligible:
            logger.info(
                "Hard distance QC failed for %s; skipping PoseBusters and Privateer",
                pose.pose_id,
            )
            _append_pose_timing(
                {
                    "step": "hard_qc.posebusters",
                    "pose_id": pose.pose_id,
                    "condition_id": "",
                    "worker_count": worker_count,
                    "wall_time_s": 0.0,
                    "status": "skipped",
                    "detail": "hard_distance_qc_failed",
                }
            )
        else:
            _append_pose_timing(
                {
                    "step": "hard_qc.posebusters",
                    "pose_id": pose.pose_id,
                    "condition_id": "",
                    "worker_count": worker_count,
                    "wall_time_s": 0.0,
                    "status": "disabled",
                    "detail": "run_posebusters_false",
                }
            )

        return _PoseEvaluation(
            pose=pose,
            proximity_result=proximity_result,
            pb_result=pb_result,
            geom_result=geom_result,
            timing_events=pose_timing_events,
        )

    if worker_count <= 1 or len(poses) <= 1:
        evaluations = [_evaluate_pose(pose) for pose in poses]
    else:
        with ThreadPoolExecutor(max_workers=min(worker_count, len(poses))) as executor:
            evaluations = list(executor.map(_evaluate_pose, poses))

    timing_events: list[dict[str, Any]] = []
    if collect_timing_events:
        for evaluation in evaluations:
            timing_events.extend(evaluation.timing_events)

    privateer_results: dict[str, PrivateerResult] = {}
    privateer_inputs = [
        PrivateerBatchInput(cif_path=evaluation.pose.privateer_cif_path, pose_id=evaluation.pose.pose_id)
        for evaluation in evaluations
        if run_privateer
        and evaluation.pose.privateer_cif_path is not None
        and evaluation.proximity_result.passed
        and evaluation.geom_result is not None
        and evaluation.geom_result.passed
    ]
    if privateer_inputs:
        privateer_start = perf_counter()
        try:
            for result in run_privateer_batch(privateer_inputs, max_workers=worker_count):
                privateer_results[result.pose_id] = result
            missing_pose_ids = [
                item.pose_id for item in privateer_inputs if item.pose_id not in privateer_results
            ]
            for pose_id in missing_pose_ids:
                privateer_results[pose_id] = _privateer_runner_error_result(
                    pose_id,
                    "privateer_missing_result",
                )
            if collect_timing_events:
                timing_events.append(
                    _timing_event(
                        step="hard_qc.privateer_batch",
                        start=privateer_start,
                        status="ok",
                        detail=f"n_inputs={len(privateer_inputs)}",
                    )
                )
        except Exception as exc:
            logger.exception(
                "Privateer batch raised an unexpected error; failing eligible poses closed"
            )
            for item in privateer_inputs:
                privateer_results[item.pose_id] = _privateer_runner_error_result(
                    item.pose_id,
                    f"privateer_batch_runner_error:{exc}",
                )
            if collect_timing_events:
                timing_events.append(
                    _timing_event(
                        step="hard_qc.privateer_batch",
                        start=privateer_start,
                        status="error",
                        detail=f"n_inputs={len(privateer_inputs)}",
                    )
                )
        if collect_timing_events and privateer_inputs:
            missing_pose_ids = [
                item.pose_id
                for item in privateer_inputs
                if privateer_results.get(item.pose_id, PrivateerResult()).runner_error
                == "privateer_missing_result"
            ]
            if missing_pose_ids:
                timing_events.append(
                    _timing_event(
                        step="hard_qc.privateer_batch_missing_results",
                        start=privateer_start,
                        status="error",
                        detail=f"pose_ids={','.join(missing_pose_ids)}",
                    )
                )
    else:
        if collect_timing_events:
            timing_events.append(
                {
                    "step": "hard_qc.privateer_batch",
                    "pose_id": "",
                    "condition_id": "",
                    "worker_count": worker_count,
                    "wall_time_s": 0.0,
                    "status": "skipped" if run_privateer else "disabled",
                    "detail": "n_inputs=0",
                }
            )

    for evaluation in evaluations:
        verdict = compute_verdict(
            pose_id=evaluation.pose.pose_id,
            pb_result=evaluation.pb_result,
            priv_result=privateer_results.get(evaluation.pose.pose_id),
            geom_result=evaluation.geom_result,
            proximity_result=evaluation.proximity_result,
        )
        verdicts.append(verdict)

        logger.info(
            "Hard QC verdict for %s: %s (drop_reasons=%s)",
            evaluation.pose.pose_id,
            verdict.status,
            verdict.drop_reasons or "none",
        )

    if collect_timing_events:
        timing_events.append(
            _timing_event(
                step="hard_qc.total",
                start=total_start,
                status="ok",
                detail=f"n_poses={len(poses)}",
            )
        )
    report = build_qc_report(run_id=run_id, verdicts=verdicts, timing_events=timing_events)
    logger.info(
        "Hard QC complete: %d total, %d passed, %d flagged, %d dropped",
        report.total,
        report.passed,
        report.flagged,
        report.dropped,
    )
    return report
