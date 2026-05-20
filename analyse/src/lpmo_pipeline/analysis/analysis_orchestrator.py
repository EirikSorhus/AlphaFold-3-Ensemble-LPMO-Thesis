"""Production analysis orchestration through ProLIF and Stage 6 clustering.

This module provides the real production control path for discovery,
normalization, hard QC, downstream geometry, ligand-resolved ProLIF IFP,
condition-wise HDBSCAN clustering, and report generation.
"""
from __future__ import annotations

import json
import logging
import re
import traceback
import csv
from time import perf_counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import yaml
from jsonschema import validate

from lpmo_pipeline.config import load_defaults_config, load_runtime_paths_config
from lpmo_pipeline.analysis.clustering_agglomerative import (
    AgglomerativeJaccardClusterer,
    AgglomerativeJaccardConfig,
)
from lpmo_pipeline.analysis.clustering_hdbscan import (
    ClusteringResult,
    HDBSCANClusterer,
    build_cluster_assignment_rows,
    build_condition_cluster_summary,
    build_medoid_rows,
)
from lpmo_pipeline.analysis.clustering_pilot import (
    DEFAULT_EXCLUDED_CLUSTERING_INTERACTION_TYPES,
    DEFAULT_INSUFFICIENT_CLUSTERABLE_SIGNAL_LABEL,
    DEFAULT_MAIN_CLUSTERING_INTERACTION_TYPES,
    DEFAULT_MINIMUM_CLUSTERABLE_N,
    PilotConditionIFP,
    build_feature_prevalence,
    build_filtered_ifp_matrix,
    build_interaction_type_prevalence,
    select_main_clustering_features,
    summarize_condition_matrices,
)
from lpmo_pipeline.analysis.cluster_signatures import (
    build_cluster_signature_tables,
    write_cluster_ifp_signature_table,
    write_cluster_residue_signature_table,
    write_cluster_signature_summary_json,
)
from lpmo_pipeline.analysis.convergence_metrics import (
    ConditionConvergenceSummary,
    ConvergencePoseInput,
    PoseConvergenceMetrics,
    compute_condition_convergence,
    write_condition_convergence_summary_tsv,
    write_pose_convergence_tsv,
)
from lpmo_pipeline.analysis.prolif_ifp import (
    ContactEligibility,
    IFPResult,
    compute_ifp_batch,
    evaluate_contact_eligibility,
    load_contact_eligibility_rule,
    load_prolif_features_config,
    write_ifp_matrix,
    write_pose_ifp_table,
)
from lpmo_pipeline.analysis.residue_contact_extraction import (
    build_pose_residue_contact_rows,
    write_pose_residue_contact_table,
)
from lpmo_pipeline.analysis.residue_importance import (
    compute_residue_importance_outputs,
    write_condition_patch_summary,
    write_protein_condition_residue_scores,
    write_protein_patch_summary,
    write_protein_residue_regio_delta,
)
from lpmo_pipeline.analysis.crystal_anchoring import (
    run_crystal_reference_screen,
    write_crystal_reference_screen_report,
)
from lpmo_pipeline.analysis.mdanalysis_metrics import (
    PoseGeometryMetrics,
    compute_pose_metrics_from_structure,
    write_pose_geometry_tsv,
)
from lpmo_pipeline.io.cif_to_pdb import convert_cif_to_pdb
from lpmo_pipeline.io.discovery import discover_work_root
from lpmo_pipeline.io.gemmi_compat import gemmi
from lpmo_pipeline.io.normalize_mmcif import NormalizeMMCIFRunner
from lpmo_pipeline.io.protonate_export import protonate_and_export
from lpmo_pipeline.qc.hard_qc_orchestrator import HardQCInput, run_hard_qc
from lpmo_pipeline.qc.privateer_runner import prepare_privateer_input
from lpmo_pipeline.qc.qc_report import write_qc_report
from lpmo_pipeline.report.build_metrics_csv import build_metrics_csv, merge_pose_record
from lpmo_pipeline.report.build_report_html import build_report_html
from lpmo_pipeline.report.build_summary_json import build_summary, write_summary_json

logger = logging.getLogger(__name__)

_RUNTIME_PATHS = load_runtime_paths_config()
_DEFAULTS_CONFIG = load_defaults_config()
_TARGET_PREFIX_TO_SUBSTRATE = {
    str(prefix).upper(): str(substrate)
    for prefix, substrate in (
        (_DEFAULTS_CONFIG.get("target_metadata") or {}).get("substrate_by_prefix")
        or {
            "NAG": "chitin",
            "CEL": "cellulose",
            "BGC": "cellulose",
            "STA": "amylose",
            "GLC": "amylose",
        }
    ).items()
}


@dataclass(frozen=True)
class ProductionRunOptions:
    run_id: str
    output_dir: Path
    del_variant: str
    work_root: Path
    af3_only: bool = True
    latest_only: bool = True
    max_cases: int | None = None
    include_targets: tuple[str, ...] = ()
    include_proteins: tuple[str, ...] = ()
    construct_type: str = "domain_only"
    run_posebusters: bool = True
    run_privateer: bool = True
    clustering_pilot: "ClusteringPilotOptions | None" = None
    n_jobs: int = 1
    collect_timing_events: bool = False


@dataclass(frozen=True)
class ClusteringPilotOptions:
    enabled: bool = False
    label: str = "default"
    output_dirname: str = "clustering_pilot"
    main_include_interaction_types: tuple[str, ...] = DEFAULT_MAIN_CLUSTERING_INTERACTION_TYPES
    extra_include_interaction_types: tuple[str, ...] = ()
    excluded_interaction_types: tuple[str, ...] = DEFAULT_EXCLUDED_CLUSTERING_INTERACTION_TYPES
    rare_feature_pose_prevalence_lt: float = 0.01
    rare_feature_condition_prevalence_lt_n_conditions: int = 2
    minimum_clusterable_n: int = DEFAULT_MINIMUM_CLUSTERABLE_N
    insufficient_clusterable_signal_label: str = DEFAULT_INSUFFICIENT_CLUSTERABLE_SIGNAL_LABEL
    agglomerative_linkage: str = "average"
    agglomerative_distance_threshold: float = 0.5
    agglomerative_min_cluster_size: int = DEFAULT_MINIMUM_CLUSTERABLE_N


@dataclass(frozen=True)
class PoseInputRecord:
    cif_path: Path
    pose_id: str
    protein_id: str
    ligand_id: str
    model: str
    source_run_id: str
    discovered_run_id: str | None = None
    confidence_json_path: Path | None = None
    run_status: str = ""
    seed: int | None = None
    sample: int | None = None


@dataclass
class PreparedPose:
    pose: PoseInputRecord
    case_dir: Path
    normalized_cif: Path
    posebusters_pdb: Path
    privateer_input_cif: Path
    structure: Any


@dataclass(frozen=True)
class AnalysisCoreResult:
    run_id: str
    output_dir: Path
    summary_path: Path
    qc_report_path: Path | None
    pose_manifest_tsv_path: Path | None
    pose_confidence_tsv_path: Path | None
    structure_index_tsv_path: Path | None
    qc_attrition_tsv_path: Path | None
    pose_geometry_tsv_path: Path | None
    metrics_csv_path: Path | None
    summary_json_path: Path | None
    report_html_path: Path | None
    n_discovered: int
    n_prepared: int
    n_analyzed: int
    success: bool
    pose_ifp_table_tsv_path: Path | None = None
    pose_residue_contact_tsv_path: Path | None = None
    pose_convergence_tsv_path: Path | None = None
    condition_convergence_summary_tsv_path: Path | None = None
    cluster_assignments_tsv_path: Path | None = None
    medoid_manifest_tsv_path: Path | None = None
    condition_cluster_summary_tsv_path: Path | None = None
    cluster_ifp_signature_tsv_path: Path | None = None
    cluster_residue_signature_tsv_path: Path | None = None
    cluster_signatures_json_path: Path | None = None
    cluster_annotation_stage_completed: bool = False
    protein_condition_residue_scores_tsv_path: Path | None = None
    protein_residue_regio_delta_tsv_path: Path | None = None
    condition_patch_summary_tsv_path: Path | None = None
    protein_patch_summary_tsv_path: Path | None = None
    crystal_anchor_tsv_path: Path | None = None
    n_cluster_conditions: int = 0
    crystal_anchoring_stage_completed: bool = False
    n_crystal_anchoring_conditions: int = 0
    n_crystal_anchoring_errors: int = 0


def _load_qc_report_schema() -> dict[str, Any]:
    return json.loads(_RUNTIME_PATHS.pipeline_assets.qc_report_schema.read_text())


def _slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def _infer_source_run_id(cif_path: Path, fallback_run_id: str) -> str:
    resolved_path = cif_path.resolve()
    for ancestor in resolved_path.parents:
        if ancestor.parent.name == "runs" and ancestor.name.isdigit():
            return ancestor.name
    return fallback_run_id


def _select_confidence_json(cif_path: Path, candidates: list[Path]) -> Path | None:
    """Choose the confidence JSON that most likely belongs to one CIF pose."""
    if not candidates:
        return None

    cif_stem = cif_path.stem
    for candidate in candidates:
        candidate_stem = candidate.stem
        if candidate_stem == cif_stem:
            return candidate
        if cif_stem in candidate_stem or candidate_stem in cif_stem:
            return candidate
    return candidates[0]


def _flatten_numeric_values(value: Any) -> list[float]:
    if isinstance(value, (int, float)):
        return [float(value)]
    if isinstance(value, list):
        values: list[float] = []
        for item in value:
            values.extend(_flatten_numeric_values(item))
        return values
    return []


def _coerce_scalar(value: Any) -> Any:
    """Convert common confidence JSON values to TSV-friendly scalars."""
    if value is None:
        return None
    if isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, list):
        numeric_values = _flatten_numeric_values(value)
        if numeric_values:
            return float(np.mean(numeric_values))
    return json.dumps(value, sort_keys=True)


def _first_present(payload: dict[str, Any], keys: tuple[str, ...]) -> Any:
    for key in keys:
        if key in payload:
            return _coerce_scalar(payload.get(key))
    return None


def _read_confidence_summary(confidence_json_path: Path | None) -> dict[str, Any]:
    """Read the AF3 confidence JSON fields used by pose_confidence.tsv."""
    if confidence_json_path is None:
        return {
            "confidence_json_status": "missing",
            "ranking_score": None,
            "iptm": None,
            "ptm": None,
            "mean_plddt": None,
            "ligand_interface_confidence": None,
            "pae_summary": None,
        }

    try:
        payload = json.loads(confidence_json_path.read_text())
    except Exception as exc:
        return {
            "confidence_json_status": f"parse_error:{exc.__class__.__name__}",
            "ranking_score": None,
            "iptm": None,
            "ptm": None,
            "mean_plddt": None,
            "ligand_interface_confidence": None,
            "pae_summary": None,
        }

    return {
        "confidence_json_status": "ok",
        "ranking_score": _first_present(
            payload,
            ("ranking_score", "ranking_confidence", "confidence_score", "score"),
        ),
        "iptm": _first_present(
            payload,
            ("iptm", "ipTM", "iptm_score", "interface_predicted_tm_score"),
        ),
        "ptm": _first_present(payload, ("ptm", "pTM", "predicted_tm_score")),
        "mean_plddt": _first_present(
            payload,
            ("mean_plddt", "mean_pLDDT", "plddt_mean", "plddt"),
        ),
        "ligand_interface_confidence": _first_present(
            payload,
            (
                "ligand_interface_confidence",
                "ligand_iptm",
                "chain_pair_iptm",
                "interface_confidence",
            ),
        ),
        "pae_summary": _first_present(
            payload,
            ("pae_summary", "mean_pae", "predicted_aligned_error"),
        ),
    }


def _normalize_del_variant(value: str) -> str:
    return value.removeprefix("del_")


def _parse_target_metadata(ligand_id: str) -> tuple[str, int]:
    match = re.match(r"^([A-Za-z]+)(\d+)$", ligand_id)
    if not match:
        return "", 0
    prefix = match.group(1).upper()
    dp = int(match.group(2))
    return _TARGET_PREFIX_TO_SUBSTRATE.get(prefix, ""), dp


def load_production_options(
    config: dict[str, Any],
    output_dir: Path,
    del_variant: str,
    n_jobs: int | None = None,
) -> ProductionRunOptions:
    production = config.get("production") or {}
    construct_type = str(production.get("construct_type") or "domain_only")
    work_roots = production.get("work_roots") or config.get("work_roots") or {}
    if work_roots and not isinstance(work_roots, dict):
        raise TypeError("Production config field 'work_roots' must be a mapping")

    work_root_value = (
        work_roots.get(construct_type)
        or production.get("work_root")
        or config.get("work_root")
    )
    if not work_root_value:
        raise ValueError(
            "Production config must define production.work_roots[construct_type] or production.work_root"
        )

    max_cases_raw = production.get("max_cases")
    max_cases = int(max_cases_raw) if max_cases_raw is not None else None
    include_targets = tuple(str(value) for value in production.get("include_targets", ()))
    include_proteins = tuple(str(value) for value in production.get("include_proteins", ()))
    run_posebusters = bool(production.get("run_posebusters", True))
    run_privateer = bool(production.get("run_privateer", True))
    production_n_jobs = int(production.get("n_jobs", 1) or 1)
    clustering_pilot_config = production.get("clustering_pilot") or {}
    clustering_pilot: ClusteringPilotOptions | None = None
    if bool(clustering_pilot_config.get("enabled", False)):
        clustering_pilot = ClusteringPilotOptions(
            enabled=True,
            label=str(clustering_pilot_config.get("label") or "default"),
            output_dirname=str(clustering_pilot_config.get("output_dirname") or "clustering_pilot"),
            main_include_interaction_types=tuple(
                str(value)
                for value in clustering_pilot_config.get(
                    "main_include_interaction_types",
                    DEFAULT_MAIN_CLUSTERING_INTERACTION_TYPES,
                )
            ),
            extra_include_interaction_types=tuple(
                str(value)
                for value in clustering_pilot_config.get("extra_include_interaction_types", ())
            ),
            excluded_interaction_types=tuple(
                str(value)
                for value in clustering_pilot_config.get(
                    "excluded_interaction_types",
                    DEFAULT_EXCLUDED_CLUSTERING_INTERACTION_TYPES,
                )
            ),
            rare_feature_pose_prevalence_lt=float(
                clustering_pilot_config.get("rare_feature_pose_prevalence_lt", 0.01)
            ),
            rare_feature_condition_prevalence_lt_n_conditions=int(
                clustering_pilot_config.get("rare_feature_condition_prevalence_lt_n_conditions", 2)
            ),
            minimum_clusterable_n=int(
                clustering_pilot_config.get("minimum_clusterable_n", DEFAULT_MINIMUM_CLUSTERABLE_N)
            ),
            insufficient_clusterable_signal_label=str(
                clustering_pilot_config.get(
                    "insufficient_clusterable_signal_label",
                    DEFAULT_INSUFFICIENT_CLUSTERABLE_SIGNAL_LABEL,
                )
            ),
            agglomerative_linkage=str(clustering_pilot_config.get("agglomerative_linkage", "average")),
            agglomerative_distance_threshold=float(
                clustering_pilot_config.get("agglomerative_distance_threshold", 0.5)
            ),
            agglomerative_min_cluster_size=int(
                clustering_pilot_config.get(
                    "agglomerative_min_cluster_size",
                    clustering_pilot_config.get("minimum_clusterable_n", DEFAULT_MINIMUM_CLUSTERABLE_N),
                )
            ),
        )

    return ProductionRunOptions(
        run_id=str(production.get("run_id") or output_dir.name),
        output_dir=output_dir.resolve(),
        del_variant=_normalize_del_variant(del_variant),
        work_root=Path(work_root_value).resolve(),
        af3_only=bool(production.get("af3_only", True)),
        latest_only=bool(production.get("latest_only", True)),
        max_cases=max_cases,
        include_targets=include_targets,
        include_proteins=include_proteins,
        construct_type=construct_type,
        run_posebusters=run_posebusters,
        run_privateer=run_privateer,
        clustering_pilot=clustering_pilot,
        n_jobs=max(1, int(n_jobs if n_jobs is not None else production_n_jobs)),
        collect_timing_events=bool(production.get("collect_timing_events", False)),
    )


def _condition_id_for_pose(pose: PoseInputRecord, construct_type: str) -> str:
    substrate_type, dp = _parse_target_metadata(pose.ligand_id)
    substrate_token = substrate_type or "unknown"
    return f"{pose.protein_id}__{construct_type}__{substrate_token}_DP{dp}"


def _make_ifp_failure_result(pose_id: str, error: str) -> IFPResult:
    return IFPResult(
        pose_id=pose_id,
        status="input_missing",
        error=error,
        interaction_counts={},
    )


def _write_tsv_rows(output_path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_ifp_matrix_csv(
    output_path: Path,
    pose_ids: list[str],
    feature_names: list[str],
    matrix: list[list[int]],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["pose_id"] + feature_names)
        for pose_id, row in zip(pose_ids, matrix, strict=True):
            writer.writerow([pose_id] + row)


def _write_clustering_pilot_outputs(
    options: ClusteringPilotOptions,
    output_dir: Path,
    pilot_conditions: list[PilotConditionIFP],
) -> dict[str, Any]:
    pilot_output_dir = output_dir / options.output_dirname / _slug(options.label)
    pilot_output_dir.mkdir(parents=True, exist_ok=True)

    interaction_type_prevalence = build_interaction_type_prevalence(pilot_conditions)
    feature_prevalence = build_feature_prevalence(pilot_conditions)
    interaction_type_rows = [row.to_row() for row in interaction_type_prevalence]
    feature_rows = [row.to_row() for row in feature_prevalence]
    selected_main_feature_names = select_main_clustering_features(
        feature_prevalence,
        default_include_interaction_types=options.main_include_interaction_types,
        extra_include_interaction_types=options.extra_include_interaction_types,
        excluded_interaction_types=options.excluded_interaction_types,
        rare_feature_pose_prevalence_lt=options.rare_feature_pose_prevalence_lt,
        rare_feature_condition_prevalence_lt_n_conditions=options.rare_feature_condition_prevalence_lt_n_conditions,
    )

    interaction_type_prevalence_path = pilot_output_dir / "pilot_ifp_interaction_type_prevalence.tsv"
    _write_tsv_rows(
        interaction_type_prevalence_path,
        [
            "interaction_type",
            "n_selected_poses_with_type",
            "pose_prevalence",
            "n_conditions_with_type",
            "condition_prevalence",
            "median_count_per_pose_when_present",
        ],
        interaction_type_rows,
    )

    feature_prevalence_path = pilot_output_dir / "pilot_ifp_feature_prevalence.tsv"
    _write_tsv_rows(
        feature_prevalence_path,
        [
            "feature_name",
            "interaction_type",
            "n_selected_poses_with_feature",
            "pose_prevalence",
            "n_conditions_with_feature",
            "condition_prevalence",
            "within_condition_max_prevalence",
            "within_condition_variance",
        ],
        feature_rows,
    )

    feature_selection_path = pilot_output_dir / "pilot_feature_selection.json"
    feature_selection_path.write_text(
        json.dumps(
            {
                "label": options.label,
                "main_include_interaction_types": list(options.main_include_interaction_types),
                "extra_include_interaction_types": list(options.extra_include_interaction_types),
                "excluded_interaction_types": list(options.excluded_interaction_types),
                "rare_feature_pose_prevalence_lt": options.rare_feature_pose_prevalence_lt,
                "rare_feature_condition_prevalence_lt_n_conditions": options.rare_feature_condition_prevalence_lt_n_conditions,
                "minimum_clusterable_n": options.minimum_clusterable_n,
                "insufficient_clusterable_signal_label": options.insufficient_clusterable_signal_label,
                "agglomerative_linkage": options.agglomerative_linkage,
                "agglomerative_distance_threshold": options.agglomerative_distance_threshold,
                "agglomerative_min_cluster_size": options.agglomerative_min_cluster_size,
                "selected_main_feature_count": len(selected_main_feature_names),
                "selected_main_feature_names": selected_main_feature_names,
            },
            indent=2,
        )
    )

    pilot_method_summary_rows: list[dict[str, Any]] = []
    condition_matrix_entries: list[dict[str, Any]] = []
    for condition in pilot_conditions:
        raw_matrix = build_filtered_ifp_matrix(
            condition.batch,
            selected_pose_ids=condition.selected_pose_ids,
        )
        main_matrix = build_filtered_ifp_matrix(
            condition.batch,
            selected_pose_ids=condition.selected_pose_ids,
            allowed_feature_names=selected_main_feature_names,
        )
        condition_output_dir = pilot_output_dir / "conditions" / _slug(condition.condition_id)
        raw_matrix_path = condition_output_dir / "raw_contact_eligible_ifp_matrix.csv"
        main_matrix_path = condition_output_dir / "main_contact_eligible_ifp_matrix.csv"
        _write_ifp_matrix_csv(raw_matrix_path, raw_matrix.pose_ids, raw_matrix.feature_names, raw_matrix.matrix)
        _write_ifp_matrix_csv(main_matrix_path, main_matrix.pose_ids, main_matrix.feature_names, main_matrix.matrix)

        matrix_summary = summarize_condition_matrices(
            condition.condition_id,
            raw_matrix,
            main_matrix,
            minimum_clusterable_n=options.minimum_clusterable_n,
            insufficient_clusterable_signal_label=options.insufficient_clusterable_signal_label,
        )

        method_entries: list[dict[str, Any]] = []
        method_specs = [
            (
                "hdbscan",
                (
                    lambda method_output_dir: HDBSCANClusterer(output_dir=method_output_dir),
                    lambda clusterer: (
                        f"min_cluster_size={clusterer.config.min_cluster_size};"
                        f"min_samples={clusterer.config.min_samples};"
                        f"selection={clusterer.config.cluster_selection_method}"
                    ),
                ),
            ),
            (
                "agglomerative_jaccard",
                (
                    lambda method_output_dir: AgglomerativeJaccardClusterer(
                        output_dir=method_output_dir,
                        config=AgglomerativeJaccardConfig(
                            linkage=options.agglomerative_linkage,
                            distance_threshold=options.agglomerative_distance_threshold,
                            min_cluster_size=options.agglomerative_min_cluster_size,
                        ),
                    ),
                    lambda clusterer: (
                        f"linkage={clusterer.config.linkage};"
                        f"distance_threshold={clusterer.config.distance_threshold};"
                        f"min_cluster_size={clusterer.config.min_cluster_size}"
                    ),
                ),
            ),
        ]

        for method_name, (build_clusterer, build_parameter_label) in method_specs:
            method_output_dir = condition_output_dir / method_name
            cluster_assignments_path = method_output_dir / "cluster_assignments.tsv"
            medoid_manifest_path = method_output_dir / "medoid_manifest.tsv"
            clusterer = build_clusterer(method_output_dir)
            parameter_label = build_parameter_label(clusterer)

            if matrix_summary.formal_clustering_allowed:
                clustering_input = np.asarray(main_matrix.matrix, dtype=np.uint8)
                clustering_result = clusterer.cluster(clustering_input, main_matrix.pose_ids)
                distance_matrix = clusterer.compute_jaccard_distances(clustering_input)
                cluster_assignment_rows = build_cluster_assignment_rows(
                    condition.condition_id,
                    main_matrix.pose_ids,
                    clustering_result.cluster_labels,
                    distance_matrix=distance_matrix,
                    medoids=clustering_result.medoids,
                )
                medoid_rows = build_medoid_rows(
                    condition.condition_id,
                    main_matrix.pose_ids,
                    clustering_result.medoids,
                    clustering_result.medoid_distance_sums,
                )
            else:
                clustering_result = _empty_clustering_result()
                cluster_assignment_rows = []
                medoid_rows = []

            _write_tsv_rows(
                cluster_assignments_path,
                [
                    "pose_id",
                    "condition_id",
                    "cluster_id",
                    "cluster_member_flag",
                    "noise_flag",
                    "distance_to_cluster_representative",
                ],
                cluster_assignment_rows,
            )
            _write_tsv_rows(
                medoid_manifest_path,
                [
                    "condition_id",
                    "cluster_id",
                    "medoid_pose_id",
                    "medoid_structure_path",
                    "medoid_ifp_distance_sum",
                ],
                medoid_rows,
            )

            condition_summary = build_condition_cluster_summary(
                condition.condition_id,
                clustering_result,
                n_qc_pass_poses=(
                    condition.n_qc_pass_poses
                    if condition.n_qc_pass_poses is not None
                    else len(condition.batch.results)
                ),
                n_ifp_success=(
                    len(condition.contact_eligibilities)
                    if condition.contact_eligibilities
                    else len(condition.batch.results)
                ),
                contact_eligibilities=list(condition.contact_eligibilities),
            )
            pilot_method_summary_rows.append(
                {
                    "method": method_name,
                    "parameter_label": parameter_label,
                    **matrix_summary.to_row(),
                    **condition_summary.to_row(),
                    "cluster_assignments_tsv": str(cluster_assignments_path),
                    "medoid_manifest_tsv": str(medoid_manifest_path),
                }
            )
            method_entries.append(
                {
                    "method": method_name,
                    "parameter_label": parameter_label,
                    "clustering_status": matrix_summary.clustering_status,
                    "formal_clustering_allowed": matrix_summary.formal_clustering_allowed,
                    "n_clusters": clustering_result.n_clusters,
                    "n_noise": clustering_result.n_outliers,
                    "cluster_assignments_tsv": str(cluster_assignments_path),
                    "medoid_manifest_tsv": str(medoid_manifest_path),
                }
            )

        condition_matrix_entries.append(
            {
                **matrix_summary.to_row(),
                "raw_matrix_path": str(raw_matrix_path),
                "main_matrix_path": str(main_matrix_path),
                "pilot_method_entries": method_entries,
            }
        )

    pilot_method_summary_path = pilot_output_dir / "pilot_clustering_method_summary.tsv"
    _write_tsv_rows(
        pilot_method_summary_path,
        [
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
        ],
        pilot_method_summary_rows,
    )

    return {
        "label": options.label,
        "output_dir": str(pilot_output_dir),
        "pilot_ifp_interaction_type_prevalence_tsv": str(interaction_type_prevalence_path),
        "pilot_ifp_feature_prevalence_tsv": str(feature_prevalence_path),
        "pilot_feature_selection_json": str(feature_selection_path),
        "pilot_clustering_method_summary_tsv": str(pilot_method_summary_path),
        "selected_main_feature_count": len(selected_main_feature_names),
        "n_conditions": len(pilot_conditions),
        "n_conditions_with_contact_eligible_signal": sum(
            1 for entry in condition_matrix_entries if entry["n_selected_poses"] > 0
        ),
        "n_conditions_formal_clustering_allowed": sum(
            1 for entry in condition_matrix_entries if entry["formal_clustering_allowed"]
        ),
        "n_conditions_insufficient_clusterable_signal": sum(
            1
            for entry in condition_matrix_entries
            if entry["clustering_status"] == options.insufficient_clusterable_signal_label
        ),
        "n_conditions_empty_main_matrix": sum(
            1 for entry in condition_matrix_entries if entry["clustering_status"] == "empty_main_matrix"
        ),
        "condition_matrix_entries": condition_matrix_entries,
    }


def _empty_clustering_result() -> ClusteringResult:
    return ClusteringResult(
        n_clusters=0,
        n_outliers=0,
        outlier_rate=0.0,
        cluster_sizes={},
        cluster_occupancy={},
        outlier_indices=[],
        cluster_labels=np.array([], dtype=int),
        medoids={},
        medoid_distance_sums={},
    )


def _discover_pose_inputs(options: ProductionRunOptions) -> tuple[list[str], list[PoseInputRecord], dict[str, Any]]:
    manifest = discover_work_root(
        options.work_root,
        af3_only=options.af3_only,
        latest_only=options.latest_only,
        include_targets=options.include_targets,
    )

    discovered: list[PoseInputRecord] = []
    for target, model, run_id, ut in manifest.iter_all_uniprot_targets():
        if options.include_proteins and ut.uniprot_id not in options.include_proteins:
            continue

        sample_inputs: list[PoseInputRecord] = []
        for sample in ut.samples:
            for cif_path in sample.cif_paths:
                resolved_cif_path = Path(cif_path).resolve()
                confidence_json_path = _select_confidence_json(
                    Path(cif_path),
                    sample.confidence_json_paths,
                )
                sample_inputs.append(
                    PoseInputRecord(
                        cif_path=resolved_cif_path,
                        pose_id=Path(cif_path).stem,
                        protein_id=ut.uniprot_id,
                        ligand_id=target,
                        model=model,
                        source_run_id=_infer_source_run_id(resolved_cif_path, run_id),
                        discovered_run_id=run_id,
                        confidence_json_path=(
                            Path(confidence_json_path).resolve()
                            if confidence_json_path is not None
                            else None
                        ),
                        run_status=sample.status.value,
                        seed=sample.seed,
                        sample=sample.sample,
                    )
                )

        if not sample_inputs and ut.model_cif_path is not None:
            resolved_cif_path = Path(ut.model_cif_path).resolve()
            confidence_json_path = _select_confidence_json(
                Path(ut.model_cif_path),
                ut.confidence_json_paths,
            )
            sample_inputs.append(
                PoseInputRecord(
                    cif_path=resolved_cif_path,
                    pose_id=Path(ut.model_cif_path).stem,
                    protein_id=ut.uniprot_id,
                    ligand_id=target,
                    model=model,
                    source_run_id=_infer_source_run_id(resolved_cif_path, run_id),
                    discovered_run_id=run_id,
                    confidence_json_path=(
                        Path(confidence_json_path).resolve()
                        if confidence_json_path is not None
                        else None
                    ),
                    run_status="model_cif",
                )
            )

        discovered.extend(sample_inputs)

    discovered.sort(
        key=lambda pose: (
            pose.ligand_id,
            pose.protein_id,
            pose.model,
            pose.source_run_id,
            -1 if pose.seed is None else pose.seed,
            -1 if pose.sample is None else pose.sample,
            str(pose.cif_path),
        )
    )
    if options.max_cases is not None:
        discovered = discovered[: options.max_cases]

    return manifest.errors, discovered, manifest.summary


def _prepare_pose_case(
    pose: PoseInputRecord,
    *,
    index: int,
    output_dir: Path,
) -> tuple[dict[str, Any], PreparedPose | None]:
    case_dir = output_dir / "cases" / f"{index:04d}_{_slug(pose.pose_id)}"
    case_dir.mkdir(parents=True, exist_ok=True)
    case: dict[str, Any] = {
        "index": index,
        "pose_id": pose.pose_id,
        "protein_id": pose.protein_id,
        "ligand_id": pose.ligand_id,
        "model": pose.model,
        "source_run_id": pose.source_run_id,
        "discovered_run_id": pose.discovered_run_id,
        "seed": pose.seed,
        "sample": pose.sample,
        "cif_path": str(pose.cif_path),
        "case_dir": str(case_dir),
        "status": "preparing",
    }

    try:
        normalize_dir = case_dir / "normalize"
        normalize_ok, normalized_path = NormalizeMMCIFRunner(pose.cif_path, normalize_dir).run()
        case["normalize_report_path"] = str(normalize_dir / "normalize_report.json")
        if not normalize_ok or normalized_path is None:
            raise RuntimeError("Normalization failed")

        normalized_path = Path(normalized_path).resolve()
        posebusters_dir = case_dir / "posebusters_input"
        posebusters_ok, posebusters_pdb = convert_cif_to_pdb(normalized_path, posebusters_dir)
        case["cif_to_pdb_report_path"] = str(posebusters_dir / "cif_to_pdb_report.json")
        if not posebusters_ok or posebusters_pdb is None:
            raise RuntimeError("cif_to_pdb failed")

        posebusters_pdb = Path(posebusters_pdb).resolve()
        privateer_input_cif = Path(
            prepare_privateer_input(normalized_path, case_dir / "privateer_input.cif")
        ).resolve()
        structure = gemmi.read_structure(str(normalized_path))

        prepared = PreparedPose(
            pose=pose,
            case_dir=case_dir,
            normalized_cif=normalized_path,
            posebusters_pdb=posebusters_pdb,
            privateer_input_cif=privateer_input_cif,
            structure=structure,
        )
        case.update(
            {
                "status": "prepared",
                "normalized_cif": str(normalized_path),
                "posebusters_pdb": str(posebusters_pdb),
                "privateer_input_cif": str(privateer_input_cif),
            }
        )
        return case, prepared
    except Exception as exc:
        case.update(
            {
                "status": "prep_error",
                "error": f"{exc.__class__.__name__}: {exc}",
                "traceback": traceback.format_exc(),
            }
        )
        return case, None


def _prepare_pose_case_worker(args: tuple[int, PoseInputRecord, Path]) -> tuple[int, dict[str, Any], dict[str, str] | None]:
    index, pose, output_dir = args
    case, prepared = _prepare_pose_case(pose, index=index, output_dir=output_dir)
    if prepared is None:
        return index, case, None
    return (
        index,
        case,
        {
            "case_dir": str(prepared.case_dir),
            "normalized_cif": str(prepared.normalized_cif),
            "posebusters_pdb": str(prepared.posebusters_pdb),
            "privateer_input_cif": str(prepared.privateer_input_cif),
        },
    )


def _rebuild_prepared_pose(
    pose: PoseInputRecord,
    payload: dict[str, str],
) -> PreparedPose:
    normalized_cif = Path(payload["normalized_cif"]).resolve()
    return PreparedPose(
        pose=pose,
        case_dir=Path(payload["case_dir"]).resolve(),
        normalized_cif=normalized_cif,
        posebusters_pdb=Path(payload["posebusters_pdb"]).resolve(),
        privateer_input_cif=Path(payload["privateer_input_cif"]).resolve(),
        structure=gemmi.read_structure(str(normalized_cif)),
    )


def _prepare_pose_cases(
    pose_inputs: list[PoseInputRecord],
    *,
    output_dir: Path,
    n_jobs: int,
) -> tuple[list[dict[str, Any]], list[PreparedPose]]:
    if n_jobs <= 1 or len(pose_inputs) <= 1:
        cases: list[dict[str, Any]] = []
        prepared_poses: list[PreparedPose] = []
        for index, pose in enumerate(pose_inputs, start=1):
            case, prepared_pose = _prepare_pose_case(pose, index=index, output_dir=output_dir)
            cases.append(case)
            if prepared_pose is not None:
                prepared_poses.append(prepared_pose)
        return cases, prepared_poses

    indexed_poses = list(enumerate(pose_inputs, start=1))
    cases_by_index: dict[int, dict[str, Any]] = {}
    prepared_payloads_by_index: dict[int, dict[str, str]] = {}
    max_workers = min(max(1, n_jobs), len(indexed_poses))

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(_prepare_pose_case_worker, (index, pose, output_dir)): (index, pose)
            for index, pose in indexed_poses
        }
        for future in as_completed(futures):
            index, pose = futures[future]
            try:
                result_index, case, prepared_payload = future.result()
            except Exception as exc:
                case_dir = output_dir / "cases" / f"{index:04d}_{_slug(pose.pose_id)}"
                case = {
                    "index": index,
                    "pose_id": pose.pose_id,
                    "protein_id": pose.protein_id,
                    "ligand_id": pose.ligand_id,
                    "model": pose.model,
                    "source_run_id": pose.source_run_id,
                    "discovered_run_id": pose.discovered_run_id,
                    "seed": pose.seed,
                    "sample": pose.sample,
                    "cif_path": str(pose.cif_path),
                    "case_dir": str(case_dir),
                    "status": "prep_error",
                    "error": f"{exc.__class__.__name__}: {exc}",
                    "traceback": traceback.format_exc(),
                }
                result_index = index
                prepared_payload = None
            cases_by_index[result_index] = case
            if prepared_payload is not None:
                prepared_payloads_by_index[result_index] = prepared_payload

    cases = [cases_by_index[index] for index, _pose in indexed_poses]
    prepared_poses = [
        _rebuild_prepared_pose(pose, prepared_payloads_by_index[index])
        for index, pose in indexed_poses
        if index in prepared_payloads_by_index
    ]
    return cases, prepared_poses


def _hard_qc_input(prepared_pose: PreparedPose) -> HardQCInput:
    return HardQCInput(
        pose_id=prepared_pose.pose.pose_id,
        mol_pred_path=prepared_pose.posebusters_pdb,
        structure=prepared_pose.structure,
        privateer_cif_path=prepared_pose.privateer_input_cif,
    )


def _write_pose_manifest_tsv(
    pose_inputs: list[PoseInputRecord],
    *,
    output_path: Path,
    options: ProductionRunOptions,
    case_by_pose_id: dict[str, dict[str, Any]],
) -> None:
    rows: list[dict[str, Any]] = []
    for pose in pose_inputs:
        case = case_by_pose_id.get(pose.pose_id, {})
        substrate_type, dp = _parse_target_metadata(pose.ligand_id)
        condition_id = _condition_id_for_pose(pose, options.construct_type)
        rows.append(
            {
                "analysis_id": options.run_id,
                "pose_id": pose.pose_id,
                "condition_id": condition_id,
                "protein_id": pose.protein_id,
                "construct_type": options.construct_type,
                "ligand_id": pose.ligand_id,
                "substrate_class": substrate_type,
                "dp": dp,
                "model": pose.model,
                "source_run_id": pose.source_run_id,
                "discovered_run_id": pose.discovered_run_id or "",
                "seed": "" if pose.seed is None else pose.seed,
                "sample_index": "" if pose.sample is None else pose.sample,
                "raw_structure_path": str(pose.cif_path),
                "confidence_json_path": str(pose.confidence_json_path or ""),
                "run_status": pose.run_status,
                "parse_status": case.get("status", "not_prepared"),
                "qc_status": (case.get("qc_verdict") or {}).get("status", ""),
                "analysis_status": case.get("analysis_status", ""),
                "case_dir": case.get("case_dir", ""),
                "normalized_cif": case.get("normalized_cif", ""),
            }
        )

    _write_tsv_rows(
        output_path,
        [
            "analysis_id",
            "pose_id",
            "condition_id",
            "protein_id",
            "construct_type",
            "ligand_id",
            "substrate_class",
            "dp",
            "model",
            "source_run_id",
            "discovered_run_id",
            "seed",
            "sample_index",
            "raw_structure_path",
            "confidence_json_path",
            "run_status",
            "parse_status",
            "qc_status",
            "analysis_status",
            "case_dir",
            "normalized_cif",
        ],
        rows,
    )


def _write_pose_confidence_tsv(
    pose_inputs: list[PoseInputRecord],
    *,
    output_path: Path,
    options: ProductionRunOptions,
) -> None:
    rows: list[dict[str, Any]] = []
    for pose in pose_inputs:
        substrate_type, dp = _parse_target_metadata(pose.ligand_id)
        confidence = _read_confidence_summary(pose.confidence_json_path)
        rows.append(
            {
                "pose_id": pose.pose_id,
                "condition_id": _condition_id_for_pose(pose, options.construct_type),
                "protein_id": pose.protein_id,
                "construct_type": options.construct_type,
                "ligand_id": pose.ligand_id,
                "substrate_class": substrate_type,
                "dp": dp,
                "seed": "" if pose.seed is None else pose.seed,
                "sample_index": "" if pose.sample is None else pose.sample,
                "confidence_json_path": str(pose.confidence_json_path or ""),
                **confidence,
            }
        )

    _write_tsv_rows(
        output_path,
        [
            "pose_id",
            "condition_id",
            "protein_id",
            "construct_type",
            "ligand_id",
            "substrate_class",
            "dp",
            "seed",
            "sample_index",
            "confidence_json_path",
            "confidence_json_status",
            "ranking_score",
            "iptm",
            "ptm",
            "mean_plddt",
            "ligand_interface_confidence",
            "pae_summary",
        ],
        rows,
    )


def _write_structure_index_tsv(
    pose_inputs: list[PoseInputRecord],
    *,
    output_path: Path,
    options: ProductionRunOptions,
    case_by_pose_id: dict[str, dict[str, Any]],
) -> None:
    rows: list[dict[str, Any]] = []
    for pose in pose_inputs:
        case = case_by_pose_id.get(pose.pose_id, {})
        rows.append(
            {
                "pose_id": pose.pose_id,
                "condition_id": _condition_id_for_pose(pose, options.construct_type),
                "raw_structure_path": str(pose.cif_path),
                "normalized_cif": case.get("normalized_cif", ""),
                "posebusters_pdb": case.get("posebusters_pdb", ""),
                "privateer_input_cif": case.get("privateer_input_cif", ""),
                "complex_h_pdb": case.get("complex_h_pdb", ""),
                "ligand_mol2": case.get("ligand_mol2", ""),
                "case_dir": case.get("case_dir", ""),
            }
        )

    _write_tsv_rows(
        output_path,
        [
            "pose_id",
            "condition_id",
            "raw_structure_path",
            "normalized_cif",
            "posebusters_pdb",
            "privateer_input_cif",
            "complex_h_pdb",
            "ligand_mol2",
            "case_dir",
        ],
        rows,
    )


def _write_qc_attrition_tsv(
    pose_inputs: list[PoseInputRecord],
    *,
    output_path: Path,
    options: ProductionRunOptions,
    case_by_pose_id: dict[str, dict[str, Any]],
) -> None:
    grouped: dict[str, dict[str, Any]] = {}
    for pose in pose_inputs:
        condition_id = _condition_id_for_pose(pose, options.construct_type)
        row = grouped.setdefault(
            condition_id,
            {
                "condition_id": condition_id,
                "protein_id": pose.protein_id,
                "construct_type": options.construct_type,
                "ligand_id": pose.ligand_id,
                "substrate_class": _parse_target_metadata(pose.ligand_id)[0],
                "dp": _parse_target_metadata(pose.ligand_id)[1],
                "n_generated": 0,
                "n_prepared": 0,
                "n_prep_error": 0,
                "n_stage1_hard_fail": 0,
                "n_stage1_soft_flag": 0,
                "n_stage1_pass": 0,
                "hard_fail_reason_counts": {},
                "soft_flag_reason_counts": {},
            },
        )
        row["n_generated"] += 1

        case = case_by_pose_id.get(pose.pose_id, {})
        if case.get("status") == "prep_error":
            row["n_prep_error"] += 1
            row["n_stage1_hard_fail"] += 1
            reason_counts = row["hard_fail_reason_counts"]
            reason_counts["prep_error"] = reason_counts.get("prep_error", 0) + 1
            continue
        if case.get("status") == "prepared":
            row["n_prepared"] += 1

        verdict = case.get("qc_verdict") or {}
        status = verdict.get("status")
        if status == "dropped":
            row["n_stage1_hard_fail"] += 1
            reason_counts = row["hard_fail_reason_counts"]
            for reason in verdict.get("drop_reasons", []) or ["unknown"]:
                reason_counts[reason] = reason_counts.get(reason, 0) + 1
        elif status == "flagged":
            row["n_stage1_soft_flag"] += 1
            row["n_stage1_pass"] += 1
            reason_counts = row["soft_flag_reason_counts"]
            for reason in verdict.get("warnings", []) or ["unknown"]:
                reason_counts[reason] = reason_counts.get(reason, 0) + 1
        elif status == "passed":
            row["n_stage1_pass"] += 1

    rows: list[dict[str, Any]] = []
    for row in grouped.values():
        n_generated = int(row["n_generated"])
        rows.append(
            {
                **row,
                "hard_fail_rate": (
                    int(row["n_stage1_hard_fail"]) / n_generated
                    if n_generated
                    else 0.0
                ),
                "hard_fail_reason_counts": json.dumps(
                    row["hard_fail_reason_counts"],
                    sort_keys=True,
                ),
                "soft_flag_reason_counts": json.dumps(
                    row["soft_flag_reason_counts"],
                    sort_keys=True,
                ),
            }
        )

    _write_tsv_rows(
        output_path,
        [
            "condition_id",
            "protein_id",
            "construct_type",
            "ligand_id",
            "substrate_class",
            "dp",
            "n_generated",
            "n_prepared",
            "n_prep_error",
            "n_stage1_hard_fail",
            "n_stage1_soft_flag",
            "n_stage1_pass",
            "hard_fail_rate",
            "hard_fail_reason_counts",
            "soft_flag_reason_counts",
        ],
        rows,
    )


def _summary_geometry_row(
    pose: PoseInputRecord,
    metrics: PoseGeometryMetrics | None,
) -> dict[str, Any]:
    base = {
        "pose_id": pose.pose_id,
        "protein_id": pose.protein_id,
        "ligand_id": pose.ligand_id,
        "model": pose.model,
    }
    if metrics is None:
        return base
    return {
        **base,
        **metrics.to_legacy_geometry_dict(),
        **metrics.to_row(),
    }


def _choose_crystal_representative_pose_id(
    batch: IFPBatch,
    clustering_result: ClusteringResult,
    cluster_pose_ids: list[str],
) -> str | None:
    if clustering_result.medoids and clustering_result.cluster_sizes and cluster_pose_ids:
        best_cluster_id = max(
            sorted(clustering_result.cluster_sizes),
            key=lambda cluster_id: clustering_result.cluster_sizes[cluster_id],
        )
        medoid_index = clustering_result.medoids.get(best_cluster_id)
        if medoid_index is not None and 0 <= medoid_index < len(cluster_pose_ids):
            return cluster_pose_ids[medoid_index]

    for result in batch.results:
        if result.status == "ok":
            return result.pose_id
    return None


def run_analysis_core(
    config: dict[str, Any],
    output_dir: Path,
    del_variant: str,
    n_jobs: int | None = None,
) -> AnalysisCoreResult:
    run_start = perf_counter()
    timing_events: list[dict[str, Any]] = []

    def _record_timing(
        *,
        step: str,
        start: float,
        pose_id: str = "",
        condition_id: str = "",
        status: str = "ok",
        detail: str = "",
    ) -> None:
        if not options.collect_timing_events:
            return
        timing_events.append(
            {
                "step": step,
                "pose_id": pose_id,
                "condition_id": condition_id,
                "worker_count": options.n_jobs,
                "wall_time_s": round(perf_counter() - start, 6),
                "status": status,
                "detail": detail,
            }
        )

    options = load_production_options(config, output_dir, del_variant, n_jobs=n_jobs)
    options.output_dir.mkdir(parents=True, exist_ok=True)

    discovery_start = perf_counter()
    discovery_errors, pose_inputs, discovery_summary = _discover_pose_inputs(options)
    _record_timing(
        step="analysis.discovery",
        start=discovery_start,
        status="ok" if not discovery_errors else "error",
        detail=f"n_discovered={len(pose_inputs)};n_errors={len(discovery_errors)}",
    )
    analysis_summary: dict[str, Any] = {
        "run_id": options.run_id,
        "output_dir": str(options.output_dir),
        "del_variant": options.del_variant,
        "work_root": str(options.work_root),
        "af3_only": options.af3_only,
        "latest_only": options.latest_only,
        "max_cases": options.max_cases,
        "include_targets": list(options.include_targets),
        "include_proteins": list(options.include_proteins),
        "run_posebusters": options.run_posebusters,
        "run_privateer": options.run_privateer,
        "n_jobs": options.n_jobs,
        "collect_timing_events": options.collect_timing_events,
        "discovery_summary": discovery_summary,
        "discovery_errors": discovery_errors,
        "n_discovered": len(pose_inputs),
        "cases": [],
    }
    if options.clustering_pilot is not None:
        analysis_summary["clustering_pilot"] = {
            "enabled": True,
            "label": options.clustering_pilot.label,
            "output_dirname": options.clustering_pilot.output_dirname,
        }

    prepare_start = perf_counter()
    cases, prepared_poses = _prepare_pose_cases(
        pose_inputs,
        output_dir=options.output_dir,
        n_jobs=options.n_jobs,
    )
    _record_timing(
        step="analysis.prepare_cases",
        start=prepare_start,
        status="ok" if prepared_poses else "error",
        detail=f"n_inputs={len(pose_inputs)};n_prepared={len(prepared_poses)}",
    )
    analysis_summary["cases"].extend(cases)
    case_by_pose_id: dict[str, dict[str, Any]] = {
        str(case["pose_id"]): case for case in cases
    }

    prepared_pose_by_id = {
        prepared_pose.pose.pose_id: prepared_pose
        for prepared_pose in prepared_poses
    }

    analysis_summary["n_prepared"] = len(prepared_poses)
    analysis_summary["n_prep_errors"] = sum(
        1 for case in analysis_summary["cases"] if case["status"] == "prep_error"
    )

    summary_path = options.output_dir / "analysis_core_summary.json"
    qc_report_path: Path | None = None
    pose_manifest_tsv_path: Path | None = None
    pose_confidence_tsv_path: Path | None = None
    structure_index_tsv_path: Path | None = None
    qc_attrition_tsv_path: Path | None = None
    pose_geometry_tsv_path: Path | None = None
    pose_ifp_table_tsv_path: Path | None = None
    pose_residue_contact_tsv_path: Path | None = None
    pose_convergence_tsv_path: Path | None = None
    condition_convergence_summary_tsv_path: Path | None = None
    cluster_assignments_tsv_path: Path | None = None
    medoid_manifest_tsv_path: Path | None = None
    condition_cluster_summary_tsv_path: Path | None = None
    cluster_ifp_signature_tsv_path: Path | None = None
    cluster_residue_signature_tsv_path: Path | None = None
    cluster_signatures_json_path: Path | None = None
    cluster_annotation_stage_completed = False
    protein_condition_residue_scores_tsv_path: Path | None = None
    protein_residue_regio_delta_tsv_path: Path | None = None
    condition_patch_summary_tsv_path: Path | None = None
    protein_patch_summary_tsv_path: Path | None = None
    crystal_anchor_tsv_path: Path | None = None
    metrics_csv_path: Path | None = None
    summary_json_path: Path | None = None
    report_html_path: Path | None = None
    n_analyzed = 0
    n_cluster_conditions = 0
    crystal_anchoring_stage_completed = False
    n_crystal_anchoring_conditions = 0
    n_crystal_anchoring_errors = 0

    if prepared_poses:
        qc_report_path = options.output_dir / "qc_report.json"
        hard_qc_start = perf_counter()
        report = run_hard_qc(
            [_hard_qc_input(prepared_pose) for prepared_pose in prepared_poses],
            run_id=options.run_id,
            run_posebusters=options.run_posebusters,
            run_privateer=options.run_privateer,
            max_workers=options.n_jobs,
            collect_timing_events=options.collect_timing_events,
        )
        _record_timing(
            step="analysis.hard_qc",
            start=hard_qc_start,
            status="ok",
            detail=f"n_poses={len(prepared_poses)}",
        )
        if options.collect_timing_events:
            timing_events.extend(
                {
                    **event,
                    "worker_count": event.get("worker_count", options.n_jobs),
                }
                for event in report.timing_events
            )
        write_qc_report(report, qc_report_path)

        qc_report_payload = json.loads(qc_report_path.read_text())
        validate(instance=qc_report_payload, schema=_load_qc_report_schema())
        analysis_summary["qc_report_path"] = str(qc_report_path)
        analysis_summary["qc_counts"] = {
            "total": report.total,
            "passed": report.passed,
            "flagged": report.flagged,
            "dropped": report.dropped,
        }

        verdict_by_pose = {
            verdict["pose_id"]: verdict
            for verdict in qc_report_payload.get("verdicts", [])
        }

        prolif_config = load_prolif_features_config()
        contact_eligibility_rule = load_contact_eligibility_rule()
        geometry_metrics: list[PoseGeometryMetrics] = []
        summary_geometry_rows: list[dict[str, Any]] = []
        metrics_record_by_pose_id: dict[str, dict[str, Any]] = {}
        all_ifp_results: list[IFPResult] = []
        residue_contact_rows: list[dict[str, Any]] = []
        contact_eligibility_by_condition: dict[str, list[ContactEligibility]] = {}
        ifp_pose_inputs_by_condition: dict[str, list[dict[str, Any]]] = {}
        n_qc_pass_by_condition: dict[str, int] = {}
        normalized_path_by_pose_id: dict[str, str] = {}

        for prepared_pose in prepared_poses:
            pose = prepared_pose.pose
            case = case_by_pose_id[pose.pose_id]
            verdict = verdict_by_pose.get(
                pose.pose_id,
                {"pose_id": pose.pose_id, "status": "unknown", "warnings": []},
            )
            case["qc_verdict"] = verdict
            condition_id = _condition_id_for_pose(pose, options.construct_type)
            case["condition_id"] = condition_id
            normalized_path_by_pose_id[pose.pose_id] = str(prepared_pose.normalized_cif)

            geometry_row: dict[str, Any] = {}
            metrics: PoseGeometryMetrics | None = None
            if verdict.get("status") != "dropped":
                n_qc_pass_by_condition[condition_id] = n_qc_pass_by_condition.get(condition_id, 0) + 1
                geometry_start = perf_counter()
                try:
                    metrics = compute_pose_metrics_from_structure(
                        structure=prepared_pose.structure,
                        pose_id=pose.pose_id,
                        model=pose.model,
                        protein_id=pose.protein_id,
                        ligand_id=pose.ligand_id,
                    )
                    geometry_metrics.append(metrics)
                    geometry_row = metrics.to_row()
                    case["analysis_status"] = "analyzed"
                    case["geometry_row"] = geometry_row
                    n_analyzed += 1
                    _record_timing(
                        step="analysis.geometry_metrics",
                        start=geometry_start,
                        pose_id=pose.pose_id,
                        condition_id=condition_id,
                        status="ok",
                    )
                except Exception as exc:
                    case["analysis_status"] = "analysis_error"
                    case["analysis_error"] = f"{exc.__class__.__name__}: {exc}"
                    case["analysis_traceback"] = traceback.format_exc()
                    _record_timing(
                        step="analysis.geometry_metrics",
                        start=geometry_start,
                        pose_id=pose.pose_id,
                        condition_id=condition_id,
                        status="error",
                        detail=f"{exc.__class__.__name__}: {exc}",
                    )

                protonation_dir = prepared_pose.case_dir / "protonated"
                protonation_start = perf_counter()
                protonate_ok, protonation_report = protonate_and_export(
                    prepared_pose.normalized_cif,
                    protonation_dir,
                )
                _record_timing(
                    step="analysis.protonation_export",
                    start=protonation_start,
                    pose_id=pose.pose_id,
                    condition_id=condition_id,
                    status="ok" if protonate_ok else "failed",
                    detail=";".join((protonation_report or {}).get("blockers", []) or []),
                )
                case["protonate_ok"] = bool(protonate_ok)
                case["protonation_report_path"] = str(protonation_dir / "protonation_report.json")
                case["protonation_blockers"] = (protonation_report or {}).get("blockers", []) if protonation_report else []
                case["protonation_warnings"] = (protonation_report or {}).get("warnings", []) if protonation_report else []
                case["complex_h_pdb"] = str(protonation_dir / "complex_H.pdb")
                case["ligand_mol2"] = str(protonation_dir / "ligand_for_prolif.mol2")
                if protonate_ok:
                    ifp_pose_inputs_by_condition.setdefault(condition_id, []).append(
                        {
                            "pose_id": pose.pose_id,
                            "complex_pdb": protonation_dir / "complex_H.pdb",
                            "ligand_mol2": protonation_dir / "ligand_for_prolif.mol2",
                            "protein_id": pose.protein_id,
                            "ligand_id": pose.ligand_id,
                            "model": pose.model,
                            "seed": pose.seed,
                            "sample": pose.sample,
                        }
                    )
                else:
                    failed_ifp = _make_ifp_failure_result(
                        pose.pose_id,
                        "Protonation blockers: " + ", ".join(case["protonation_blockers"] or ["unknown"]),
                    )
                    all_ifp_results.append(failed_ifp)
                    case["ifp_generation_status"] = failed_ifp.status
                    case["ifp_error"] = failed_ifp.error
            else:
                case["analysis_status"] = "skipped_dropped"
                if options.collect_timing_events:
                    timing_events.append(
                        {
                            "step": "analysis.geometry_metrics",
                            "pose_id": pose.pose_id,
                            "condition_id": condition_id,
                            "worker_count": options.n_jobs,
                            "wall_time_s": 0.0,
                            "status": "skipped",
                            "detail": "qc_dropped",
                        }
                    )
                    timing_events.append(
                        {
                            "step": "analysis.protonation_export",
                            "pose_id": pose.pose_id,
                            "condition_id": condition_id,
                            "worker_count": options.n_jobs,
                            "wall_time_s": 0.0,
                            "status": "skipped",
                            "detail": "qc_dropped",
                        }
                    )

            summary_geometry_rows.append(_summary_geometry_row(pose, metrics))
            substrate_type, dp = _parse_target_metadata(pose.ligand_id)
            metrics_record_by_pose_id[pose.pose_id] = merge_pose_record(
                pose_id=pose.pose_id,
                run_id=options.run_id,
                protein_id=pose.protein_id,
                ligand_id=pose.ligand_id,
                model=pose.model,
                seed=0 if pose.seed is None else pose.seed,
                qc_verdict=verdict,
                geometry=geometry_row,
                del_variant=options.del_variant,
                substrate_type=substrate_type,
                dp=dp,
                cbm_present=options.del_variant == "b",
            )

        cluster_assignment_rows: list[dict[str, Any]] = []
        medoid_rows: list[dict[str, Any]] = []
        condition_summary_rows: list[dict[str, Any]] = []
        cluster_results_for_summary: list[dict[str, Any]] = []
        crystal_reports_for_summary: list[dict[str, Any]] = []
        crystal_report_entries: list[dict[str, Any]] = []
        crystal_report_errors: list[dict[str, Any]] = []
        crystal_anchor_rows: list[dict[str, Any]] = []
        residue_contact_rows: list[dict[str, Any]] = []
        pose_convergence_metrics: list[PoseConvergenceMetrics] = []
        condition_convergence_summaries: list[ConditionConvergenceSummary] = []
        pilot_conditions: list[PilotConditionIFP] = []

        clusterer = HDBSCANClusterer(output_dir=options.output_dir)
        ifp_matrix_root = options.output_dir / "ifp_matrices"
        for condition_id in sorted(n_qc_pass_by_condition):
            condition_pose_inputs = ifp_pose_inputs_by_condition.get(condition_id, [])
            clustering_result = _empty_clustering_result()
            if condition_pose_inputs:
                convergence_start = perf_counter()
                convergence_metrics, convergence_summary = compute_condition_convergence(
                    [
                        ConvergencePoseInput(
                            pose_id=str(entry["pose_id"]),
                            condition_id=condition_id,
                            complex_pdb=Path(entry["complex_pdb"]),
                            seed=int(entry["seed"]) if entry.get("seed") is not None else None,
                            sample=int(entry["sample"]) if entry.get("sample") is not None else None,
                        )
                        for entry in condition_pose_inputs
                    ]
                )
                _record_timing(
                    step="analysis.convergence",
                    start=convergence_start,
                    condition_id=condition_id,
                    status="ok",
                    detail=f"n_poses={len(condition_pose_inputs)}",
                )
                pose_convergence_metrics.extend(convergence_metrics)
                condition_convergence_summaries.append(convergence_summary)
                for metric in convergence_metrics:
                    case_by_pose_id[metric.pose_id]["ligand_rmsd_to_reference"] = metric.ligand_rmsd_to_reference
                    case_by_pose_id[metric.pose_id]["convergent_flag"] = metric.convergent_flag
                    case_by_pose_id[metric.pose_id]["convergence_reference_pose_id"] = metric.reference_pose_id

                first_input = condition_pose_inputs[0]
                ifp_start = perf_counter()
                batch = compute_ifp_batch(
                    [
                        {
                            "pose_id": entry["pose_id"],
                            "complex_pdb": entry["complex_pdb"],
                            "ligand_mol2": entry["ligand_mol2"],
                        }
                        for entry in condition_pose_inputs
                    ],
                    protein_id=str(first_input["protein_id"]),
                    ligand_id=str(first_input["ligand_id"]),
                    model=str(first_input["model"]),
                    max_workers=options.n_jobs,
                )
                _record_timing(
                    step="analysis.prolif_ifp_batch",
                    start=ifp_start,
                    condition_id=condition_id,
                    status="ok",
                    detail=f"n_poses={len(condition_pose_inputs)}",
                )
                all_ifp_results.extend(batch.results)

                condition_ifp_dir = ifp_matrix_root / _slug(condition_id)
                write_ifp_matrix(batch, condition_ifp_dir / "ifp_matrix.csv")
                cluster_pose_ids: list[str] = []

                eligibility_by_index: dict[int, ContactEligibility] = {}
                for result in batch.results:
                    case = case_by_pose_id[result.pose_id]
                    case["ifp_generation_status"] = result.status
                    case["ifp_error"] = result.error
                    case["n_ifp_contacts"] = result.n_total_contacts
                    metrics_record_by_pose_id[result.pose_id]["n_ifp_contacts"] = result.n_total_contacts

                for index, result in enumerate(batch.results):
                    if result.status not in {"ok", "zero_contacts"}:
                        continue

                    eligibility = evaluate_contact_eligibility(result, contact_eligibility_rule)
                    eligibility_by_index[index] = eligibility
                    contact_eligibility_by_condition.setdefault(condition_id, []).append(eligibility)

                    case = case_by_pose_id[result.pose_id]
                    case["contact_eligible"] = eligibility.eligible
                    case["contact_exclusion_class"] = eligibility.exclusion_class
                    case["n_vdw_interactions"] = eligibility.n_vdw_interactions
                    case["n_non_vdw_interactions"] = eligibility.n_non_vdw_interactions
                    case["n_non_vdw_contact_residues"] = eligibility.n_non_vdw_contact_residues

                    metrics_record_by_pose_id[result.pose_id]["contact_eligible"] = eligibility.eligible
                    metrics_record_by_pose_id[result.pose_id]["contact_exclusion_class"] = eligibility.exclusion_class
                    metrics_record_by_pose_id[result.pose_id]["n_vdw_interactions"] = eligibility.n_vdw_interactions
                    metrics_record_by_pose_id[result.pose_id]["n_non_vdw_interactions"] = eligibility.n_non_vdw_interactions
                    metrics_record_by_pose_id[result.pose_id]["n_non_vdw_contact_residues"] = eligibility.n_non_vdw_contact_residues

                clusterable_indices = [
                    index
                    for index, result in enumerate(batch.results)
                    if result.status == "ok" and eligibility_by_index.get(index) is not None and eligibility_by_index[index].eligible
                ]
                if options.clustering_pilot is not None:
                    pilot_conditions.append(
                        PilotConditionIFP(
                            condition_id=condition_id,
                            batch=batch,
                            selected_pose_ids=tuple(
                                batch.results[index].pose_id for index in clusterable_indices
                            ),
                            contact_eligibilities=tuple(contact_eligibility_by_condition.get(condition_id, [])),
                            n_qc_pass_poses=n_qc_pass_by_condition.get(condition_id, 0),
                        )
                    )
                if clusterable_indices:
                    clustering_start = perf_counter()
                    cluster_matrix = np.asarray([batch.matrix[index] for index in clusterable_indices], dtype=np.uint8)
                    cluster_pose_ids = [batch.results[index].pose_id for index in clusterable_indices]
                    clustering_result = clusterer.cluster(cluster_matrix, cluster_pose_ids)
                    distance_matrix = clusterer.compute_jaccard_distances(cluster_matrix)
                    cluster_assignment_rows.extend(
                        build_cluster_assignment_rows(
                            condition_id,
                            cluster_pose_ids,
                            clustering_result.cluster_labels,
                            distance_matrix=distance_matrix,
                            medoids=clustering_result.medoids,
                        )
                    )
                    medoid_rows.extend(
                        build_medoid_rows(
                            condition_id,
                            cluster_pose_ids,
                            clustering_result.medoids,
                            clustering_result.medoid_distance_sums,
                            structure_paths={
                                pose_id: normalized_path_by_pose_id.get(pose_id, "")
                                for pose_id in cluster_pose_ids
                            },
                        )
                    )
                    for pose_id, label in zip(cluster_pose_ids, clustering_result.cluster_labels, strict=True):
                        metrics_record_by_pose_id[pose_id]["cluster_id"] = int(label)
                        case_by_pose_id[pose_id]["cluster_id"] = int(label)
                    _record_timing(
                        step="analysis.production_clustering",
                        start=clustering_start,
                        condition_id=condition_id,
                        status="ok",
                        detail=f"n_clusterable={len(clusterable_indices)}",
                    )
                else:
                    if options.collect_timing_events:
                        timing_events.append(
                            {
                                "step": "analysis.production_clustering",
                                "pose_id": "",
                                "condition_id": condition_id,
                                "worker_count": options.n_jobs,
                                "wall_time_s": 0.0,
                                "status": "skipped",
                                "detail": "no_clusterable_poses",
                            }
                        )

                representative_pose_id = _choose_crystal_representative_pose_id(
                    batch,
                    clustering_result,
                    cluster_pose_ids,
                )
                representative_pose = prepared_pose_by_id.get(str(representative_pose_id or ""))
                if representative_pose is not None:
                    crystal_output_dir = options.output_dir / "crystal_anchoring" / _slug(condition_id)
                    crystal_report_path = crystal_output_dir / "crystal_reference_screen.json"
                    try:
                        crystal_start = perf_counter()
                        crystal_report = run_crystal_reference_screen(
                            representative_pose.pose.cif_path,
                            protein_id=representative_pose.pose.protein_id,
                            ligand_id=representative_pose.pose.ligand_id,
                            representative_pose_id=representative_pose.pose.pose_id,
                            output_dir=crystal_output_dir,
                        )
                        _record_timing(
                            step="analysis.crystal_anchoring",
                            start=crystal_start,
                            pose_id=representative_pose.pose.pose_id,
                            condition_id=condition_id,
                            status="ok",
                            detail=f"n_comparisons={len(crystal_report.comparisons)}",
                        )
                        write_crystal_reference_screen_report(crystal_report, crystal_report_path)
                        best_tanimoto = max(
                            (
                                comparison.ifp_tanimoto
                                for comparison in crystal_report.comparisons
                                if comparison.ifp_tanimoto is not None
                            ),
                            default=None,
                        )
                        best_pocket_rmsd = min(
                            (
                                comparison.pocket_rmsd
                                for comparison in crystal_report.comparisons
                                if comparison.pocket_rmsd is not None
                            ),
                            default=None,
                        )
                        if best_tanimoto is not None:
                            metrics_record_by_pose_id[representative_pose.pose.pose_id]["ifp_similarity_crystal"] = best_tanimoto
                        if best_pocket_rmsd is not None:
                            metrics_record_by_pose_id[representative_pose.pose.pose_id]["pocket_rmsd_vs_crystal"] = best_pocket_rmsd
                        case_by_pose_id[representative_pose.pose.pose_id]["crystal_anchoring_representative_pose"] = True
                        case_by_pose_id[representative_pose.pose.pose_id]["crystal_anchoring_report_path"] = str(crystal_report_path)
                        case_by_pose_id[representative_pose.pose.pose_id]["crystal_anchoring_best_tanimoto"] = best_tanimoto
                        case_by_pose_id[representative_pose.pose.pose_id]["crystal_anchoring_best_pocket_rmsd"] = best_pocket_rmsd
                        case_by_pose_id[representative_pose.pose.pose_id]["crystal_anchoring_n_comparisons"] = len(crystal_report.comparisons)
                        report_entry = {
                            "condition_id": condition_id,
                            "protein_id": representative_pose.pose.protein_id,
                            "ligand_id": representative_pose.pose.ligand_id,
                            "representative_pose_id": representative_pose.pose.pose_id,
                            "report_path": str(crystal_report_path),
                            "comparison_count": len(crystal_report.comparisons),
                            "best_tanimoto": best_tanimoto,
                            "best_pocket_rmsd": best_pocket_rmsd,
                        }
                        crystal_report_entries.append(report_entry)
                        crystal_reports_for_summary.append(report_entry)
                        representative_cluster_id = case_by_pose_id[
                            representative_pose.pose.pose_id
                        ].get("cluster_id", "")
                        for comparison in crystal_report.comparisons:
                            crystal_anchor_rows.append(
                                {
                                    "condition_id": condition_id,
                                    "protein_id": representative_pose.pose.protein_id,
                                    "ligand_id": representative_pose.pose.ligand_id,
                                    "cluster_id": representative_cluster_id,
                                    "medoid_pose_id": representative_pose.pose.pose_id,
                                    "crystal_reference_id": comparison.pdb_code,
                                    "source_cif": comparison.source_cif,
                                    "prepared_subset_cif": comparison.prepared_subset_cif,
                                    "selected_protein_chain": comparison.selected_protein_chain,
                                    "local_pocket_rmsd": comparison.pocket_rmsd,
                                    "ligand_rmsd_if_comparable": "",
                                    "proximal_sugar_rmsd": "",
                                    "contact_overlap_score": comparison.ifp_tanimoto,
                                    "same_binding_region_flag": comparison.pocket_rmsd_below_threshold,
                                    "same_general_orientation_flag": "",
                                    "ifp_tanimoto": comparison.ifp_tanimoto,
                                    "pocket_residues": json.dumps(comparison.pocket_residues),
                                    "comparison_status": comparison.status,
                                    "notes": comparison.error or "",
                                }
                            )
                    except Exception as exc:
                        if "crystal_start" in locals():
                            _record_timing(
                                step="analysis.crystal_anchoring",
                                start=crystal_start,
                                pose_id=representative_pose.pose.pose_id,
                                condition_id=condition_id,
                                status="error",
                                detail=f"{exc.__class__.__name__}: {exc}",
                            )
                        error_entry = {
                            "condition_id": condition_id,
                            "protein_id": representative_pose.pose.protein_id,
                            "ligand_id": representative_pose.pose.ligand_id,
                            "representative_pose_id": representative_pose.pose.pose_id,
                            "error": f"{exc.__class__.__name__}: {exc}",
                        }
                        crystal_report_errors.append(error_entry)
                        case_by_pose_id[representative_pose.pose.pose_id]["crystal_anchoring_error"] = error_entry["error"]
                        crystal_anchor_rows.append(
                            {
                                "condition_id": condition_id,
                                "protein_id": representative_pose.pose.protein_id,
                                "ligand_id": representative_pose.pose.ligand_id,
                                "cluster_id": case_by_pose_id[
                                    representative_pose.pose.pose_id
                                ].get("cluster_id", ""),
                                "medoid_pose_id": representative_pose.pose.pose_id,
                                "crystal_reference_id": "",
                                "source_cif": "",
                                "prepared_subset_cif": "",
                                "selected_protein_chain": "",
                                "local_pocket_rmsd": "",
                                "ligand_rmsd_if_comparable": "",
                                "proximal_sugar_rmsd": "",
                                "contact_overlap_score": "",
                                "same_binding_region_flag": "",
                                "same_general_orientation_flag": "",
                                "ifp_tanimoto": "",
                                "pocket_residues": "[]",
                                "comparison_status": "error",
                                "notes": error_entry["error"],
                            }
                        )

            condition_summary = build_condition_cluster_summary(
                condition_id,
                clustering_result,
                n_qc_pass_poses=n_qc_pass_by_condition[condition_id],
                n_ifp_success=len(contact_eligibility_by_condition.get(condition_id, [])),
                contact_eligibilities=contact_eligibility_by_condition.get(condition_id, []),
            )
            condition_summary_rows.append(condition_summary.to_row())
            cluster_results_for_summary.append(
                {
                    **condition_summary.to_row(),
                    "outlier_rate": condition_summary.noise_fraction,
                }
            )

        if options.clustering_pilot is not None and pilot_conditions:
            pilot_outputs_start = perf_counter()
            analysis_summary["clustering_pilot"].update(
                _write_clustering_pilot_outputs(
                    options.clustering_pilot,
                    options.output_dir,
                    pilot_conditions,
                )
            )
            _record_timing(
                step="analysis.clustering_pilot_outputs",
                start=pilot_outputs_start,
                status="ok",
                detail=f"n_conditions={len(pilot_conditions)}",
            )

        if all_ifp_results:
            pose_ifp_table_tsv_path = options.output_dir / "pose_ifp_table.tsv"
            write_pose_ifp_table(all_ifp_results, pose_ifp_table_tsv_path)

            pose_metadata_by_id = {
                prepared_pose.pose.pose_id: {
                    "protein_id": prepared_pose.pose.protein_id,
                    "ligand_id": prepared_pose.pose.ligand_id,
                    "condition_id": str(
                        case_by_pose_id[prepared_pose.pose.pose_id].get("condition_id", "")
                    ),
                }
                for prepared_pose in prepared_poses
            }
            residue_contact_rows = build_pose_residue_contact_rows(
                all_ifp_results,
                pose_metadata_by_id,
            )
            pose_residue_contact_tsv_path = options.output_dir / "pose_residue_contact_table.tsv"
            write_pose_residue_contact_table(
                residue_contact_rows,
                pose_residue_contact_tsv_path,
            )

        if pose_convergence_metrics:
            pose_convergence_tsv_path = options.output_dir / "pose_convergence.tsv"
            write_pose_convergence_tsv(pose_convergence_metrics, pose_convergence_tsv_path)

            condition_convergence_summary_tsv_path = options.output_dir / "condition_convergence_summary.tsv"
            write_condition_convergence_summary_tsv(
                condition_convergence_summaries,
                condition_convergence_summary_tsv_path,
            )

        if n_qc_pass_by_condition:
            cluster_assignments_tsv_path = options.output_dir / "cluster_assignments.tsv"
            _write_tsv_rows(
                cluster_assignments_tsv_path,
                [
                    "pose_id",
                    "condition_id",
                    "cluster_id",
                    "cluster_member_flag",
                    "noise_flag",
                    "distance_to_cluster_representative",
                ],
                cluster_assignment_rows,
            )

            medoid_manifest_tsv_path = options.output_dir / "medoid_manifest.tsv"
            _write_tsv_rows(
                medoid_manifest_tsv_path,
                [
                    "condition_id",
                    "cluster_id",
                    "medoid_pose_id",
                    "medoid_structure_path",
                    "medoid_ifp_distance_sum",
                ],
                medoid_rows,
            )

            condition_cluster_summary_tsv_path = options.output_dir / "condition_cluster_summary.tsv"
            _write_tsv_rows(
                condition_cluster_summary_tsv_path,
                [
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
                ],
                condition_summary_rows,
            )
            n_cluster_conditions = len(condition_summary_rows)

            pose_metadata_by_id = {
                prepared_pose.pose.pose_id: {
                    "protein_id": prepared_pose.pose.protein_id,
                    "ligand_id": prepared_pose.pose.ligand_id,
                    "condition_id": str(
                        case_by_pose_id[prepared_pose.pose.pose_id].get("condition_id", "")
                    ),
                }
                for prepared_pose in prepared_poses
            }
            cluster_annotation_start = perf_counter()
            cluster_signature_tables = build_cluster_signature_tables(
                cluster_assignment_rows=cluster_assignment_rows,
                medoid_rows=medoid_rows,
                geometry_rows_by_pose_id={
                    metric.pose_id: metric.to_row()
                    for metric in geometry_metrics
                },
                ifp_results=all_ifp_results,
                residue_contact_rows=residue_contact_rows,
                pose_metadata_by_id=pose_metadata_by_id,
            )
            cluster_ifp_signature_tsv_path = options.output_dir / "cluster_ifp_signature.tsv"
            write_cluster_ifp_signature_table(
                cluster_signature_tables.ifp_signature_rows,
                cluster_ifp_signature_tsv_path,
            )
            cluster_residue_signature_tsv_path = options.output_dir / "cluster_residue_signature.tsv"
            write_cluster_residue_signature_table(
                cluster_signature_tables.residue_signature_rows,
                cluster_residue_signature_tsv_path,
            )
            cluster_signatures_json_path = options.output_dir / "cluster_signatures.json"
            write_cluster_signature_summary_json(
                cluster_signature_tables.cluster_summaries,
                cluster_signatures_json_path,
            )
            cluster_annotation_stage_completed = True
            _record_timing(
                step="analysis.cluster_annotation",
                start=cluster_annotation_start,
                status="ok",
                detail=f"n_cluster_summaries={len(cluster_signature_tables.cluster_summaries)}",
            )

            condition_metadata_by_id: dict[str, dict[str, Any]] = {}
            for prepared_pose in prepared_poses:
                pose_id = prepared_pose.pose.pose_id
                case = case_by_pose_id.get(pose_id, {})
                condition_id = str(case.get("condition_id", ""))
                if not condition_id or condition_id in condition_metadata_by_id:
                    continue
                substrate_class, dp = _parse_target_metadata(prepared_pose.pose.ligand_id)
                condition_metadata_by_id[condition_id] = {
                    "protein_id": prepared_pose.pose.protein_id,
                    "construct_type": options.construct_type,
                    "substrate_class": str(case.get("substrate_class", "")) or substrate_class,
                    "dp": dp,
                }

            residue_importance_start = perf_counter()
            residue_importance_outputs = compute_residue_importance_outputs(
                cluster_signature_tables.residue_signature_rows,
                cluster_signature_tables.cluster_summaries,
                cluster_ifp_signature_rows=cluster_signature_tables.ifp_signature_rows,
                observed_residue_contact_rows=residue_contact_rows,
                condition_metadata_by_id=condition_metadata_by_id,
            )
            protein_condition_residue_scores_tsv_path = (
                options.output_dir / "protein_condition_residue_scores.tsv"
            )
            write_protein_condition_residue_scores(
                residue_importance_outputs.protein_condition_residue_scores,
                protein_condition_residue_scores_tsv_path,
            )
            protein_residue_regio_delta_tsv_path = options.output_dir / "protein_residue_regio_delta.tsv"
            write_protein_residue_regio_delta(
                residue_importance_outputs.protein_residue_regio_delta,
                protein_residue_regio_delta_tsv_path,
            )
            condition_patch_summary_tsv_path = options.output_dir / "condition_patch_summary.tsv"
            write_condition_patch_summary(
                residue_importance_outputs.condition_patch_summary,
                condition_patch_summary_tsv_path,
            )
            protein_patch_summary_tsv_path = options.output_dir / "protein_patch_summary.tsv"
            write_protein_patch_summary(
                residue_importance_outputs.protein_patch_summary,
                protein_patch_summary_tsv_path,
            )
            _record_timing(
                step="analysis.residue_importance",
                start=residue_importance_start,
                status="ok",
                detail=(
                    "n_residue_rows="
                    f"{len(residue_importance_outputs.protein_condition_residue_scores)}"
                ),
            )

        n_crystal_anchoring_conditions = len(crystal_report_entries) + len(crystal_report_errors)
        n_crystal_anchoring_errors = len(crystal_report_errors)
        crystal_anchoring_stage_completed = (
            n_crystal_anchoring_conditions > 0 and n_crystal_anchoring_errors == 0
        )

        output_write_start = perf_counter()
        crystal_anchor_tsv_path = options.output_dir / "crystal_anchor_table.tsv"
        _write_tsv_rows(
            crystal_anchor_tsv_path,
            [
                "condition_id",
                "protein_id",
                "ligand_id",
                "cluster_id",
                "medoid_pose_id",
                "crystal_reference_id",
                "source_cif",
                "prepared_subset_cif",
                "selected_protein_chain",
                "local_pocket_rmsd",
                "ligand_rmsd_if_comparable",
                "proximal_sugar_rmsd",
                "contact_overlap_score",
                "same_binding_region_flag",
                "same_general_orientation_flag",
                "ifp_tanimoto",
                "pocket_residues",
                "comparison_status",
                "notes",
            ],
            crystal_anchor_rows,
        )

        pose_geometry_tsv_path = options.output_dir / "pose_geometry.tsv"
        write_pose_geometry_tsv(geometry_metrics, pose_geometry_tsv_path)

        metrics_csv_path = options.output_dir / "metrics.csv"
        build_metrics_csv(list(metrics_record_by_pose_id.values()), metrics_csv_path)

        summary_json_path = options.output_dir / "summary.json"
        summary_json = build_summary(
            run_id=options.run_id,
            mode="production",
            del_variant=options.del_variant,
            qc_reports=[qc_report_payload],
            cluster_results=cluster_results_for_summary,
            geometry_stats=summary_geometry_rows,
            crystal_reports=crystal_reports_for_summary,
            pipeline_version=str(config.get("pipeline_version", "")),
        )
        write_summary_json(summary_json, summary_json_path)

        report_html_path = options.output_dir / "report.html"
        build_report_html(summary_json, metrics_csv_path, report_html_path)
        _record_timing(
            step="analysis.report_outputs",
            start=output_write_start,
            status="ok",
        )

        analysis_summary.update(
            {
                "pose_geometry_tsv": str(pose_geometry_tsv_path),
                "pose_ifp_table_tsv": str(pose_ifp_table_tsv_path) if pose_ifp_table_tsv_path else None,
                "pose_residue_contact_tsv": str(pose_residue_contact_tsv_path) if pose_residue_contact_tsv_path else None,
                "pose_convergence_tsv": str(pose_convergence_tsv_path) if pose_convergence_tsv_path else None,
                "condition_convergence_summary_tsv": str(condition_convergence_summary_tsv_path) if condition_convergence_summary_tsv_path else None,
                "cluster_assignments_tsv": str(cluster_assignments_tsv_path) if cluster_assignments_tsv_path else None,
                "medoid_manifest_tsv": str(medoid_manifest_tsv_path) if medoid_manifest_tsv_path else None,
                "condition_cluster_summary_tsv": str(condition_cluster_summary_tsv_path) if condition_cluster_summary_tsv_path else None,
                "cluster_ifp_signature_tsv": str(cluster_ifp_signature_tsv_path) if cluster_ifp_signature_tsv_path else None,
                "cluster_residue_signature_tsv": str(cluster_residue_signature_tsv_path) if cluster_residue_signature_tsv_path else None,
                "cluster_signatures_json": str(cluster_signatures_json_path) if cluster_signatures_json_path else None,
                "cluster_annotation_stage_completed": cluster_annotation_stage_completed,
                "protein_condition_residue_scores_tsv": str(protein_condition_residue_scores_tsv_path) if protein_condition_residue_scores_tsv_path else None,
                "protein_residue_regio_delta_tsv": str(protein_residue_regio_delta_tsv_path) if protein_residue_regio_delta_tsv_path else None,
                "condition_patch_summary_tsv": str(condition_patch_summary_tsv_path) if condition_patch_summary_tsv_path else None,
                "protein_patch_summary_tsv": str(protein_patch_summary_tsv_path) if protein_patch_summary_tsv_path else None,
                "crystal_anchor_tsv": str(crystal_anchor_tsv_path),
                "crystal_anchoring_reports": crystal_report_entries,
                "crystal_anchoring_errors": crystal_report_errors,
                "crystal_anchoring_stage_completed": crystal_anchoring_stage_completed,
                "n_crystal_anchoring_conditions": n_crystal_anchoring_conditions,
                "n_crystal_anchoring_errors": n_crystal_anchoring_errors,
                "metrics_csv": str(metrics_csv_path),
                "summary_json": str(summary_json_path),
                "report_html": str(report_html_path),
                "n_analyzed": n_analyzed,
                "n_cluster_conditions": n_cluster_conditions,
            }
        )
    else:
        analysis_summary["qc_report_error"] = "No poses were prepared successfully for analysis"
        analysis_summary["n_analyzed"] = 0

    pose_manifest_tsv_path = options.output_dir / "pose_manifest.tsv"
    _write_pose_manifest_tsv(
        pose_inputs,
        output_path=pose_manifest_tsv_path,
        options=options,
        case_by_pose_id=case_by_pose_id,
    )
    pose_confidence_tsv_path = options.output_dir / "pose_confidence.tsv"
    _write_pose_confidence_tsv(
        pose_inputs,
        output_path=pose_confidence_tsv_path,
        options=options,
    )
    structure_index_tsv_path = options.output_dir / "structure_index.tsv"
    _write_structure_index_tsv(
        pose_inputs,
        output_path=structure_index_tsv_path,
        options=options,
        case_by_pose_id=case_by_pose_id,
    )
    qc_attrition_tsv_path = options.output_dir / "qc_attrition_table.tsv"
    _write_qc_attrition_tsv(
        pose_inputs,
        output_path=qc_attrition_tsv_path,
        options=options,
        case_by_pose_id=case_by_pose_id,
    )
    analysis_summary.update(
        {
            "pose_manifest_tsv": str(pose_manifest_tsv_path),
            "pose_confidence_tsv": str(pose_confidence_tsv_path),
            "structure_index_tsv": str(structure_index_tsv_path),
            "qc_attrition_tsv": str(qc_attrition_tsv_path),
        }
    )

    if options.collect_timing_events:
        analysis_summary["timing_events"] = timing_events
    analysis_summary["total_wall_time_s"] = round(perf_counter() - run_start, 6)
    summary_path.write_text(json.dumps(analysis_summary, indent=2))
    success = bool(prepared_poses)
    return AnalysisCoreResult(
        run_id=options.run_id,
        output_dir=options.output_dir,
        summary_path=summary_path,
        qc_report_path=qc_report_path,
        pose_manifest_tsv_path=pose_manifest_tsv_path,
        pose_confidence_tsv_path=pose_confidence_tsv_path,
        structure_index_tsv_path=structure_index_tsv_path,
        qc_attrition_tsv_path=qc_attrition_tsv_path,
        pose_geometry_tsv_path=pose_geometry_tsv_path,
        metrics_csv_path=metrics_csv_path,
        summary_json_path=summary_json_path,
        report_html_path=report_html_path,
        n_discovered=len(pose_inputs),
        n_prepared=len(prepared_poses),
        n_analyzed=n_analyzed,
        success=success,
        pose_ifp_table_tsv_path=pose_ifp_table_tsv_path,
        pose_residue_contact_tsv_path=pose_residue_contact_tsv_path,
        pose_convergence_tsv_path=pose_convergence_tsv_path,
        condition_convergence_summary_tsv_path=condition_convergence_summary_tsv_path,
        cluster_assignments_tsv_path=cluster_assignments_tsv_path,
        medoid_manifest_tsv_path=medoid_manifest_tsv_path,
        condition_cluster_summary_tsv_path=condition_cluster_summary_tsv_path,
        cluster_ifp_signature_tsv_path=cluster_ifp_signature_tsv_path,
        cluster_residue_signature_tsv_path=cluster_residue_signature_tsv_path,
        cluster_signatures_json_path=cluster_signatures_json_path,
        cluster_annotation_stage_completed=cluster_annotation_stage_completed,
        protein_condition_residue_scores_tsv_path=protein_condition_residue_scores_tsv_path,
        protein_residue_regio_delta_tsv_path=protein_residue_regio_delta_tsv_path,
        condition_patch_summary_tsv_path=condition_patch_summary_tsv_path,
        protein_patch_summary_tsv_path=protein_patch_summary_tsv_path,
        crystal_anchor_tsv_path=crystal_anchor_tsv_path,
        n_cluster_conditions=n_cluster_conditions,
        crystal_anchoring_stage_completed=crystal_anchoring_stage_completed,
        n_crystal_anchoring_conditions=n_crystal_anchoring_conditions,
        n_crystal_anchoring_errors=n_crystal_anchoring_errors,
    )
