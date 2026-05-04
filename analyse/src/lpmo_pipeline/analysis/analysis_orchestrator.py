"""Production analysis orchestration up to downstream geometry.

This module intentionally stops before ProLIF/IFP generation and clustering.
It provides a real production control path for discovery, normalization,
hard QC, downstream geometry, and report generation.
"""
from __future__ import annotations

import json
import logging
import re
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from jsonschema import validate

from lpmo_pipeline.analysis.mdanalysis_metrics import (
    PoseGeometryMetrics,
    compute_pose_metrics_from_structure,
    write_pose_geometry_tsv,
)
from lpmo_pipeline.io.cif_to_pdb import convert_cif_to_pdb
from lpmo_pipeline.io.discovery import discover_work_root
from lpmo_pipeline.io.gemmi_compat import gemmi
from lpmo_pipeline.io.normalize_mmcif import NormalizeMMCIFRunner
from lpmo_pipeline.qc.hard_qc_orchestrator import HardQCInput, run_hard_qc
from lpmo_pipeline.qc.privateer_runner import prepare_privateer_input
from lpmo_pipeline.qc.qc_report import write_qc_report
from lpmo_pipeline.report.build_metrics_csv import build_metrics_csv, merge_pose_record
from lpmo_pipeline.report.build_report_html import build_report_html
from lpmo_pipeline.report.build_summary_json import build_summary, write_summary_json

logger = logging.getLogger(__name__)

_QC_REPORT_SCHEMA_PATH = Path(__file__).resolve().parents[3] / "schemas" / "qc_report_schema.json"
_TARGET_PREFIX_TO_SUBSTRATE = {
    "NAG": "chitin",
    "CEL": "cellulose",
    "BGC": "cellulose",
    "STA": "amylose",
    "GLC": "amylose",
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


@dataclass(frozen=True)
class PoseInputRecord:
    cif_path: Path
    pose_id: str
    protein_id: str
    ligand_id: str
    model: str
    source_run_id: str
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
    pose_geometry_tsv_path: Path | None
    metrics_csv_path: Path | None
    summary_json_path: Path | None
    report_html_path: Path | None
    n_discovered: int
    n_prepared: int
    n_analyzed: int
    success: bool


def _load_qc_report_schema() -> dict[str, Any]:
    return json.loads(_QC_REPORT_SCHEMA_PATH.read_text())


def _slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


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
) -> ProductionRunOptions:
    production = config.get("production") or {}
    work_root_value = production.get("work_root") or config.get("work_root")
    if not work_root_value:
        raise ValueError("Production config must define production.work_root or work_root")

    max_cases_raw = production.get("max_cases")
    max_cases = int(max_cases_raw) if max_cases_raw is not None else None
    include_targets = tuple(str(value) for value in production.get("include_targets", ()))

    return ProductionRunOptions(
        run_id=str(production.get("run_id") or output_dir.name),
        output_dir=output_dir.resolve(),
        del_variant=_normalize_del_variant(del_variant),
        work_root=Path(work_root_value).resolve(),
        af3_only=bool(production.get("af3_only", True)),
        latest_only=bool(production.get("latest_only", True)),
        max_cases=max_cases,
        include_targets=include_targets,
    )


def _discover_pose_inputs(options: ProductionRunOptions) -> tuple[list[str], list[PoseInputRecord], dict[str, Any]]:
    manifest = discover_work_root(
        options.work_root,
        af3_only=options.af3_only,
        latest_only=options.latest_only,
    )

    discovered: list[PoseInputRecord] = []
    for target, model, run_id, ut in manifest.iter_all_uniprot_targets():
        if options.include_targets and target not in options.include_targets:
            continue

        sample_inputs: list[PoseInputRecord] = []
        for sample in ut.samples:
            for cif_path in sample.cif_paths:
                sample_inputs.append(
                    PoseInputRecord(
                        cif_path=Path(cif_path).resolve(),
                        pose_id=Path(cif_path).stem,
                        protein_id=ut.uniprot_id,
                        ligand_id=target,
                        model=model,
                        source_run_id=run_id,
                        seed=sample.seed,
                        sample=sample.sample,
                    )
                )

        if not sample_inputs and ut.model_cif_path is not None:
            sample_inputs.append(
                PoseInputRecord(
                    cif_path=Path(ut.model_cif_path).resolve(),
                    pose_id=Path(ut.model_cif_path).stem,
                    protein_id=ut.uniprot_id,
                    ligand_id=target,
                    model=model,
                    source_run_id=run_id,
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


def _hard_qc_input(prepared_pose: PreparedPose) -> HardQCInput:
    return HardQCInput(
        pose_id=prepared_pose.pose.pose_id,
        mol_pred_path=prepared_pose.posebusters_pdb,
        structure=prepared_pose.structure,
        privateer_cif_path=prepared_pose.privateer_input_cif,
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


def run_analysis_core(
    config: dict[str, Any],
    output_dir: Path,
    del_variant: str,
) -> AnalysisCoreResult:
    options = load_production_options(config, output_dir, del_variant)
    options.output_dir.mkdir(parents=True, exist_ok=True)

    discovery_errors, pose_inputs, discovery_summary = _discover_pose_inputs(options)
    analysis_summary: dict[str, Any] = {
        "run_id": options.run_id,
        "output_dir": str(options.output_dir),
        "del_variant": options.del_variant,
        "work_root": str(options.work_root),
        "af3_only": options.af3_only,
        "latest_only": options.latest_only,
        "max_cases": options.max_cases,
        "include_targets": list(options.include_targets),
        "discovery_summary": discovery_summary,
        "discovery_errors": discovery_errors,
        "n_discovered": len(pose_inputs),
        "cases": [],
    }

    prepared_poses: list[PreparedPose] = []
    case_by_pose_id: dict[str, dict[str, Any]] = {}
    for index, pose in enumerate(pose_inputs, start=1):
        case, prepared_pose = _prepare_pose_case(pose, index=index, output_dir=options.output_dir)
        analysis_summary["cases"].append(case)
        case_by_pose_id[pose.pose_id] = case
        if prepared_pose is not None:
            prepared_poses.append(prepared_pose)

    analysis_summary["n_prepared"] = len(prepared_poses)
    analysis_summary["n_prep_errors"] = sum(
        1 for case in analysis_summary["cases"] if case["status"] == "prep_error"
    )

    summary_path = options.output_dir / "analysis_core_summary.json"
    qc_report_path: Path | None = None
    pose_geometry_tsv_path: Path | None = None
    metrics_csv_path: Path | None = None
    summary_json_path: Path | None = None
    report_html_path: Path | None = None
    n_analyzed = 0

    if prepared_poses:
        qc_report_path = options.output_dir / "qc_report.json"
        report = run_hard_qc([_hard_qc_input(prepared_pose) for prepared_pose in prepared_poses], run_id=options.run_id)
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

        geometry_metrics: list[PoseGeometryMetrics] = []
        summary_geometry_rows: list[dict[str, Any]] = []
        metrics_records: list[dict[str, Any]] = []

        for prepared_pose in prepared_poses:
            pose = prepared_pose.pose
            case = case_by_pose_id[pose.pose_id]
            verdict = verdict_by_pose.get(
                pose.pose_id,
                {"pose_id": pose.pose_id, "status": "unknown", "warnings": []},
            )
            case["qc_verdict"] = verdict

            geometry_row: dict[str, Any] = {}
            metrics: PoseGeometryMetrics | None = None
            if verdict.get("status") != "dropped":
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
                except Exception as exc:
                    case["analysis_status"] = "analysis_error"
                    case["analysis_error"] = f"{exc.__class__.__name__}: {exc}"
                    case["analysis_traceback"] = traceback.format_exc()
            else:
                case["analysis_status"] = "skipped_dropped"

            summary_geometry_rows.append(_summary_geometry_row(pose, metrics))
            substrate_type, dp = _parse_target_metadata(pose.ligand_id)
            metrics_records.append(
                merge_pose_record(
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
            )

        pose_geometry_tsv_path = options.output_dir / "pose_geometry.tsv"
        write_pose_geometry_tsv(geometry_metrics, pose_geometry_tsv_path)

        metrics_csv_path = options.output_dir / "metrics.csv"
        build_metrics_csv(metrics_records, metrics_csv_path)

        summary_json_path = options.output_dir / "summary.json"
        summary_json = build_summary(
            run_id=options.run_id,
            mode="production",
            del_variant=options.del_variant,
            qc_reports=[qc_report_payload],
            cluster_results=[],
            geometry_stats=summary_geometry_rows,
            pipeline_version=str(config.get("pipeline_version", "")),
        )
        write_summary_json(summary_json, summary_json_path)

        report_html_path = options.output_dir / "report.html"
        build_report_html(summary_json, metrics_csv_path, report_html_path)

        analysis_summary.update(
            {
                "pose_geometry_tsv": str(pose_geometry_tsv_path),
                "metrics_csv": str(metrics_csv_path),
                "summary_json": str(summary_json_path),
                "report_html": str(report_html_path),
                "n_analyzed": n_analyzed,
            }
        )
    else:
        analysis_summary["qc_report_error"] = "No poses were prepared successfully for analysis"
        analysis_summary["n_analyzed"] = 0

    summary_path.write_text(json.dumps(analysis_summary, indent=2))
    success = bool(prepared_poses)
    return AnalysisCoreResult(
        run_id=options.run_id,
        output_dir=options.output_dir,
        summary_path=summary_path,
        qc_report_path=qc_report_path,
        pose_geometry_tsv_path=pose_geometry_tsv_path,
        metrics_csv_path=metrics_csv_path,
        summary_json_path=summary_json_path,
        report_html_path=report_html_path,
        n_discovered=len(pose_inputs),
        n_prepared=len(prepared_poses),
        n_analyzed=n_analyzed,
        success=success,
    )