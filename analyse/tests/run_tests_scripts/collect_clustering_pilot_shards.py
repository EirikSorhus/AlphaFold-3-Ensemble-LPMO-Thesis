#!/usr/bin/env python3
"""Collect and merge production outputs from staged Slurm shard runs."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

from lpmo_pipeline.analysis.cbm_comparison import run_cbm_paired_analysis
from lpmo_pipeline.analysis.family_enrichment_postprocess import run_family_enrichment_postprocess
from lpmo_pipeline.analysis.predictive_postprocess import run_predictive_postprocess
from lpmo_pipeline.report.build_report_html import build_report_html


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize and merge staged production shard outputs.")
    parser.add_argument("--run-root", required=True, help="Root used by submit_clustering_pilot_staged.sh.")
    parser.add_argument("--output-json", required=True, help="Aggregate JSON summary.")
    parser.add_argument("--output-tsv", required=True, help="Per-shard TSV status table.")
    parser.add_argument(
        "--merge-output-root",
        default="",
        help="Merged production output root. Defaults to <run-root>/merged_production_output.",
    )
    parser.add_argument(
        "--run-global-postprocess",
        action="store_true",
        help="Run predictive/CBM/family postprocess after merge when required inputs exist.",
    )
    parser.add_argument(
        "--protein-metadata",
        default="/cluster/work/projects/nn1003k/eirik/Masteroppgave_clean/analyse/input_data/metadata_final_ec_fixed.tsv",
        help="Protein metadata TSV for predictive/CBM/family postprocess.",
    )
    parser.add_argument(
        "--core-fasta",
        default=(
            "/cluster/work/projects/nn1003k/eirik/Masteroppgave_clean/analyse/input_data/"
            "lpmo_core_domain_2026-03-14_06-52-15_deduplicated.fasta"
        ),
        help="Core-domain FASTA used by family enrichment postprocess.",
    )
    parser.add_argument(
        "--family-alignment-dir",
        default="",
        help="Optional precomputed AA9/AA10 alignment directory for family enrichment.",
    )
    parser.add_argument(
        "--mafft-executable",
        default="mafft",
        help="MAFFT executable when family alignments must be generated.",
    )
    return parser.parse_args()


def _read_index(index_path: Path) -> list[dict[str, str]]:
    if not index_path.exists():
        return []
    with index_path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _read_json(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text())
    except Exception as exc:
        return {"summary_read_error": f"{exc.__class__.__name__}: {exc}"}


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key in seen:
                continue
            seen.add(key)
            fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def _safe_int(value: Any) -> int:
    try:
        return int(value)
    except Exception:
        return 0


MERGE_SURFACES: tuple[str, ...] = (
    "pose_manifest.tsv",
    "pose_confidence.tsv",
    "qc_attrition_table.tsv",
    "pose_geometry.tsv",
    "pose_ifp_table.tsv",
    "pose_residue_contact_table.tsv",
    "cluster_assignments.tsv",
    "medoid_manifest.tsv",
    "condition_table.tsv",
    "cluster_table.tsv",
    "protein_summary_table.tsv",
    "protein_condition_residue_scores.tsv",
    "protein_residue_regio_delta.tsv",
    "condition_patch_summary.tsv",
    "protein_patch_summary.tsv",
    "crystal_anchor_table.tsv",
    "crystal_geometry_table.tsv",
    "crystal_ifp_diagnostic_summary.tsv",
    "metrics.csv",
)


def _shard_rows(run_root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for construct_type in ("domain_only", "full_length"):
        index_path = run_root / "manifests" / "shards" / construct_type / "shard_index.tsv"
        for index_row in _read_index(index_path):
            shard_id = str(index_row["shard_id"])
            shard_run_dir = run_root / f"{construct_type}_shards" / shard_id
            output_root = str(index_row.get("output_root") or "").strip()
            if output_root:
                production_output = Path(output_root)
            else:
                production_output = shard_run_dir / "production_output"
            summary_path = shard_run_dir / "clustering_pilot_real_case_summary.json"
            summary = _read_json(summary_path) if summary_path.exists() else {}
            analysis_summary_path = production_output / "analysis_core_summary.json"
            analysis_summary = _read_json(analysis_summary_path) if analysis_summary_path.exists() else {}

            run_step_status = str(summary.get("run_step_status") or "")
            if not run_step_status and analysis_summary_path.exists():
                run_step_status = "executed"
            if not run_step_status:
                run_step_status = "missing_summary"

            pose_manifest_path = str(analysis_summary.get("pose_manifest_tsv") or "")
            cluster_assignments_path = str(analysis_summary.get("cluster_assignments_tsv") or "")
            medoid_manifest_path = str(analysis_summary.get("medoid_manifest_tsv") or "")
            crystal_anchor_path = str(analysis_summary.get("crystal_anchor_tsv") or "")
            metrics_csv_path = str(analysis_summary.get("metrics_csv") or "")
            summary_json_path = str(analysis_summary.get("summary_json") or "")
            report_html_path = str(analysis_summary.get("report_html") or "")
            row = {
                "construct_type": construct_type,
                "array_task_id": index_row.get("array_task_id", ""),
                "shard_id": shard_id,
                "n_selections": _safe_int(index_row.get("n_selections")),
                "estimated_pose_count": _safe_int(index_row.get("estimated_pose_count")),
                "protein_ids": index_row.get("protein_ids", ""),
                "targets": index_row.get("targets", ""),
                "manifest_path": index_row.get("manifest_path", ""),
                "run_dir": str(shard_run_dir),
                "production_output": str(production_output),
                "summary_path": str(summary_path) if summary_path.exists() else "",
                "analysis_core_summary_path": str(analysis_summary_path) if analysis_summary_path.exists() else "",
                "run_step_status": run_step_status,
                "summary_step_status": summary.get("summary_step_status", ""),
                "cli_exit_code": summary.get("cli_exit_code", ""),
                "production_outputs_ready": bool(analysis_summary_path.exists()),
                "n_discovered_poses": _safe_int(summary.get("n_discovered_poses")),
                "analysis_core_n_discovered": _safe_int(
                    analysis_summary.get("n_discovered", summary.get("analysis_core_n_discovered"))
                ),
                "analysis_core_n_prepared": _safe_int(
                    analysis_summary.get("n_prepared", summary.get("analysis_core_n_prepared"))
                ),
                "analysis_core_n_analyzed": _safe_int(
                    analysis_summary.get("n_analyzed", summary.get("analysis_core_n_analyzed"))
                ),
                "pose_manifest_tsv": pose_manifest_path,
                "cluster_assignments_tsv": cluster_assignments_path,
                "medoid_manifest_tsv": medoid_manifest_path,
                "crystal_anchor_tsv": crystal_anchor_path,
                "metrics_csv": metrics_csv_path,
                "summary_json": summary_json_path,
                "report_html": report_html_path,
                "error": summary.get("error", summary.get("summary_read_error", analysis_summary.get("summary_read_error", ""))),
            }
            rows.append(row)
    return rows


def _merge_surface_rows(
    rows: list[dict[str, Any]],
    *,
    merge_output_root: Path,
) -> dict[str, Any]:
    merge_output_root.mkdir(parents=True, exist_ok=True)
    merged_outputs: dict[str, Any] = {}

    for filename in MERGE_SURFACES:
        merged_rows: list[dict[str, Any]] = []
        source_files: list[str] = []
        missing_for_surface: list[str] = []
        for shard_row in rows:
            shard_output = Path(str(shard_row["production_output"]))
            source_path = shard_output / filename
            if not source_path.exists():
                missing_for_surface.append(str(source_path))
                continue
            source_files.append(str(source_path))
            shard_records = _read_tsv(source_path)
            for record in shard_records:
                merged_rows.append(
                    {
                        "_construct_type": shard_row["construct_type"],
                        "_shard_id": shard_row["shard_id"],
                        **record,
                    }
                )

        if not source_files:
            merged_outputs[filename] = {
                "merged_path": "",
                "n_rows": 0,
                "n_source_files": 0,
                "missing_source_files": missing_for_surface,
            }
            continue

        merged_path = merge_output_root / filename
        _write_tsv(merged_path, merged_rows)
        merged_outputs[filename] = {
            "merged_path": str(merged_path),
            "n_rows": len(merged_rows),
            "n_source_files": len(source_files),
            "source_files": source_files,
            "missing_source_files": missing_for_surface,
        }

    summary_payload = {
        "run_root": str(merge_output_root.parent),
        "n_shards": len(rows),
        "n_shards_with_analysis_core_summary": sum(1 for row in rows if row["analysis_core_summary_path"]),
        "n_rows": {
            filename: int(details["n_rows"])
            for filename, details in sorted(merged_outputs.items())
        },
        "sources": {
            filename: details.get("source_files", [])
            for filename, details in sorted(merged_outputs.items())
            if details.get("source_files")
        },
    }
    summary_json_path = merge_output_root / "summary.json"
    summary_json_path.write_text(json.dumps(summary_payload, indent=2))
    merged_outputs["summary.json"] = {
        "merged_path": str(summary_json_path),
        "n_rows": 1,
        "n_source_files": 0,
    }

    metrics_rows = _read_tsv(merge_output_root / "metrics.csv")
    report_html_path = merge_output_root / "report.html"
    build_report_html(summary_payload, merge_output_root / "metrics.csv", report_html_path)
    merged_outputs["report.html"] = {
        "merged_path": str(report_html_path),
        "n_rows": len(metrics_rows),
        "n_source_files": 0,
    }
    return merged_outputs


def _run_global_postprocesses(
    *,
    merge_output_root: Path,
    protein_metadata_path: Path,
    core_fasta_path: Path,
    family_alignment_dir: Path | None,
    mafft_executable: str,
) -> dict[str, Any]:
    postprocess: dict[str, Any] = {}
    condition_table_path = merge_output_root / "condition_table.tsv"
    cluster_table_path = merge_output_root / "cluster_table.tsv"
    residue_scores_path = merge_output_root / "protein_condition_residue_scores.tsv"
    residue_delta_path = merge_output_root / "protein_residue_regio_delta.tsv"

    if condition_table_path.exists() and protein_metadata_path.exists():
        predictive = run_predictive_postprocess(
            condition_table_path=condition_table_path,
            protein_metadata_path=protein_metadata_path,
            output_dir=merge_output_root,
        )
        postprocess["predictive"] = {
            "summary_path": str(predictive.summary_path),
            "modeling_tables": {name: str(path) for name, path in sorted(predictive.modeling_table_paths.items())},
            "metrics": {name: str(path) for name, path in sorted(predictive.metrics_paths.items())},
            "predictions": {name: str(path) for name, path in sorted(predictive.predictions_paths.items())},
        }
    else:
        postprocess["predictive"] = {
            "skipped": True,
            "reason": "missing_condition_table_or_protein_metadata",
        }

    if condition_table_path.exists() and cluster_table_path.exists() and protein_metadata_path.exists():
        cbm = run_cbm_paired_analysis(
            condition_table_path=condition_table_path,
            cluster_table_path=cluster_table_path,
            protein_metadata_path=protein_metadata_path,
            output_dir=merge_output_root,
        )
        postprocess["cbm_paired"] = {
            "summary_path": str(cbm.summary_path),
            "table_paths": {name: str(path) for name, path in sorted(cbm.table_paths.items())},
        }
    else:
        postprocess["cbm_paired"] = {
            "skipped": True,
            "reason": "missing_condition_table_or_cluster_table_or_protein_metadata",
        }

    if (
        residue_scores_path.exists()
        and residue_delta_path.exists()
        and protein_metadata_path.exists()
        and core_fasta_path.exists()
    ):
        family = run_family_enrichment_postprocess(
            protein_condition_residue_scores_path=residue_scores_path,
            protein_residue_regio_delta_path=residue_delta_path,
            protein_metadata_path=protein_metadata_path,
            core_fasta_path=core_fasta_path,
            output_dir=merge_output_root,
            alignment_dir=family_alignment_dir,
            mafft_executable=mafft_executable,
        )
        postprocess["family_enrichment"] = {
            "summary_path": str(family.summary_path),
            "processed_families": family.processed_families,
            "skipped_families": family.skipped_families,
            "family_aligned_residue_table_path": str(family.family_aligned_residue_table_path),
            "family_residue_enrichment_path": str(family.family_residue_enrichment_path),
            "family_substrate_residue_enrichment_path": str(
                family.family_substrate_residue_enrichment_path
            ),
            "family_wrong_ligand_residue_enrichment_path": str(
                family.family_wrong_ligand_residue_enrichment_path
            ),
        }
    else:
        postprocess["family_enrichment"] = {
            "skipped": True,
            "reason": "missing_residue_tables_or_metadata_or_core_fasta",
        }

    return postprocess


def main() -> int:
    args = _parse_args()
    run_root = Path(args.run_root).resolve()
    output_json = Path(args.output_json).resolve()
    output_tsv = Path(args.output_tsv).resolve()
    merge_output_root = (
        Path(args.merge_output_root).resolve()
        if args.merge_output_root
        else (run_root / "merged_production_output").resolve()
    )
    protein_metadata_path = Path(args.protein_metadata).resolve()
    core_fasta_path = Path(args.core_fasta).resolve()
    family_alignment_dir = Path(args.family_alignment_dir).resolve() if args.family_alignment_dir else None
    rows = _shard_rows(run_root)

    fieldnames = [
        "construct_type",
        "array_task_id",
        "shard_id",
        "n_selections",
        "estimated_pose_count",
        "protein_ids",
        "targets",
        "manifest_path",
        "run_dir",
        "production_output",
        "summary_path",
        "analysis_core_summary_path",
        "run_step_status",
        "summary_step_status",
        "cli_exit_code",
        "production_outputs_ready",
        "n_discovered_poses",
        "analysis_core_n_discovered",
        "analysis_core_n_prepared",
        "analysis_core_n_analyzed",
        "pose_manifest_tsv",
        "cluster_assignments_tsv",
        "medoid_manifest_tsv",
        "crystal_anchor_tsv",
        "metrics_csv",
        "summary_json",
        "report_html",
        "error",
    ]
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    failed_rows = [
        row
        for row in rows
        if row["run_step_status"] not in {"executed", "reused_existing_output"}
        or not row["production_outputs_ready"]
    ]

    merged_outputs = _merge_surface_rows(rows, merge_output_root=merge_output_root) if rows else {}
    postprocess_summary: dict[str, Any] = {}
    if args.run_global_postprocess and rows:
        postprocess_summary = _run_global_postprocesses(
            merge_output_root=merge_output_root,
            protein_metadata_path=protein_metadata_path,
            core_fasta_path=core_fasta_path,
            family_alignment_dir=family_alignment_dir,
            mafft_executable=args.mafft_executable,
        )

    summary = {
        "run_root": str(run_root),
        "output_tsv": str(output_tsv),
        "merge_output_root": str(merge_output_root),
        "n_shards": len(rows),
        "n_failed_or_incomplete_shards": len(failed_rows),
        "all_shards_complete": bool(rows) and not failed_rows,
        "totals": {
            "n_selections": sum(_safe_int(row["n_selections"]) for row in rows),
            "estimated_pose_count": sum(_safe_int(row["estimated_pose_count"]) for row in rows),
            "n_discovered_poses": sum(_safe_int(row["n_discovered_poses"]) for row in rows),
            "analysis_core_n_discovered": sum(_safe_int(row["analysis_core_n_discovered"]) for row in rows),
            "analysis_core_n_prepared": sum(_safe_int(row["analysis_core_n_prepared"]) for row in rows),
            "analysis_core_n_analyzed": sum(_safe_int(row["analysis_core_n_analyzed"]) for row in rows),
        },
        "merged_outputs": merged_outputs,
        "postprocess": postprocess_summary,
        "failed_or_incomplete_shards": failed_rows,
        "shards": rows,
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0 if summary["all_shards_complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
