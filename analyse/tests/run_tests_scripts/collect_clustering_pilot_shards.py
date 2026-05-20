#!/usr/bin/env python3
"""Collect compact status from clustering-pilot Slurm shard runs."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize staged clustering-pilot shard outputs.")
    parser.add_argument("--run-root", required=True, help="Root used by submit_clustering_pilot_staged.sh.")
    parser.add_argument("--output-json", required=True, help="Compact aggregate JSON summary.")
    parser.add_argument("--output-tsv", required=True, help="Per-shard TSV summary.")
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


def _shard_rows(run_root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for construct_type in ("domain_only", "full_length"):
        index_path = run_root / "manifests" / "shards" / construct_type / "shard_index.tsv"
        for index_row in _read_index(index_path):
            shard_id = str(index_row["shard_id"])
            shard_run_dir = run_root / f"{construct_type}_shards" / shard_id
            summary_path = shard_run_dir / "clustering_pilot_real_case_summary.json"
            summary = _read_json(summary_path) if summary_path.exists() else {}
            row = {
                "construct_type": construct_type,
                "array_task_id": index_row.get("array_task_id", ""),
                "shard_id": shard_id,
                "n_selections": int(index_row.get("n_selections") or 0),
                "protein_ids": index_row.get("protein_ids", ""),
                "targets": index_row.get("targets", ""),
                "manifest_path": index_row.get("manifest_path", ""),
                "run_dir": str(shard_run_dir),
                "summary_path": str(summary_path) if summary_path.exists() else "",
                "run_step_status": summary.get("run_step_status", "missing_summary"),
                "summary_step_status": summary.get("summary_step_status", ""),
                "cli_exit_code": summary.get("cli_exit_code", ""),
                "pilot_outputs_ready": bool(summary.get("pilot_outputs_ready", False)),
                "n_discovered_poses": int(summary.get("n_discovered_poses") or 0),
                "analysis_core_n_discovered": int(summary.get("analysis_core_n_discovered") or 0),
                "analysis_core_n_prepared": int(summary.get("analysis_core_n_prepared") or 0),
                "analysis_core_n_analyzed": int(summary.get("analysis_core_n_analyzed") or 0),
                "pilot_n_conditions": int(summary.get("pilot_n_conditions") or 0),
                "pilot_n_conditions_formal_clustering_allowed": int(
                    summary.get("pilot_n_conditions_formal_clustering_allowed") or 0
                ),
                "pilot_method_summary_tsv": summary.get("pilot_method_summary_tsv", ""),
                "error": summary.get("error", summary.get("summary_read_error", "")),
            }
            rows.append(row)
    return rows


def main() -> int:
    args = _parse_args()
    run_root = Path(args.run_root).resolve()
    output_json = Path(args.output_json).resolve()
    output_tsv = Path(args.output_tsv).resolve()
    rows = _shard_rows(run_root)

    fieldnames = [
        "construct_type",
        "array_task_id",
        "shard_id",
        "n_selections",
        "protein_ids",
        "targets",
        "manifest_path",
        "run_dir",
        "summary_path",
        "run_step_status",
        "summary_step_status",
        "cli_exit_code",
        "pilot_outputs_ready",
        "n_discovered_poses",
        "analysis_core_n_discovered",
        "analysis_core_n_prepared",
        "analysis_core_n_analyzed",
        "pilot_n_conditions",
        "pilot_n_conditions_formal_clustering_allowed",
        "pilot_method_summary_tsv",
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
        or not row["pilot_outputs_ready"]
    ]
    summary = {
        "run_root": str(run_root),
        "output_tsv": str(output_tsv),
        "n_shards": len(rows),
        "n_failed_or_incomplete_shards": len(failed_rows),
        "all_shards_complete": bool(rows) and not failed_rows,
        "totals": {
            "n_selections": sum(int(row["n_selections"]) for row in rows),
            "n_discovered_poses": sum(int(row["n_discovered_poses"]) for row in rows),
            "analysis_core_n_discovered": sum(int(row["analysis_core_n_discovered"]) for row in rows),
            "analysis_core_n_prepared": sum(int(row["analysis_core_n_prepared"]) for row in rows),
            "analysis_core_n_analyzed": sum(int(row["analysis_core_n_analyzed"]) for row in rows),
            "pilot_n_conditions": sum(int(row["pilot_n_conditions"]) for row in rows),
            "pilot_n_conditions_formal_clustering_allowed": sum(
                int(row["pilot_n_conditions_formal_clustering_allowed"]) for row in rows
            ),
        },
        "failed_or_incomplete_shards": failed_rows,
        "shards": rows,
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0 if summary["all_shards_complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
