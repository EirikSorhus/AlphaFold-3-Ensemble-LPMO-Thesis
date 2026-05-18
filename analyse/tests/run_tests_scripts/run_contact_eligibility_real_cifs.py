#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import traceback
from pathlib import Path
from typing import Any

import yaml

from lpmo_pipeline.cli import cmd_run


DEFAULT_WORK_ROOT = Path(
    "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core"
)
DEFAULT_TARGET = "CEL8"
DEFAULT_PROTEIN_ID = "A0A0S2GKZ1"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the production path on the live work_root for one protein-ligand selection and summarize contact-eligibility outputs.",
    )
    parser.add_argument("--run-dir", required=True, help="Output directory for staged work root and artifacts.")
    parser.add_argument("--run-id", default="", help="Optional run identifier written into the config/output.")
    parser.add_argument(
        "--del-branch",
        default="del_a",
        choices=["del_a", "del_b"],
        help="Analysis branch to run through the production CLI.",
    )
    parser.add_argument(
        "--work-root",
        default=str(DEFAULT_WORK_ROOT),
        help="Existing structure_pipeline work root to analyze directly.",
    )
    parser.add_argument(
        "--target",
        default=DEFAULT_TARGET,
        help="Ligand target to analyze, e.g. STA8.",
    )
    parser.add_argument(
        "--protein-id",
        default=DEFAULT_PROTEIN_ID,
        help="UniProt protein identifier to analyze, e.g. Q7SCE9.",
    )
    parser.add_argument(
        "--all-runs",
        action="store_true",
        help="Disable latest_only and scan all numeric runs under the selected target/model.",
    )
    return parser.parse_args()


def _write_config(
    run_dir: Path,
    work_root: Path,
    run_id: str,
    *,
    target: str,
    protein_id: str,
    latest_only: bool,
) -> Path:
    config = {
        "pipeline_version": "2.1",
        "production": {
            "run_id": run_id,
            "work_root": str(work_root),
            "af3_only": True,
            "latest_only": latest_only,
            "include_targets": [target],
            "include_proteins": [protein_id],
            "run_posebusters": False,
            "run_privateer": False,
        },
    }
    config_path = run_dir / "production.contact_eligibility.real_cifs.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    return config_path


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _count_by_key(rows: list[dict[str, Any]], key: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        value = str(row.get(key, ""))
        counts[value] = counts.get(value, 0) + 1
    return counts


def main() -> int:
    args = _parse_args()
    run_dir = Path(args.run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    run_id = args.run_id or run_dir.name
    work_root = Path(args.work_root).resolve()
    latest_only = not args.all_runs

    summary: dict[str, Any] = {
        "run_dir": str(run_dir),
        "run_id": run_id,
        "work_root": str(work_root),
        "selected_target": args.target,
        "selected_protein_id": args.protein_id,
        "latest_only": latest_only,
    }

    if not work_root.is_dir():
        summary["error"] = f"Missing work_root: {work_root}"
        print(json.dumps(summary, indent=2))
        return 1

    try:
        config_path = _write_config(
            run_dir,
            work_root,
            run_id,
            target=args.target,
            protein_id=args.protein_id,
            latest_only=latest_only,
        )

        production_output = run_dir / "production_output"
        exit_code = cmd_run(
            argparse.Namespace(
                mode="production",
                config=config_path,
                output=production_output,
                del_branch=args.del_branch,
                n_jobs=1,
            )
        )

        analysis_summary_path = production_output / "analysis_core_summary.json"
        condition_cluster_summary_path = production_output / "condition_cluster_summary.tsv"
        pose_ifp_table_path = production_output / "pose_ifp_table.tsv"
        cluster_assignments_path = production_output / "cluster_assignments.tsv"

        analysis_summary = json.loads(analysis_summary_path.read_text()) if analysis_summary_path.exists() else {}
        condition_rows = _read_tsv_rows(condition_cluster_summary_path)
        pose_ifp_rows = _read_tsv_rows(pose_ifp_table_path)
        cluster_assignment_rows = _read_tsv_rows(cluster_assignments_path)

        selected_cases = [
            case
            for case in analysis_summary.get("cases", [])
            if case.get("protein_id") == args.protein_id and case.get("ligand_id") == args.target
        ]
        analyzed_cases = [case for case in selected_cases if case.get("analysis_status") == "analyzed"]
        selected_condition_ids = sorted(
            {
                str(case.get("condition_id"))
                for case in selected_cases
                if case.get("condition_id")
            }
        )
        selected_condition_rows = [
            row for row in condition_rows if row.get("condition_id") in set(selected_condition_ids)
        ]
        eligibility_cases = [
            {
                "pose_id": case.get("pose_id"),
                "condition_id": case.get("condition_id"),
                "contact_eligible": case.get("contact_eligible"),
                "contact_exclusion_class": case.get("contact_exclusion_class"),
                "n_vdw_interactions": case.get("n_vdw_interactions"),
                "n_non_vdw_interactions": case.get("n_non_vdw_interactions"),
                "n_non_vdw_contact_residues": case.get("n_non_vdw_contact_residues"),
            }
            for case in analyzed_cases
        ]

        required_condition_columns = [
            "n_ifp_success",
            "n_contact_eligible",
            "contact_eligible_fraction",
            "null_ifp_fraction",
            "vdw_only_fraction",
            "low_specific_contact_fraction",
            "median_n_non_vdw_interactions",
            "median_n_non_vdw_contact_residues",
        ]
        missing_condition_columns = [
            column_name
            for column_name in required_condition_columns
            if selected_condition_rows and column_name not in selected_condition_rows[0]
        ]

        summary.update(
            {
                "config_path": str(config_path),
                "production_output": str(production_output),
                "cli_exit_code": exit_code,
                "analysis_core_summary_path": str(analysis_summary_path),
                "condition_cluster_summary_tsv": str(condition_cluster_summary_path),
                "pose_ifp_table_tsv": str(pose_ifp_table_path),
                "cluster_assignments_tsv": str(cluster_assignments_path),
                "n_discovered": analysis_summary.get("n_discovered", 0),
                "n_prepared": analysis_summary.get("n_prepared", 0),
                "n_analyzed": analysis_summary.get("n_analyzed", 0),
                "selected_case_count": len(selected_cases),
                "selected_case_status_counts": _count_by_key(selected_cases, "status"),
                "selected_analysis_status_counts": _count_by_key(selected_cases, "analysis_status"),
                "selected_condition_ids": selected_condition_ids,
                "selected_condition_rows": selected_condition_rows,
                "pose_ifp_row_count": len(pose_ifp_rows),
                "cluster_assignment_row_count": len(cluster_assignment_rows),
                "eligibility_cases": eligibility_cases,
                "clustering_invoked": any(
                    int(float(row.get("n_contact_eligible", "0"))) > 0
                    for row in selected_condition_rows
                ),
                "meaningful_clustering_possible": any(
                    int(float(row.get("n_contact_eligible", "0"))) > 1
                    for row in selected_condition_rows
                ),
                "required_condition_columns": required_condition_columns,
                "missing_condition_columns": missing_condition_columns,
            }
        )

        summary_path = run_dir / "contact_eligibility_real_cifs_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))

        if exit_code != 0:
            return exit_code
        if not analysis_summary_path.exists() or not condition_cluster_summary_path.exists():
            return 1
        if analysis_summary.get("n_discovered", 0) <= 1:
            return 1
        if len(selected_cases) <= 1:
            return 1
        if not eligibility_cases:
            return 1
        if missing_condition_columns:
            return 1
        if not selected_condition_rows:
            return 1
        if any(case.get("protein_id") != args.protein_id or case.get("ligand_id") != args.target for case in selected_cases):
            return 1
        if any("contact_eligible" not in case or "contact_exclusion_class" not in case for case in analyzed_cases):
            return 1
        return 0
    except Exception as exc:
        summary["error"] = f"{exc.__class__.__name__}: {exc}"
        summary["traceback"] = traceback.format_exc()
        summary_path = run_dir / "contact_eligibility_real_cifs_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())