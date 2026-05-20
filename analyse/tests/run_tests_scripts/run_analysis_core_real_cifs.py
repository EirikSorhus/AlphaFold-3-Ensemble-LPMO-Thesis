#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import traceback
from pathlib import Path
from typing import Any

import yaml

from lpmo_pipeline.cli import cmd_run


DEFAULT_CIF_PATHS: tuple[Path, ...] = (
    Path(
        "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core/"
        "NAG4/af3/runs/408002/Q7SCE9_NAG4/seed-4_sample-1/Q7SCE9_NAG4_seed-4_sample-1_model.cif"
    ),
    Path(
        "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core/"
        "STA6/af3/runs/408006/Q59930_STA6/seed-9_sample-0/Q59930_STA6_seed-9_sample-0_model.cif"
    ),
    Path(
        "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core/"
        "STA4/af3/runs/408005/A0A0S2GKZ1_STA4/seed-9_sample-0/A0A0S2GKZ1_STA4_seed-9_sample-0_model.cif"
    ),
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the production analysis-core path on a staged set of real AF3 CIFs.",
    )
    parser.add_argument(
        "--run-dir",
        required=True,
        help="Output directory for the staged work root and production artifacts.",
    )
    parser.add_argument(
        "--run-id",
        default="",
        help="Optional run identifier written into the generated config/output.",
    )
    parser.add_argument(
        "--del-branch",
        default="del_a",
        choices=["del_a", "del_b"],
        help="Analysis branch to run through the production CLI.",
    )
    parser.add_argument(
        "cifs",
        nargs="*",
        help="Optional explicit CIF paths. Defaults to the three standard real-CIF cases.",
    )
    return parser.parse_args()


def _parse_pose_id(pose_id: str) -> tuple[str, str, int, int]:
    match = re.match(r"^(?P<protein>[^_]+)_(?P<ligand>[^_]+)_seed-(?P<seed>\d+)_sample-(?P<sample>\d+)_model$", pose_id)
    if not match:
        raise ValueError(f"Unexpected pose_id format: {pose_id}")
    return (
        match.group("protein"),
        match.group("ligand"),
        int(match.group("seed")),
        int(match.group("sample")),
    )


def _stage_work_root(cif_paths: list[Path], work_root: Path) -> list[dict[str, Any]]:
    staged_cases: list[dict[str, Any]] = []
    for index, cif_path in enumerate(cif_paths, start=1):
        pose_id = cif_path.stem
        protein_id, ligand_id, seed, sample = _parse_pose_id(pose_id)
        run_id = "000001"
        model_dir = work_root / ligand_id / "af3"
        runs_dir = model_dir / "runs"
        run_dir = runs_dir / run_id
        ut_dir = run_dir / f"{protein_id}_{ligand_id}"
        sample_dir = ut_dir / f"seed-{seed}_sample-{sample}"
        sample_dir.mkdir(parents=True, exist_ok=True)

        staged_cif = sample_dir / cif_path.name
        if staged_cif.exists() or staged_cif.is_symlink():
            staged_cif.unlink()
        staged_cif.symlink_to(cif_path)

        latest_link = model_dir / "latest"
        if latest_link.exists() or latest_link.is_symlink():
            latest_link.unlink()
        latest_link.symlink_to(Path("runs") / run_id)

        staged_cases.append(
            {
                "index": index,
                "pose_id": pose_id,
                "protein_id": protein_id,
                "ligand_id": ligand_id,
                "seed": seed,
                "sample": sample,
                "source_cif": str(cif_path),
                "staged_cif": str(staged_cif),
            }
        )
    return staged_cases


def _write_config(run_dir: Path, work_root: Path, run_id: str, include_targets: list[str]) -> Path:
    config = {
        "pipeline_version": "2.1",
        "production": {
            "run_id": run_id,
            "work_root": str(work_root),
            "af3_only": True,
            "latest_only": True,
            "include_targets": include_targets,
        },
    }
    config_path = run_dir / "production.analysis_core.real_cifs.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    return config_path


def _count_table_rows(path: Path) -> int:
    if not path.exists():
        return 0
    lines = path.read_text().splitlines()
    return max(len(lines) - 1, 0)


def main() -> int:
    args = _parse_args()
    run_dir = Path(args.run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    cif_paths = [Path(path).resolve() for path in args.cifs] if args.cifs else list(DEFAULT_CIF_PATHS)
    run_id = args.run_id or run_dir.name

    summary: dict[str, Any] = {
        "run_dir": str(run_dir),
        "run_id": run_id,
        "requested_cifs": [str(path) for path in cif_paths],
        "n_requested": len(cif_paths),
    }

    for cif_path in cif_paths:
        if not cif_path.exists():
            summary["error"] = f"Missing CIF: {cif_path}"
            print(json.dumps(summary, indent=2))
            return 1

    try:
        staged_work_root = run_dir / "staged_work_root"
        staged_cases = _stage_work_root(cif_paths, staged_work_root)
        include_targets = sorted({case["ligand_id"] for case in staged_cases})
        config_path = _write_config(run_dir, staged_work_root, run_id, include_targets)

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

        qc_report_path = production_output / "qc_report.json"
        pose_manifest_path = production_output / "pose_manifest.tsv"
        pose_confidence_path = production_output / "pose_confidence.tsv"
        structure_index_path = production_output / "structure_index.tsv"
        qc_attrition_path = production_output / "qc_attrition_table.tsv"
        pose_geometry_path = production_output / "pose_geometry.tsv"
        pose_ifp_table_path = production_output / "pose_ifp_table.tsv"
        cluster_assignments_path = production_output / "cluster_assignments.tsv"
        medoid_manifest_path = production_output / "medoid_manifest.tsv"
        condition_cluster_summary_path = production_output / "condition_cluster_summary.tsv"
        cluster_ifp_signature_path = production_output / "cluster_ifp_signature.tsv"
        cluster_residue_signature_path = production_output / "cluster_residue_signature.tsv"
        cluster_signatures_json_path = production_output / "cluster_signatures.json"
        protein_condition_residue_scores_path = production_output / "protein_condition_residue_scores.tsv"
        protein_residue_regio_delta_path = production_output / "protein_residue_regio_delta.tsv"
        condition_patch_summary_path = production_output / "condition_patch_summary.tsv"
        protein_patch_summary_path = production_output / "protein_patch_summary.tsv"
        crystal_anchor_path = production_output / "crystal_anchor_table.tsv"
        metrics_csv_path = production_output / "metrics.csv"
        summary_json_path = production_output / "summary.json"
        report_html_path = production_output / "report.html"
        manifest_path = production_output / "run_manifest.json"
        analysis_summary_path = production_output / "analysis_core_summary.json"
        debug_pdb_paths = [str(path) for path in production_output.rglob("geometry_debug.pdb")]

        summary.update(
            {
                "config_path": str(config_path),
                "staged_work_root": str(staged_work_root),
                "staged_cases": staged_cases,
                "production_output": str(production_output),
                "cli_exit_code": exit_code,
                "qc_report_path": str(qc_report_path),
                "pose_manifest_tsv": str(pose_manifest_path),
                "pose_confidence_tsv": str(pose_confidence_path),
                "structure_index_tsv": str(structure_index_path),
                "qc_attrition_tsv": str(qc_attrition_path),
                "pose_geometry_tsv": str(pose_geometry_path),
                "pose_ifp_table_tsv": str(pose_ifp_table_path),
                "cluster_assignments_tsv": str(cluster_assignments_path),
                "medoid_manifest_tsv": str(medoid_manifest_path),
                "condition_cluster_summary_tsv": str(condition_cluster_summary_path),
                "cluster_ifp_signature_tsv": str(cluster_ifp_signature_path),
                "cluster_residue_signature_tsv": str(cluster_residue_signature_path),
                "cluster_signatures_json": str(cluster_signatures_json_path),
                "protein_condition_residue_scores_tsv": str(protein_condition_residue_scores_path),
                "protein_residue_regio_delta_tsv": str(protein_residue_regio_delta_path),
                "condition_patch_summary_tsv": str(condition_patch_summary_path),
                "protein_patch_summary_tsv": str(protein_patch_summary_path),
                "crystal_anchor_tsv": str(crystal_anchor_path),
                "metrics_csv": str(metrics_csv_path),
                "summary_json": str(summary_json_path),
                "report_html": str(report_html_path),
                "manifest_path": str(manifest_path),
                "analysis_core_summary_path": str(analysis_summary_path),
                "pose_manifest_rows": _count_table_rows(pose_manifest_path),
                "pose_confidence_rows": _count_table_rows(pose_confidence_path),
                "structure_index_rows": _count_table_rows(structure_index_path),
                "qc_attrition_rows": _count_table_rows(qc_attrition_path),
                "pose_geometry_rows": _count_table_rows(pose_geometry_path),
                "pose_ifp_rows": _count_table_rows(pose_ifp_table_path),
                "cluster_assignment_rows": _count_table_rows(cluster_assignments_path),
                "medoid_manifest_rows": _count_table_rows(medoid_manifest_path),
                "condition_cluster_summary_rows": _count_table_rows(condition_cluster_summary_path),
                "cluster_ifp_signature_rows": _count_table_rows(cluster_ifp_signature_path),
                "cluster_residue_signature_rows": _count_table_rows(cluster_residue_signature_path),
                "protein_condition_residue_scores_rows": _count_table_rows(protein_condition_residue_scores_path),
                "protein_residue_regio_delta_rows": _count_table_rows(protein_residue_regio_delta_path),
                "condition_patch_summary_rows": _count_table_rows(condition_patch_summary_path),
                "protein_patch_summary_rows": _count_table_rows(protein_patch_summary_path),
                "crystal_anchor_rows": _count_table_rows(crystal_anchor_path),
                "metrics_rows": _count_table_rows(metrics_csv_path),
                "debug_pdb_paths": debug_pdb_paths,
            }
        )

        if analysis_summary_path.exists():
            summary["analysis_core_summary"] = json.loads(analysis_summary_path.read_text())
        if manifest_path.exists():
            summary["manifest"] = json.loads(manifest_path.read_text())
        if cluster_signatures_json_path.exists():
            summary["cluster_signatures"] = json.loads(cluster_signatures_json_path.read_text())

        summary_path = run_dir / "analysis_core_real_cifs_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))

        if exit_code != 0:
            return exit_code
        if debug_pdb_paths:
            return 1
        expected_paths = [
            qc_report_path,
            pose_manifest_path,
            pose_confidence_path,
            structure_index_path,
            qc_attrition_path,
            pose_geometry_path,
            crystal_anchor_path,
            summary_json_path,
            report_html_path,
            manifest_path,
            analysis_summary_path,
        ]
        if summary.get("analysis_core_summary", {}).get("n_analyzed", 0) > 0:
            expected_paths.extend(
                [
                    pose_ifp_table_path,
                    cluster_assignments_path,
                    medoid_manifest_path,
                    condition_cluster_summary_path,
                    cluster_ifp_signature_path,
                    cluster_residue_signature_path,
                    cluster_signatures_json_path,
                    protein_condition_residue_scores_path,
                    protein_residue_regio_delta_path,
                    condition_patch_summary_path,
                    protein_patch_summary_path,
                    metrics_csv_path,
                ]
            )
        else:
            expected_paths.append(metrics_csv_path)
        if any(not path.exists() for path in expected_paths):
            return 1
        if summary.get("analysis_core_summary", {}).get("n_analyzed", 0) > 0:
            if (
                summary.get("manifest", {})
                .get("gates_passed", {})
                .get("cluster_annotation_stage_completed")
                is not True
            ):
                return 1
            if summary.get("analysis_core_summary", {}).get("cluster_annotation_stage_completed") is not True:
                return 1
        n_clusters = len(summary.get("cluster_signatures", {}).get("clusters", []))
        if n_clusters > 0 and (
            summary["cluster_ifp_signature_rows"] == 0
            or summary["cluster_residue_signature_rows"] == 0
            or summary["protein_condition_residue_scores_rows"] == 0
            or summary["protein_residue_regio_delta_rows"] == 0
            or summary["condition_patch_summary_rows"] == 0
            or summary["protein_patch_summary_rows"] == 0
        ):
            return 1
        return 0
    except Exception as exc:
        summary["error"] = f"{exc.__class__.__name__}: {exc}"
        summary["traceback"] = traceback.format_exc()
        summary_path = run_dir / "analysis_core_real_cifs_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
