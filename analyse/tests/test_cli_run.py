from __future__ import annotations

import argparse
import json
from pathlib import Path

from lpmo_pipeline.analysis.analysis_orchestrator import AnalysisCoreResult
from lpmo_pipeline.cli import cmd_run


def test_cmd_run_calls_analysis_core_and_writes_manifest(tmp_path, monkeypatch) -> None:
    output_dir = tmp_path / "results"
    config_path = tmp_path / "production.yaml"
    config_path.write_text("pipeline_version: '2.1'\nproduction:\n  work_root: /tmp/work\n")

    called = {}

    def _fake_run_analysis_core(config, output_dir: Path, del_variant: str) -> AnalysisCoreResult:
        called["config"] = config
        called["output_dir"] = output_dir
        called["del_variant"] = del_variant
        output_dir.mkdir(parents=True, exist_ok=True)

        qc_report = output_dir / "qc_report.json"
        pose_manifest = output_dir / "pose_manifest.tsv"
        pose_confidence = output_dir / "pose_confidence.tsv"
        structure_index = output_dir / "structure_index.tsv"
        qc_attrition = output_dir / "qc_attrition_table.tsv"
        pose_geometry = output_dir / "pose_geometry.tsv"
        pose_ifp_table = output_dir / "pose_ifp_table.tsv"
        cluster_assignments = output_dir / "cluster_assignments.tsv"
        medoid_manifest = output_dir / "medoid_manifest.tsv"
        condition_cluster_summary = output_dir / "condition_cluster_summary.tsv"
        crystal_anchor = output_dir / "crystal_anchor_table.tsv"
        metrics_csv = output_dir / "metrics.csv"
        summary_json = output_dir / "summary.json"
        report_html = output_dir / "report.html"
        run_summary = output_dir / "analysis_core_summary.json"
        for path in (
            qc_report,
            pose_manifest,
            pose_confidence,
            structure_index,
            qc_attrition,
            pose_geometry,
            pose_ifp_table,
            cluster_assignments,
            medoid_manifest,
            condition_cluster_summary,
            crystal_anchor,
            metrics_csv,
            summary_json,
            report_html,
            run_summary,
        ):
            path.write_text("{}\n")

        return AnalysisCoreResult(
            run_id="results",
            output_dir=output_dir,
            summary_path=run_summary,
            qc_report_path=qc_report,
            pose_manifest_tsv_path=pose_manifest,
            pose_confidence_tsv_path=pose_confidence,
            structure_index_tsv_path=structure_index,
            qc_attrition_tsv_path=qc_attrition,
            pose_geometry_tsv_path=pose_geometry,
            pose_ifp_table_tsv_path=pose_ifp_table,
            cluster_assignments_tsv_path=cluster_assignments,
            medoid_manifest_tsv_path=medoid_manifest,
            condition_cluster_summary_tsv_path=condition_cluster_summary,
            crystal_anchor_tsv_path=crystal_anchor,
            metrics_csv_path=metrics_csv,
            summary_json_path=summary_json,
            report_html_path=report_html,
            n_discovered=1,
            n_prepared=1,
            n_analyzed=1,
            success=True,
            crystal_anchoring_stage_completed=True,
            n_crystal_anchoring_conditions=1,
        )

    monkeypatch.setattr("lpmo_pipeline.cli.run_analysis_core", _fake_run_analysis_core)
    monkeypatch.setattr("lpmo_pipeline.cli.ToolVersionFetcher.get_versions", lambda: {})

    args = argparse.Namespace(
        mode="production",
        config=config_path,
        output=output_dir,
        del_branch="del_a",
        n_jobs=1,
    )

    exit_code = cmd_run(args)

    assert exit_code == 0
    assert called["del_variant"] == "del_a"
    assert called["output_dir"] == output_dir
    manifest = json.loads((output_dir / "run_manifest.json").read_text())
    assert manifest["gates_passed"]["production_pipeline_started"] is True
    assert manifest["gates_passed"]["analysis_core_completed"] is True
    assert manifest["gates_passed"]["hard_qc_completed"] is True
    assert manifest["gates_passed"]["geometry_stage_completed"] is True
    assert manifest["gates_passed"]["ifp_stage_completed"] is True
    assert manifest["gates_passed"]["clustering_stage_completed"] is True
    assert manifest["gates_passed"]["crystal_anchoring_stage_completed"] is True
