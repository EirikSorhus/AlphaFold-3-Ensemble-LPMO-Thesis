from __future__ import annotations

import argparse
import json
from pathlib import Path

from lpmo_pipeline.analysis.analysis_orchestrator import AnalysisCoreResult
from lpmo_pipeline.cli import cmd_cbm_paired, cmd_family_enrichment, cmd_predictive, cmd_run, cmd_tune


def test_cmd_tune_calls_run_tuning_and_writes_manifest(tmp_path, monkeypatch) -> None:
    output_dir = tmp_path / "tune_results"
    config_path = tmp_path / "tuning.yaml"
    config_path.write_text("model: AF3\n")

    called = {}

    def _fake_run_tuning(
        *,
        model: str,
        tuning_config_path: Path,
        test_cases,
        output_dir: Path,
        pipeline_config,
        max_parallel: int,
    ):
        called["model"] = model
        called["tuning_config_path"] = tuning_config_path
        called["test_cases"] = test_cases
        called["output_dir"] = output_dir
        called["pipeline_config"] = pipeline_config
        called["max_parallel"] = max_parallel
        return {"model": model, "best_params": {}}

    monkeypatch.setattr("lpmo_pipeline.cli.run_tuning", _fake_run_tuning)
    monkeypatch.setattr("lpmo_pipeline.cli.ToolVersionFetcher.get_versions", lambda: {})

    args = argparse.Namespace(
        model="AF3",
        config=config_path,
        output=output_dir,
        n_jobs=2,
        run_config=None,
    )

    exit_code = cmd_tune(args)

    assert exit_code == 0
    assert called["model"] == "AF3"
    assert called["tuning_config_path"] == config_path
    assert called["test_cases"] == []
    assert called["output_dir"] == output_dir
    assert called["pipeline_config"] == {"model": "AF3"}
    assert called["max_parallel"] == 2

    manifest = json.loads((output_dir / "run_manifest.json").read_text())
    assert manifest["gates_passed"]["tuning_completed"] is True


def test_cmd_run_calls_analysis_core_and_writes_manifest(tmp_path, monkeypatch) -> None:
    output_dir = tmp_path / "results"
    config_path = tmp_path / "production.yaml"
    config_path.write_text("pipeline_version: '2.1'\nproduction:\n  work_root: /tmp/work\n")

    called = {}

    def _fake_run_analysis_core(
        config,
        output_dir: Path,
        del_variant: str,
        n_jobs: int | None = None,
    ) -> AnalysisCoreResult:
        called["config"] = config
        called["output_dir"] = output_dir
        called["del_variant"] = del_variant
        called["n_jobs"] = n_jobs
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
        cluster_ifp_signature = output_dir / "cluster_ifp_signature.tsv"
        cluster_residue_signature = output_dir / "cluster_residue_signature.tsv"
        cluster_signatures_json = output_dir / "cluster_signatures.json"
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
            cluster_ifp_signature,
            cluster_residue_signature,
            cluster_signatures_json,
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
            cluster_ifp_signature_tsv_path=cluster_ifp_signature,
            cluster_residue_signature_tsv_path=cluster_residue_signature,
            cluster_signatures_json_path=cluster_signatures_json,
            cluster_annotation_stage_completed=True,
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
        n_jobs=3,
    )

    exit_code = cmd_run(args)

    assert exit_code == 0
    assert called["del_variant"] == "del_a"
    assert called["output_dir"] == output_dir
    assert called["n_jobs"] == 3
    manifest = json.loads((output_dir / "run_manifest.json").read_text())
    assert manifest["gates_passed"]["production_pipeline_started"] is True
    assert manifest["gates_passed"]["analysis_core_completed"] is True
    assert manifest["gates_passed"]["hard_qc_completed"] is True
    assert manifest["gates_passed"]["geometry_stage_completed"] is True
    assert manifest["gates_passed"]["ifp_stage_completed"] is True
    assert manifest["gates_passed"]["clustering_stage_completed"] is True
    assert manifest["gates_passed"]["cluster_annotation_stage_completed"] is True
    assert manifest["gates_passed"]["crystal_anchoring_stage_completed"] is True


def test_cmd_predictive_calls_postprocess(tmp_path, monkeypatch) -> None:
    condition_table = tmp_path / "condition_table.tsv"
    protein_metadata = tmp_path / "protein_metadata.tsv"
    output_dir = tmp_path / "predictive_results"
    condition_table.write_text("condition_id\tprotein_id\n")
    protein_metadata.write_text("protein_id\texperimental_regio_label\n")

    called = {}

    def _fake_run_predictive_postprocess(
        *,
        condition_table_path: Path,
        protein_metadata_path: Path,
        output_dir: Path,
        task: str,
        n_folds: int,
        random_state: int,
    ):
        called["condition_table_path"] = condition_table_path
        called["protein_metadata_path"] = protein_metadata_path
        called["output_dir"] = output_dir
        called["task"] = task
        called["n_folds"] = n_folds
        called["random_state"] = random_state
        output_dir.mkdir(parents=True, exist_ok=True)
        summary_path = output_dir / "10_predictive" / "predictive_summary.json"
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text("{}\n")

        class _Result:
            def __init__(self) -> None:
                self.summary_path = summary_path
                self.modeling_table_paths = {"c1_c4": output_dir / "10_predictive" / "modeling_tables" / "c1.tsv"}
                self.metrics_paths = {"c1_activity": output_dir / "10_predictive" / "cv_results" / "c1_metrics.tsv"}
                self.predictions_paths = {"c1_activity": output_dir / "10_predictive" / "cv_results" / "c1_predictions.tsv"}

        return _Result()

    monkeypatch.setattr("lpmo_pipeline.cli.run_predictive_postprocess", _fake_run_predictive_postprocess)

    args = argparse.Namespace(
        condition_table=condition_table,
        protein_metadata=protein_metadata,
        output=output_dir,
        task="c1_c4",
        n_folds=3,
        random_state=11,
    )

    exit_code = cmd_predictive(args)

    assert exit_code == 0
    assert called["condition_table_path"] == condition_table
    assert called["protein_metadata_path"] == protein_metadata
    assert called["output_dir"] == output_dir
    assert called["task"] == "c1_c4"
    assert called["n_folds"] == 3
    assert called["random_state"] == 11


def test_cmd_cbm_paired_calls_postprocess(tmp_path, monkeypatch) -> None:
    condition_table = tmp_path / "condition_table.tsv"
    cluster_table = tmp_path / "cluster_table.tsv"
    protein_metadata = tmp_path / "protein_metadata.tsv"
    output_dir = tmp_path / "cbm_results"
    condition_table.write_text("condition_id\tprotein_id\n")
    cluster_table.write_text("condition_id\tcluster_type\n")
    protein_metadata.write_text("protein_id\tcbm_type\n")

    called = {}

    def _fake_run_cbm_paired_analysis(
        *,
        condition_table_path: Path,
        output_dir: Path,
        cluster_table_path: Path | None,
        cluster_residue_signature_table_path: Path | None,
        protein_metadata_path: Path | None,
        random_state: int,
    ):
        called["condition_table_path"] = condition_table_path
        called["cluster_table_path"] = cluster_table_path
        called["cluster_residue_signature_table_path"] = cluster_residue_signature_table_path
        called["protein_metadata_path"] = protein_metadata_path
        called["output_dir"] = output_dir
        called["random_state"] = random_state
        summary_path = output_dir / "15_cbm_paired_analysis" / "cbm_paired_analysis_summary.json"
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text("{}\n")

        class _Result:
            def __init__(self) -> None:
                self.summary_path = summary_path
                self.table_paths = {
                    "paired_comparison": output_dir / "15_cbm_paired_analysis" / "cbm_paired_comparison_table.tsv"
                }

        return _Result()

    monkeypatch.setattr("lpmo_pipeline.cli.run_cbm_paired_analysis", _fake_run_cbm_paired_analysis)

    args = argparse.Namespace(
        condition_table=condition_table,
        cluster_table=cluster_table,
        cluster_residue_signature_table=None,
        protein_metadata=protein_metadata,
        output=output_dir,
        random_state=13,
    )

    exit_code = cmd_cbm_paired(args)

    assert exit_code == 0
    assert called["condition_table_path"] == condition_table
    assert called["cluster_table_path"] == cluster_table
    assert called["cluster_residue_signature_table_path"] is None
    assert called["protein_metadata_path"] == protein_metadata
    assert called["output_dir"] == output_dir
    assert called["random_state"] == 13


def test_cmd_family_enrichment_calls_postprocess(tmp_path, monkeypatch) -> None:
    scores_path = tmp_path / "protein_condition_residue_scores.tsv"
    delta_path = tmp_path / "protein_residue_regio_delta.tsv"
    metadata_path = tmp_path / "protein_metadata.tsv"
    core_fasta_path = tmp_path / "core.fasta"
    output_dir = tmp_path / "family_results"
    scores_path.write_text("protein_id\tcondition_id\n")
    delta_path.write_text("protein_id\tresidue_label\n")
    metadata_path.write_text("UniProt_ID\tCAZy_family\n")
    core_fasta_path.write_text(">UniProtIDs|P1|Example|AA10 enzyme\nMNAS\n")

    called = {}

    def _fake_run_family_enrichment_postprocess(
        *,
        protein_condition_residue_scores_path: Path,
        protein_residue_regio_delta_path: Path,
        protein_metadata_path: Path,
        core_fasta_path: Path,
        output_dir: Path,
        alignment_dir: Path | None,
        families: tuple[str, ...],
        substrate_classes: tuple[str, ...],
        mafft_executable: str,
    ):
        called["protein_condition_residue_scores_path"] = protein_condition_residue_scores_path
        called["protein_residue_regio_delta_path"] = protein_residue_regio_delta_path
        called["protein_metadata_path"] = protein_metadata_path
        called["core_fasta_path"] = core_fasta_path
        called["output_dir"] = output_dir
        called["alignment_dir"] = alignment_dir
        called["families"] = families
        called["substrate_classes"] = substrate_classes
        called["mafft_executable"] = mafft_executable
        family_root = output_dir / "08_family_residue_enrichment"
        family_root.mkdir(parents=True, exist_ok=True)

        class _Result:
            def __init__(self) -> None:
                self.family_aligned_residue_table_path = family_root / "family_aligned_residue_table.tsv"
                self.family_residue_enrichment_path = family_root / "family_residue_enrichment.tsv"
                self.family_substrate_residue_enrichment_path = family_root / "family_substrate_residue_enrichment.tsv"
                self.family_wrong_ligand_residue_enrichment_path = family_root / "family_wrong_ligand_residue_enrichment.tsv"
                self.alignment_manifest_path = family_root / "family_alignment_manifest.tsv"
                self.summary_path = family_root / "family_enrichment_summary.json"

        result = _Result()
        result.family_aligned_residue_table_path.write_text("\n")
        result.family_residue_enrichment_path.write_text("\n")
        result.family_substrate_residue_enrichment_path.write_text("\n")
        result.family_wrong_ligand_residue_enrichment_path.write_text("\n")
        result.alignment_manifest_path.write_text("\n")
        result.summary_path.write_text("{}\n")
        return result

    monkeypatch.setattr(
        "lpmo_pipeline.cli.run_family_enrichment_postprocess",
        _fake_run_family_enrichment_postprocess,
    )

    args = argparse.Namespace(
        protein_condition_residue_scores=scores_path,
        protein_residue_regio_delta=delta_path,
        protein_metadata=metadata_path,
        core_fasta=core_fasta_path,
        output=output_dir,
        alignment_dir=None,
        families=["AA9", "AA10"],
        substrates=["cellulose", "chitin"],
        mafft_executable="mafft",
    )

    exit_code = cmd_family_enrichment(args)

    assert exit_code == 0
    assert called["protein_condition_residue_scores_path"] == scores_path
    assert called["protein_residue_regio_delta_path"] == delta_path
    assert called["protein_metadata_path"] == metadata_path
    assert called["core_fasta_path"] == core_fasta_path
    assert called["output_dir"] == output_dir
    assert called["alignment_dir"] is None
    assert called["families"] == ("AA9", "AA10")
    assert called["substrate_classes"] == ("cellulose", "chitin")
    assert called["mafft_executable"] == "mafft"


def test_cmd_predictive_resolves_values_from_run_config(tmp_path, monkeypatch) -> None:
    condition_table = tmp_path / "condition_table.tsv"
    protein_metadata = tmp_path / "protein_metadata.tsv"
    output_dir = tmp_path / "predictive_results"
    run_config = tmp_path / "pipeline_run.yaml"
    condition_table.write_text("condition_id\tprotein_id\n")
    protein_metadata.write_text("protein_id\texperimental_regio_label\n")
    run_config.write_text(
        f"""
commands:
  predictive:
    condition_table: {condition_table}
    protein_metadata: {protein_metadata}
    output: {output_dir}
    task: substrate
    n_folds: 7
    random_state: 99
""".lstrip()
    )

    called = {}

    def _fake_run_predictive_postprocess(
        *,
        condition_table_path: Path,
        protein_metadata_path: Path,
        output_dir: Path,
        task: str,
        n_folds: int,
        random_state: int,
    ):
        called["condition_table_path"] = condition_table_path
        called["protein_metadata_path"] = protein_metadata_path
        called["output_dir"] = output_dir
        called["task"] = task
        called["n_folds"] = n_folds
        called["random_state"] = random_state
        summary_path = output_dir / "10_predictive" / "predictive_summary.json"
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text("{}\n")

        class _Result:
            def __init__(self) -> None:
                self.summary_path = summary_path
                self.modeling_table_paths = {}
                self.metrics_paths = {}
                self.predictions_paths = {}

        return _Result()

    monkeypatch.setattr("lpmo_pipeline.cli.run_predictive_postprocess", _fake_run_predictive_postprocess)

    args = argparse.Namespace(
        run_config=run_config,
        condition_table=None,
        protein_metadata=None,
        output=None,
        task=None,
        n_folds=None,
        random_state=None,
    )

    exit_code = cmd_predictive(args)
    assert exit_code == 0
    assert called["condition_table_path"] == condition_table
    assert called["protein_metadata_path"] == protein_metadata
    assert called["output_dir"] == output_dir
    assert called["task"] == "substrate"
    assert called["n_folds"] == 7
    assert called["random_state"] == 99


def test_cmd_run_resolves_values_from_run_config(tmp_path, monkeypatch) -> None:
    output_dir = tmp_path / "results"
    production_config = tmp_path / "production.yaml"
    run_config = tmp_path / "pipeline_run.yaml"
    production_config.write_text("pipeline_version: '2.1'\nproduction:\n  work_root: /tmp/work\n")
    run_config.write_text(
        f"""
commands:
  run:
    mode: production
    config: {production_config}
    output: {output_dir}
    del: del_b
    n_jobs: 2
""".lstrip()
    )

    called = {}

    def _fake_run_analysis_core(
        config,
        output_dir: Path,
        del_variant: str,
        n_jobs: int | None = None,
    ) -> AnalysisCoreResult:
        called["output_dir"] = output_dir
        called["del_variant"] = del_variant
        called["n_jobs"] = n_jobs
        output_dir.mkdir(parents=True, exist_ok=True)
        summary = output_dir / "analysis_core_summary.json"
        summary.write_text("{}\n")
        return AnalysisCoreResult(
            run_id="cfg",
            output_dir=output_dir,
            summary_path=summary,
            qc_report_path=None,
            pose_manifest_tsv_path=None,
            pose_confidence_tsv_path=None,
            structure_index_tsv_path=None,
            qc_attrition_tsv_path=None,
            pose_geometry_tsv_path=None,
            metrics_csv_path=None,
            summary_json_path=None,
            report_html_path=None,
            n_discovered=0,
            n_prepared=0,
            n_analyzed=0,
            success=False,
        )

    monkeypatch.setattr("lpmo_pipeline.cli.run_analysis_core", _fake_run_analysis_core)
    monkeypatch.setattr("lpmo_pipeline.cli.ToolVersionFetcher.get_versions", lambda: {})

    args = argparse.Namespace(
        run_config=run_config,
        mode=None,
        config=None,
        output=None,
        del_branch=None,
        n_jobs=None,
    )

    exit_code = cmd_run(args)
    assert exit_code == 1
    assert called["output_dir"] == output_dir
    assert called["del_variant"] == "del_b"
    assert called["n_jobs"] == 2
