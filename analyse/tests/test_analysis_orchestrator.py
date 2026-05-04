from __future__ import annotations

import csv
import json
from pathlib import Path

from lpmo_pipeline.analysis.analysis_orchestrator import run_analysis_core
from lpmo_pipeline.analysis.mdanalysis_metrics import PoseGeometryMetrics
from lpmo_pipeline.qc.qc_report import PoseQCVerdict, build_qc_report


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("test\n")
    return path


def _make_work_root(tmp_path: Path) -> Path:
    work_root = tmp_path / "work"
    ut_dir = work_root / "NAG4" / "af3" / "runs" / "100001" / "Q7SCE9_NAG4"
    _touch(ut_dir / "Q7SCE9_NAG4_model.cif")
    _touch(ut_dir / "seed-1_sample-0" / "Q7SCE9_NAG4_seed-1_sample-0_model.cif")
    _touch(ut_dir / "seed-1_sample-1" / "Q7SCE9_NAG4_seed-1_sample-1_model.cif")
    return work_root


def test_run_analysis_core_writes_qc_geometry_and_reports(tmp_path, monkeypatch) -> None:
    work_root = _make_work_root(tmp_path)
    output_dir = tmp_path / "output"
    config = {
        "pipeline_version": "2.1",
        "production": {
            "work_root": str(work_root),
            "af3_only": True,
            "latest_only": False,
        },
    }

    class _FakeNormalizeRunner:
        def __init__(self, input_cif: Path, output_dir: Path, ccd_client=None) -> None:
            self.output_dir = output_dir

        def run(self):
            self.output_dir.mkdir(parents=True, exist_ok=True)
            normalized = self.output_dir / "normalized.cif"
            normalized.write_text("data_test\n")
            (self.output_dir / "normalize_report.json").write_text("{}\n")
            return True, normalized

    def _fake_convert_cif_to_pdb(normalized_path: Path, output_dir: Path):
        output_dir.mkdir(parents=True, exist_ok=True)
        pdb_path = output_dir / "for_posebusters.pdb"
        pdb_path.write_text("ATOM      1  CA  ALA A   1       0.0     0.0     0.0\nEND\n")
        (output_dir / "cif_to_pdb_report.json").write_text("{}\n")
        return True, pdb_path

    def _fake_prepare_privateer_input(normalized_path: Path, output_path: Path):
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("data_privateer\n")
        return output_path

    class _FakeGemmi:
        @staticmethod
        def read_structure(path: str):
            return {"path": path}

    def _fake_run_hard_qc(poses, run_id: str = "", config=None, config_path=None):
        verdicts = []
        for pose in poses:
            if pose.pose_id.endswith("sample-0_model"):
                verdicts.append(PoseQCVerdict(pose_id=pose.pose_id, status="passed"))
            else:
                verdicts.append(
                    PoseQCVerdict(
                        pose_id=pose.pose_id,
                        status="dropped",
                        drop_reasons=["test_drop"],
                    )
                )
        return build_qc_report(run_id=run_id, verdicts=verdicts)

    geometry_calls: list[str] = []

    def _fake_compute_pose_metrics_from_structure(
        structure,
        pose_id: str = "",
        protein_id: str = "",
        ligand_id: str = "",
        model: str = "",
        **kwargs,
    ) -> PoseGeometryMetrics:
        geometry_calls.append(pose_id)
        return PoseGeometryMetrics(
            pose_id=pose_id,
            model=model,
            protein_id=protein_id,
            ligand_id=ligand_id,
            proximal_sugar_id="B:NAG4",
            brace_integrity_flag=True,
            cu_c1_distance=4.1,
            cu_c4_distance=5.2,
            oxyl_h_c1_distance=2.1,
            oxyl_h_c4_distance=2.7,
            attack_angle_c1=118.0,
            attack_angle_c4=111.0,
            sugar_face_orientation="ambiguous",
            ring_normal_vs_brace_normal=90.0,
            oxyl_h_score_c1=0.95,
            oxyl_h_score_c4=0.68,
            geometry_status_c1="geometry_plausible",
            geometry_status_c4="geometry_plausible",
            his_brace_angle_deg=175.0,
        )

    monkeypatch.setattr(
        "lpmo_pipeline.analysis.analysis_orchestrator.NormalizeMMCIFRunner",
        _FakeNormalizeRunner,
    )
    monkeypatch.setattr(
        "lpmo_pipeline.analysis.analysis_orchestrator.convert_cif_to_pdb",
        _fake_convert_cif_to_pdb,
    )
    monkeypatch.setattr(
        "lpmo_pipeline.analysis.analysis_orchestrator.prepare_privateer_input",
        _fake_prepare_privateer_input,
    )
    monkeypatch.setattr("lpmo_pipeline.analysis.analysis_orchestrator.gemmi", _FakeGemmi())
    monkeypatch.setattr(
        "lpmo_pipeline.analysis.analysis_orchestrator.run_hard_qc",
        _fake_run_hard_qc,
    )
    monkeypatch.setattr(
        "lpmo_pipeline.analysis.analysis_orchestrator.compute_pose_metrics_from_structure",
        _fake_compute_pose_metrics_from_structure,
    )
    monkeypatch.setattr(
        "lpmo_pipeline.analysis.analysis_orchestrator.validate",
        lambda instance, schema: None,
    )

    result = run_analysis_core(config, output_dir, del_variant="del_a")

    assert result.success is True
    assert result.n_discovered == 2
    assert result.n_prepared == 2
    assert result.n_analyzed == 1
    assert geometry_calls == ["Q7SCE9_NAG4_seed-1_sample-0_model"]

    assert result.qc_report_path is not None and result.qc_report_path.exists()
    assert result.pose_geometry_tsv_path is not None and result.pose_geometry_tsv_path.exists()
    assert result.metrics_csv_path is not None and result.metrics_csv_path.exists()
    assert result.summary_json_path is not None and result.summary_json_path.exists()
    assert result.report_html_path is not None and result.report_html_path.exists()
    assert result.summary_path.exists()
    assert list(output_dir.rglob("geometry_debug.pdb")) == []

    pose_geometry_lines = result.pose_geometry_tsv_path.read_text().splitlines()
    assert len(pose_geometry_lines) == 2
    assert "Q7SCE9_NAG4_seed-1_sample-0_model" in pose_geometry_lines[1]

    with open(result.metrics_csv_path, newline="") as handle:
        metrics_rows = list(csv.DictReader(handle))
    assert len(metrics_rows) == 1
    assert metrics_rows[0]["pose_id"] == "Q7SCE9_NAG4_seed-1_sample-0_model"

    summary_json = json.loads(result.summary_json_path.read_text())
    assert summary_json["dataset_stats"]["n_total_poses"] == 2
    assert summary_json["geometry_stats"]["cu_c1_range"] == [4.1, 4.1]

    analysis_summary = json.loads(result.summary_path.read_text())
    assert analysis_summary["qc_counts"]["dropped"] == 1
    dropped_case = next(
        case
        for case in analysis_summary["cases"]
        if case["pose_id"].endswith("sample-1_model")
    )
    assert dropped_case["analysis_status"] == "skipped_dropped"