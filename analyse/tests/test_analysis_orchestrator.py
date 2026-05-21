from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from lpmo_pipeline.analysis.analysis_orchestrator import (
    PreparedPose,
    ProductionRunOptions,
    _choose_crystal_representatives,
    _discover_pose_inputs,
    _discover_top_model_fallback_inputs,
    _prepare_hard_qc_passing_top_model_fallback,
    load_production_options,
    run_analysis_core,
)
from lpmo_pipeline.analysis.clustering_hdbscan import ClusteringResult
from lpmo_pipeline.analysis.crystal_anchoring import (
    CrystalReferencePoseComparison,
    CrystalReferenceScreenReport,
)
from lpmo_pipeline.analysis.convergence_metrics import (
    ConditionConvergenceSummary,
    PoseConvergenceMetrics,
)
from lpmo_pipeline.analysis.prolif_ifp import IFPBatch, IFPResult
from lpmo_pipeline.analysis.mdanalysis_metrics import PoseGeometryMetrics
from lpmo_pipeline.qc.qc_report import PoseQCVerdict, build_qc_report


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("test\n")
    return path


def _make_work_root(tmp_path: Path, *, include_other_protein: bool = False) -> Path:
    work_root = tmp_path / "work"
    ut_dir = work_root / "NAG4" / "af3" / "runs" / "100001" / "Q7SCE9_NAG4"
    _touch(ut_dir / "Q7SCE9_NAG4_model.cif")
    _touch(ut_dir / "seed-1_sample-0" / "Q7SCE9_NAG4_seed-1_sample-0_model.cif")
    _touch(ut_dir / "seed-1_sample-1" / "Q7SCE9_NAG4_seed-1_sample-1_model.cif")
    if include_other_protein:
        other_ut_dir = work_root / "NAG4" / "af3" / "runs" / "100001" / "Q59930_NAG4"
        _touch(other_ut_dir / "Q59930_NAG4_model.cif")
        _touch(other_ut_dir / "seed-2_sample-0" / "Q59930_NAG4_seed-2_sample-0_model.cif")
    return work_root


def test_load_production_options_selects_work_root_from_construct_type(tmp_path: Path) -> None:
    output_dir = tmp_path / "output"
    config = {
        "production": {
            "construct_type": "full_length",
            "work_roots": {
                "domain_only": str(tmp_path / "work_core"),
                "full_length": str(tmp_path / "work_full_length"),
            },
            "include_targets": ["NAG4"],
            "include_proteins": ["Q7SCE9"],
            "run_posebusters": False,
            "run_privateer": False,
            "clustering_pilot": {
                "enabled": True,
                "label": "pilot-a",
                "extra_include_interaction_types": ["Hydrophobic"],
                "minimum_clusterable_n": 12,
                "agglomerative_distance_threshold": 0.35,
            },
        }
    }

    options = load_production_options(config, output_dir, "del_a")

    assert options.construct_type == "full_length"
    assert options.work_root == (tmp_path / "work_full_length").resolve()
    assert options.include_targets == ("NAG4",)
    assert options.include_proteins == ("Q7SCE9",)
    assert options.run_posebusters is False
    assert options.run_privateer is False
    assert options.clustering_pilot is not None
    assert options.clustering_pilot.enabled is True
    assert options.clustering_pilot.label == "pilot-a"
    assert options.clustering_pilot.extra_include_interaction_types == ("Hydrophobic",)
    assert options.clustering_pilot.minimum_clusterable_n == 12
    assert options.clustering_pilot.agglomerative_distance_threshold == 0.35


def test_discover_pose_inputs_filters_by_include_proteins(tmp_path: Path) -> None:
    work_root = _make_work_root(tmp_path, include_other_protein=True)
    options = ProductionRunOptions(
        run_id="test-run",
        output_dir=tmp_path / "output",
        del_variant="a",
        work_root=work_root,
        af3_only=True,
        latest_only=False,
        include_targets=("NAG4",),
        include_proteins=("Q7SCE9",),
    )

    errors, poses, _summary = _discover_pose_inputs(options)

    assert errors == []
    assert len(poses) == 2
    assert {pose.protein_id for pose in poses} == {"Q7SCE9"}
    assert {pose.pose_id for pose in poses} == {
        "Q7SCE9_NAG4_seed-1_sample-0_model",
        "Q7SCE9_NAG4_seed-1_sample-1_model",
    }


def test_discover_top_model_fallback_inputs_keeps_best_af3_model_separate(tmp_path: Path) -> None:
    work_root = _make_work_root(tmp_path)
    options = ProductionRunOptions(
        run_id="test-run",
        output_dir=tmp_path / "output",
        del_variant="a",
        work_root=work_root,
        af3_only=True,
        latest_only=False,
        include_targets=("NAG4",),
    )

    errors, poses, _summary = _discover_pose_inputs(options)
    fallback_errors, fallback_by_condition = _discover_top_model_fallback_inputs(options)

    assert errors == []
    assert fallback_errors == []
    assert {pose.pose_id for pose in poses} == {
        "Q7SCE9_NAG4_seed-1_sample-0_model",
        "Q7SCE9_NAG4_seed-1_sample-1_model",
    }
    fallback = fallback_by_condition["Q7SCE9__domain_only__chitin_DP4"]
    assert fallback.pose_id == "Q7SCE9_NAG4_model"
    assert fallback.run_status == "af3_top_model_fallback"


def test_choose_crystal_representatives_returns_all_cluster_medoids(tmp_path: Path) -> None:
    work_root = _make_work_root(tmp_path)
    options = ProductionRunOptions(
        run_id="test-run",
        output_dir=tmp_path / "output",
        del_variant="a",
        work_root=work_root,
        af3_only=True,
        latest_only=False,
        include_targets=("NAG4",),
    )
    _errors, poses, _summary = _discover_pose_inputs(options)
    prepared_by_id = {
        pose.pose_id: PreparedPose(
            pose=pose,
            case_dir=tmp_path / pose.pose_id,
            normalized_cif=pose.cif_path,
            posebusters_pdb=pose.cif_path,
            privateer_input_cif=pose.cif_path,
            structure={},
        )
        for pose in poses
    }
    clustering_result = ClusteringResult(
        n_clusters=2,
        n_outliers=0,
        outlier_rate=0.0,
        cluster_sizes={0: 1, 1: 1},
        cluster_occupancy={0: 0.5, 1: 0.5},
        outlier_indices=[],
        cluster_labels=np.asarray([0, 1], dtype=int),
        medoids={0: 0, 1: 1},
        medoid_distance_sums={0: 0.0, 1: 0.0},
    )

    representatives = _choose_crystal_representatives(
        clustering_result,
        [pose.pose_id for pose in poses],
        prepared_by_id,
    )

    assert [representative.role for representative in representatives] == [
        "cluster_medoid",
        "cluster_medoid",
    ]
    assert [representative.medoid_pose_id for representative in representatives] == [
        "Q7SCE9_NAG4_seed-1_sample-0_model",
        "Q7SCE9_NAG4_seed-1_sample-1_model",
    ]


def test_top_model_fallback_requires_non_dropped_hard_qc(tmp_path: Path, monkeypatch) -> None:
    work_root = _make_work_root(tmp_path)
    options = ProductionRunOptions(
        run_id="test-run",
        output_dir=tmp_path / "output",
        del_variant="a",
        work_root=work_root,
        af3_only=True,
        latest_only=False,
        include_targets=("NAG4",),
    )
    _errors, fallback_by_condition = _discover_top_model_fallback_inputs(options)
    fallback_pose = fallback_by_condition["Q7SCE9__domain_only__chitin_DP4"]

    def _fake_prepare_pose_case(pose, *, index: int, output_dir: Path):
        del index
        prepared = PreparedPose(
            pose=pose,
            case_dir=output_dir / "case",
            normalized_cif=pose.cif_path,
            posebusters_pdb=pose.cif_path,
            privateer_input_cif=pose.cif_path,
            structure={},
        )
        return {"pose_id": pose.pose_id, "status": "prepared"}, prepared

    monkeypatch.setattr(
        "lpmo_pipeline.analysis.analysis_orchestrator._prepare_pose_case",
        _fake_prepare_pose_case,
    )
    monkeypatch.setattr(
        "lpmo_pipeline.analysis.analysis_orchestrator.run_hard_qc",
        lambda poses, **kwargs: build_qc_report(
            run_id="fallback",
            verdicts=[PoseQCVerdict(pose_id=poses[0].pose_id, status="dropped")],
        ),
    )

    representative, summary = _prepare_hard_qc_passing_top_model_fallback(
        fallback_pose,
        output_dir=tmp_path / "fallback",
        options=options,
    )

    assert representative is None
    assert summary["selected"] is False
    assert summary["skip_reason"] == "fallback_hard_qc_failed"

    monkeypatch.setattr(
        "lpmo_pipeline.analysis.analysis_orchestrator.run_hard_qc",
        lambda poses, **kwargs: build_qc_report(
            run_id="fallback",
            verdicts=[PoseQCVerdict(pose_id=poses[0].pose_id, status="passed")],
        ),
    )

    representative, summary = _prepare_hard_qc_passing_top_model_fallback(
        fallback_pose,
        output_dir=tmp_path / "fallback",
        options=options,
    )

    assert representative is not None
    assert representative.role == "af3_top_model_fallback"
    assert representative.medoid_pose_id == ""
    assert summary["selected"] is True


def test_discover_pose_inputs_preserves_upstream_run_id_for_symlinked_stage(tmp_path: Path) -> None:
    source_cif = _touch(
        tmp_path
        / "source_work"
        / "NAG4"
        / "af3"
        / "runs"
        / "408002"
        / "Q7SCE9_NAG4"
        / "seed-1_sample-0"
        / "Q7SCE9_NAG4_seed-1_sample-0_model.cif"
    )
    staged_cif = (
        tmp_path
        / "staged_work"
        / "NAG4"
        / "af3"
        / "runs"
        / "000001"
        / "Q7SCE9_NAG4"
        / "seed-1_sample-0"
        / "Q7SCE9_NAG4_seed-1_sample-0_model.cif"
    )
    staged_cif.parent.mkdir(parents=True, exist_ok=True)
    staged_cif.symlink_to(source_cif)

    options = ProductionRunOptions(
        run_id="test-run",
        output_dir=tmp_path / "output",
        del_variant="a",
        work_root=tmp_path / "staged_work",
        af3_only=True,
        latest_only=False,
        include_targets=("NAG4",),
    )

    errors, poses, _summary = _discover_pose_inputs(options)

    assert errors == []
    assert len(poses) == 1
    assert poses[0].source_run_id == "408002"
    assert poses[0].discovered_run_id == "000001"
    assert poses[0].cif_path == source_cif.resolve()


def test_run_analysis_core_writes_qc_geometry_and_reports(tmp_path, monkeypatch) -> None:
    work_root = _make_work_root(tmp_path)
    output_dir = tmp_path / "output"
    config = {
        "pipeline_version": "2.1",
        "production": {
            "work_root": str(work_root),
            "af3_only": True,
            "latest_only": False,
            "run_posebusters": False,
            "run_privateer": False,
            "clustering_pilot": {
                "enabled": True,
                "label": "test-pilot",
            },
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

    def _fake_protonate_and_export(normalized_path: Path, output_dir: Path):
        output_dir.mkdir(parents=True, exist_ok=True)
        complex_h = output_dir / "complex_H.pdb"
        ligand_mol2 = output_dir / "ligand_for_prolif.mol2"
        complex_h.write_text("ATOM      1  CA  ALA A   1       0.0     0.0     0.0\nEND\n")
        ligand_mol2.write_text("@<TRIPOS>ATOM\n      1 C1 0.0 0.0 0.0 C.3 1 NAG1 0.0\n@<TRIPOS>BOND\n")
        (output_dir / "protonation_report.json").write_text("{}\n")
        return True, {
            "blockers": [],
            "warnings": [],
            "complex_h_pdb": str(complex_h),
            "ligand_mol2": str(ligand_mol2),
        }

    class _FakeGemmi:
        @staticmethod
        def read_structure(path: str):
            return {"path": path}

    hard_qc_kwargs: dict[str, bool] = {}

    def _fake_run_hard_qc(
        poses,
        run_id: str = "",
        config=None,
        config_path=None,
        *,
        run_posebusters: bool = True,
        run_privateer: bool = True,
        max_workers: int | None = None,
        collect_timing_events: bool = False,
    ):
        hard_qc_kwargs["run_posebusters"] = run_posebusters
        hard_qc_kwargs["run_privateer"] = run_privateer
        hard_qc_kwargs["max_workers"] = max_workers
        hard_qc_kwargs["collect_timing_events"] = collect_timing_events
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

    def _fake_compute_ifp_batch(
        pose_data,
        protein_id: str = "",
        ligand_id: str = "",
        model: str = "",
        protein_chain: str = "A",
        max_workers: int | None = None,
    ) -> IFPBatch:
        del protein_chain, max_workers
        results = []
        matrix = []
        for pose_entry in pose_data:
            results.append(
                IFPResult(
                    pose_id=pose_entry["pose_id"],
                    status="ok",
                    n_residues=1,
                    n_interaction_types=2,
                    residue_names=["NAG1.B|ASN10.A"],
                    interaction_types=["HBDonor", "HBAcceptor"],
                    feature_names=[
                        "NAG1.B|ASN10.A|HBDonor",
                        "NAG1.B|ASN10.A|HBAcceptor",
                    ],
                    fingerprint=[[1, 1]],
                    flat_bitvector=[1, 1],
                    n_total_contacts=2,
                    interaction_counts={"HBDonor": 1, "HBAcceptor": 1},
                )
            )
            matrix.append([1, 1])
        return IFPBatch(
            protein_id=protein_id,
            ligand_id=ligand_id,
            model=model,
            results=results,
            matrix=matrix,
            feature_names=["NAG1.B|ASN10.A|HBDonor", "NAG1.B|ASN10.A|HBAcceptor"],
        )

    class _FakeClusterer:
        def __init__(self, output_dir: Path | None = None, config=None):
            del output_dir, config
            self.config = type(
                "FakeClusterConfig",
                (),
                {
                    "min_cluster_size": 1,
                    "min_samples": None,
                    "cluster_selection_method": "eom",
                    "linkage": "average",
                    "distance_threshold": 0.55,
                },
            )()

        def cluster(self, ifp_matrix, pose_ids):
            del ifp_matrix, pose_ids
            return ClusteringResult(
                n_clusters=1,
                n_outliers=0,
                outlier_rate=0.0,
                cluster_sizes={0: 1},
                cluster_occupancy={0: 1.0},
                outlier_indices=[],
                cluster_labels=np.asarray([0], dtype=int),
                medoids={0: 0},
                medoid_distance_sums={0: 0.0},
            )

        def compute_jaccard_distances(self, ifp_matrix):
            return np.zeros((len(ifp_matrix), len(ifp_matrix)), dtype=float)

    def _fake_compute_condition_convergence(pose_inputs, *, config=None):
        del pose_inputs, config
        return (
            [
                PoseConvergenceMetrics(
                    pose_id="Q7SCE9_NAG4_seed-1_sample-0_model",
                    condition_id="Q7SCE9__domain_only__chitin_DP4",
                    reference_pose_id="Q7SCE9_NAG4_seed-1_sample-0_model",
                    ligand_rmsd_to_reference=0.0,
                    convergent_flag=True,
                )
            ],
            ConditionConvergenceSummary(
                condition_id="Q7SCE9__domain_only__chitin_DP4",
                reference_pose_id="Q7SCE9_NAG4_seed-1_sample-0_model",
                n_qc_pass_poses=1,
                convergence_fraction=1.0,
                median_ligand_rmsd=0.0,
                iqr_ligand_rmsd=0.0,
                low_convergence_flag=False,
            ),
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
    monkeypatch.setattr(
        "lpmo_pipeline.analysis.analysis_orchestrator.protonate_and_export",
        _fake_protonate_and_export,
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
        "lpmo_pipeline.analysis.analysis_orchestrator.compute_ifp_batch",
        _fake_compute_ifp_batch,
    )
    monkeypatch.setattr(
        "lpmo_pipeline.analysis.analysis_orchestrator.compute_condition_convergence",
        _fake_compute_condition_convergence,
    )
    monkeypatch.setattr(
        "lpmo_pipeline.analysis.analysis_orchestrator.AgglomerativeJaccardClusterer",
        _FakeClusterer,
    )
    monkeypatch.setattr(
        "lpmo_pipeline.analysis.analysis_orchestrator.validate",
        lambda instance, schema: None,
    )

    def _fake_run_crystal_reference_screen(
        representative_pose_cif: Path,
        *,
        protein_id: str,
        output_dir: Path,
        representative_pose_id: str = "",
        ligand_id: str = "",
        reference_index_csv=None,
        crystal_root=None,
    ) -> CrystalReferenceScreenReport:
        del representative_pose_cif, output_dir, reference_index_csv, crystal_root
        return CrystalReferenceScreenReport(
            protein_id=protein_id,
            representative_pose_id=representative_pose_id,
            representative_pose_cif="/tmp/representative.cif",
            ligand_id=ligand_id,
            representative_pose_ifp_result_json="/tmp/representative_ifp.json",
            representative_pose_pose_ifp_table_tsv="/tmp/representative_pose_ifp.tsv",
            representative_pose_ifp_matrix_csv="/tmp/representative_ifp_matrix.csv",
            comparisons=[
                CrystalReferencePoseComparison(
                    pdb_code="5ACI",
                    source_cif="/tmp/5ACI.cif",
                    prepared_subset_cif="/tmp/5ACI_subset.cif",
                    representative_pose_id=representative_pose_id,
                    representative_pose_ifp_status="ok",
                    crystal_ifp_status="ok",
                    crystal_ifp_result_json="/tmp/5ACI_ifp.json",
                    crystal_pose_ifp_table_tsv="/tmp/5ACI_ifp.tsv",
                    crystal_ifp_matrix_csv="/tmp/5ACI_ifp.csv",
                    crystal_ifp_contact_eligible=True,
                    crystal_n_non_vdw_interactions=2,
                    crystal_n_non_vdw_contact_residues=1,
                    ifp_comparison_eligible=True,
                    crystal_geometry={
                        "pose_id": "5ACI",
                        "model": "crystal",
                        "protein_id": protein_id,
                        "ligand_id": ligand_id,
                        "Cu_C1_distance": 3.2,
                        "Cu_C4_distance": 4.8,
                        "geometry_status_C1": "geometry_plausible",
                        "geometry_status_C4": "geometry_computable_implausible",
                    },
                    pocket_residues=[10, 11, 12],
                    pocket_rmsd=1.5,
                    pocket_rmsd_below_threshold=True,
                    ifp_tanimoto=0.25,
                    status="ok",
                )
            ],
        )

    def _fake_write_crystal_reference_screen_report(report, output_path: Path) -> None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps({"protein_id": report.protein_id}, indent=2))

    monkeypatch.setattr(
        "lpmo_pipeline.analysis.analysis_orchestrator.run_crystal_reference_screen",
        _fake_run_crystal_reference_screen,
    )
    monkeypatch.setattr(
        "lpmo_pipeline.analysis.analysis_orchestrator.write_crystal_reference_screen_report",
        _fake_write_crystal_reference_screen_report,
    )

    result = run_analysis_core(config, output_dir, del_variant="del_a")

    assert result.success is True
    assert result.n_discovered == 2
    assert result.n_prepared == 2
    assert result.n_analyzed == 1
    assert result.crystal_anchoring_stage_completed is True
    assert result.n_crystal_anchoring_conditions == 1
    assert result.n_crystal_anchoring_errors == 0
    assert hard_qc_kwargs == {
        "run_posebusters": False,
        "run_privateer": False,
        "max_workers": 1,
        "collect_timing_events": False,
    }
    assert geometry_calls == ["Q7SCE9_NAG4_seed-1_sample-0_model"]

    assert result.qc_report_path is not None and result.qc_report_path.exists()
    assert result.pose_manifest_tsv_path is not None and result.pose_manifest_tsv_path.exists()
    assert result.pose_confidence_tsv_path is not None and result.pose_confidence_tsv_path.exists()
    assert result.structure_index_tsv_path is not None and result.structure_index_tsv_path.exists()
    assert result.qc_attrition_tsv_path is not None and result.qc_attrition_tsv_path.exists()
    assert result.pose_geometry_tsv_path is not None and result.pose_geometry_tsv_path.exists()
    assert result.pose_ifp_table_tsv_path is not None and result.pose_ifp_table_tsv_path.exists()
    assert result.pose_residue_contact_tsv_path is not None and result.pose_residue_contact_tsv_path.exists()
    assert result.pose_convergence_tsv_path is not None and result.pose_convergence_tsv_path.exists()
    assert result.condition_convergence_summary_tsv_path is not None and result.condition_convergence_summary_tsv_path.exists()
    assert result.cluster_assignments_tsv_path is not None and result.cluster_assignments_tsv_path.exists()
    assert result.medoid_manifest_tsv_path is not None and result.medoid_manifest_tsv_path.exists()
    assert result.condition_cluster_summary_tsv_path is not None and result.condition_cluster_summary_tsv_path.exists()
    assert result.cluster_ifp_signature_tsv_path is not None and result.cluster_ifp_signature_tsv_path.exists()
    assert result.cluster_residue_signature_tsv_path is not None and result.cluster_residue_signature_tsv_path.exists()
    assert result.cluster_signatures_json_path is not None and result.cluster_signatures_json_path.exists()
    assert result.cluster_annotation_stage_completed is True
    assert result.protein_condition_residue_scores_tsv_path is not None and result.protein_condition_residue_scores_tsv_path.exists()
    assert result.protein_residue_regio_delta_tsv_path is not None and result.protein_residue_regio_delta_tsv_path.exists()
    assert result.condition_patch_summary_tsv_path is not None and result.condition_patch_summary_tsv_path.exists()
    assert result.protein_patch_summary_tsv_path is not None and result.protein_patch_summary_tsv_path.exists()
    assert result.condition_table_tsv_path is not None and result.condition_table_tsv_path.exists()
    assert result.protein_summary_table_tsv_path is not None and result.protein_summary_table_tsv_path.exists()
    assert result.crystal_anchor_tsv_path is not None and result.crystal_anchor_tsv_path.exists()
    assert result.crystal_geometry_tsv_path is not None and result.crystal_geometry_tsv_path.exists()
    assert (
        result.crystal_ifp_diagnostic_summary_tsv_path is not None
        and result.crystal_ifp_diagnostic_summary_tsv_path.exists()
    )
    assert result.metrics_csv_path is not None and result.metrics_csv_path.exists()
    assert result.summary_json_path is not None and result.summary_json_path.exists()
    assert result.report_html_path is not None and result.report_html_path.exists()
    assert result.summary_path.exists()
    assert list(output_dir.rglob("geometry_debug.pdb")) == []

    pose_geometry_lines = result.pose_geometry_tsv_path.read_text().splitlines()
    assert len(pose_geometry_lines) == 2
    assert "Q7SCE9_NAG4_seed-1_sample-0_model" in pose_geometry_lines[1]

    with open(result.pose_manifest_tsv_path, newline="") as handle:
        manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(manifest_rows) == 2
    assert manifest_rows[0]["pose_id"] == "Q7SCE9_NAG4_seed-1_sample-0_model"
    assert manifest_rows[0]["condition_id"] == "Q7SCE9__domain_only__chitin_DP4"
    assert manifest_rows[0]["parse_status"] == "prepared"
    assert manifest_rows[0]["qc_status"] == "passed"
    assert manifest_rows[1]["qc_status"] == "dropped"
    assert manifest_rows[1]["analysis_status"] == "skipped_dropped"

    with open(result.pose_confidence_tsv_path, newline="") as handle:
        confidence_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(confidence_rows) == 2
    assert confidence_rows[0]["confidence_json_status"] == "missing"

    with open(result.structure_index_tsv_path, newline="") as handle:
        structure_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(structure_rows) == 2
    assert structure_rows[0]["normalized_cif"].endswith("normalized.cif")
    assert structure_rows[0]["complex_h_pdb"].endswith("complex_H.pdb")

    with open(result.qc_attrition_tsv_path, newline="") as handle:
        attrition_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert attrition_rows == [
        {
            "condition_id": "Q7SCE9__domain_only__chitin_DP4",
            "protein_id": "Q7SCE9",
            "construct_type": "domain_only",
            "ligand_id": "NAG4",
            "substrate_class": "chitin",
            "dp": "4",
            "n_generated": "2",
            "n_prepared": "2",
            "n_prep_error": "0",
            "n_stage1_hard_fail": "1",
            "n_stage1_soft_flag": "0",
            "n_stage1_pass": "1",
            "hard_fail_rate": "0.5",
            "hard_fail_reason_counts": '{"test_drop": 1}',
            "soft_flag_reason_counts": "{}",
        }
    ]

    with open(result.metrics_csv_path, newline="") as handle:
        metrics_rows = list(csv.DictReader(handle))
    assert len(metrics_rows) == 1
    assert metrics_rows[0]["pose_id"] == "Q7SCE9_NAG4_seed-1_sample-0_model"
    assert metrics_rows[0]["pocket_rmsd_vs_crystal"] == "1.5"
    assert metrics_rows[0]["ifp_similarity_crystal"] == "0.25"
    assert metrics_rows[0]["n_ifp_contacts"] == "2"
    assert metrics_rows[0]["cluster_id"] == "0"

    pose_ifp_lines = result.pose_ifp_table_tsv_path.read_text().splitlines()
    assert len(pose_ifp_lines) == 2
    assert "Q7SCE9_NAG4_seed-1_sample-0_model" in pose_ifp_lines[1]

    with open(result.pose_residue_contact_tsv_path, newline="") as handle:
        residue_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(residue_rows) == 2
    assert residue_rows[0]["pose_id"] == "Q7SCE9_NAG4_seed-1_sample-0_model"
    assert residue_rows[0]["protein_id"] == "Q7SCE9"
    assert residue_rows[0]["condition_id"] == "Q7SCE9__domain_only__chitin_DP4"
    assert residue_rows[0]["ligand_residue_label"] == "NAG1.B"
    assert residue_rows[0]["residue_chain"] == "A"
    assert residue_rows[0]["residue_number"] == "10"
    assert residue_rows[0]["residue_name"] == "ASN"

    with open(result.pose_convergence_tsv_path, newline="") as handle:
        convergence_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert convergence_rows == [
        {
            "pose_id": "Q7SCE9_NAG4_seed-1_sample-0_model",
            "condition_id": "Q7SCE9__domain_only__chitin_DP4",
            "reference_pose_id": "Q7SCE9_NAG4_seed-1_sample-0_model",
            "ligand_rmsd_to_reference": "0.0",
            "convergent_flag": "True",
        }
    ]

    with open(result.condition_convergence_summary_tsv_path, newline="") as handle:
        condition_convergence_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert condition_convergence_rows == [
        {
            "condition_id": "Q7SCE9__domain_only__chitin_DP4",
            "reference_pose_id": "Q7SCE9_NAG4_seed-1_sample-0_model",
            "n_qc_pass_poses": "1",
            "convergence_fraction": "1.0",
            "median_ligand_rmsd": "0.0",
            "iqr_ligand_rmsd": "0.0",
            "low_convergence_flag": "False",
        }
    ]

    with open(result.cluster_assignments_tsv_path, newline="") as handle:
        cluster_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(cluster_rows) == 1
    assert cluster_rows[0]["pose_id"] == "Q7SCE9_NAG4_seed-1_sample-0_model"
    assert cluster_rows[0]["cluster_id"] == "0"

    with open(result.medoid_manifest_tsv_path, newline="") as handle:
        medoid_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(medoid_rows) == 1
    assert medoid_rows[0]["medoid_pose_id"] == "Q7SCE9_NAG4_seed-1_sample-0_model"

    with open(result.condition_cluster_summary_tsv_path, newline="") as handle:
        condition_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(condition_rows) == 1
    assert condition_rows[0]["condition_id"] == "Q7SCE9__domain_only__chitin_DP4"
    assert condition_rows[0]["n_qc_pass_poses"] == "1"
    assert condition_rows[0]["n_ifp_success"] == "1"
    assert condition_rows[0]["n_contact_eligible"] == "1"
    assert condition_rows[0]["contact_eligible_fraction"] == "1.0"
    assert condition_rows[0]["low_specific_contact_fraction"] == "0.0"
    assert condition_rows[0]["n_ifp_clustered"] == "1"
    assert condition_rows[0]["n_noise"] == "0"

    with open(result.cluster_ifp_signature_tsv_path, newline="") as handle:
        cluster_ifp_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(cluster_ifp_rows) == 2
    assert cluster_ifp_rows[0]["condition_id"] == "Q7SCE9__domain_only__chitin_DP4"
    assert {row["feature_name"] for row in cluster_ifp_rows} == {
        "NAG1.B|ASN10.A|HBDonor",
        "NAG1.B|ASN10.A|HBAcceptor",
    }

    with open(result.cluster_residue_signature_tsv_path, newline="") as handle:
        cluster_residue_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(cluster_residue_rows) == 1
    assert cluster_residue_rows[0]["residue_chain"] == "A"
    assert cluster_residue_rows[0]["residue_number"] == "10"
    assert cluster_residue_rows[0]["interaction_type"] == "HBAcceptor,HBDonor"
    assert cluster_residue_rows[0]["contact_frequency"] == "1.0"

    with open(result.protein_condition_residue_scores_tsv_path, newline="") as handle:
        protein_condition_residue_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert protein_condition_residue_rows == [
        {
            "protein_id": "Q7SCE9",
            "condition_id": "Q7SCE9__domain_only__chitin_DP4",
            "construct_type": "domain_only",
            "substrate_class": "chitin",
            "dp": "4",
            "residue_chain": "A",
            "residue_number": "10",
            "residue_name": "ASN",
            "residue_label": "ASN10.A",
            "residue_contact_score": "1.0",
            "c1_weighted_residue_score": "1.0",
            "c4_weighted_residue_score": "1.0",
            "c1_minus_c4_weighted_delta": "0.0",
            "cluster_support_count": "1",
            "total_cluster_occupancy_with_contact": "1.0",
            "max_cluster_residue_frequency": "1.0",
            "is_catalytic_surface_region": "False",
            "is_cbm_region": "False",
            "is_linker_region": "False",
        }
    ]

    with open(result.protein_residue_regio_delta_tsv_path, newline="") as handle:
        protein_regio_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert protein_regio_rows == [
        {
            "protein_id": "Q7SCE9",
            "residue_chain": "A",
            "residue_number": "10",
            "residue_name": "ASN",
            "residue_label": "ASN10.A",
            "n_conditions_with_valid_clusters": "1",
            "n_conditions_with_contact": "1",
            "mean_c1_weighted_residue_score": "1.0",
            "mean_c4_weighted_residue_score": "1.0",
            "c1_minus_c4_weighted_delta": "0.0",
            "is_catalytic_surface_region": "False",
            "is_cbm_region": "False",
            "is_linker_region": "False",
        }
    ]

    with open(result.condition_patch_summary_tsv_path, newline="") as handle:
        condition_patch_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert condition_patch_rows == [
        {
            "protein_id": "Q7SCE9",
            "condition_id": "Q7SCE9__domain_only__chitin_DP4",
            "construct_type": "domain_only",
            "substrate_class": "chitin",
            "dp": "4",
            "any_valid_cluster": "True",
            "n_clusters_considered": "1",
            "total_nonnoise_cluster_occupancy": "1.0",
            "weighted_contact_feature_mass": "2.0",
            "aromatic_contact_fraction": "0.0",
            "polar_contact_fraction": "1.0",
            "charged_contact_fraction": "0.0",
            "hydrophobic_contact_fraction": "0.0",
            "hbond_contact_fraction": "1.0",
            "catalytic_surface_contact_fraction": "0.0",
            "cbm_contact_fraction": "0.0",
            "linker_contact_fraction": "0.0",
        }
    ]

    with open(result.protein_patch_summary_tsv_path, newline="") as handle:
        protein_patch_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert protein_patch_rows == [
        {
            "protein_id": "Q7SCE9",
            "n_conditions_total": "1",
            "n_conditions_with_valid_clusters": "1",
            "mean_aromatic_contact_fraction": "0.0",
            "mean_polar_contact_fraction": "1.0",
            "mean_charged_contact_fraction": "0.0",
            "mean_hydrophobic_contact_fraction": "0.0",
            "mean_hbond_contact_fraction": "1.0",
            "mean_catalytic_surface_contact_fraction": "0.0",
            "mean_cbm_contact_fraction": "0.0",
            "mean_linker_contact_fraction": "0.0",
        }
    ]

    with open(result.condition_table_tsv_path, newline="") as handle:
        condition_table_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(condition_table_rows) == 1
    assert condition_table_rows[0]["condition_id"] == "Q7SCE9__domain_only__chitin_DP4"
    assert condition_table_rows[0]["n_generated"] == "2"
    assert condition_table_rows[0]["n_stage1_hard_fail"] == "1"
    assert condition_table_rows[0]["n_clusters"] == "1"
    assert condition_table_rows[0]["convergence_fraction"] == "1.0"
    assert condition_table_rows[0]["any_valid_cluster"] == "True"
    assert condition_table_rows[0]["n_confidence_rows"] == "2"
    assert condition_table_rows[0]["n_cluster_rows"] == "1"
    assert condition_table_rows[0]["cluster_total_occupancy"] == "1.0"
    assert condition_table_rows[0]["occupancy_weighted_c1_plausible_fraction"] == "1.0"

    with open(result.protein_summary_table_tsv_path, newline="") as handle:
        protein_summary_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(protein_summary_rows) == 1
    assert protein_summary_rows[0]["protein_id"] == "Q7SCE9"
    assert protein_summary_rows[0]["n_conditions"] == "1"
    assert protein_summary_rows[0]["n_generated"] == "2"
    assert protein_summary_rows[0]["n_conditions_with_clusters"] == "1"

    cluster_signatures = json.loads(result.cluster_signatures_json_path.read_text())
    assert cluster_signatures["clusters"][0]["cluster_id"] == 0
    assert cluster_signatures["clusters"][0]["medoid_pose_id"] == "Q7SCE9_NAG4_seed-1_sample-0_model"

    with open(result.crystal_anchor_tsv_path, newline="") as handle:
        crystal_anchor_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(crystal_anchor_rows) == 1
    assert crystal_anchor_rows[0]["condition_id"] == "Q7SCE9__domain_only__chitin_DP4"
    assert crystal_anchor_rows[0]["representative_pose_id"] == "Q7SCE9_NAG4_seed-1_sample-0_model"
    assert crystal_anchor_rows[0]["representative_role"] == "cluster_medoid"
    assert crystal_anchor_rows[0]["medoid_pose_id"] == "Q7SCE9_NAG4_seed-1_sample-0_model"
    assert crystal_anchor_rows[0]["crystal_reference_id"] == "5ACI"
    assert crystal_anchor_rows[0]["local_pocket_rmsd"] == "1.5"
    assert crystal_anchor_rows[0]["contact_overlap_score"] == "0.25"
    assert crystal_anchor_rows[0]["crystal_ifp_contact_eligible"] == "True"
    assert crystal_anchor_rows[0]["ifp_comparison_eligible"] == "True"
    assert crystal_anchor_rows[0]["crystal_Cu_C1_distance"] == "3.2"
    assert crystal_anchor_rows[0]["crystal_geometry_status_C1"] == "geometry_plausible"
    assert crystal_anchor_rows[0]["comparison_status"] == "ok"

    with open(result.crystal_geometry_tsv_path, newline="") as handle:
        crystal_geometry_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(crystal_geometry_rows) == 1
    assert crystal_geometry_rows[0]["pdb_code"] == "5ACI"
    assert crystal_geometry_rows[0]["Cu_C1_distance"] == "3.2"
    assert crystal_geometry_rows[0]["geometry_status_C4"] == "geometry_computable_implausible"

    with open(result.crystal_ifp_diagnostic_summary_tsv_path, newline="") as handle:
        crystal_ifp_summary_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert crystal_ifp_summary_rows[0]["scope"] == "unique_crystal_reference"
    assert crystal_ifp_summary_rows[0]["n_contact_eligible"] == "1"
    assert crystal_ifp_summary_rows[1]["scope"] == "medoid_or_fallback_comparison"
    assert crystal_ifp_summary_rows[1]["contact_eligible_fraction"] == "1.0"

    summary_json = json.loads(result.summary_json_path.read_text())
    assert summary_json["dataset_stats"]["n_total_poses"] == 2
    assert summary_json["geometry_stats"]["cu_c1_range"] == [4.1, 4.1]
    assert summary_json["cluster_stats"]["n_protein_ligand_combinations"] == 1
    assert summary_json["crystal_stats"]["n_comparisons"] == 1
    assert summary_json["crystal_stats"]["mean_best_tanimoto"] == 0.25

    analysis_summary = json.loads(result.summary_path.read_text())
    assert analysis_summary["qc_counts"]["dropped"] == 1
    assert analysis_summary["n_cluster_conditions"] == 1
    assert analysis_summary["primary_clustering"] == {
        "method": "agglomerative_jaccard",
        "metric": "jaccard",
        "linkage": "average",
        "distance_threshold": 0.55,
        "min_cluster_size": 3,
    }
    assert analysis_summary["run_posebusters"] is False
    assert analysis_summary["run_privateer"] is False
    assert analysis_summary["pose_ifp_table_tsv"].endswith("pose_ifp_table.tsv")
    assert analysis_summary["pose_manifest_tsv"].endswith("pose_manifest.tsv")
    assert analysis_summary["pose_confidence_tsv"].endswith("pose_confidence.tsv")
    assert analysis_summary["structure_index_tsv"].endswith("structure_index.tsv")
    assert analysis_summary["qc_attrition_tsv"].endswith("qc_attrition_table.tsv")
    assert analysis_summary["pose_residue_contact_tsv"].endswith("pose_residue_contact_table.tsv")
    assert analysis_summary["pose_convergence_tsv"].endswith("pose_convergence.tsv")
    assert analysis_summary["condition_convergence_summary_tsv"].endswith("condition_convergence_summary.tsv")
    assert analysis_summary["cluster_ifp_signature_tsv"].endswith("cluster_ifp_signature.tsv")
    assert analysis_summary["cluster_residue_signature_tsv"].endswith("cluster_residue_signature.tsv")
    assert analysis_summary["cluster_signatures_json"].endswith("cluster_signatures.json")
    assert analysis_summary["cluster_annotation_stage_completed"] is True
    assert analysis_summary["protein_condition_residue_scores_tsv"].endswith("protein_condition_residue_scores.tsv")
    assert analysis_summary["protein_residue_regio_delta_tsv"].endswith("protein_residue_regio_delta.tsv")
    assert analysis_summary["condition_patch_summary_tsv"].endswith("condition_patch_summary.tsv")
    assert analysis_summary["protein_patch_summary_tsv"].endswith("protein_patch_summary.tsv")
    assert analysis_summary["condition_table_tsv"].endswith("condition_table.tsv")
    assert analysis_summary["protein_summary_table_tsv"].endswith("protein_summary_table.tsv")
    assert analysis_summary["clustering_pilot"]["label"] == "test-pilot"
    assert analysis_summary["clustering_pilot"]["selected_main_feature_count"] == 2
    assert analysis_summary["clustering_pilot"]["n_conditions"] == 1
    assert analysis_summary["clustering_pilot"]["n_conditions_with_contact_eligible_signal"] == 1
    assert analysis_summary["clustering_pilot"]["n_conditions_formal_clustering_allowed"] == 0
    assert analysis_summary["clustering_pilot"]["n_conditions_insufficient_clusterable_signal"] == 1
    assert Path(analysis_summary["clustering_pilot"]["pilot_ifp_interaction_type_prevalence_tsv"]).exists()
    assert Path(analysis_summary["clustering_pilot"]["pilot_ifp_feature_prevalence_tsv"]).exists()
    assert Path(analysis_summary["clustering_pilot"]["pilot_feature_selection_json"]).exists()
    assert Path(analysis_summary["clustering_pilot"]["pilot_clustering_method_summary_tsv"]).exists()
    assert len(analysis_summary["clustering_pilot"]["condition_matrix_entries"]) == 1
    assert Path(analysis_summary["clustering_pilot"]["condition_matrix_entries"][0]["raw_matrix_path"]).exists()
    assert Path(analysis_summary["clustering_pilot"]["condition_matrix_entries"][0]["main_matrix_path"]).exists()
    assert (
        analysis_summary["clustering_pilot"]["condition_matrix_entries"][0]["clustering_status"]
        == "insufficient_clusterable_signal"
    )
    assert analysis_summary["clustering_pilot"]["condition_matrix_entries"][0]["formal_clustering_allowed"] is False
    assert len(analysis_summary["clustering_pilot"]["condition_matrix_entries"][0]["pilot_method_entries"]) == 2
    assert Path(
        analysis_summary["clustering_pilot"]["condition_matrix_entries"][0]["pilot_method_entries"][0]["cluster_assignments_tsv"]
    ).exists()
    assert Path(
        analysis_summary["clustering_pilot"]["condition_matrix_entries"][0]["pilot_method_entries"][1]["medoid_manifest_tsv"]
    ).exists()
    assert analysis_summary["crystal_anchoring_stage_completed"] is True
    assert analysis_summary["crystal_anchor_tsv"].endswith("crystal_anchor_table.tsv")
    assert analysis_summary["crystal_geometry_tsv"].endswith("crystal_geometry_table.tsv")
    assert analysis_summary["crystal_ifp_diagnostic_summary_tsv"].endswith("crystal_ifp_diagnostic_summary.tsv")
    assert analysis_summary["crystal_anchoring_reports"][0]["report_path"].endswith("crystal_reference_screen.json")
    assert analysis_summary["crystal_anchoring_reports"][0]["representative_role"] == "cluster_medoid"
    assert analysis_summary["crystal_anchoring_reports"][0]["best_tanimoto"] == 0.25
    assert analysis_summary["crystal_anchoring_reports"][0]["best_pocket_rmsd"] == 1.5
    dropped_case = next(
        case
        for case in analysis_summary["cases"]
        if case["pose_id"].endswith("sample-1_model")
    )
    assert dropped_case["analysis_status"] == "skipped_dropped"
    analyzed_case = next(
        case
        for case in analysis_summary["cases"]
        if case["pose_id"] == "Q7SCE9_NAG4_seed-1_sample-0_model"
    )
    assert analyzed_case["source_run_id"] == "100001"
    assert analyzed_case["discovered_run_id"] == "100001"
    assert analyzed_case["crystal_anchoring_representative_pose"] is True
    assert analyzed_case["crystal_anchoring_report_path"].endswith("crystal_reference_screen.json")
    assert analyzed_case["crystal_anchoring_best_tanimoto"] == 0.25
    assert analyzed_case["crystal_anchoring_best_pocket_rmsd"] == 1.5
    assert analyzed_case["contact_eligible"] is True
    assert analyzed_case["contact_exclusion_class"] is None
    assert analyzed_case["n_non_vdw_interactions"] == 2
