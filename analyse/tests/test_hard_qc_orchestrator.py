# tests/test_hard_qc_orchestrator.py
"""
Tests for the hard QC orchestrator.
Verifies end-to-end orchestration logic, failure policy, and the
"Ikke-slett regel" (numeric metrics preserved even for dropped poses).
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from lpmo_pipeline.qc.active_site_proximity import ActiveSiteProximityResult
from lpmo_pipeline.qc.custom_geometry_checks import CuHisMeasurement, GeometryResult
from lpmo_pipeline.qc.hard_qc_orchestrator import HardQCInput, run_hard_qc
from lpmo_pipeline.qc.privateer_runner import PrivateerResult
from lpmo_pipeline.qc.posebusters_runner import PoseBustersSingleResult


def _make_input(pose_id: str = "pose_001") -> HardQCInput:
    return HardQCInput(
        pose_id=pose_id,
        mol_pred_path=Path("/fake/pose.pdb"),
        structure=MagicMock(),
    )


@pytest.fixture(autouse=True)
def _mock_proximity_gate():
    """Default Stage 8 behavior in tests: proximity passes unless overridden."""
    with patch("lpmo_pipeline.qc.hard_qc_orchestrator.check_active_site_proximity") as mock_prox:
        mock_prox.return_value = ActiveSiteProximityResult(
            pose_id="pose_default",
            cu_found=True,
            min_cu_ligand_distance=5.0,
            min_cu_c1=4.2,
            min_cu_c4=4.7,
            passed=True,
        )
        yield mock_prox


class TestHardQCOrchestrator:
    """Orchestration logic for hard QC pipeline."""

    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.run_posebusters_single")
    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.check_geometry")
    def test_pose_passes_all_gates(self, mock_geom, mock_pb) -> None:
        """Happy path: PB pass + geometry pass -> status='passed'."""
        mock_pb.return_value = PoseBustersSingleResult(
            pose_id="pose_001", passed=True,
        )
        mock_geom.return_value = GeometryResult(
            pose_id="pose_001", cu_found=True, cu_his_all_in_range=True, passed=True,
        )
        report = run_hard_qc([_make_input()], run_id="test_run")
        assert report.total == 1
        assert report.passed == 1
        assert report.dropped == 0

    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.run_posebusters_single")
    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.check_geometry")
    def test_pose_dropped_by_posebusters(self, mock_geom, mock_pb) -> None:
        """PB critical error -> status='dropped'."""
        mock_pb.return_value = PoseBustersSingleResult(
            pose_id="pose_001", passed=False,
            critical_errors=["internal_steric_clash"],
        )
        mock_geom.return_value = GeometryResult(
            pose_id="pose_001", cu_found=True, cu_his_all_in_range=True, passed=True,
        )
        report = run_hard_qc([_make_input()], run_id="test_run")
        assert report.total == 1
        assert report.dropped == 1

    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.run_posebusters_single")
    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.check_geometry")
    def test_pose_dropped_by_geometry(self, mock_geom, mock_pb) -> None:
        """Cu-His out of range -> status='dropped'."""
        mock_geom.return_value = GeometryResult(
            pose_id="pose_001", cu_found=True, cu_his_all_in_range=False, passed=False,
            failure_reasons=["Cu-His distance out of range"],
        )
        report = run_hard_qc([_make_input()], run_id="test_run")
        assert report.total == 1
        assert report.dropped == 1
        mock_pb.assert_not_called()

    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.run_posebusters_single")
    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.check_geometry")
    def test_metrics_preserved_on_drop(self, mock_geom, mock_pb) -> None:
        """Ikke-slett regel: even dropped poses retain numeric metrics."""
        mock_pb.return_value = PoseBustersSingleResult(
            pose_id="pose_001", passed=False,
            critical_errors=["sanitization"],
        )
        mock_geom.return_value = GeometryResult(
            pose_id="pose_001", cu_found=True, cu_his_all_in_range=True,
            passed=True, min_cu_c1=5.3, min_cu_c4=6.1,
        )
        report = run_hard_qc([_make_input()], run_id="test_run")
        verdict = report.verdicts[0]
        assert verdict.status == "dropped"
        # Canonical min_cu_* keys remain bound to the pre-QC proximity gate.
        assert "min_cu_c1" in verdict.metrics
        assert verdict.metrics["min_cu_c1"] == pytest.approx(4.2)
        assert "min_cu_c4" in verdict.metrics
        assert verdict.metrics["min_cu_c4"] == pytest.approx(4.7)
        # Geometry metrics are still preserved explicitly for downstream use.
        assert verdict.metrics["geometry_min_cu_c1"] == pytest.approx(5.3)
        assert verdict.metrics["geometry_min_cu_c4"] == pytest.approx(6.1)
        assert verdict.cu_geometry["cu_c1_dist_a"] == pytest.approx(5.3)
        assert verdict.cu_geometry["cu_c4_dist_a"] == pytest.approx(6.1)

    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.run_posebusters_single")
    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.check_geometry")
    def test_privateer_skipped_gracefully(self, mock_geom, mock_pb) -> None:
        """Privateer=None should not cause errors; privateer_passed defaults True."""
        mock_pb.return_value = PoseBustersSingleResult(
            pose_id="pose_001", passed=True,
        )
        mock_geom.return_value = GeometryResult(
            pose_id="pose_001", cu_found=True, cu_his_all_in_range=True, passed=True,
        )
        report = run_hard_qc([_make_input()], run_id="test_run")
        verdict = report.verdicts[0]
        assert verdict.privateer_passed is True
        assert verdict.status == "passed"

    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.run_privateer_batch")
    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.run_posebusters_single")
    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.check_geometry")
    def test_privateer_runs_as_batch_for_multiple_eligible_poses(
        self,
        mock_geom,
        mock_pb,
        mock_privateer_batch,
    ) -> None:
        mock_pb.side_effect = [
            PoseBustersSingleResult(pose_id="pose_1", passed=True),
            PoseBustersSingleResult(pose_id="pose_2", passed=True),
        ]
        mock_geom.side_effect = [
            GeometryResult(pose_id="pose_1", cu_found=True, cu_his_all_in_range=True, passed=True),
            GeometryResult(pose_id="pose_2", cu_found=True, cu_his_all_in_range=True, passed=True),
        ]
        mock_privateer_batch.return_value = [
            PrivateerResult(pose_id="pose_1", all_pass=True),
            PrivateerResult(pose_id="pose_2", all_pass=True),
        ]

        inputs = [
            HardQCInput(
                pose_id="pose_1",
                mol_pred_path=Path("/fake/pose1.pdb"),
                structure=MagicMock(),
                privateer_cif_path=Path("/fake/privateer1.cif"),
            ),
            HardQCInput(
                pose_id="pose_2",
                mol_pred_path=Path("/fake/pose2.pdb"),
                structure=MagicMock(),
                privateer_cif_path=Path("/fake/privateer2.cif"),
            ),
        ]

        report = run_hard_qc(inputs, run_id="batch_privateer")

        assert report.total == 2
        assert mock_privateer_batch.call_count == 1
        batch_inputs = mock_privateer_batch.call_args.args[0]
        assert [item.pose_id for item in batch_inputs] == ["pose_1", "pose_2"]
        assert [item.cif_path for item in batch_inputs] == [
            Path("/fake/privateer1.cif"),
            Path("/fake/privateer2.cif"),
        ]

    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.run_posebusters_single")
    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.check_geometry")
    def test_pb_soft_warning_flags_pose(self, mock_geom, mock_pb) -> None:
        """PB soft warning only -> status='flagged', not 'dropped'."""
        mock_pb.return_value = PoseBustersSingleResult(
            pose_id="pose_001", passed=True,
            warnings=["aromatic_ring_flatness"],
        )
        mock_geom.return_value = GeometryResult(
            pose_id="pose_001", cu_found=True, cu_his_all_in_range=True, passed=True,
        )
        report = run_hard_qc([_make_input()], run_id="test_run")
        assert report.flagged == 1
        assert report.dropped == 0

    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.run_posebusters_single")
    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.check_geometry")
    def test_geometry_soft_warning_flags_pose(self, mock_geom, mock_pb) -> None:
        """Preferred QC band miss should flag, not drop."""
        mock_pb.return_value = PoseBustersSingleResult(
            pose_id="pose_001", passed=True,
        )
        mock_geom.return_value = GeometryResult(
            pose_id="pose_001",
            cu_found=True,
            cu_his_all_in_range=True,
            cu_his_all_in_soft_range=False,
            passed=True,
            warnings=["Cu-His distance outside preferred QC range"],
        )

        report = run_hard_qc([_make_input()], run_id="test_run")

        assert report.flagged == 1
        assert report.dropped == 0
        assert any(
            warning.startswith("geometry_soft:")
            for warning in report.verdicts[0].warnings
        )

    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.run_posebusters_single")
    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.check_geometry")
    def test_threshold_config_is_passed_into_qc_checks(
        self,
        mock_geom,
        mock_pb,
        _mock_proximity_gate,
        tmp_path: Path,
    ) -> None:
        config_path = tmp_path / "thresholds.yaml"
        config_path.write_text(
            """
hard_gates:
    active_site_proximity_max_a: 9.5
    cu_his_distance_min_a: 1.5
    cu_his_distance_max_a: 3.0
soft_thresholds:
    cu_c_proximity_max_a: 6.2
geometry_rules:
    his_brace:
        max_search_dist_a: 3.4
qc:
    pre_qc_active_site_max_a: 9.0
    cu_his_min_a: 1.8
    cu_his_max_a: 2.6
""".lstrip()
        )

        mock_pb.return_value = PoseBustersSingleResult(pose_id="pose_001", passed=True)
        mock_geom.return_value = GeometryResult(
            pose_id="pose_001", cu_found=True, cu_his_all_in_range=True, passed=True,
        )

        run_hard_qc([_make_input()], run_id="test_run", config_path=config_path)

        proximity_kwargs = _mock_proximity_gate.call_args.kwargs
        assert proximity_kwargs["hard_cutoff_a"] == pytest.approx(9.5)
        assert proximity_kwargs["cu_c_soft_flag_a"] == pytest.approx(6.2)

        geometry_kwargs = mock_geom.call_args.kwargs
        assert geometry_kwargs["hard_cu_his_min_a"] == pytest.approx(1.5)
        assert geometry_kwargs["hard_cu_his_max_a"] == pytest.approx(3.0)
        assert geometry_kwargs["soft_cu_his_min_a"] == pytest.approx(1.8)
        assert geometry_kwargs["soft_cu_his_max_a"] == pytest.approx(2.6)
        assert geometry_kwargs["cu_c_soft_flag_a"] == pytest.approx(6.2)
        assert geometry_kwargs["his_brace_max_search_a"] == pytest.approx(3.4)

    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.run_posebusters_single")
    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.check_geometry")
    def test_multiple_poses(self, mock_geom, mock_pb) -> None:
        """Batch of 3 poses: 1 pass, 1 flagged, 1 dropped."""
        poses = [_make_input(f"pose_{i}") for i in range(3)]
        mock_pb.side_effect = [
            PoseBustersSingleResult(pose_id="pose_0", passed=True),
            PoseBustersSingleResult(pose_id="pose_1", passed=True,
                                    warnings=["double_bond_flatness"]),
            PoseBustersSingleResult(pose_id="pose_2", passed=False,
                                    critical_errors=["bond_lengths"]),
        ]
        mock_geom.side_effect = [
            GeometryResult(pose_id="pose_0", cu_found=True,
                           cu_his_all_in_range=True, passed=True),
            GeometryResult(pose_id="pose_1", cu_found=True,
                           cu_his_all_in_range=True, passed=True),
            GeometryResult(pose_id="pose_2", cu_found=True,
                           cu_his_all_in_range=True, passed=True),
        ]
        report = run_hard_qc(poses, run_id="batch_test")
        assert report.total == 3
        assert report.passed == 1
        assert report.flagged == 1
        assert report.dropped == 1

    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.run_posebusters_single")
    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.check_geometry")
    def test_pb_infrastructure_failure_drops_pose(self, mock_geom, mock_pb) -> None:
        """If PB throws an unexpected exception, the pose fails closed."""
        mock_pb.side_effect = RuntimeError("unexpected PB crash")
        mock_geom.return_value = GeometryResult(
            pose_id="pose_001", cu_found=True, cu_his_all_in_range=True, passed=True,
        )
        report = run_hard_qc([_make_input()], run_id="test_run")
        assert report.total == 1
        verdict = report.verdicts[0]
        assert verdict.status == "dropped"
        assert verdict.posebusters_passed is False
        assert "posebusters_critical:posebusters_runner_error" in verdict.drop_reasons

    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.run_privateer_batch")
    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.run_posebusters_single")
    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.check_geometry")
    def test_privateer_batch_failure_drops_eligible_pose(
        self,
        mock_geom,
        mock_pb,
        mock_privateer_batch,
    ) -> None:
        """If Privateer batch crashes, eligible poses fail closed."""
        mock_pb.return_value = PoseBustersSingleResult(pose_id="pose_001", passed=True)
        mock_geom.return_value = GeometryResult(
            pose_id="pose_001", cu_found=True, cu_his_all_in_range=True, passed=True,
        )
        mock_privateer_batch.side_effect = RuntimeError("privateer crashed")
        pose = HardQCInput(
            pose_id="pose_001",
            mol_pred_path=Path("/fake/pose.pdb"),
            structure=MagicMock(),
            privateer_cif_path=Path("/fake/privateer.cif"),
        )

        report = run_hard_qc([pose], run_id="test_run")

        verdict = report.verdicts[0]
        assert verdict.status == "dropped"
        assert verdict.privateer_passed is False
        assert any(
            reason.startswith("privateer_runner_error:privateer_batch_runner_error:")
            for reason in verdict.drop_reasons
        )

    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.run_posebusters_single")
    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.check_geometry")
    def test_pre_qc_failure_stops_remaining_qc(self, mock_geom, mock_pb, _mock_proximity_gate) -> None:
        """Stage 8 fail should stop remaining QC, and the pose is dropped."""
        _mock_proximity_gate.return_value = ActiveSiteProximityResult(
            pose_id="pose_001",
            cu_found=True,
            min_cu_ligand_distance=12.1,
            min_cu_c1=11.5,
            min_cu_c4=10.9,
            passed=False,
            failure_reasons=["Ligand too far from active site"],
        )
        mock_pb.return_value = PoseBustersSingleResult(
            pose_id="pose_001",
            passed=True,
        )
        mock_geom.return_value = GeometryResult(
            pose_id="pose_001",
            cu_found=True,
            cu_his_all_in_range=True,
            passed=True,
        )

        report = run_hard_qc([_make_input()], run_id="test_run")
        verdict = report.verdicts[0]

        assert verdict.status == "dropped"
        assert any(r.startswith("active_site_proximity:") for r in verdict.drop_reasons)
        mock_pb.assert_not_called()
        mock_geom.assert_not_called()

    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.run_posebusters_single")
    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.check_geometry")
    def test_pre_qc_failure_preserves_proximity_metrics(self, mock_geom, mock_pb, _mock_proximity_gate) -> None:
        """Ikke-slett regel also applies to Stage 8 pre-QC metrics."""
        _mock_proximity_gate.return_value = ActiveSiteProximityResult(
            pose_id="pose_001",
            cu_found=True,
            min_cu_ligand_distance=13.0,
            min_cu_c1=12.0,
            min_cu_c4=9.5,
            nearest_ligand_atom="B:NAG1:C1",
            passed=False,
            failure_reasons=["Ligand too far from active site"],
        )
        report = run_hard_qc([_make_input()], run_id="test_run")
        verdict = report.verdicts[0]

        assert verdict.status == "dropped"
        assert verdict.metrics["min_cu_ligand_distance"] == pytest.approx(13.0)
        assert verdict.metrics["min_cu_c1"] == pytest.approx(12.0)
        assert verdict.metrics["min_cu_c4"] == pytest.approx(9.5)
        assert verdict.metrics["nearest_ligand_atom"] == "B:NAG1:C1"
        assert verdict.cu_geometry == {}
        mock_pb.assert_not_called()
        mock_geom.assert_not_called()
