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
        mock_pb.return_value = PoseBustersSingleResult(
            pose_id="pose_001", passed=True,
        )
        mock_geom.return_value = GeometryResult(
            pose_id="pose_001", cu_found=True, cu_his_all_in_range=False, passed=False,
            failure_reasons=["Cu-His distance out of range"],
        )
        report = run_hard_qc([_make_input()], run_id="test_run")
        assert report.total == 1
        assert report.dropped == 1

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
        # Metrics from geometry check must be preserved
        assert "min_cu_c1" in verdict.metrics
        assert verdict.metrics["min_cu_c1"] == pytest.approx(5.3)
        assert "min_cu_c4" in verdict.metrics
        assert verdict.metrics["min_cu_c4"] == pytest.approx(6.1)

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
    def test_pb_infrastructure_failure_returns_none(self, mock_geom, mock_pb) -> None:
        """If PB throws an unexpected exception, pose is still processed."""
        mock_pb.side_effect = RuntimeError("unexpected PB crash")
        mock_geom.return_value = GeometryResult(
            pose_id="pose_001", cu_found=True, cu_his_all_in_range=True, passed=True,
        )
        report = run_hard_qc([_make_input()], run_id="test_run")
        # PB returned None, so posebusters_passed defaults True in verdict
        assert report.total == 1
        verdict = report.verdicts[0]
        assert verdict.posebusters_passed is True  # None result -> default

    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.run_posebusters_single")
    @patch("lpmo_pipeline.qc.hard_qc_orchestrator.check_geometry")
    def test_pre_qc_failure_skips_pb_and_geometry(self, mock_geom, mock_pb, _mock_proximity_gate) -> None:
        """Stage 8 fail should drop pose before PB/Privateer/geometry."""
        _mock_proximity_gate.return_value = ActiveSiteProximityResult(
            pose_id="pose_001",
            cu_found=True,
            min_cu_ligand_distance=12.1,
            min_cu_c1=11.5,
            min_cu_c4=10.9,
            passed=False,
            failure_reasons=["Ligand too far from active site"],
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
        mock_pb.assert_not_called()
        mock_geom.assert_not_called()
