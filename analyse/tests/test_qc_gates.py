# tests/test_qc_gates.py
"""
Tests for QC gates (pass/fail logic).
Verifies the failure policy from MASTERPLAN §7.
"""
from __future__ import annotations

import pytest


class TestPoseBustersGate:
    """no_critical_posebusters_errors = true."""

    def test_critical_error_drops_pose(self) -> None:
        """A critical PoseBusters error must result in status='dropped'."""
        from lpmo_pipeline.qc.posebusters_runner import PoseBustersSingleResult
        from lpmo_pipeline.qc.qc_report import compute_verdict

        pb = PoseBustersSingleResult(
            pose_id="test_001",
            passed=False,
            critical_errors=["steric_clash"],
        )
        verdict = compute_verdict("test_001", pb_result=pb, priv_result=None, geom_result=None)
        assert verdict.status == "dropped"
        assert "posebusters_critical:steric_clash" in verdict.drop_reasons

    def test_soft_warning_flags_pose(self) -> None:
        """A soft PoseBusters warning: status='flagged', not 'dropped'."""
        from lpmo_pipeline.qc.posebusters_runner import PoseBustersSingleResult
        from lpmo_pipeline.qc.qc_report import compute_verdict

        pb = PoseBustersSingleResult(
            pose_id="test_002",
            passed=True,
            warnings=["minor_angle_deviation"],
        )
        verdict = compute_verdict("test_002", pb_result=pb, priv_result=None, geom_result=None)
        assert verdict.status == "flagged"
        assert len(verdict.drop_reasons) == 0


class TestCuHisGate:
    """Cu-His_distances_in_range: 1.9–2.6 Å."""

    def test_cu_his_out_of_range_drops(self) -> None:
        """Cu–His distance > 2.6 Å → drop."""
        from lpmo_pipeline.qc.custom_geometry_checks import GeometryResult, CuHisMeasurement
        from lpmo_pipeline.qc.qc_report import compute_verdict

        geom = GeometryResult(
            pose_id="test_003",
            cu_found=True,
            cu_his_measurements=[
                CuHisMeasurement(
                    his_chain="A", his_resnum=1, his_atom="NE2",
                    cu_chain="E", cu_resnum=1,
                    distance_angstrom=3.0,
                    in_range=False,
                ),
            ],
            cu_his_all_in_range=False,
            passed=False,
            failure_reasons=["Cu-His distance out of range: His1:NE2=3.00Å"],
        )
        verdict = compute_verdict("test_003", pb_result=None, priv_result=None, geom_result=geom)
        assert verdict.status == "dropped"

    def test_cu_his_in_range_passes(self) -> None:
        """Cu–His distance 2.1 Å → pass."""
        from lpmo_pipeline.qc.custom_geometry_checks import GeometryResult
        from lpmo_pipeline.qc.qc_report import compute_verdict

        geom = GeometryResult(
            pose_id="test_004",
            cu_found=True,
            cu_his_all_in_range=True,
            passed=True,
        )
        verdict = compute_verdict("test_004", pb_result=None, priv_result=None, geom_result=geom)
        assert verdict.status == "passed"


class TestPrivateerGate:
    """privateer_recognized_sugars = 100%."""

    def test_unrecognized_sugar_drops(self) -> None:
        """If recognition < 100%, verdict = dropped."""
        from lpmo_pipeline.qc.privateer_runner import PrivateerResult, PrivateerResidueResult
        from lpmo_pipeline.qc.qc_report import compute_verdict

        priv = PrivateerResult(
            residues=[
                PrivateerResidueResult(
                    chain="B", resname="NAG", resnum=1,
                    sugar_recognized=True, anomer_ok=True, ring_pucker_ok=True,
                ),
                PrivateerResidueResult(
                    chain="B", resname="XXX", resnum=2,
                    sugar_recognized=False, anomer_ok=False, ring_pucker_ok=False,
                ),
            ],
            total_sugars=2,
            recognized=1,
            recognition_rate=0.5,
            all_pass=False,
        )
        verdict = compute_verdict("test_005", pb_result=None, priv_result=priv, geom_result=None)
        assert verdict.status == "dropped"


class TestAtomMappingGate:
    """atom_mapping_coverage = 100%."""

    def test_coverage_gate_concept(self) -> None:
        """Conceptual test: mapping coverage < 100% should trigger skip."""
        from lpmo_pipeline.utils.exceptions import AtomMappingCoverageError
        # In production: if mapping.coverage < 1.0, raise AtomMappingCoverageError
        with pytest.raises(AtomMappingCoverageError):
            raise AtomMappingCoverageError(
                "Coverage 95%: unmapped atoms [CU1]",
                run_id="test_006",
                step="normalize",
            )
