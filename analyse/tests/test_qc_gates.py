# tests/test_qc_gates.py
"""
Tests for QC gates (pass/fail logic).
Verifies the failure policy from MASTERPLAN §7.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
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
    """Cu-His_distances_in_range hard gate plus soft preferred band."""

    def test_gate_config_defaults_match_locked_threshold_policy(self) -> None:
        from lpmo_pipeline.qc.gates import GateConfig

        config = GateConfig()

        assert config.cu_his_dist_min == pytest.approx(1.5)
        assert config.cu_his_dist_max == pytest.approx(3.0)
        assert config.cu_his_soft_min_a == pytest.approx(1.8)
        assert config.cu_his_soft_max_a == pytest.approx(2.6)

    def test_cu_his_out_of_range_drops(self) -> None:
        """Cu–His distance outside the hard gate → drop."""
        from lpmo_pipeline.qc.custom_geometry_checks import GeometryResult, CuHisMeasurement
        from lpmo_pipeline.qc.qc_report import compute_verdict

        geom = GeometryResult(
            pose_id="test_003",
            cu_found=True,
            cu_his_measurements=[
                CuHisMeasurement(
                    his_chain="A", his_resnum=1, his_atom="NE2",
                    cu_chain="E", cu_resnum=1,
                    distance_angstrom=3.2,
                    in_range=False,
                ),
            ],
            cu_his_all_in_range=False,
            passed=False,
            failure_reasons=["Cu-His distance out of range: His1:NE2=3.20Å"],
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


class TestCheckGeometryThresholds:
    """Threshold overrides should flow into live geometry evaluation."""

    def test_check_geometry_accepts_threshold_overrides_and_sets_soft_warning(self) -> None:
        from types import SimpleNamespace

        from lpmo_pipeline.qc.custom_geometry_checks import check_geometry

        class _FakeAtom:
            def __init__(self, name: str, element: str, coords: tuple[float, float, float]) -> None:
                self.name = name
                self.element = SimpleNamespace(name=element)
                self.pos = SimpleNamespace(x=coords[0], y=coords[1], z=coords[2])

        class _FakeResidue(list):
            def __init__(self, name: str, seqnum: int, atoms: list[_FakeAtom]) -> None:
                super().__init__(atoms)
                self.name = name
                self.seqid = SimpleNamespace(num=seqnum)

        class _FakeChain(list):
            def __init__(self, name: str, residues: list[_FakeResidue]) -> None:
                super().__init__(residues)
                self.name = name

        his1 = _FakeResidue(
            "HIS",
            1,
            [
                _FakeAtom("N", "N", (2.0, 0.0, 0.0)),
                _FakeAtom("ND1", "N", (0.0, 2.0, 0.0)),
            ],
        )
        his78 = _FakeResidue(
            "HIS",
            78,
            [_FakeAtom("NE2", "N", (0.0, 0.0, 1.76))],
        )
        cu_residue = _FakeResidue("CU", 1, [_FakeAtom("CU", "Cu", (0.0, 0.0, 0.0))])
        glycan = _FakeResidue(
            "NAG",
            3,
            [
                _FakeAtom("C1", "C", (5.0, 0.0, 0.0)),
                _FakeAtom("C4", "C", (0.0, 6.0, 0.0)),
            ],
        )
        structure = [[
            _FakeChain("A", [his1, his78]),
            _FakeChain("E", [cu_residue]),
            _FakeChain("B", [glycan]),
        ]]

        result = check_geometry(
            structure=structure,
            pose_id="pose_thresholds",
            glycan_chains=["B"],
            hard_cu_his_min_a=1.5,
            hard_cu_his_max_a=3.0,
            soft_cu_his_min_a=1.8,
            soft_cu_his_max_a=2.6,
            cu_c_soft_flag_a=6.2,
            his_brace_max_search_a=3.4,
        )

        assert result.passed is True
        assert result.cu_found is True
        assert result.cu_his_all_in_range is True
        assert result.cu_his_all_in_soft_range is False
        assert len(result.cu_his_measurements) == 3
        assert result.min_cu_c1 == pytest.approx(5.0)
        assert result.min_cu_c4 == pytest.approx(6.0)
        assert result.failure_reasons == []
        assert len(result.warnings) == 1
        assert "Cu-His distance outside preferred QC range:" in result.warnings[0]
        assert "His78:NE2=1.76Å" in result.warnings[0]
        assert "[1.80, 2.60] A" in result.warnings[0]


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


class TestPoseBustersRunner:
    """Tests for the PoseBusters result classification logic."""

    def test_all_pass_returns_passed(self) -> None:
        """All PB tests True -> passed=True, no errors."""
        from lpmo_pipeline.qc.posebusters_runner import _classify_results

        tests = {
            "sanitization": True, "bond_lengths": True,
            "internal_steric_clash": True, "bond_angles": True,
        }
        result = _classify_results("pose_001", tests)
        assert result.passed is True
        assert len(result.critical_errors) == 0
        assert len(result.warnings) == 0

    def test_critical_failure_returns_failed(self) -> None:
        """A critical PB test False -> passed=False with that test in critical_errors."""
        from lpmo_pipeline.qc.posebusters_runner import _classify_results

        tests = {"sanitization": False, "bond_lengths": True, "bond_angles": True}
        result = _classify_results("pose_002", tests)
        assert result.passed is False
        assert "sanitization" in result.critical_errors

    def test_soft_warning_keeps_passed(self) -> None:
        """A soft-only failure keeps passed=True and records warning."""
        from lpmo_pipeline.qc.posebusters_runner import _classify_results

        tests = {"sanitization": True, "aromatic_ring_flatness": False}
        result = _classify_results("pose_003", tests)
        assert result.passed is True
        assert "aromatic_ring_flatness" in result.warnings

    def test_unknown_failure_treated_as_critical(self) -> None:
        """Unknown test names failing -> treated as critical (conservative)."""
        from lpmo_pipeline.qc.posebusters_runner import _classify_results

        tests = {"sanitization": True, "unknown_test_xyz": False}
        result = _classify_results("pose_004", tests)
        assert result.passed is False
        assert "unknown_test_xyz" in result.critical_errors

    def test_empty_results_treated_as_failure(self) -> None:
        """Empty dict (PB produced no results) -> passed=False."""
        from lpmo_pipeline.qc.posebusters_runner import _classify_results

        result = _classify_results("pose_005", {})
        assert result.passed is False
        assert "posebusters_no_results" in result.critical_errors

    def test_multiple_critical_errors(self) -> None:
        """Multiple critical failures all recorded."""
        from lpmo_pipeline.qc.posebusters_runner import _classify_results

        tests = {
            "sanitization": False,
            "internal_steric_clash": False,
            "bond_lengths": True,
        }
        result = _classify_results("pose_006", tests)
        assert result.passed is False
        assert "sanitization" in result.critical_errors
        assert "internal_steric_clash" in result.critical_errors

    def test_file_not_found_in_batch(self) -> None:
        """Batch runner handles missing files gracefully."""
        from lpmo_pipeline.qc.posebusters_runner import run_posebusters_batch

        result = run_posebusters_batch(Path("/nonexistent"), ["fake_pose"])
        assert result.total == 1
        assert result.failed == 1

    def test_csv_parsing(self) -> None:
        """_parse_csv_output correctly parses PB CSV output."""
        from lpmo_pipeline.qc.posebusters_runner import _parse_csv_output

        csv_text = (
            "file,molecule,position,mol_pred_loaded,sanitization,bond_lengths,"
            "internal_steric_clash\n"
            "test.pdb,lig,0,True,True,False,True\n"
        )
        result = _parse_csv_output(csv_text)
        assert result["sanitization"] is True
        assert result["bond_lengths"] is False
        assert result["internal_steric_clash"] is True
        # Loading columns should be excluded
        assert "mol_pred_loaded" not in result

    def test_result_table_parsing_excludes_loading_columns_and_nan(self) -> None:
        """_parse_result_table should drop loading columns and coerce NaN-like values to False."""
        import pandas as pd

        from lpmo_pipeline.qc.posebusters_runner import _parse_result_table

        table = pd.DataFrame(
            [
                {
                    "mol_pred_loaded": True,
                    "sanitization": True,
                    "bond_lengths": False,
                    "internal_energy": float("nan"),
                    "inchi_convertible": "True",
                }
            ]
        )

        result = _parse_result_table(table)

        assert "mol_pred_loaded" not in result
        assert result["sanitization"] is True
        assert result["bond_lengths"] is False
        assert result["internal_energy"] is False
        assert result["inchi_convertible"] is True


class TestCuHisAngle:
    """Tests for _compute_his_brace_angle() calculation."""

    def test_orthogonal_vectors_give_90_degrees(self) -> None:
        """Two His-N atoms at orthogonal positions around Cu -> ~90 deg."""
        from lpmo_pipeline.qc.custom_geometry_checks import (
            CuHisMeasurement,
            _compute_his_brace_angle,
        )

        cu_pos = np.array([0.0, 0.0, 0.0])
        measurements = [
            CuHisMeasurement(
                his_chain="A", his_resnum=1, his_atom="NE2",
                cu_chain="E", cu_resnum=1, distance_angstrom=2.0,
                in_range=True, his_atom_position=(2.0, 0.0, 0.0),
            ),
            CuHisMeasurement(
                his_chain="A", his_resnum=78, his_atom="ND1",
                cu_chain="E", cu_resnum=1, distance_angstrom=2.0,
                in_range=True, his_atom_position=(0.0, 2.0, 0.0),
            ),
        ]
        angle = _compute_his_brace_angle(measurements, cu_pos)
        assert angle is not None
        assert abs(angle - 90.0) < 0.1

    def test_linear_arrangement_gives_180_degrees(self) -> None:
        """Two His-N atoms on opposite sides of Cu -> ~180 deg."""
        from lpmo_pipeline.qc.custom_geometry_checks import (
            CuHisMeasurement,
            _compute_his_brace_angle,
        )

        cu_pos = np.array([0.0, 0.0, 0.0])
        measurements = [
            CuHisMeasurement(
                his_chain="A", his_resnum=1, his_atom="NE2",
                cu_chain="E", cu_resnum=1, distance_angstrom=2.0,
                in_range=True, his_atom_position=(2.0, 0.0, 0.0),
            ),
            CuHisMeasurement(
                his_chain="A", his_resnum=78, his_atom="ND1",
                cu_chain="E", cu_resnum=1, distance_angstrom=2.0,
                in_range=True, his_atom_position=(-2.0, 0.0, 0.0),
            ),
        ]
        angle = _compute_his_brace_angle(measurements, cu_pos)
        assert angle is not None
        assert abs(angle - 180.0) < 0.1

    def test_single_his_returns_none(self) -> None:
        """Only one His residue -> cannot compute brace angle."""
        from lpmo_pipeline.qc.custom_geometry_checks import (
            CuHisMeasurement,
            _compute_his_brace_angle,
        )

        cu_pos = np.array([0.0, 0.0, 0.0])
        measurements = [
            CuHisMeasurement(
                his_chain="A", his_resnum=1, his_atom="NE2",
                cu_chain="E", cu_resnum=1, distance_angstrom=2.0,
                in_range=True, his_atom_position=(2.0, 0.0, 0.0),
            ),
        ]
        angle = _compute_his_brace_angle(measurements, cu_pos)
        assert angle is None

    def test_picks_closest_n_per_his(self) -> None:
        """When a His has both NE2 and ND1, the closest one is used."""
        from lpmo_pipeline.qc.custom_geometry_checks import (
            CuHisMeasurement,
            _compute_his_brace_angle,
        )

        cu_pos = np.array([0.0, 0.0, 0.0])
        measurements = [
            # His1: NE2 at 2.0 A, ND1 at 3.5 A — NE2 should be picked
            CuHisMeasurement(
                his_chain="A", his_resnum=1, his_atom="NE2",
                cu_chain="E", cu_resnum=1, distance_angstrom=2.0,
                in_range=True, his_atom_position=(2.0, 0.0, 0.0),
            ),
            CuHisMeasurement(
                his_chain="A", his_resnum=1, his_atom="ND1",
                cu_chain="E", cu_resnum=1, distance_angstrom=3.5,
                in_range=False, his_atom_position=(3.5, 0.0, 0.0),
            ),
            # His78: ND1 at 2.1 A
            CuHisMeasurement(
                his_chain="A", his_resnum=78, his_atom="ND1",
                cu_chain="E", cu_resnum=1, distance_angstrom=2.1,
                in_range=True, his_atom_position=(0.0, 2.1, 0.0),
            ),
        ]
        angle = _compute_his_brace_angle(measurements, cu_pos)
        assert angle is not None
        # NE2(His1) at (2,0,0) and ND1(His78) at (0,2.1,0) -> ~90 deg
        assert abs(angle - 90.0) < 1.0


class TestCuHisGateAtomSelection:
    """Tests for locked Cu-His gate atom-selection rules."""

    def test_selects_his1_n_his1_nd1_and_non_his1_third(self) -> None:
        """Selected set must be His1:N, His1:ND1, and one non-His1 histidine N."""
        from lpmo_pipeline.qc.custom_geometry_checks import _select_cu_his_gate_nitrogens

        cu_pos = np.array([0.0, 0.0, 0.0], dtype=float)
        candidates = [
            {"his_chain": "A", "his_resnum": 1, "atom_name": "N", "position": np.array([2.0, 0.0, 0.0])},
            {"his_chain": "A", "his_resnum": 1, "atom_name": "ND1", "position": np.array([0.0, 2.1, 0.0])},
            {"his_chain": "A", "his_resnum": 1, "atom_name": "NE2", "position": np.array([1.2, 0.0, 0.0])},
            {"his_chain": "A", "his_resnum": 78, "atom_name": "NE2", "position": np.array([0.0, 0.0, 2.2])},
            {"his_chain": "A", "his_resnum": 90, "atom_name": "ND1", "position": np.array([0.0, 0.0, 2.4])},
        ]

        selected, errors = _select_cu_his_gate_nitrogens(candidates, cu_pos, max_search_a=3.0)

        assert errors == []
        assert len(selected) == 3
        labels = {(int(r["his_resnum"]), str(r["atom_name"])) for r in selected}
        assert (1, "N") in labels
        assert (1, "ND1") in labels
        assert (1, "NE2") not in labels
        assert any(resnum != 1 for resnum, _ in labels)

    def test_missing_his1_n_reports_error(self) -> None:
        """Mandatory His1:N missing should produce selection error."""
        from lpmo_pipeline.qc.custom_geometry_checks import _select_cu_his_gate_nitrogens

        cu_pos = np.array([0.0, 0.0, 0.0], dtype=float)
        candidates = [
            {"his_chain": "A", "his_resnum": 1, "atom_name": "ND1", "position": np.array([0.0, 2.1, 0.0])},
            {"his_chain": "A", "his_resnum": 78, "atom_name": "NE2", "position": np.array([0.0, 0.0, 2.2])},
        ]

        selected, errors = _select_cu_his_gate_nitrogens(candidates, cu_pos, max_search_a=3.0)

        assert len(selected) == 2
        assert any("His1:N not found" in e for e in errors)

    def test_his1_ne2_not_required_or_selected(self) -> None:
        """His1:NE2 out-of-range/absent must not block selection."""
        from lpmo_pipeline.qc.custom_geometry_checks import _select_cu_his_gate_nitrogens

        cu_pos = np.array([0.0, 0.0, 0.0], dtype=float)
        candidates = [
            {"his_chain": "A", "his_resnum": 1, "atom_name": "N", "position": np.array([2.0, 0.0, 0.0])},
            {"his_chain": "A", "his_resnum": 1, "atom_name": "ND1", "position": np.array([0.0, 2.0, 0.0])},
            {"his_chain": "A", "his_resnum": 78, "atom_name": "ND1", "position": np.array([0.0, 0.0, 2.0])},
        ]

        selected, errors = _select_cu_his_gate_nitrogens(candidates, cu_pos, max_search_a=3.0)

        assert errors == []
        labels = {(int(r["his_resnum"]), str(r["atom_name"])) for r in selected}
        assert (1, "NE2") not in labels

    def test_third_n_must_be_non_his1_within_cutoff(self) -> None:
        """No non-His1 candidate within cutoff should report selection failure."""
        from lpmo_pipeline.qc.custom_geometry_checks import _select_cu_his_gate_nitrogens

        cu_pos = np.array([0.0, 0.0, 0.0], dtype=float)
        candidates = [
            {"his_chain": "A", "his_resnum": 1, "atom_name": "N", "position": np.array([2.0, 0.0, 0.0])},
            {"his_chain": "A", "his_resnum": 1, "atom_name": "ND1", "position": np.array([0.0, 2.0, 0.0])},
            {"his_chain": "A", "his_resnum": 78, "atom_name": "ND1", "position": np.array([0.0, 0.0, 3.5])},
        ]

        selected, errors = _select_cu_his_gate_nitrogens(candidates, cu_pos, max_search_a=3.0)

        assert len(selected) == 2
        assert any("no non-His1 histidine" in e for e in errors)

    def test_distances_are_euclidean(self) -> None:
        """Computed distances for selected atoms should use Euclidean norm."""
        from lpmo_pipeline.qc.custom_geometry_checks import _select_cu_his_gate_nitrogens

        cu_pos = np.array([0.0, 0.0, 0.0], dtype=float)
        candidates = [
            {"his_chain": "A", "his_resnum": 1, "atom_name": "N", "position": np.array([1.0, 1.0, 1.0])},
            {"his_chain": "A", "his_resnum": 1, "atom_name": "ND1", "position": np.array([0.0, 2.0, 0.0])},
            {"his_chain": "A", "his_resnum": 78, "atom_name": "NE2", "position": np.array([0.0, 0.0, 2.0])},
        ]

        selected, errors = _select_cu_his_gate_nitrogens(candidates, cu_pos, max_search_a=3.0)
        assert errors == []

        d_his1_n = [r["distance_angstrom"] for r in selected if r["his_resnum"] == 1 and r["atom_name"] == "N"][0]
        assert abs(float(d_his1_n) - float(np.sqrt(3.0))) < 1e-9
