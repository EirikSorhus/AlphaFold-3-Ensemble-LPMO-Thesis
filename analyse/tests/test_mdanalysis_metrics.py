from __future__ import annotations

import math

import pytest

from lpmo_pipeline.analysis.cluster_signatures import compute_cluster_signatures
from lpmo_pipeline.analysis.mdanalysis_metrics import (
    GEOMETRY_STATUS_COMPUTABLE_IMPLAUSIBLE,
    GEOMETRY_STATUS_HIGHLY_PLAUSIBLE,
    GEOMETRY_STATUS_NOT_COMPUTABLE,
    GeometryThresholds,
    POSE_GEOMETRY_COLUMNS,
    _append_virtual_geometry_atoms_to_pdb_text,
    _assign_geometry_status,
    _infer_bonded_heavy_neighbors,
    _score_oxyl_h_distance,
    compute_pose_metrics_from_structure,
    write_pose_geometry_tsv,
)
from lpmo_pipeline.report.build_metrics_csv import merge_pose_record


class _FakeAtom:
    def __init__(self, name: str, element: str, coords: tuple[float, float, float]) -> None:
        self.name = name
        self.element = type("Element", (), {"name": element})()
        self.pos = type("Position", (), {"x": coords[0], "y": coords[1], "z": coords[2]})()


class _FakeResidue(list):
    def __init__(self, name: str, seqnum: int, atoms: list[_FakeAtom]) -> None:
        super().__init__(atoms)
        self.name = name
        self.seqid = type("SeqId", (), {"num": seqnum})()


class _FakeChain(list):
    def __init__(self, name: str, residues: list[_FakeResidue]) -> None:
        super().__init__(residues)
        self.name = name


class _FakeStructure(list):
    def __init__(self, models: list[list[_FakeChain]], connections: list[object] | None = None) -> None:
        super().__init__(models)
        self.connections = connections or []


class _FakeAtomAddress:
    def __init__(self, chain_name: str, resname: str, resnum: int, atom_name: str) -> None:
        self.chain_name = chain_name
        self.atom_name = atom_name
        self.res_id = type(
            "ResId",
            (),
            {
                "name": resname,
                "seqid": type("SeqId", (), {"num": resnum})(),
            },
        )()


class _FakeConnection:
    def __init__(self, partner1: _FakeAtomAddress, partner2: _FakeAtomAddress) -> None:
        self.partner1 = partner1
        self.partner2 = partner2


def _thresholds() -> GeometryThresholds:
    return GeometryThresholds(
        oxyl_h_min_a=1.5,
        oxyl_h_max_a=4.0,
        oxyl_h_optimum_a=2.1,
        oxyl_h_score_normalization=1.9,
        oxyl_h_highly_plausible_min_a=1.8,
        oxyl_h_highly_plausible_max_a=2.5,
        his_brace_max_search_a=3.0,
        cu_oxyl_bond_length_a=1.9,
        c_h_bond_length_a=1.1,
    )


def test_default_geometry_thresholds_load_operational_oxyl_bond_length() -> None:
    from lpmo_pipeline.analysis.mdanalysis_metrics import load_geometry_thresholds

    load_geometry_thresholds.cache_clear()
    thresholds = load_geometry_thresholds()

    assert thresholds.cu_oxyl_bond_length_a == pytest.approx(1.9)


def _build_structure(include_c4: bool = True) -> _FakeStructure:
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
        [_FakeAtom("NE2", "N", (0.0, 0.0, 1.9))],
    )
    cu = _FakeResidue("CU", 1, [_FakeAtom("CU", "Cu", (0.0, 0.0, 0.0))])

    glycan_atoms = [
        _FakeAtom("O5", "O", (-1.05, 0.95, 0.0)),
        _FakeAtom("C1", "C", (-2.20, 1.40, 0.0)),
        _FakeAtom("O1", "O", (-3.35, 2.15, 0.0)),
        _FakeAtom("C2", "C", (-3.05, 0.15, 0.0)),
        _FakeAtom("C3", "C", (-4.45, 0.55, 0.0)),
        _FakeAtom("C5", "C", (-3.65, 2.75, 0.0)),
        _FakeAtom("O4", "O", (-5.95, 1.45, 0.0)),
    ]
    if include_c4:
        glycan_atoms.append(_FakeAtom("C4", "C", (-4.80, 1.95, 0.0)))

    glycan = _FakeResidue("NAG", 3, glycan_atoms)
    return _FakeStructure(
        [[_FakeChain("A", [his1, his78]), _FakeChain("E", [cu]), _FakeChain("B", [glycan])]]
    )


def _build_linked_glycan_structure() -> tuple[_FakeStructure, _FakeResidue, _FakeResidue, _FakeAtom]:
    his1 = _FakeResidue(
        "HIS",
        1,
        [
            _FakeAtom("N", "N", (2.0, 0.0, 0.0)),
            _FakeAtom("ND1", "N", (0.0, 2.0, 0.0)),
        ],
    )
    his78 = _FakeResidue("HIS", 78, [_FakeAtom("NE2", "N", (0.0, 0.0, 1.9))])
    cu = _FakeResidue("CU", 1, [_FakeAtom("CU", "Cu", (0.0, 0.0, 0.0))])
    glycan3 = _FakeResidue(
        "GLC",
        3,
        [
            _FakeAtom("O5", "O", (-1.05, 0.95, 0.0)),
            _FakeAtom("C1", "C", (-2.20, 1.40, 0.0)),
            _FakeAtom("C2", "C", (-3.05, 0.15, 0.0)),
            _FakeAtom("C3", "C", (-4.45, 0.55, 0.0)),
            _FakeAtom("C4", "C", (-4.80, 1.95, 0.0)),
            _FakeAtom("C5", "C", (-3.65, 2.75, 0.0)),
        ],
    )
    glycan4 = _FakeResidue(
        "GLC",
        4,
        [
            _FakeAtom("O4", "O", (-3.35, 2.15, 0.0)),
            _FakeAtom("C4", "C", (-2.40, 2.95, 0.8)),
        ],
    )
    structure = _FakeStructure(
        [[_FakeChain("A", [his1, his78]), _FakeChain("E", [cu]), _FakeChain("B", [glycan3, glycan4])]],
        connections=[
            _FakeConnection(
                _FakeAtomAddress("B", "GLC", 3, "C1"),
                _FakeAtomAddress("B", "GLC", 4, "O4"),
            )
        ],
    )
    c1_atom = next(atom for atom in glycan3 if atom.name == "C1")
    return structure, glycan3, glycan4, c1_atom


def test_compute_pose_metrics_from_structure_populates_pose_geometry_fields() -> None:
    metrics = compute_pose_metrics_from_structure(
        structure=_build_structure(),
        pose_id="pose_geo",
        thresholds=_thresholds(),
        protein_id="P12345",
        ligand_id="NAG4",
        model="af3",
    )

    assert metrics.pose_id == "pose_geo"
    assert metrics.proximal_sugar_id == "B:NAG3"
    assert metrics.brace_integrity_flag is True
    assert metrics.cu_c1_distance is not None
    assert metrics.cu_c4_distance is not None
    assert metrics.repositioned_cu_coordinates is not None
    assert metrics.virtual_oxyl_coordinates is not None
    assert metrics.virtual_h_c1_coordinates is not None
    assert metrics.virtual_h_c4_coordinates is not None
    assert metrics.geometry_status_c1 in {
        GEOMETRY_STATUS_COMPUTABLE_IMPLAUSIBLE,
        GEOMETRY_STATUS_HIGHLY_PLAUSIBLE,
        "geometry_plausible",
    }
    assert metrics.geometry_status_c4 in {
        GEOMETRY_STATUS_COMPUTABLE_IMPLAUSIBLE,
        GEOMETRY_STATUS_HIGHLY_PLAUSIBLE,
        "geometry_plausible",
    }
    assert metrics.ring_normal_vs_brace_normal is not None
    assert metrics.sugar_face_orientation in {
        "ambiguous",
        "ring_normal_toward_cu",
        "ring_normal_away_from_cu",
    }
    assert metrics.cu_oxyl_h_c1_angle is not None
    assert metrics.to_row()["Cu_oxyl_H_C1_angle"] == metrics.cu_oxyl_h_c1_angle


def test_cu_oxyl_h_angle_uses_oxyl_vertex_not_target_carbon() -> None:
    structure = _build_structure()
    metrics = compute_pose_metrics_from_structure(
        structure=structure,
        pose_id="pose_geo",
        thresholds=_thresholds(),
    )
    glycan = structure[0][2][0]
    c1_atom = next(atom for atom in glycan if atom.name == "C1")

    def _angle(a: tuple[float, float, float], b: tuple[float, float, float], c: tuple[float, float, float]) -> float:
        ba = [a[index] - b[index] for index in range(3)]
        bc = [c[index] - b[index] for index in range(3)]
        dot = sum(ba[index] * bc[index] for index in range(3))
        norm_ba = math.sqrt(sum(value * value for value in ba))
        norm_bc = math.sqrt(sum(value * value for value in bc))
        return math.degrees(math.acos(max(-1.0, min(1.0, dot / (norm_ba * norm_bc)))))

    assert metrics.repositioned_cu_coordinates is not None
    assert metrics.virtual_oxyl_coordinates is not None
    assert metrics.virtual_h_c1_coordinates is not None

    expected_cu_oxyl_h = _angle(
        metrics.repositioned_cu_coordinates,
        metrics.virtual_oxyl_coordinates,
        metrics.virtual_h_c1_coordinates,
    )
    old_oxyl_c_h = _angle(
        metrics.virtual_oxyl_coordinates,
        (c1_atom.pos.x, c1_atom.pos.y, c1_atom.pos.z),
        metrics.virtual_h_c1_coordinates,
    )
    assert metrics.cu_oxyl_h_c1_angle == pytest.approx(expected_cu_oxyl_h)
    assert metrics.cu_oxyl_h_c1_angle != pytest.approx(old_oxyl_c_h)


def test_compute_pose_metrics_marks_missing_target_not_computable() -> None:
    metrics = compute_pose_metrics_from_structure(
        structure=_build_structure(include_c4=False),
        pose_id="pose_missing_c4",
        thresholds=_thresholds(),
    )

    assert metrics.cu_c4_distance is None
    assert metrics.oxyl_h_c4_distance is None
    assert metrics.geometry_status_c4 == GEOMETRY_STATUS_NOT_COMPUTABLE


def test_score_and_status_logic_follow_locked_thresholds() -> None:
    thresholds = _thresholds()

    best_score = _score_oxyl_h_distance(2.1, thresholds)
    assert best_score == pytest.approx(1.0)
    assert _assign_geometry_status(2.1, best_score, thresholds) == GEOMETRY_STATUS_HIGHLY_PLAUSIBLE

    outside_score = _score_oxyl_h_distance(4.2, thresholds)
    assert math.isnan(outside_score)
    assert (
        _assign_geometry_status(4.2, outside_score, thresholds)
        == GEOMETRY_STATUS_COMPUTABLE_IMPLAUSIBLE
    )


def test_write_pose_geometry_tsv_writes_expected_columns(tmp_path) -> None:
    metrics = compute_pose_metrics_from_structure(
        structure=_build_structure(),
        pose_id="pose_geo",
        thresholds=_thresholds(),
    )
    output_path = tmp_path / "pose_geometry.tsv"

    write_pose_geometry_tsv([metrics], output_path)

    lines = output_path.read_text().splitlines()
    assert lines[0].split("\t") == POSE_GEOMETRY_COLUMNS
    assert lines[1].split("\t")[0] == "pose_geo"


def test_infer_bonded_heavy_neighbors_includes_glycosidic_cross_residue_neighbor() -> None:
    structure, glycan3, glycan4, c1_atom = _build_linked_glycan_structure()

    neighbors = _infer_bonded_heavy_neighbors(
        structure=structure,
        chain_name="B",
        residue=glycan3,
        target_atom=c1_atom,
    )

    neighbor_names = [atom.name for atom in neighbors]
    assert "O4" in neighbor_names
    linked_o4 = next(atom for atom in glycan4 if atom.name == "O4")
    assert linked_o4 in neighbors


def test_append_virtual_geometry_atoms_to_pdb_text_includes_debug_atoms() -> None:
    metrics = compute_pose_metrics_from_structure(
        structure=_build_structure(),
        pose_id="pose_geo",
        thresholds=_thresholds(),
    )

    base_pdb = (
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )
    debug_pdb = _append_virtual_geometry_atoms_to_pdb_text(base_pdb, metrics)

    assert "REMARK 700 Virtual geometry atoms below are test-only debug output" in debug_pdb
    assert " CUX " in debug_pdb
    assert " HC1 " in debug_pdb
    assert " HC4 " in debug_pdb
    assert " OX " in debug_pdb
    assert debug_pdb.rstrip().endswith("END")


def test_merge_pose_record_accepts_pose_geometry_row_keys() -> None:
    metrics = compute_pose_metrics_from_structure(
        structure=_build_structure(),
        pose_id="pose_geo",
        thresholds=_thresholds(),
    )

    record = merge_pose_record(
        pose_id="pose_geo",
        run_id="run_1",
        protein_id="P12345",
        ligand_id="NAG4",
        model="af3",
        seed=1,
        qc_verdict={"status": "passed", "warnings": []},
        geometry=metrics.to_row(),
    )

    assert record["min_cu_c1"] == metrics.cu_c1_distance
    assert record["min_cu_c4"] == metrics.cu_c4_distance


def test_cluster_signatures_accept_pose_geometry_row_keys() -> None:
    metrics_a = compute_pose_metrics_from_structure(
        structure=_build_structure(),
        pose_id="pose_a",
        thresholds=_thresholds(),
    )
    metrics_b = compute_pose_metrics_from_structure(
        structure=_build_structure(),
        pose_id="pose_b",
        thresholds=_thresholds(),
    )

    signatures = compute_cluster_signatures(
        cluster_labels={"pose_a": 0, "pose_b": 0},
        geometry_metrics={
            "pose_a": metrics_a.to_row(),
            "pose_b": metrics_b.to_row(),
        },
        ifp_matrix={"pose_a": [1, 0], "pose_b": [1, 1]},
        feature_names=["His1_HBond", "His78_Hydrophobic"],
        total_poses=2,
    )

    assert len(signatures) == 1
    assert signatures[0].mean_cu_c1 == pytest.approx(metrics_a.cu_c1_distance)
    assert signatures[0].mean_cu_c4 == pytest.approx(metrics_a.cu_c4_distance)
