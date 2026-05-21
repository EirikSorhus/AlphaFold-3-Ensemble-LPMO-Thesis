from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
import pytest

from lpmo_pipeline.analysis import prolif_ifp
from lpmo_pipeline.analysis.prolif_ifp import (
    ContactEligibilityRule,
    IFPResult,
    compute_tanimoto_similarity,
)


def test_compute_ifp_single_missing_inputs_returns_status() -> None:
    result = prolif_ifp.compute_ifp_single(
        Path("/tmp/does-not-exist-complex.pdb"),
        Path("/tmp/does-not-exist-ligand.mol2"),
        pose_id="missing",
    )

    assert result.status == "input_missing"
    assert result.n_total_contacts == 0
    assert result.flat_bitvector == []
    assert result.error is not None


def test_compute_ifp_batch_aligns_union_of_features(monkeypatch: pytest.MonkeyPatch) -> None:
    fake_results = [
        IFPResult(
            pose_id="pose-a",
            status="ok",
            interaction_types=["HBDonor", "HBAcceptor"],
            feature_names=[
                "NAG4.B|GLU26.A|HBDonor",
                "NAG4.B|GLU26.A|HBAcceptor",
            ],
            flat_bitvector=[1, 0],
        ),
        IFPResult(
            pose_id="pose-b",
            status="ok",
            interaction_types=["HBDonor", "HBAcceptor"],
            feature_names=["NAG2.B|SER10.A|HBAcceptor"],
            flat_bitvector=[1],
        ),
    ]

    def fake_compute_ifp_single(**_: object) -> IFPResult:
        return fake_results.pop(0)

    monkeypatch.setattr(prolif_ifp, "compute_ifp_single", fake_compute_ifp_single)

    batch = prolif_ifp.compute_ifp_batch(
        [
            {"pose_id": "pose-a", "complex_pdb": Path("a.pdb"), "ligand_mol2": Path("a.mol2")},
            {"pose_id": "pose-b", "complex_pdb": Path("b.pdb"), "ligand_mol2": Path("b.mol2")},
        ]
    )

    assert batch.feature_names == [
        "NAG2.B|SER10.A|HBAcceptor",
        "NAG4.B|GLU26.A|HBAcceptor",
        "NAG4.B|GLU26.A|HBDonor",
    ]
    assert batch.matrix == [[0, 0, 1], [1, 0, 0]]


def test_result_from_dataframe_keeps_ligand_residues_separate() -> None:
    dataframe = pd.DataFrame(
        [[1, 1, 0]],
        columns=pd.MultiIndex.from_tuples(
            [
                ("NAG1.B", "ASN88.A", "HBDonor"),
                ("NAG4.B", "ASN88.A", "HBDonor"),
                ("NAG4.B", "ASN88.A", "Hydrophobic"),
            ]
        ),
    )

    result = prolif_ifp._result_from_dataframe(
        pose_id="pose",
        interaction_types=["HBDonor", "Hydrophobic"],
        dataframe=dataframe,
    )

    assert result.status == "ok"
    assert result.n_residues == 2
    assert result.feature_names == [
        "NAG1.B|ASN88.A|HBDonor",
        "NAG1.B|ASN88.A|Hydrophobic",
        "NAG4.B|ASN88.A|HBDonor",
        "NAG4.B|ASN88.A|Hydrophobic",
    ]
    assert result.flat_bitvector == [1, 0, 1, 0]


def test_ligand_annotation_uses_companion_pdb_chains_for_duplicate_residue_numbers(
    tmp_path: Path,
) -> None:
    ligand_mol2 = tmp_path / "ligand_for_prolif.mol2"
    ligand_mol2.write_text("@<TRIPOS>MOLECULE\nfixture\n")
    (tmp_path / "ligand_only_for_prolif.pdb").write_text(
        "\n".join(
            [
                "HETATM    1  O1  BGC B   1       0.000   0.000   0.000  1.00  0.00           O  ",
                "HETATM    2  O1  BGC C   1       5.000   0.000   0.000  1.00  0.00           O  ",
                "END",
            ]
        )
    )

    class FakeAtom:
        def __init__(
            self,
            index: int,
            atomic_number: int,
            symbol: str,
            neighbors: list["FakeAtom"] | None = None,
        ) -> None:
            self._index = index
            self._atomic_number = atomic_number
            self._symbol = symbol
            self._neighbors = neighbors or []

        def GetAtomicNum(self) -> int:
            return self._atomic_number

        def GetSymbol(self) -> str:
            return self._symbol

        def GetIdx(self) -> int:
            return self._index

        def GetNeighbors(self) -> list["FakeAtom"]:
            return self._neighbors

    first = FakeAtom(0, 8, "O")
    second = FakeAtom(1, 8, "O")
    hydrogen = FakeAtom(2, 1, "H", neighbors=[second])

    class FakeLigand:
        def GetAtoms(self) -> list[FakeAtom]:
            return [first, second, hydrogen]

    atom_records = [
        {
            "atom_name": "O1",
            "x": 0.0,
            "y": 0.0,
            "z": 0.0,
            "substructure_id": 1,
            "substructure_name": "BGC1",
        },
        {
            "atom_name": "O1",
            "x": 5.0,
            "y": 0.0,
            "z": 0.0,
            "substructure_id": 1,
            "substructure_name": "BGC1",
        },
        {
            "atom_name": "H1",
            "x": 5.8,
            "y": 0.0,
            "z": 0.0,
            "substructure_id": 1,
            "substructure_name": "BGC1",
        },
    ]

    annotations = prolif_ifp._resolve_ligand_annotations(
        FakeLigand(),
        atom_records,
        ligand_mol2,
        "B",
    )

    assert annotations == [("BGC", 1, "B"), ("BGC", 1, "C"), ("BGC", 1, "C")]


def test_compute_ifp_single_real_probe_artifacts_if_available() -> None:
    probe_dir = (
        Path(__file__).resolve().parent
        / "tests_results"
        / "protonation_contracts_573050"
        / "protonated"
    )
    complex_pdb = probe_dir / "complex_H.pdb"
    ligand_mol2 = probe_dir / "ligand_for_prolif.mol2"

    if not complex_pdb.exists() or not ligand_mol2.exists():
        pytest.skip("Real protonation probe artifacts are not available in this workspace")

    result = prolif_ifp.compute_ifp_single(
        complex_pdb=complex_pdb,
        ligand_mol2=ligand_mol2,
        pose_id="probe",
    )

    assert result.status == "ok"
    assert result.n_total_contacts >= 1
    assert result.feature_names
    assert len(result.feature_names) == len(result.flat_bitvector)
    assert result.interaction_counts["HBDonor"] >= 1
    assert any(feature_name.startswith("NAG") for feature_name in result.feature_names)


def test_compute_ifp_single_real_6ydc_crystal_uses_multiple_glycan_residues_if_available() -> None:
    probe_dir = (
        Path(__file__).resolve().parent
        / "tests_results"
        / "crystal_anchoring_real_cifs_1144620"
        / "crystal_anchoring_output"
        / "references"
        / "6YDC_A0A223GEC9"
        / "protonated"
    )
    complex_pdb = probe_dir / "complex_H.pdb"
    ligand_mol2 = probe_dir / "ligand_for_prolif.mol2"

    if not complex_pdb.exists() or not ligand_mol2.exists():
        pytest.skip("Real 6YDC crystal anchoring artifacts are not available in this workspace")

    result = prolif_ifp.compute_ifp_single(
        complex_pdb=complex_pdb,
        ligand_mol2=ligand_mol2,
        pose_id="6YDC",
    )

    active_ligand_residues = {
        feature_name.split("|", maxsplit=1)[0]
        for feature_name, value in zip(result.feature_names, result.flat_bitvector, strict=True)
        if value
    }

    assert result.status == "ok"
    assert len(active_ligand_residues) > 1
    assert "BGC4.C" in active_ligand_residues


def test_compute_tanimoto_similarity_handles_empty_union() -> None:
    assert compute_tanimoto_similarity([], []) == 0.0
    assert compute_tanimoto_similarity([1, 0, 1], [1, 1, 0]) == pytest.approx(1 / 3)


def test_write_pose_ifp_table_emits_expected_columns(tmp_path: Path) -> None:
    output_path = tmp_path / "pose_ifp_table.tsv"
    result = IFPResult(
        pose_id="pose-1",
        status="ok",
        feature_names=["NAG4.B|GLU26.A|HBDonor", "NAG4.B|GLU26.A|HBAcceptor"],
        flat_bitvector=[1, 0],
        n_total_contacts=1,
        interaction_counts={
            "HBDonor": 1,
            "HBAcceptor": 0,
            "Hydrophobic": 0,
            "PiStacking": 0,
            "Anionic": 0,
            "Cationic": 0,
            "CationPi": 0,
            "PiCation": 0,
            "VdWContact": 2,
        },
    )

    prolif_ifp.write_pose_ifp_table([result], output_path)

    with output_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert list(rows[0].keys()) == [
        "pose_id",
        "ifp_generation_status",
        "ifp_vector",
        "ifp_feature_names",
        "n_total_contacts",
        "n_hbond_donor",
        "n_hbond_acceptor",
        "n_hydrophobic",
        "n_aromatic",
        "n_anionic",
        "n_cationic",
        "n_cation_pi",
        "n_pi_cation",
        "n_vdw_contact",
        "ifp_interaction_counts",
        "ifp_error",
    ]
    assert rows[0]["pose_id"] == "pose-1"
    assert rows[0]["ifp_vector"] == "[1, 0]"
    assert rows[0]["ifp_feature_names"] == "[\"NAG4.B|GLU26.A|HBDonor\", \"NAG4.B|GLU26.A|HBAcceptor\"]"
    assert rows[0]["n_vdw_contact"] == "2"


def test_load_contact_eligibility_rule_defaults_to_main_rule() -> None:
    rule = prolif_ifp.load_contact_eligibility_rule()

    assert rule.min_non_vdw_interactions == 2
    assert rule.min_non_vdw_contact_residues == 1


def test_evaluate_contact_eligibility_marks_null_ifp() -> None:
    result = IFPResult(
        pose_id="pose-null",
        status="zero_contacts",
        feature_names=[],
        flat_bitvector=[],
        n_total_contacts=0,
        interaction_counts={"VdWContact": 0},
    )

    eligibility = prolif_ifp.evaluate_contact_eligibility(
        result,
        ContactEligibilityRule(min_non_vdw_interactions=2, min_non_vdw_contact_residues=1),
    )

    assert eligibility.eligible is False
    assert eligibility.exclusion_class == "null_ifp"
    assert eligibility.n_non_vdw_interactions == 0


def test_evaluate_contact_eligibility_marks_vdw_only() -> None:
    result = IFPResult(
        pose_id="pose-vdw",
        status="ok",
        feature_names=["NAG1.B|ASN10.A|VdWContact"],
        flat_bitvector=[1],
        n_total_contacts=1,
        interaction_counts={"VdWContact": 1},
    )

    eligibility = prolif_ifp.evaluate_contact_eligibility(
        result,
        ContactEligibilityRule(min_non_vdw_interactions=2, min_non_vdw_contact_residues=1),
    )

    assert eligibility.eligible is False
    assert eligibility.exclusion_class == "vdw_only"
    assert eligibility.n_vdw_interactions == 1
    assert eligibility.n_non_vdw_contact_residues == 0


def test_evaluate_contact_eligibility_marks_low_specific_contact() -> None:
    result = IFPResult(
        pose_id="pose-low-specific",
        status="ok",
        feature_names=["NAG1.B|ASN10.A|HBDonor"],
        flat_bitvector=[1],
        n_total_contacts=1,
        interaction_counts={"HBDonor": 1, "VdWContact": 0},
    )

    eligibility = prolif_ifp.evaluate_contact_eligibility(
        result,
        ContactEligibilityRule(min_non_vdw_interactions=2, min_non_vdw_contact_residues=1),
    )

    assert eligibility.eligible is False
    assert eligibility.exclusion_class == "low_specific_contact"
    assert eligibility.n_non_vdw_interactions == 1
    assert eligibility.n_non_vdw_contact_residues == 1


def test_evaluate_contact_eligibility_accepts_two_non_vdw_interactions() -> None:
    result = IFPResult(
        pose_id="pose-eligible",
        status="ok",
        feature_names=[
            "NAG1.B|ASN10.A|HBDonor",
            "NAG1.B|ASN10.A|HBAcceptor",
            "NAG1.B|ASN10.A|VdWContact",
        ],
        flat_bitvector=[1, 1, 1],
        n_total_contacts=3,
        interaction_counts={"HBDonor": 1, "HBAcceptor": 1, "VdWContact": 1},
    )

    eligibility = prolif_ifp.evaluate_contact_eligibility(
        result,
        ContactEligibilityRule(min_non_vdw_interactions=2, min_non_vdw_contact_residues=1),
    )

    assert eligibility.eligible is True
    assert eligibility.exclusion_class is None
    assert eligibility.n_non_vdw_interactions == 2
    assert eligibility.n_non_vdw_contact_residues == 1
