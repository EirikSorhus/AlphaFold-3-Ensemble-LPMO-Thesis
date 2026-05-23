from __future__ import annotations

import csv
from pathlib import Path
import tempfile

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
            interaction_types=["ImplicitHBDonor", "ImplicitHBAcceptor"],
            feature_names=[
                "NAG4.B|GLU26.A|ImplicitHBDonor",
                "NAG4.B|GLU26.A|ImplicitHBAcceptor",
            ],
            flat_bitvector=[1, 0],
        ),
        IFPResult(
            pose_id="pose-b",
            status="ok",
            interaction_types=["ImplicitHBDonor", "ImplicitHBAcceptor"],
            feature_names=["NAG2.B|SER10.A|ImplicitHBAcceptor"],
            flat_bitvector=[1],
        ),
    ]

    def fake_compute_ifp_single(**_: object) -> IFPResult:
        return fake_results.pop(0)

    monkeypatch.setattr(prolif_ifp, "compute_ifp_single", fake_compute_ifp_single)

    batch = prolif_ifp.compute_ifp_batch(
        [
            {"pose_id": "pose-a", "complex_pdb": Path("a.pdb"), "ligand_pdb": Path("a.pdb")},
            {"pose_id": "pose-b", "complex_pdb": Path("b.pdb"), "ligand_pdb": Path("b.pdb")},
        ]
    )

    assert batch.feature_names == [
        "NAG2.B|SER10.A|ImplicitHBAcceptor",
        "NAG4.B|GLU26.A|ImplicitHBAcceptor",
        "NAG4.B|GLU26.A|ImplicitHBDonor",
    ]
    assert batch.matrix == [[0, 0, 1], [1, 0, 0]]


def test_result_from_dataframe_keeps_ligand_residues_separate() -> None:
    dataframe = pd.DataFrame(
        [[1, 1, 0]],
        columns=pd.MultiIndex.from_tuples(
            [
                ("NAG1.B", "ASN88.A", "ImplicitHBDonor"),
                ("NAG4.B", "ASN88.A", "ImplicitHBDonor"),
                ("NAG4.B", "ASN88.A", "VdWContact"),
            ]
        ),
    )

    result = prolif_ifp._result_from_dataframe(
        pose_id="pose",
        interaction_types=["ImplicitHBDonor", "VdWContact"],
        dataframe=dataframe,
    )

    assert result.status == "ok"
    assert result.n_residues == 2
    assert result.feature_names == [
        "NAG1.B|ASN88.A|ImplicitHBDonor",
        "NAG1.B|ASN88.A|VdWContact",
        "NAG4.B|ASN88.A|ImplicitHBDonor",
        "NAG4.B|ASN88.A|VdWContact",
    ]
    assert result.flat_bitvector == [1, 0, 1, 0]


def test_compute_ifp_single_real_probe_artifacts_if_available() -> None:
    probe_dir = (
        Path(__file__).resolve().parent
        / "tests_results"
        / "protonation_contracts_573050"
        / "analysis_export"
    )
    complex_pdb = probe_dir / "complex_for_prolif.pdb"
    ligand_pdb = probe_dir / "ligand_only_for_prolif.pdb"

    if not complex_pdb.exists() or not ligand_pdb.exists():
        pytest.skip("Real analysis-export probe artifacts are not available in this workspace")

    result = prolif_ifp.compute_ifp_single(
        complex_pdb=complex_pdb,
        ligand_pdb=ligand_pdb,
        pose_id="probe",
    )

    assert result.status == "ok"
    assert result.n_total_contacts >= 1
    assert result.feature_names
    assert len(result.feature_names) == len(result.flat_bitvector)
    assert (
        result.interaction_counts["ImplicitHBAcceptor"]
        + result.interaction_counts["ImplicitHBDonor"]
    ) >= 1
    assert any(feature_name.startswith("NAG") for feature_name in result.feature_names)


def test_compute_ifp_single_real_6ydc_crystal_uses_multiple_glycan_residues_if_available() -> None:
    probe_dir = (
        Path(__file__).resolve().parent
        / "tests_results"
        / "crystal_anchoring_real_cifs_1144620"
        / "crystal_anchoring_output"
        / "references"
        / "6YDC_A0A223GEC9"
        / "analysis_export"
    )
    complex_pdb = probe_dir / "complex_for_prolif.pdb"
    ligand_pdb = probe_dir / "ligand_only_for_prolif.pdb"

    if not complex_pdb.exists() or not ligand_pdb.exists():
        pytest.skip("Real 6YDC crystal analysis-export artifacts are not available in this workspace")

    result = prolif_ifp.compute_ifp_single(
        complex_pdb=complex_pdb,
        ligand_pdb=ligand_pdb,
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
        feature_names=[
            "NAG4.B|GLU26.A|ImplicitHBDonor",
            "NAG4.B|GLU26.A|ImplicitHBAcceptor",
        ],
        flat_bitvector=[1, 0],
        n_total_contacts=1,
        interaction_counts={
            "ImplicitHBDonor": 1,
            "ImplicitHBAcceptor": 0,
            "VdWContact": 2,
        },
        interaction_occurrence_counts={
            "ImplicitHBDonor": 1,
            "ImplicitHBAcceptor": 0,
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
        "n_implicit_hbond_acceptor",
        "n_implicit_hbond_donor",
        "n_vdw_contact",
        "ifp_interaction_counts",
        "ifp_interaction_occurrence_counts",
        "ifp_error",
    ]
    assert rows[0]["pose_id"] == "pose-1"
    assert rows[0]["ifp_vector"] == "[1, 0]"
    assert rows[0]["ifp_feature_names"] == (
        "[\"NAG4.B|GLU26.A|ImplicitHBDonor\", \"NAG4.B|GLU26.A|ImplicitHBAcceptor\"]"
    )
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
        feature_names=["NAG1.B|ASN10.A|ImplicitHBDonor"],
        flat_bitvector=[1],
        n_total_contacts=1,
        interaction_counts={"ImplicitHBDonor": 1, "VdWContact": 0},
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
            "NAG1.B|ASN10.A|ImplicitHBDonor",
            "NAG1.B|ASN10.A|ImplicitHBAcceptor",
            "NAG1.B|ASN10.A|VdWContact",
        ],
        flat_bitvector=[1, 1, 1],
        n_total_contacts=3,
        interaction_counts={"ImplicitHBDonor": 1, "ImplicitHBAcceptor": 1, "VdWContact": 1},
    )

    eligibility = prolif_ifp.evaluate_contact_eligibility(
        result,
        ContactEligibilityRule(min_non_vdw_interactions=2, min_non_vdw_contact_residues=1),
    )

    assert eligibility.eligible is True
    assert eligibility.exclusion_class is None
    assert eligibility.n_non_vdw_interactions == 2
    assert eligibility.n_non_vdw_contact_residues == 1


def _run_implicit_hbond_probe(complex_pdb: Path, ligand_pdb: Path) -> tuple[int, set[str]]:
    import MDAnalysis as mda
    import prolif as plf
    from prolif import Molecule

    try:
        from prolif.io.protein_helper import ProteinHelper
    except ImportError as exc:
        pytest.skip(f"Installed ProLIF does not expose implicit H-bond helpers: {exc}")

    from rdkit import Chem
    from rdkit.Chem.rdDetermineBonds import DetermineConnectivity

    try:
        fingerprint = plf.Fingerprint(
            interactions=["ImplicitHBAcceptor", "ImplicitHBDonor"],
            count=True,
        )
    except Exception as exc:
        pytest.skip(f"Installed ProLIF does not expose implicit H-bond interactions: {exc}")

    with tempfile.TemporaryDirectory() as tmp_dir:
        protein_pdb = Path(tmp_dir) / "protein_only.pdb"
        universe = mda.Universe(str(complex_pdb))
        protein_atoms = universe.select_atoms("protein and chainID A")

        assert protein_atoms.n_atoms > 0
        protein_atoms.write(str(protein_pdb))

        protein_rdkit = Chem.MolFromPDBFile(str(protein_pdb), removeHs=False)
        assert protein_rdkit is not None

        protein_helper = ProteinHelper()
        try:
            standardized_protein = protein_helper.standardize_protein(protein_rdkit)
            protein = Molecule.from_rdkit(standardized_protein)
        except Exception:
            protein = Molecule.from_rdkit(protein_rdkit)
            protein = protein_helper.standardize_protein(protein)

        ligand_rdkit = Chem.MolFromPDBFile(str(ligand_pdb), removeHs=False)
        assert ligand_rdkit is not None
        try:
            DetermineConnectivity(ligand_rdkit, useHueckel=True)
        except TypeError:
            DetermineConnectivity(ligand_rdkit)
        for atom in ligand_rdkit.GetAtoms():
            atom.SetNoImplicit(False)
        ligand = Molecule.from_rdkit(ligand_rdkit)

        fingerprint.run_from_iterable([ligand], protein)
        dataframe = fingerprint.to_dataframe()

    active_contacts = int(dataframe.to_numpy().sum()) if not dataframe.empty else 0
    interaction_types = (
        {str(value) for value in dataframe.columns.get_level_values(-1)}
        if getattr(dataframe.columns, "nlevels", 0) >= 3
        else set()
    )
    return active_contacts, interaction_types


def test_implicit_hbond_real_crystal_artifacts_if_supported() -> None:
    probe_dirs = {
        "7PXW": (
            Path(__file__).resolve().parent
            / "tests_results"
            / "crystal_fil_test_1156504"
            / "cases"
            / "0001_A0A0S2GKZ1_CEL4"
            / "crystal_anchoring"
            / "references"
            / "7PXW_A0A0S2GKZ1"
            / "protonated"
        ),
        "5ACI": (
            Path(__file__).resolve().parent
            / "tests_results"
            / "crystal_fil_test_1156504"
            / "cases"
            / "0002_A0A0S2GKZ1_CEL6"
            / "crystal_anchoring"
            / "references"
            / "5ACI_A0A0S2GKZ1"
            / "protonated"
        ),
    }

    required_paths = [
        artifact_path
        for probe_dir in probe_dirs.values()
        for artifact_path in (
            probe_dir / "for_posebusters.pdb",
            probe_dir / "ligand_only_for_prolif.pdb",
        )
    ]
    missing_paths = [path for path in required_paths if not path.exists()]
    if missing_paths:
        pytest.skip(
            "Real crystal implicit-H probe artifacts are not available in this workspace: "
            + ", ".join(str(path) for path in missing_paths)
        )

    total_active_contacts = 0
    observed_interactions: set[str] = set()
    for label, probe_dir in probe_dirs.items():
        active_contacts, interaction_types = _run_implicit_hbond_probe(
            probe_dir / "for_posebusters.pdb",
            probe_dir / "ligand_only_for_prolif.pdb",
        )

        assert interaction_types <= {"ImplicitHBAcceptor", "ImplicitHBDonor"}
        assert active_contacts >= 1, f"No implicit H-bond contacts detected for {label}"

        total_active_contacts += active_contacts
        observed_interactions.update(interaction_types)

    assert total_active_contacts >= len(probe_dirs)
    assert observed_interactions
