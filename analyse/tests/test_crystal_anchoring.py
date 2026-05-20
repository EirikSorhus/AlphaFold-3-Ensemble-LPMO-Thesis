import json
from pathlib import Path

import pytest

from lpmo_pipeline.analysis import crystal_anchoring


ANALYSE_ROOT = Path(__file__).resolve().parents[1]
REFERENCE_INDEX_CSV = ANALYSE_ROOT / "input_data" / "pdb_structure_data.csv"
CRYSTAL_ROOT = ANALYSE_ROOT / "crystal_structures"
WORK_CORE_ROOT = ANALYSE_ROOT.parent / "structure_pipeline" / "work_core"


def _write_pocket_test_pdb(
    pdb_path: Path,
    *,
    translation: tuple[float, float, float] = (0.0, 0.0, 0.0),
    rotate_z_deg: float = 0.0,
    perturb_residue_2: tuple[float, float, float] | None = None,
    residue_numbers: tuple[int, int, int] = (1, 2, 3),
    residue_names: tuple[str, str, str] = ("ALA", "GLY", "SER"),
    ligand_point: tuple[float, float, float] | None = (1.0, 0.6, 0.0),
    copper_point: tuple[float, float, float] | None = None,
) -> None:
    import math

    tx, ty, tz = translation
    angle = math.radians(rotate_z_deg)
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)

    def _transform(point: tuple[float, float, float]) -> tuple[float, float, float]:
        x, y, z = point
        rx = x * cos_a - y * sin_a
        ry = x * sin_a + y * cos_a
        return (rx + tx, ry + ty, z + tz)

    residue_positions = ((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0))
    residues = {
        residue_number: _transform(position)
        for residue_number, position in zip(residue_numbers, residue_positions, strict=True)
    }
    if perturb_residue_2 is not None:
        px, py, pz = perturb_residue_2
        residue_number = residue_numbers[1]
        x, y, z = residues[residue_number]
        residues[residue_number] = (x + px, y + py, z + pz)

    ligand = _transform(ligand_point) if ligand_point is not None else None
    copper = _transform(copper_point) if copper_point is not None else None

    pdb_lines: list[str] = []
    for atom_serial, (residue_number, residue_name) in enumerate(
        zip(residue_numbers, residue_names, strict=True),
        start=1,
    ):
        x, y, z = residues[residue_number]
        pdb_lines.append(
            f"ATOM  {atom_serial:5d}  CA  {residue_name:>3} A{residue_number:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 10.00           C"
        )
    atom_serial = len(pdb_lines) + 1
    if ligand is not None:
        pdb_lines.append(
            f"HETATM{atom_serial:5d}  C1  BGC B   1    "
            f"{ligand[0]:8.3f}{ligand[1]:8.3f}{ligand[2]:8.3f}  1.00 10.00           C"
        )
        atom_serial += 1
    if copper is not None:
        pdb_lines.append(
            f"HETATM{atom_serial:5d} CU   CU  C   1    "
            f"{copper[0]:8.3f}{copper[1]:8.3f}{copper[2]:8.3f}  1.00 10.00          CU"
        )

    pdb_path.write_text(
        "\n".join(
            [*pdb_lines, "END", ""]
        )
    )


def _write_fallback_reference_cif(cif_path: Path) -> None:
    cif_path.write_text(
        """data_fallback_reference
_entry.id fallback_reference
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA  A 1 1 ?   0.0 0.0 0.0 1.00 10.0 1 A 1
ATOM 2 C CA . GLY  A 1 2 ?   5.0 0.0 0.0 1.00 10.0 2 A 1
ATOM 3 C CA . SER  A 1 3 ?  10.0 0.0 0.0 1.00 10.0 3 A 1
ATOM 4 C CA . ALA  B 2 1 ? 100.0 0.0 0.0 1.00 10.0 1 B 1
ATOM 5 C CA . GLY  B 2 2 ? 105.0 0.0 0.0 1.00 10.0 2 B 1
ATOM 6 C CA . SER  B 2 3 ? 110.0 0.0 0.0 1.00 10.0 3 B 1
HETATM 7 C C1 . BGC C 3 1 ? 101.5 0.0 0.0 1.00 10.0 1 C 1
#
loop_
_pdbx_branch_scheme.asym_id
_pdbx_branch_scheme.entity_id
_pdbx_branch_scheme.mon_id
_pdbx_branch_scheme.num
_pdbx_branch_scheme.pdb_asym_id
_pdbx_branch_scheme.pdb_mon_id
_pdbx_branch_scheme.pdb_seq_num
_pdbx_branch_scheme.auth_asym_id
_pdbx_branch_scheme.auth_mon_id
_pdbx_branch_scheme.auth_seq_num
_pdbx_branch_scheme.hetero
C 3 BGC 1 C BGC 1 C BGC 1 n
#
loop_
_struct_asym.id
_struct_asym.pdbx_blank_PDB_chainid_flag
_struct_asym.pdbx_modified
_struct_asym.entity_id
_struct_asym.details
A N N 1 ?
B N N 2 ?
C N N 3 ?
#
"""
    )


def _best_ranked_af3_pose_cif(
    protein_id: str,
    ligand_id: str,
    *,
    run_id: str,
) -> tuple[Path, float]:
    case_dir = WORK_CORE_ROOT / ligand_id / "af3" / "runs" / run_id / f"{protein_id}_{ligand_id}"
    best_score: float | None = None
    best_cif: Path | None = None

    for conf_path in sorted(case_dir.glob("seed-*_sample-*/*_summary_confidences.json")):
        confidence_data = json.loads(conf_path.read_text())
        ranking_score = float(confidence_data.get("ranking_score", 0.0))
        cif_path = conf_path.with_name(
            f"{conf_path.name.removesuffix('_summary_confidences.json')}_model.cif"
        )
        if not cif_path.exists():
            continue
        if best_score is None or ranking_score > best_score:
            best_score = ranking_score
            best_cif = cif_path

    if best_cif is None or best_score is None:
        raise FileNotFoundError(
            f"No ranked AF3 pose found for {protein_id} {ligand_id} under {case_dir}"
        )

    return best_cif.resolve(), best_score


def _comparison_by_pdb(
    report: crystal_anchoring.CrystalReferenceScreenReport,
    pdb_code: str,
) -> crystal_anchoring.CrystalReferencePoseComparison:
    for comparison in report.comparisons:
        if comparison.pdb_code == pdb_code:
            return comparison
    raise AssertionError(f"Missing comparison for crystal {pdb_code}: {report.comparisons}")


def test_load_crystal_reference_records_returns_both_a0a0s2gkz1_references() -> None:
    records = crystal_anchoring.load_crystal_reference_records(
        "A0A0S2GKZ1",
        reference_index_csv=REFERENCE_INDEX_CSV,
        crystal_root=CRYSTAL_ROOT,
    )

    assert [record.pdb_code for record in records] == ["5ACI", "7PXW"]
    assert {record.dp for record in records} == {4, 6}
    assert all(record.expected_ligand for record in records)
    assert all(record.source_cif.exists() for record in records)


def test_select_crystal_reference_site_prefers_chain_a_for_real_5aci() -> None:
    source_cif = CRYSTAL_ROOT / "A0A0S2GKZ1" / "5ACI_A0A0S2GKZ1_ligand.cif"
    selection = crystal_anchoring.select_crystal_reference_site(
        source_cif,
        preferred_protein_chain="A",
        expected_ligand=True,
    )

    assert selection.selected_protein_chain == "A"
    assert selection.ligand_chain_ids == ("B",)
    assert selection.preferred_chain_has_ligand is True
    assert selection.selected_chain_has_ligand is True
    assert selection.used_fallback_protein_chain is False
    assert "C" in selection.copper_chain_ids


def test_write_selected_reference_cif_keeps_real_selected_chains(tmp_path: Path) -> None:
    source_cif = CRYSTAL_ROOT / "A0A0S2GKZ1" / "5ACI_A0A0S2GKZ1_ligand.cif"
    selection = crystal_anchoring.select_crystal_reference_site(
        source_cif,
        preferred_protein_chain="A",
        expected_ligand=True,
    )

    subset_cif = crystal_anchoring._write_selected_reference_cif(
        source_cif,
        selection,
        tmp_path / "selected_reference.cif",
    )

    subset_text = subset_cif.read_text()
    subset_structure = crystal_anchoring.gemmi.read_structure(str(subset_cif))

    assert subset_cif.exists()
    assert {chain.name for chain in subset_structure[0]} == {"A", "B"}
    assert "_pdbx_branch_scheme.asym_id" in subset_text
    assert "_pdbx_nonpoly_scheme.asym_id" in subset_text
    assert " CU " in subset_text


def test_select_crystal_reference_site_resolves_real_7pxw_branch_ligand() -> None:
    source_cif = CRYSTAL_ROOT / "A0A0S2GKZ1" / "7PXW_A0A0S2GKZ1_ligand.cif"
    selection = crystal_anchoring.select_crystal_reference_site(
        source_cif,
        preferred_protein_chain="A",
        expected_ligand=True,
    )

    assert selection.selected_protein_chain == "A"
    assert selection.ligand_chain_ids == ("B",)
    assert selection.preferred_chain_has_ligand is True
    assert selection.selected_chain_has_ligand is True
    assert selection.used_fallback_protein_chain is False
    assert "E" in selection.copper_chain_ids


def test_compute_feature_aligned_tanimoto_aligns_union_feature_space() -> None:
    tanimoto = crystal_anchoring.compute_feature_aligned_tanimoto(
        ["feat_a", "feat_b"],
        [1, 1],
        ["feat_b", "feat_c"],
        [1, 1],
    )

    assert tanimoto == pytest.approx(1.0 / 3.0)


def test_crystal_ifp_contact_eligibility_uses_non_vdw_rule() -> None:
    vdw_only = crystal_anchoring.IFPResult(
        pose_id="crystal_vdw",
        status="ok",
        n_residues=1,
        n_interaction_types=1,
        residue_names=["BGC1|ASN10"],
        interaction_types=["VdWContact"],
        feature_names=["BGC1|ASN10|VdWContact"],
        fingerprint=[[1]],
        flat_bitvector=[1],
        n_total_contacts=1,
        interaction_counts={"VdWContact": 1},
    )
    eligible = crystal_anchoring.IFPResult(
        pose_id="crystal_specific",
        status="ok",
        n_residues=1,
        n_interaction_types=2,
        residue_names=["BGC1|ASN10"],
        interaction_types=["HBDonor", "HBAcceptor"],
        feature_names=["BGC1|ASN10|HBDonor", "BGC1|ASN10|HBAcceptor"],
        fingerprint=[[1, 1]],
        flat_bitvector=[1, 1],
        n_total_contacts=2,
        interaction_counts={"HBDonor": 1, "HBAcceptor": 1},
    )

    vdw_eligibility = crystal_anchoring._contact_eligibility_for_ifp(vdw_only)
    specific_eligibility = crystal_anchoring._contact_eligibility_for_ifp(eligible)

    assert vdw_eligibility.eligible is False
    assert vdw_eligibility.exclusion_class == "vdw_only"
    assert specific_eligibility.eligible is True
    assert specific_eligibility.exclusion_class is None


def test_select_crystal_reference_site_falls_back_when_chain_a_is_not_ligand_bound(tmp_path: Path) -> None:
    cif_path = tmp_path / "fallback_reference.cif"
    _write_fallback_reference_cif(cif_path)

    selection = crystal_anchoring.select_crystal_reference_site(
        cif_path,
        preferred_protein_chain="A",
        expected_ligand=True,
        ligand_distance_cutoff_a=6.0,
    )

    assert selection.selected_protein_chain == "B"
    assert selection.ligand_chain_ids == ("C",)
    assert selection.preferred_chain_has_ligand is False
    assert selection.selected_chain_has_ligand is True
    assert selection.used_fallback_protein_chain is True


def test_write_selected_reference_cif_remaps_fallback_site_to_default_chains(tmp_path: Path) -> None:
    source_cif = tmp_path / "fallback_reference.cif"
    _write_fallback_reference_cif(source_cif)
    selection = crystal_anchoring.select_crystal_reference_site(
        source_cif,
        preferred_protein_chain="A",
        expected_ligand=True,
        ligand_distance_cutoff_a=6.0,
    )

    subset_cif = crystal_anchoring._write_selected_reference_cif(
        source_cif,
        selection,
        tmp_path / "selected_reference.cif",
    )

    subset_text = subset_cif.read_text()
    subset_structure = crystal_anchoring.gemmi.read_structure(str(subset_cif))

    assert selection.selected_protein_chain == "B"
    assert {chain.name for chain in subset_structure[0]} == {"A", "B"}
    assert "C 3 BGC 1 C BGC 1 C BGC 1 n" not in subset_text
    assert "B 3 BGC 1 B BGC 1 B BGC 1 n" in subset_text


def test_write_ifp_artifacts_persists_full_ifp_result(tmp_path: Path) -> None:
    result = crystal_anchoring.IFPResult(
        pose_id="5ACI",
        status="ok",
        n_residues=2,
        n_interaction_types=2,
        residue_names=["ASP10", "TYR24"],
        interaction_types=["HBDonor", "Hydrophobic"],
        feature_names=[
            "BGC1|ASP10|HBDonor",
            "BGC1|TYR24|Hydrophobic",
        ],
        fingerprint=[[1, 0], [0, 1]],
        flat_bitvector=[1, 0],
        n_total_contacts=1,
        interaction_counts={"HBDonor": 1, "Hydrophobic": 0},
    )

    output_dir = tmp_path / "ifp"
    artifact_dir, ifp_result_json_path, pose_ifp_table_tsv_path, ifp_matrix_csv_path = crystal_anchoring._write_ifp_artifacts(result, output_dir)

    result_json = json.loads(ifp_result_json_path.read_text())
    matrix_lines = ifp_matrix_csv_path.read_text().splitlines()
    table_text = pose_ifp_table_tsv_path.read_text()

    assert (output_dir / "pose_ifp_table.tsv").exists()
    assert (output_dir / "ifp_matrix.csv").exists()
    assert artifact_dir == output_dir
    assert ifp_result_json_path == output_dir / "ifp_result.json"
    assert pose_ifp_table_tsv_path == output_dir / "pose_ifp_table.tsv"
    assert ifp_matrix_csv_path == output_dir / "ifp_matrix.csv"
    assert result_json["pose_id"] == "5ACI"
    assert result_json["active_feature_names"] == ["BGC1|ASP10|HBDonor"]
    assert set(result_json) == {
        "pose_id",
        "status",
        "error",
        "n_residues",
        "n_interaction_types",
        "residue_names",
        "n_total_contacts",
        "interaction_counts",
        "active_feature_names",
    }
    assert matrix_lines[0] == "pose_id,BGC1|ASP10|HBDonor,BGC1|TYR24|Hydrophobic"
    assert "ifp_feature_names" in table_text


def test_identify_pocket_residues_by_proximity_combines_ligand_and_cu_nearby_residues(tmp_path: Path) -> None:
    complex_pdb = tmp_path / "complex.pdb"
    _write_pocket_test_pdb(
        complex_pdb,
        ligand_point=(1.0, 0.0, 0.0),
        copper_point=(0.0, 1.6, 0.0),
    )

    pocket_residues = crystal_anchoring.identify_pocket_residues_by_proximity(
        complex_pdb,
        ligand_chain="B",
        protein_chain="A",
        cutoff_a=1.0,
    )

    assert pocket_residues == [1, 2, 3]


def test_map_residue_number_pairs_by_sequence_normalizes_histidine_and_remaps_shifted_numbers(
    tmp_path: Path,
) -> None:
    representative_pdb = tmp_path / "representative.pdb"
    apo_crystal_pdb = tmp_path / "apo_crystal.pdb"
    _write_pocket_test_pdb(
        representative_pdb,
        residue_names=("HIS", "GLY", "SER"),
        ligand_point=(1.0, 0.0, 0.0),
        copper_point=(0.0, 1.6, 0.0),
    )
    _write_pocket_test_pdb(
        apo_crystal_pdb,
        translation=(10.0, -4.0, 3.0),
        rotate_z_deg=90.0,
        residue_numbers=(101, 102, 103),
        residue_names=("HIC", "GLY", "SER"),
        ligand_point=None,
        copper_point=(0.0, 1.6, 0.0),
    )

    representative_pocket_residues = crystal_anchoring.identify_pocket_residues_by_proximity(
        representative_pdb,
        ligand_chain="B",
        protein_chain="A",
        cutoff_a=1.0,
    )
    residue_number_pairs = crystal_anchoring._map_residue_number_pairs_by_sequence(
        representative_pdb,
        apo_crystal_pdb,
        representative_pocket_residues,
        chain_name="A",
    )
    rmsd = crystal_anchoring._compute_pocket_rmsd_from_residue_pairs(
        representative_pdb,
        apo_crystal_pdb,
        residue_number_pairs,
        protein_chain="A",
    )

    assert representative_pocket_residues == [1, 2, 3]
    assert residue_number_pairs == [(1, 101), (2, 102), (3, 103)]
    assert rmsd == pytest.approx(0.0, abs=1e-6)


def test_compute_pocket_rmsd_removes_rigid_transform(tmp_path: Path) -> None:
    crystal_pdb = tmp_path / "crystal.pdb"
    pred_pdb = tmp_path / "predicted.pdb"
    _write_pocket_test_pdb(crystal_pdb)
    _write_pocket_test_pdb(pred_pdb, translation=(10.0, -4.0, 3.0), rotate_z_deg=90.0)

    rmsd = crystal_anchoring._compute_pocket_rmsd(
        pred_pdb,
        crystal_pdb,
        [1, 2, 3],
        protein_chain="A",
    )

    assert rmsd == pytest.approx(0.0, abs=1e-6)


def test_compute_pocket_rmsd_detects_local_perturbation(tmp_path: Path) -> None:
    crystal_pdb = tmp_path / "crystal.pdb"
    pred_pdb = tmp_path / "predicted.pdb"
    _write_pocket_test_pdb(crystal_pdb)
    _write_pocket_test_pdb(
        pred_pdb,
        translation=(10.0, -4.0, 3.0),
        rotate_z_deg=90.0,
        perturb_residue_2=(0.4, 0.0, 0.0),
    )

    rmsd = crystal_anchoring._compute_pocket_rmsd(
        pred_pdb,
        crystal_pdb,
        [1, 2, 3],
        protein_chain="A",
    )

    assert rmsd is not None
    assert rmsd > 0.1


@pytest.mark.slow
def test_real_pocket_rmsd_for_best_ranked_a0a0s2gkz1_cel4_pose(tmp_path: Path) -> None:
    pose_cif, ranking_score = _best_ranked_af3_pose_cif(
        "A0A0S2GKZ1",
        "CEL4",
        run_id="407999",
    )

    report = crystal_anchoring.run_crystal_reference_screen(
        pose_cif,
        protein_id="A0A0S2GKZ1",
        ligand_id="CEL4",
        representative_pose_id=pose_cif.stem,
        output_dir=tmp_path / "a0a0s2gkz1_cel4",
    )
    comparison = _comparison_by_pdb(report, "7PXW")

    assert ranking_score == pytest.approx(0.87, abs=1e-6)
    assert comparison.status == "crystal_ifp_not_contact_eligible"
    assert comparison.crystal_ifp_contact_eligible is False
    assert comparison.crystal_ifp_exclusion_class in {"vdw_only", "low_specific_contact"}
    assert comparison.ifp_tanimoto is None
    assert comparison.ifp_comparison_eligible is False
    assert comparison.pocket_residues
    assert 1 in comparison.pocket_residues
    assert comparison.pocket_rmsd is not None
    assert comparison.pocket_rmsd == pytest.approx(0.14, abs=0.05)
    assert comparison.pocket_rmsd_below_threshold is True


@pytest.mark.slow
def test_real_pocket_rmsd_is_unavailable_for_best_ranked_q7s439_cel4_pose_against_apo_crystal(
    tmp_path: Path,
) -> None:
    pose_cif, ranking_score = _best_ranked_af3_pose_cif(
        "Q7S439",
        "CEL4",
        run_id="407999",
    )

    report = crystal_anchoring.run_crystal_reference_screen(
        pose_cif,
        protein_id="Q7S439",
        ligand_id="CEL4",
        representative_pose_id=pose_cif.stem,
        output_dir=tmp_path / "q7s439_cel4",
    )
    comparison = _comparison_by_pdb(report, "5FOH")

    assert ranking_score == pytest.approx(0.9, abs=1e-6)
    assert comparison.status == "prepared_no_ligand"
    assert comparison.pocket_residues
    assert 1 in comparison.pocket_residues
    assert comparison.pocket_rmsd is not None
    assert comparison.pocket_rmsd < crystal_anchoring.POCKET_RMSD_THRESHOLD
    assert comparison.pocket_rmsd_below_threshold is True
