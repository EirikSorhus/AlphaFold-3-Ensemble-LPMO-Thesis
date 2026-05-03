# tests/test_normalize.py
"""
Tests for lpmo_pipeline.io.normalize_mmcif — NormalizeMMCIFRunner.

Covers:
  - Chain remapping (AF3 A/B/C → canonical A/E/B)
  - Atom mapping coverage (100% for AF3 identity)
  - Struct_conn detection
  - Confidence derivation from B_iso
  - Normalized CIF output file
  - Normalize report artifact
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from lpmo_pipeline.io.gemmi_compat import gemmi
from lpmo_pipeline.io.normalize_mmcif import NormalizeMMCIFRunner


# ---------------------------------------------------------------------------
# Fixture: AF3-format mmCIF with entity/struct_asym/struct_conn
# ---------------------------------------------------------------------------
AF3_NORMALIZE_CIF = """\
data_TEST_CEL6
#
_entry.id TEST_CEL6
#
loop_
_atom_type.symbol
C
CU
N
O
#
loop_
_chem_comp.formula
_chem_comp.formula_weight
_chem_comp.id
_chem_comp.mon_nstd_flag
_chem_comp.name
_chem_comp.pdbx_smiles
_chem_comp.pdbx_synonyms
_chem_comp.type
"C3 H7 N O2" 89.093 ALA y ALANINE C[C@@H](C(=O)O)N ? "L-PEPTIDE LINKING"
"C6 H10 N3 O2" 156.162 HIS y HISTIDINE "c1c([nH+]c[nH]1)C[C@@H](C(=O)O)N" ? "L-PEPTIDE LINKING"
"C6 H12 O6" 180.156 BGC y beta-D-glucopyranose "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O" ? "D-saccharide, beta linking"
Cu 63.546 CU . "COPPER (II) ION" "[Cu+2]" ? NON-POLYMER
#
loop_
_entity.id
_entity.pdbx_description
_entity.type
1 . polymer
2 . non-polymer
3 . branched
#
_entity_poly.entity_id      1
_entity_poly.pdbx_strand_id A
_entity_poly.type           polypeptide(L)
#
loop_
_struct_asym.entity_id
_struct_asym.id
1 A
2 B
3 C
#
loop_
_struct_conn.conn_type_id
_struct_conn.id
_struct_conn.pdbx_ptnr1_PDB_ins_code
_struct_conn.pdbx_ptnr1_label_alt_id
_struct_conn.pdbx_ptnr2_PDB_ins_code
_struct_conn.pdbx_ptnr2_label_alt_id
_struct_conn.pdbx_role
_struct_conn.pdbx_value_order
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_symmetry
covale covale1 ? ? ? ? ? ? C 1 C C1 BGC . 1_555 C 2 C O4 BGC . 1_555
#
loop_
_pdbx_branch_scheme.asym_id
_pdbx_branch_scheme.auth_asym_id
_pdbx_branch_scheme.auth_seq_num
_pdbx_branch_scheme.entity_id
_pdbx_branch_scheme.hetero
_pdbx_branch_scheme.mon_id
_pdbx_branch_scheme.num
_pdbx_branch_scheme.pdb_asym_id
_pdbx_branch_scheme.pdb_ins_code
_pdbx_branch_scheme.pdb_seq_num
C C 1 3 n BGC 1 C . 1
C C 2 3 n BGC 2 C . 2
#
_pdbx_nonpoly_scheme.asym_id B
_pdbx_nonpoly_scheme.auth_seq_num 1
_pdbx_nonpoly_scheme.entity_id 2
_pdbx_nonpoly_scheme.mon_id CU
_pdbx_nonpoly_scheme.pdb_ins_code .
_pdbx_nonpoly_scheme.pdb_seq_num 1
_pdbx_nonpoly_scheme.pdb_strand_id B
#
_struct_conn_type.criteria  ?
_struct_conn_type.id        covale
_struct_conn_type.reference ?
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
ATOM   1  N  N   . ALA A 1 1  ? 1.000  2.000  3.000  1.00 85.00 1   A 1
ATOM   2  C  CA  . ALA A 1 1  ? 2.000  3.000  4.000  1.00 86.00 1   A 1
ATOM   3  C  C   . ALA A 1 1  ? 3.000  4.000  5.000  1.00 84.00 1   A 1
ATOM   4  O  O   . ALA A 1 1  ? 4.000  5.000  6.000  1.00 83.00 1   A 1
ATOM   5  C  CB  . ALA A 1 1  ? 2.500  2.500  4.500  1.00 80.00 1   A 1
ATOM   6  N  N   . HIS A 1 2  ? 3.500  3.500  4.500  1.00 90.00 2   A 1
ATOM   7  C  CA  . HIS A 1 2  ? 4.000  4.000  5.000  1.00 91.00 2   A 1
ATOM   8  C  C   . HIS A 1 2  ? 5.000  5.000  6.000  1.00 89.00 2   A 1
ATOM   9  O  O   . HIS A 1 2  ? 6.000  6.000  7.000  1.00 88.00 2   A 1
ATOM   10 C  CB  . HIS A 1 2  ? 4.500  3.500  5.500  1.00 87.00 2   A 1
HETATM 11 CU CU  . CU  B 2 .  ? 0.000  0.000  0.000  1.00 97.00 1   B 1
HETATM 12 C  C1  . BGC C 3 .  ? 10.000 10.000 10.000 1.00 65.00 1   C 1
HETATM 13 C  C2  . BGC C 3 .  ? 11.000 11.000 11.000 1.00 63.00 1   C 1
HETATM 14 O  O5  . BGC C 3 .  ? 10.500 10.500 9.500  1.00 64.00 1   C 1
HETATM 15 C  C1  . BGC C 3 .  ? 12.000 12.000 12.000 1.00 62.00 2   C 1
HETATM 16 C  C2  . BGC C 3 .  ? 13.000 13.000 13.000 1.00 60.00 2   C 1
HETATM 17 O  O5  . BGC C 3 .  ? 12.500 12.500 11.500 1.00 61.00 2   C 1
#
"""

AF3_NO_CONN_CIF = """\
data_TEST_NOCONN
#
_entry.id TEST_NOCONN
#
loop_
_chem_comp.id
_chem_comp.type
ALA "L-PEPTIDE LINKING"
CU  NON-POLYMER
#
loop_
_entity.id
_entity.pdbx_description
_entity.type
1 . polymer
2 . non-polymer
#
loop_
_struct_asym.entity_id
_struct_asym.id
1 A
2 B
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
ATOM   1  N  N   . ALA A 1 1  ? 1.0 2.0 3.0 1.00 85.0 1 A 1
ATOM   2  C  CA  . ALA A 1 1  ? 2.0 3.0 4.0 1.00 86.0 1 A 1
HETATM 3  CU CU  . CU  B 2 .  ? 0.0 0.0 0.0 1.00 97.0 1 B 1
#
"""

AF3_INVALID_GLYCAN_CIF = AF3_NORMALIZE_CIF.replace("BGC", "CEL6")


@pytest.fixture
def af3_cif(tmp_path: Path) -> Path:
    p = tmp_path / "test_model.cif"
    p.write_text(AF3_NORMALIZE_CIF)
    return p


@pytest.fixture
def no_conn_cif(tmp_path: Path) -> Path:
    p = tmp_path / "no_conn.cif"
    p.write_text(AF3_NO_CONN_CIF)
    return p


@pytest.fixture
def invalid_glycan_cif(tmp_path: Path) -> Path:
    p = tmp_path / "invalid_glycan.cif"
    p.write_text(AF3_INVALID_GLYCAN_CIF)
    return p


# ---------------------------------------------------------------------------
# Chain remapping
# ---------------------------------------------------------------------------
class TestChainRemapping:
    def test_normalizes_successfully(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = NormalizeMMCIFRunner(af3_cif, out)
        ok, path = runner.run()
        assert ok is True
        assert path is not None
        assert path.exists()

    def test_chain_mapping_in_report(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = NormalizeMMCIFRunner(af3_cif, out)
        runner.run()
        report = json.loads((out / "normalize_report.json").read_text())
        mappings = report["chain_mappings"]
        # A→A (protein), B(Cu)→E, C(glycan)→B
        orig_to_norm = {m["original"]: m["normalized"] for m in mappings}
        assert orig_to_norm["A"] == "A"  # protein stays
        assert orig_to_norm["B"] == "E"  # Cu → E
        assert orig_to_norm["C"] == "B"  # glycan → B

    def test_normalized_cif_exists(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = NormalizeMMCIFRunner(af3_cif, out)
        runner.run()
        assert (out / "normalized.cif").exists()
        content = (out / "normalized.cif").read_text()
        assert len(content) > 100

    def test_written_normalized_cif_uses_canonical_chain_ids(
        self, af3_cif: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "out"
        runner = NormalizeMMCIFRunner(af3_cif, out)
        ok, path = runner.run()

        assert ok is True
        assert path is not None

        block = gemmi.cif.read(str(path)).sole_block()

        atom_label_asym_ids = set(block.find_values("_atom_site.label_asym_id"))
        atom_auth_asym_ids = set(block.find_values("_atom_site.auth_asym_id"))

        assert atom_label_asym_ids == {"A", "B", "E"}
        assert atom_auth_asym_ids == {"A", "B", "E"}
        assert list(block.find_values("_pdbx_nonpoly_scheme.asym_id")) == ["E"]
        assert list(block.find_values("_pdbx_nonpoly_scheme.pdb_strand_id")) == ["E"]
        assert set(block.find_values("_pdbx_branch_scheme.asym_id")) == {"B"}
        assert set(block.find_values("_pdbx_branch_scheme.auth_asym_id")) == {"B"}
        assert set(block.find_values("_pdbx_branch_scheme.pdb_asym_id")) == {"B"}


# ---------------------------------------------------------------------------
# Atom mapping
# ---------------------------------------------------------------------------
class TestAtomMapping:
    def test_full_coverage(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = NormalizeMMCIFRunner(af3_cif, out)
        runner.run()
        report = json.loads((out / "normalize_report.json").read_text())
        assert report["atom_mapping"]["coverage"] == 1.0

    def test_atom_count_correct(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = NormalizeMMCIFRunner(af3_cif, out)
        runner.run()
        report = json.loads((out / "normalize_report.json").read_text())
        # 10 protein + 1 Cu + 6 BGC = 17
        assert report["atom_mapping"]["total_atoms"] == 17

    def test_atom_map_tsv_written(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = NormalizeMMCIFRunner(af3_cif, out)
        runner.run()
        tsv = out / "atom_map.tsv"
        assert tsv.exists()
        lines = tsv.read_text().strip().split("\n")
        assert lines[0].startswith("old_atom_name")  # header
        assert len(lines) == 18  # header + 17 atoms


# ---------------------------------------------------------------------------
# Connectivity
# ---------------------------------------------------------------------------
class TestConnectivity:
    def test_struct_conn_detected(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = NormalizeMMCIFRunner(af3_cif, out)
        runner.run()
        report = json.loads((out / "normalize_report.json").read_text())
        assert report["has_struct_conn"] is True

    def test_no_struct_conn_flagged(self, no_conn_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = NormalizeMMCIFRunner(no_conn_cif, out)
        runner.run()
        report = json.loads((out / "normalize_report.json").read_text())
        assert report["has_struct_conn"] is False


# ---------------------------------------------------------------------------
# Confidence
# ---------------------------------------------------------------------------
class TestConfidence:
    def test_confidence_derived(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = NormalizeMMCIFRunner(af3_cif, out)
        runner.run()
        report = json.loads((out / "normalize_report.json").read_text())
        assert report["confidence_derived"] is True
        assert report["n_residues_with_confidence"] > 0


# ---------------------------------------------------------------------------
# CCD validation wiring
# ---------------------------------------------------------------------------
class TestCCDValidation:
    def test_ccd_validation_report_written(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = NormalizeMMCIFRunner(af3_cif, out)
        ok, _ = runner.run()
        report = json.loads((out / "normalize_report.json").read_text())

        assert ok is True
        assert report["ccd_validation"]["all_valid"] is True
        assert report["ccd_validation"]["observed_comp_ids"] == ["BGC"]
        assert report["ccd_validation"]["invalid_comp_ids"] == []

    def test_invalid_glycan_fails_normalization(
        self, invalid_glycan_cif: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "out"
        runner = NormalizeMMCIFRunner(invalid_glycan_cif, out)
        ok, path = runner.run()
        report = json.loads((out / "normalize_report.json").read_text())
        failures = json.loads((out / "normalize_failures.json").read_text())

        assert ok is False
        assert path is None
        assert report["ccd_validation"]["all_valid"] is False
        assert report["ccd_validation"]["invalid_comp_ids"] == ["CEL6"]
        assert failures[-1]["reason"] == "glykan_not_ccd"
        assert not (out / "normalized.cif").exists()


# ---------------------------------------------------------------------------
# Artifacts
# ---------------------------------------------------------------------------
class TestNormalizeArtifacts:
    def test_all_artifacts_written(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = NormalizeMMCIFRunner(af3_cif, out)
        runner.run()
        assert (out / "normalized.cif").exists()
        assert (out / "atom_map.tsv").exists()
        assert (out / "rename_log.json").exists()
        assert (out / "normalize_report.json").exists()
        assert (out / "normalize_failures.json").exists()

    def test_rename_log_json(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = NormalizeMMCIFRunner(af3_cif, out)
        runner.run()
        log = json.loads((out / "rename_log.json").read_text())
        assert log["mapping_coverage"] == 1.0
        assert log["reason"] == "af3_identity_mapping"
