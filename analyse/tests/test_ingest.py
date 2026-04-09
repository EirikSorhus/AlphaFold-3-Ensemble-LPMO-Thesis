# tests/test_ingest.py
"""
Tests for lpmo_pipeline.io.mmcif_ingest — IngestQCRunner.

Covers:
  - Parsing a minimal AF3-format mmCIF
  - Required category checks
  - Atom count validation
  - Entity type validation
  - Chain summary output
  - Missing file handling
  - Invalid CIF handling
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from lpmo_pipeline.io.mmcif_ingest import IngestQCRunner


# ---------------------------------------------------------------------------
# Fixture: minimal AF3-format mmCIF
# ---------------------------------------------------------------------------
AF3_MINI_CIF = """\
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
S
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
"C6 H12 O6" 180.156 BGC y beta-D-glucopyranose "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O" ? "D-saccharide, beta linking"
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
#
"""


@pytest.fixture
def af3_cif(tmp_path: Path) -> Path:
    """Write the minimal AF3-format CIF to a temp file."""
    p = tmp_path / "test_model.cif"
    p.write_text(AF3_MINI_CIF)
    return p


@pytest.fixture
def empty_cif(tmp_path: Path) -> Path:
    """Write a CIF with no atoms."""
    p = tmp_path / "empty.cif"
    p.write_text(
        "data_empty\n_entry.id empty\n"
        "loop_\n_entity.id\n_entity.type\n1 polymer\n#\n"
        "loop_\n_chem_comp.id\nALA\n#\n"
    )
    return p


# ---------------------------------------------------------------------------
# Parse + QC checks
# ---------------------------------------------------------------------------
class TestIngestQCParse:
    def test_af3_cif_passes(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = IngestQCRunner(af3_cif, out)
        passed, path = runner.run()
        assert passed is True
        assert path == str(af3_cif)

    def test_missing_file_fails(self, tmp_path: Path) -> None:
        runner = IngestQCRunner(tmp_path / "missing.cif", tmp_path / "out")
        passed, path = runner.run()
        assert passed is False
        assert path is None

    def test_invalid_cif_fails(self, tmp_path: Path) -> None:
        bad = tmp_path / "bad.cif"
        bad.write_text("this is not a valid CIF file\n{{{garbage}}}\n")
        runner = IngestQCRunner(bad, tmp_path / "out")
        passed, path = runner.run()
        assert passed is False


class TestIngestQCCategories:
    def test_required_categories_present(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = IngestQCRunner(af3_cif, out)
        runner.run()
        cats = [f for f in runner.qc_flags if f.check_name == "required_categories"]
        assert len(cats) == 1
        assert cats[0].status.value == "pass"

    def test_entity_types_present(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = IngestQCRunner(af3_cif, out)
        runner.run()
        ent = [f for f in runner.qc_flags if f.check_name == "entity_types"]
        assert len(ent) == 1
        assert ent[0].status.value == "pass"


class TestIngestQCAtomCount:
    def test_atom_count_positive(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = IngestQCRunner(af3_cif, out)
        runner.run()
        ac = [f for f in runner.qc_flags if f.check_name == "atom_count_positive"]
        assert len(ac) == 1
        assert "14" in ac[0].message  # 10 protein + 1 Cu + 3 BGC


class TestIngestQCChainSummary:
    def test_chain_summary_contains_abc(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = IngestQCRunner(af3_cif, out)
        runner.run()
        report_path = out / "ingest_qc_report.json"
        report = json.loads(report_path.read_text())
        cs = report["chain_summary"]
        assert "A" in cs  # protein
        assert "B" in cs  # Cu
        assert "C" in cs  # glycan

    def test_chain_a_is_protein(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = IngestQCRunner(af3_cif, out)
        runner.run()
        report = json.loads((out / "ingest_qc_report.json").read_text())
        chain_a = report["chain_summary"]["A"]
        assert chain_a["n_residues"] == 2  # ALA + HIS
        assert "ALA" in chain_a["comp_ids"]
        assert "HIS" in chain_a["comp_ids"]

    def test_chain_b_is_cu(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = IngestQCRunner(af3_cif, out)
        runner.run()
        report = json.loads((out / "ingest_qc_report.json").read_text())
        chain_b = report["chain_summary"]["B"]
        assert chain_b["n_atoms"] == 1
        assert "CU" in chain_b["comp_ids"]


class TestIngestQCArtifacts:
    def test_qc_report_written(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = IngestQCRunner(af3_cif, out)
        runner.run()
        assert (out / "ingest_qc_report.json").exists()

    def test_log_written(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = IngestQCRunner(af3_cif, out)
        runner.run()
        assert (out / "ingest_qc.log").exists()

    def test_failure_log_written(self, af3_cif: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        runner = IngestQCRunner(af3_cif, out)
        runner.run()
        assert (out / "ingest_failures.json").exists()
