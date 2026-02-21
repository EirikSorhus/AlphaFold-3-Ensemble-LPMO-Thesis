"""Tests for case generation."""

import pytest

from structure_pipeline.cases import (
    Case,
    CaseStatus,
    generate_cases,
    group_cases_by_protein_model,
    group_cases_by_ligand_model,
    _generate_case_id,
)
from structure_pipeline.manifest.proteins import ProteinRecord
from structure_pipeline.manifest.ligands import LigandRecord
from structure_pipeline.manifest.msa import MSARecord


class TestCaseIdGeneration:
    """Tests for case ID generation."""

    def test_deterministic(self):
        """Test that case IDs are deterministic."""
        id1 = _generate_case_id("P12345", "L001", "af3", "native")
        id2 = _generate_case_id("P12345", "L001", "af3", "native")
        assert id1 == id2

    def test_different_inputs_different_ids(self):
        """Test that different inputs produce different IDs."""
        id1 = _generate_case_id("P12345", "L001", "af3", "native")
        id2 = _generate_case_id("P12345", "L001", "boltz", "mmseqs")
        assert id1 != id2


class TestGenerateCases:
    """Tests for case generation."""

    @pytest.fixture
    def proteins(self):
        return [
            ProteinRecord(
                protein_id="P001",
                uniprot_ids_all="P001",
                organism="Test",
                annotation="Test",
                sequence="AAAA",
                sequence_sha256="abc",
                sequence_length=4,
                fasta_header_raw=">P001",
            ),
            ProteinRecord(
                protein_id="P002",
                uniprot_ids_all="P002",
                organism="Test",
                annotation="Test",
                sequence="CCCC",
                sequence_sha256="def",
                sequence_length=4,
                fasta_header_raw=">P002",
            ),
        ]

    @pytest.fixture
    def ligands(self):
        return [
            LigandRecord(
                ligand_id="L001",
                category="amylose",
                ccd_code="AMY",
                filename="AMY.cif",
                cif_path="/path/AMY.cif",
                cif_sha256="123",
            ),
            LigandRecord(
                ligand_id="L002",
                category="cellulose",
                ccd_code="CEL",
                filename="CEL.cif",
                cif_path="/path/CEL.cif",
                cif_sha256="456",
            ),
        ]

    @pytest.fixture
    def msa_records(self):
        return [
            MSARecord(
                protein_id="P001",
                msa_path="/path/P001.a3m",
                msa_source="mmseqs",
            ),
            # P002 has no MSA
        ]

    def test_generate_all_combinations(self, proteins, ligands, msa_records):
        """Test generating all protein x ligand x model combinations."""
        cases = list(
            generate_cases(
                proteins=proteins,
                ligands=ligands,
                msa_records=msa_records,
                models=["af3", "boltz", "rf3"],
                missing_msa_proteins=["P002"],
            )
        )

        # 2 proteins x 2 ligands x 3 models = 12 cases
        assert len(cases) == 12

        # Check AF3 cases (should all be pending, uses native MSA)
        af3_cases = [c for c in cases if c.model == "af3"]
        assert len(af3_cases) == 4  # 2 proteins x 2 ligands
        assert all(c.msa_source == "native" for c in af3_cases)
        assert all(c.status == CaseStatus.PENDING.value for c in af3_cases)

        # Check Boltz cases for P001 (has MSA)
        boltz_p001 = [
            c for c in cases if c.model == "boltz" and c.protein_id == "P001"
        ]
        assert len(boltz_p001) == 2
        assert all(c.status == CaseStatus.PENDING.value for c in boltz_p001)
        assert all(c.msa_path == "/path/P001.a3m" for c in boltz_p001)

        # Check Boltz cases for P002 (no MSA - should be skipped)
        boltz_p002 = [
            c for c in cases if c.model == "boltz" and c.protein_id == "P002"
        ]
        assert len(boltz_p002) == 2
        assert all(c.status == CaseStatus.SKIPPED.value for c in boltz_p002)

    def test_output_dir_name(self, proteins, ligands, msa_records):
        """Test that output_dir_name uses protein_ccd format."""
        cases = list(
            generate_cases(
                proteins=proteins,
                ligands=ligands,
                msa_records=msa_records,
                models=["rf3"],
            )
        )
        case = [c for c in cases if c.protein_id == "P001" and c.ligand_ccd_code == "AMY"][0]
        assert case.output_dir_name == "P001_AMY"


class TestGroupCases:
    """Tests for grouping cases."""

    @pytest.fixture
    def sample_cases(self):
        return [
            Case(
                case_id="1",
                protein_id="P001",
                ligand_id="L001",
                ligand_ccd_code="AMY",
                model="af3",
                msa_source="native",
                msa_path="",
            ),
            Case(
                case_id="2",
                protein_id="P002",
                ligand_id="L001",
                ligand_ccd_code="AMY",
                model="af3",
                msa_source="native",
                msa_path="",
            ),
            Case(
                case_id="3",
                protein_id="P001",
                ligand_id="L001",
                ligand_ccd_code="AMY",
                model="boltz",
                msa_source="mmseqs",
                msa_path="/path/msa.a3m",
            ),
            Case(
                case_id="4",
                protein_id="P001",
                ligand_id="L002",
                ligand_ccd_code="CEL",
                model="af3",
                msa_source="native",
                msa_path="",
            ),
        ]

    def test_group_by_protein_model(self, sample_cases):
        """Test grouping cases by protein and model."""
        grouped = group_cases_by_protein_model(sample_cases)

        assert len(grouped) == 3
        assert ("P001", "af3") in grouped
        assert ("P002", "af3") in grouped
        assert ("P001", "boltz") in grouped

        assert len(grouped[("P001", "af3")]) == 2  # AMY + CEL
        assert len(grouped[("P002", "af3")]) == 1
        assert len(grouped[("P001", "boltz")]) == 1

    def test_group_by_ligand_model(self, sample_cases):
        """Test grouping cases by ligand and model."""
        grouped = group_cases_by_ligand_model(sample_cases)

        assert len(grouped) == 3
        assert ("AMY", "af3") in grouped
        assert ("AMY", "boltz") in grouped
        assert ("CEL", "af3") in grouped

        # AMY/af3 should have P001 + P002
        assert len(grouped[("AMY", "af3")]) == 2
        assert {c.protein_id for c in grouped[("AMY", "af3")]} == {"P001", "P002"}

        # AMY/boltz should have P001 only
        assert len(grouped[("AMY", "boltz")]) == 1

        # CEL/af3 should have P001
        assert len(grouped[("CEL", "af3")]) == 1

    def test_skip_skipped_cases(self):
        """Test that skipped cases are not grouped."""
        cases = [
            Case(
                case_id="1",
                protein_id="P001",
                ligand_id="L001",
                ligand_ccd_code="AMY",
                model="boltz",
                msa_source="mmseqs",
                msa_path="",
                status=CaseStatus.SKIPPED.value,
            ),
        ]

        grouped_protein = group_cases_by_protein_model(cases)
        assert len(grouped_protein) == 0

        grouped_ligand = group_cases_by_ligand_model(cases)
        assert len(grouped_ligand) == 0
