"""Tests for FASTA parsing and protein manifest generation."""

import tempfile
from pathlib import Path

import pytest

from structure_pipeline.manifest.proteins import (
    ProteinRecord,
    _parse_header,
    parse_fasta,
    write_proteins_manifest,
    load_proteins_manifest,
)


class TestParseHeader:
    """Tests for FASTA header parsing."""

    def test_uniprot_ids_format(self):
        """Test parsing UniProtIDs_ID1;ID2_Organism_Annotation format."""
        header = "UniProtIDs_Q1K4Q1;Q873G1_Neurospora_crassa_AA9_family_LPMO"
        result = _parse_header(header)

        assert result["protein_id"] == "Q1K4Q1"
        assert result["uniprot_ids_all"] == "Q1K4Q1;Q873G1"
        assert "Neurospora" in result["organism"]
        assert "AA9" in result["annotation"] or "LPMO" in result["annotation"]

    def test_single_uniprot_id(self):
        """Test parsing with single UniProt ID."""
        header = "UniProtIDs_A0A223GEC9_Organism_name_Some_annotation"
        result = _parse_header(header)

        assert result["protein_id"] == "A0A223GEC9"
        assert result["uniprot_ids_all"] == "A0A223GEC9"

    def test_standard_uniprot_header(self):
        """Test parsing standard UniProt header format."""
        header = "sp|P12345|PROT_HUMAN Some protein description"
        result = _parse_header(header)

        assert result["protein_id"] == "P12345"

    def test_simple_header(self):
        """Test fallback for simple headers."""
        header = "protein1 description text"
        result = _parse_header(header)

        assert result["protein_id"] == "protein1"


class TestParseFasta:
    """Tests for FASTA file parsing."""

    def test_parse_single_sequence(self):
        """Test parsing FASTA with single sequence."""
        fasta_content = """>UniProtIDs_Q1K4Q1_Organism_Annotation
MVTPETRGRLKGFAAPPPARSGRGQ
LLLLLLL
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(fasta_content)
            fasta_path = Path(f.name)

        try:
            proteins = list(parse_fasta(fasta_path))
            assert len(proteins) == 1

            protein = proteins[0]
            assert protein.protein_id == "Q1K4Q1"
            assert protein.sequence == "MVTPETRGRLKGFAAPPPARSGRGQLLLLLLL"
            assert protein.sequence_length == 32
            assert len(protein.sequence_sha256) == 64
        finally:
            fasta_path.unlink()

    def test_parse_multiple_sequences(self):
        """Test parsing FASTA with multiple sequences."""
        fasta_content = """>UniProtIDs_Q1K4Q1_Org1_Ann1
AAAAAAAAAA
>UniProtIDs_P12345_Org2_Ann2
CCCCCCCCCC
>UniProtIDs_A0A123_Org3_Ann3
GGGGGGGGGG
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(fasta_content)
            fasta_path = Path(f.name)

        try:
            proteins = list(parse_fasta(fasta_path))
            assert len(proteins) == 3

            ids = [p.protein_id for p in proteins]
            assert ids == ["Q1K4Q1", "P12345", "A0A123"]
        finally:
            fasta_path.unlink()


class TestManifestIO:
    """Tests for manifest reading/writing."""

    def test_roundtrip(self):
        """Test writing and reading manifest."""
        proteins = [
            ProteinRecord(
                protein_id="Q1K4Q1",
                uniprot_ids_all="Q1K4Q1;Q873G1",
                organism="Test organism",
                annotation="Test annotation",
                sequence="AAAAAAAAAA",
                sequence_sha256="abc123",
                sequence_length=10,
                fasta_header_raw="original header",
            ),
            ProteinRecord(
                protein_id="P12345",
                uniprot_ids_all="P12345",
                organism="Another organism",
                annotation="Another annotation",
                sequence="CCCCCCCCCC",
                sequence_sha256="def456",
                sequence_length=10,
                fasta_header_raw="another header",
            ),
        ]

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False
        ) as f:
            manifest_path = Path(f.name)

        try:
            write_proteins_manifest(proteins, manifest_path)
            loaded = load_proteins_manifest(manifest_path)

            assert len(loaded) == 2
            assert loaded[0].protein_id == "Q1K4Q1"
            assert loaded[1].protein_id == "P12345"
            assert loaded[0].sequence == "AAAAAAAAAA"
        finally:
            manifest_path.unlink()
