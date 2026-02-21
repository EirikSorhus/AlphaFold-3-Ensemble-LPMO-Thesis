"""Tests for ligand scanning and manifest generation."""

import tempfile
from pathlib import Path

import pytest

from structure_pipeline.manifest.ligands import (
    LigandRecord,
    scan_ligands,
    write_ligands_manifest,
    load_ligands_manifest,
    _extract_ccd_code,
)


class TestExtractCCDCode:
    """Tests for CCD code extraction."""

    def test_extract_from_filename(self):
        """Test extracting CCD code from filename."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".cif", delete=False
        ) as f:
            f.write("# Empty CIF\n")
            cif_path = Path(f.name)

        try:
            # Filename is something like tmpXXXXXX.cif
            code = _extract_ccd_code(cif_path)
            assert code == cif_path.stem.upper()
        finally:
            cif_path.unlink()

    def test_extract_from_cif_content(self):
        """Test extracting CCD code from CIF content."""
        cif_content = """data_GLC
_chem_comp.id                                    GLC
_chem_comp.name                                  "alpha-D-glucose"
"""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".cif", delete=False
        ) as f:
            f.write(cif_content)
            cif_path = Path(f.name)

        try:
            code = _extract_ccd_code(cif_path)
            assert code == "GLC"
        finally:
            cif_path.unlink()


class TestScanLigands:
    """Tests for ligand directory scanning."""

    def test_scan_with_categories(self):
        """Test scanning ligands with category subdirectories."""
        with tempfile.TemporaryDirectory() as tmpdir:
            ligands_dir = Path(tmpdir)

            # Create category directories
            (ligands_dir / "amylose").mkdir()
            (ligands_dir / "cellulose").mkdir()

            # Create CIF files
            (ligands_dir / "amylose" / "AMY1.cif").write_text("data_AMY1\n")
            (ligands_dir / "amylose" / "AMY2.cif").write_text("data_AMY2\n")
            (ligands_dir / "cellulose" / "CEL1.cif").write_text("data_CEL1\n")

            ligands = list(scan_ligands(ligands_dir))

            assert len(ligands) == 3

            # Check categories
            categories = {lig.category for lig in ligands}
            assert categories == {"amylose", "cellulose"}

            # Check sequential IDs
            ids = [lig.ligand_id for lig in ligands]
            assert ids == ["L001", "L002", "L003"]

    def test_scan_empty_directory(self):
        """Test scanning empty directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            ligands = list(scan_ligands(Path(tmpdir)))
            assert len(ligands) == 0

    def test_scan_nonexistent_directory(self):
        """Test scanning non-existent directory raises error."""
        with pytest.raises(FileNotFoundError):
            list(scan_ligands(Path("/nonexistent/path")))


class TestManifestIO:
    """Tests for ligand manifest I/O."""

    def test_roundtrip(self):
        """Test writing and reading ligand manifest."""
        ligands = [
            LigandRecord(
                ligand_id="L001",
                category="amylose",
                ccd_code="AMY",
                filename="AMY.cif",
                cif_path="/path/to/AMY.cif",
                cif_sha256="abc123",
            ),
            LigandRecord(
                ligand_id="L002",
                category="cellulose",
                ccd_code="CEL",
                filename="CEL.cif",
                cif_path="/path/to/CEL.cif",
                cif_sha256="def456",
            ),
        ]

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False
        ) as f:
            manifest_path = Path(f.name)

        try:
            write_ligands_manifest(ligands, manifest_path)
            loaded = load_ligands_manifest(manifest_path)

            assert len(loaded) == 2
            assert loaded[0].ligand_id == "L001"
            assert loaded[0].category == "amylose"
            assert loaded[1].ccd_code == "CEL"
        finally:
            manifest_path.unlink()
