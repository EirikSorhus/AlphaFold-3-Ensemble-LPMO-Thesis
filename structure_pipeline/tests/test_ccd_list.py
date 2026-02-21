"""Tests for CCD-list ligand loading."""

import tempfile
from pathlib import Path

import pytest

from structure_pipeline.manifest.ligands import load_ccd_list, LigandRecord


class TestLoadCCDList:
    """Tests for load_ccd_list."""

    def test_basic_list(self):
        """Simple list of CCD codes."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("CEL6\nSTA6\nNAG6\n")
            path = Path(f.name)

        try:
            records = load_ccd_list(path)
            assert len(records) == 3
            assert records[0].ccd_code == "CEL6"
            assert records[1].ccd_code == "STA6"
            assert records[2].ccd_code == "NAG6"

            # CCD-list mode: no CIF path
            for rec in records:
                assert rec.cif_path == ""
                assert rec.filename == ""
                assert rec.cif_sha256 == ""
                assert rec.category == "ccd_list"
        finally:
            path.unlink()

    def test_skip_comments_and_blanks(self):
        """Comments and blank lines should be skipped."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("# This is a comment\nCEL6\n\n# Another comment\nCU\n\n")
            path = Path(f.name)

        try:
            records = load_ccd_list(path)
            assert len(records) == 2
            assert records[0].ccd_code == "CEL6"
            assert records[1].ccd_code == "CU"
        finally:
            path.unlink()

    def test_skip_duplicates(self):
        """Duplicate CCD codes should be deduplicated."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("CEL6\nSTA6\nCEL6\n")
            path = Path(f.name)

        try:
            records = load_ccd_list(path)
            assert len(records) == 2
            codes = [r.ccd_code for r in records]
            assert codes == ["CEL6", "STA6"]
        finally:
            path.unlink()

    def test_sequential_ids(self):
        """IDs should be sequential L001, L002, ..."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("A\nB\nC\n")
            path = Path(f.name)

        try:
            records = load_ccd_list(path)
            ids = [r.ligand_id for r in records]
            assert ids == ["L001", "L002", "L003"]
        finally:
            path.unlink()

    def test_empty_file_raises(self):
        """Empty file (no valid codes) should raise ValueError."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("# Only comments\n\n")
            path = Path(f.name)

        try:
            with pytest.raises(ValueError, match="No CCD codes"):
                load_ccd_list(path)
        finally:
            path.unlink()

    def test_missing_file_raises(self):
        """Non-existent file should raise FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            load_ccd_list(Path("/nonexistent/ccd_list.txt"))
