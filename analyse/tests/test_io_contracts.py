# tests/test_io_contracts.py
"""
Contract tests for I/O modules.
Verifies that normalization, protonation, and CCD lookup
produce outputs conforming to specified schemas.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest


# ---------------------------------------------------------------------------
# CCD lookup: valid monosaccharides accepted, custom oligomers rejected
# ---------------------------------------------------------------------------
class TestCCDLookup:
    """Verify CCD monosaccharide validation (HARD RULE)."""

    def test_valid_nag_accepted(self) -> None:
        """NAG (N-acetylglucosamine) must be recognized."""
        from lpmo_pipeline.io.ccd_lookup import validate_comp_id
        result = validate_comp_id("NAG")
        assert result.is_valid_ccd_mono is True

    def test_valid_bgc_accepted(self) -> None:
        """BGC (beta-D-glucose) must be recognized."""
        from lpmo_pipeline.io.ccd_lookup import validate_comp_id
        result = validate_comp_id("BGC")
        assert result.is_valid_ccd_mono is True

    def test_cel6_rejected(self) -> None:
        """CEL6 (custom oligomer alias) MUST be rejected.
        HARD RULE: Privateer input uses valid CCD monosaccharides only.
        """
        from lpmo_pipeline.io.ccd_lookup import validate_comp_id
        result = validate_comp_id("CEL6")
        assert result.is_valid_ccd_mono is False
        assert "oligomer" in result.rejection_reason.lower()

    def test_batch_validation_all_valid(self) -> None:
        """100% recognition for valid monosaccharide list."""
        from lpmo_pipeline.io.ccd_lookup import validate_all_glycan_residues
        all_valid, results = validate_all_glycan_residues(["NAG", "NAG", "BGC"])
        assert all_valid is True
        assert len(results) == 3

    def test_batch_validation_one_invalid_fails(self) -> None:
        """Gate: any unrecognized comp_id → all_valid=False."""
        from lpmo_pipeline.io.ccd_lookup import validate_all_glycan_residues
        all_valid, results = validate_all_glycan_residues(["NAG", "CEL6"])
        assert all_valid is False


# ---------------------------------------------------------------------------
# Normalized CIF contract
# ---------------------------------------------------------------------------
class TestNormalizedCIFContract:
    """Verify normalized mmCIF structure meets the contract.

    Contract (from schemas/normalized_cif_contract.txt):
      - Chain A = protein
      - Chains B..D = glycans
      - Chain E = metal (Cu)
      - _chem_comp_bond complete for all comp_id
      - _struct_conn for glycosidic + Cu-coord
      - _atom_site.B_iso_or_equiv present for all atoms
    """

    @pytest.fixture
    def sample_cif(self, tmp_path: Path) -> Path:
        """Create a minimal test CIF (placeholder)."""
        cif_path = tmp_path / "test_normalized.cif"
        cif_path.write_text("# Placeholder test CIF\n")
        return cif_path

    def test_chain_ids_follow_schema(self, sample_cif: Path) -> None:
        """Protein=A, glycans=B/C/D, metal=E."""
        # PSEUDOCODE: parse with gemmi, check chain IDs
        # import gemmi
        # doc = gemmi.cif.read(str(sample_cif))
        # structure = gemmi.make_structure_from_block(doc[0])
        # chain_names = [chain.name for chain in structure[0]]
        # assert "A" in chain_names  # protein
        pass  # Placeholder

    def test_confidence_field_present(self, sample_cif: Path) -> None:
        """_atom_site.B_iso_or_equiv must be present for all atoms."""
        # PSEUDOCODE: parse CIF, check B_iso column
        pass  # Placeholder


# ---------------------------------------------------------------------------
# Protonation export contract
# ---------------------------------------------------------------------------
class TestProtonationContract:
    """Verify protonation outputs have required features."""

    def test_cif_to_pdb_writes_posebusters_pdb(self, tmp_path: Path) -> None:
        """Step 7b writes a PDB artifact from normalized mmCIF."""
        from lpmo_pipeline.io.cif_to_pdb import convert_cif_to_pdb

        fixture = Path(__file__).parent / "fixtures" / "test_normalized.cif"
        success, output_path = convert_cif_to_pdb(fixture, tmp_path)

        assert success is True
        assert output_path is not None
        assert output_path.name == "for_posebusters.pdb"
        pdb_text = output_path.read_text()
        assert "ATOM" in pdb_text or "HETATM" in pdb_text

        report_path = tmp_path / "cif_to_pdb_report.json"
        assert report_path.exists()

    def test_cif_to_pdb_report_has_required_fields(self, tmp_path: Path) -> None:
        """cif_to_pdb_report.json must record backend and fallback reason."""
        import json
        from lpmo_pipeline.io.cif_to_pdb import convert_cif_to_pdb

        fixture = Path(__file__).parent / "fixtures" / "test_normalized.cif"
        success, _ = convert_cif_to_pdb(fixture, tmp_path)
        assert success is True

        report = json.loads((tmp_path / "cif_to_pdb_report.json").read_text())
        assert "backend" in report
        assert report["backend"] in ("pdbfixer", "gemmi", "gemmi_atom_site_fallback")
        assert "backend_fallback_reason" in report
        assert "atom_count" in report
        assert report["atom_count"] > 0

    def test_cif_to_pdb_preferred_backend_is_pdbfixer(self, tmp_path: Path) -> None:
        """When PDBFixer is available the reported backend must be 'pdbfixer'."""
        import json
        pytest.importorskip("pdbfixer")
        pytest.importorskip("openmm")

        from lpmo_pipeline.io.cif_to_pdb import convert_cif_to_pdb

        fixture = Path(__file__).parent / "fixtures" / "test_normalized.cif"
        success, _ = convert_cif_to_pdb(fixture, tmp_path)
        assert success is True

        report = json.loads((tmp_path / "cif_to_pdb_report.json").read_text())
        assert report["backend"] == "pdbfixer", (
            f"Expected backend=pdbfixer, got {report['backend']}. "
            f"fallback_reason={report.get('backend_fallback_reason')}"
        )

    def test_protonation_report_schema(self, tmp_path: Path) -> None:
        """protonation_report.json must have all required schema fields."""
        import json
        from lpmo_pipeline.io.protonate_export import protonate_and_export
        from lpmo_pipeline.io.normalize_mmcif import NormalizeMMCIFRunner

        fixture = Path(__file__).parent / "fixtures" / "test_normalized.cif"
        # protonate_and_export accepts normalized CIF directly
        ok, report = protonate_and_export(fixture, tmp_path)

        report_path = tmp_path / "protonation_report.json"
        assert report_path.exists(), "protonation_report.json must be written"

        data = json.loads(report_path.read_text())
        required_fields = {
            "input_cif", "output_dir", "for_posebusters_pdb",
            "complex_h_pdb", "ligand_mol2", "complex_h_backend",
            "used_obabel", "blockers", "warnings",
        }
        missing = required_fields - data.keys()
        assert not missing, f"protonation_report.json missing fields: {missing}"

    @pytest.mark.slow
    def test_complex_h_has_hydrogens(self, tmp_path: Path) -> None:
        """complex_H.pdb must have more atom lines than for_posebusters.pdb."""
        pytest.importorskip("pdbfixer")
        pytest.importorskip("openmm")

        import json
        from lpmo_pipeline.io.protonate_export import protonate_and_export, _count_atom_lines

        fixture = Path(__file__).parent / "fixtures" / "test_normalized.cif"
        ok, report = protonate_and_export(fixture, tmp_path)
        assert report is not None

        data = json.loads((tmp_path / "protonation_report.json").read_text())
        assert data["complex_h_backend"] != "none", (
            "complex_H.pdb backend must not be 'none'; "
            f"warnings={data.get('warnings')}"
        )

        posebusters_pdb = Path(data["for_posebusters_pdb"])
        complex_h_pdb = Path(data["complex_h_pdb"])
        assert posebusters_pdb.exists()
        assert complex_h_pdb.exists()

        n_base = _count_atom_lines(posebusters_pdb)
        n_h = _count_atom_lines(complex_h_pdb)
        assert n_h > n_base, (
            f"complex_H.pdb must have more atoms than for_posebusters.pdb "
            f"(base={n_base}, complex_h={n_h})"
        )

    @pytest.mark.slow
    def test_mol2_has_tripos_sections(self, tmp_path: Path) -> None:
        """ligand_for_prolif.mol2 must have @<TRIPOS>ATOM and @<TRIPOS>BOND sections."""
        import json
        from lpmo_pipeline.io.protonate_export import protonate_and_export

        fixture = Path(__file__).parent / "fixtures" / "test_normalized.cif"
        ok, report = protonate_and_export(fixture, tmp_path)
        assert report is not None

        data = json.loads((tmp_path / "protonation_report.json").read_text())
        mol2_path = Path(data["ligand_mol2"])
        if not mol2_path.exists():
            pytest.skip("MOL2 not produced (obabel unavailable)")

        mol2_text = mol2_path.read_text(errors="ignore")
        assert "@<TRIPOS>ATOM" in mol2_text, "MOL2 missing @<TRIPOS>ATOM section"
        assert "@<TRIPOS>BOND" in mol2_text, "MOL2 missing @<TRIPOS>BOND section"

    @pytest.mark.slow
    def test_mol2_is_ligand_only(self, tmp_path: Path) -> None:
        """ligand_for_prolif.mol2 must not contain protein residue labels."""
        import json
        from lpmo_pipeline.io.protonate_export import protonate_and_export, _PROTEIN_RESIDUE_NAMES

        fixture = Path(__file__).parent / "fixtures" / "test_normalized.cif"
        ok, report = protonate_and_export(fixture, tmp_path)
        assert report is not None

        data = json.loads((tmp_path / "protonation_report.json").read_text())
        mol2_path = Path(data["ligand_mol2"])
        if not mol2_path.exists():
            pytest.skip("MOL2 not produced (obabel unavailable)")

        mol2_text = mol2_path.read_text(errors="ignore")
        protein_hits = [r for r in _PROTEIN_RESIDUE_NAMES if r in mol2_text]
        assert not protein_hits, (
            f"ligand_for_prolif.mol2 contains protein residue labels: {protein_hits}"
        )
