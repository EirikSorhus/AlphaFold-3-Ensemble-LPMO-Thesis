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

    def test_ligand_only_export_skips_cu_residue_in_glycan_chain(self, tmp_path: Path) -> None:
        """Cu in a copied glycan chain must not enter the ProLIF ligand PDB."""
        from lpmo_pipeline.io.protonate_export import _write_ligand_only_pdb

        cif_path = tmp_path / "glycan_with_cu.cif"
        cif_path.write_text(
            """data_glycan_with_cu
_entry.id glycan_with_cu
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
HETATM 1 C C1 . BGC B 1 1 ? 0.0 0.0 0.0 1.00 10.0 1 B 1
HETATM 2 O O5 . BGC B 1 1 ? 1.0 0.0 0.0 1.00 10.0 1 B 1
HETATM 3 CU CU . CU B 2 . ? 5.0 0.0 0.0 1.00 10.0 301 B 1
#
"""
        )
        warnings: list[str] = []
        ligand_pdb = tmp_path / "ligand_only_for_prolif.pdb"

        residue_count = _write_ligand_only_pdb(cif_path, ligand_pdb, warnings=warnings)
        ligand_text = ligand_pdb.read_text()

        assert residue_count == 1
        assert "BGC" in ligand_text
        assert " CU " not in ligand_text
        assert any("ligand_only_skipped_non_ligand_residues" in warning for warning in warnings)

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
        assert not any(
            line.startswith(("ATOM  ", "HETATM")) and line[76:78].strip().upper() == "CU"
            for line in pdb_text.splitlines()
        )

        report_path = tmp_path / "cif_to_pdb_report.json"
        assert report_path.exists()

    def test_posebusters_pdb_strips_metals_but_keeps_nonmetal_conect(self, tmp_path: Path) -> None:
        """Metal cleanup must not discard glycan connectivity between non-metal atoms."""
        from lpmo_pipeline.io.cif_to_pdb import _write_posebusters_ready_pdb

        input_pdb = tmp_path / "input_with_metal.pdb"
        output_pdb = tmp_path / "for_posebusters.pdb"
        input_pdb.write_text(
            "\n".join(
                [
                    "HETATM    1 CU    CU E   1       0.255   8.022   9.115  1.00  0.00          CU",
                    "TER       2      CU E   1",
                    "HETATM    3 C1   NAG B   1      12.500  10.000  10.000  1.00  0.00           C",
                    "HETATM    4 O4   NAG C   1      15.000  10.000  10.000  1.00  0.00           O",
                    "CONECT    1    3",
                    "CONECT    3    1    4",
                    "CONECT    4    3",
                    "END",
                ]
            )
            + "\n"
        )

        sanitization = _write_posebusters_ready_pdb(input_pdb, output_pdb)
        pdb_text = output_pdb.read_text()

        assert sanitization.stripped_metal_atom_count == 1
        assert sanitization.stripped_metal_elements == ("CU",)
        assert "HETATM    1 CU    CU E   1" not in pdb_text
        assert "TER       2      CU E   1" not in pdb_text
        assert "CONECT    3    4" in pdb_text
        assert "CONECT    4    3" in pdb_text

    def test_rewrite_conect_from_topology_adds_explicit_glycan_bonds(self, tmp_path: Path) -> None:
        """Sparse topology exports should be repaired with explicit glycan CONECT records."""
        from lpmo_pipeline.io.cif_to_pdb import _rewrite_conect_from_topology

        class _FakeAtom:
            def __init__(self, index: int) -> None:
                self.index = index

        class _FakeTopology:
            def __init__(self, atom_count: int, bond_pairs: list[tuple[int, int]] | None = None) -> None:
                self._atoms = [_FakeAtom(index) for index in range(atom_count)]
                self._bond_pairs = bond_pairs or []

            def atoms(self):
                return iter(self._atoms)

            def bonds(self):
                return iter(
                    (self._atoms[left], self._atoms[right])
                    for left, right in self._bond_pairs
                )

        pdb_path = tmp_path / "for_posebusters.pdb"
        pdb_path.write_text(
            "\n".join(
                [
                    "HETATM    1  C1  NAG B   1       0.000   0.000   0.000  1.00  0.00           C",
                    "HETATM    2  C2  NAG B   1       1.500   0.000   0.000  1.00  0.00           C",
                    "HETATM    3  C3  NAG B   1       2.500   1.000   0.000  1.00  0.00           C",
                    "HETATM    4  C4  NAG B   1       2.500   2.400   0.000  1.00  0.00           C",
                    "HETATM    5  C5  NAG B   1       1.300   3.100   0.000  1.00  0.00           C",
                    "HETATM    6  C6  NAG B   1       1.300   4.600   0.000  1.00  0.00           C",
                    "HETATM    7  O4  NAG B   1       3.700   3.100   0.000  1.00  0.00           O",
                    "HETATM    8  O5  NAG B   1       0.000   2.400   0.000  1.00  0.00           O",
                    "HETATM    9  C1  NAG B   2       5.000   3.100   0.000  1.00  0.00           C",
                    "HETATM   10  C2  NAG B   2       6.500   3.100   0.000  1.00  0.00           C",
                    "HETATM   11  C3  NAG B   2       7.500   4.100   0.000  1.00  0.00           C",
                    "HETATM   12  C4  NAG B   2       7.500   5.500   0.000  1.00  0.00           C",
                    "HETATM   13  C5  NAG B   2       6.300   6.200   0.000  1.00  0.00           C",
                    "HETATM   14  C6  NAG B   2       6.300   7.700   0.000  1.00  0.00           C",
                    "HETATM   15  O4  NAG B   2       8.700   6.200   0.000  1.00  0.00           O",
                    "HETATM   16  O5  NAG B   2       5.000   5.500   0.000  1.00  0.00           O",
                    "END",
                ]
            )
            + "\n"
        )

        _rewrite_conect_from_topology(pdb_path, _FakeTopology(atom_count=16, bond_pairs=[(0, 6)]))

        conect_map: dict[int, set[int]] = {}
        for line in pdb_path.read_text().splitlines():
            if not line.startswith("CONECT"):
                continue
            serials = [int(token) for token in line[6:].split()]
            if len(serials) < 2:
                continue
            conect_map.setdefault(serials[0], set()).update(serials[1:])

        assert 1 in conect_map
        assert 2 in conect_map[1]
        assert 8 in conect_map[1]
        assert 15 in conect_map[1]
        assert 7 not in conect_map[1]
        assert 4 in conect_map[7]

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
        assert "CU301" not in mol2_text
        assert " CU " not in mol2_text
