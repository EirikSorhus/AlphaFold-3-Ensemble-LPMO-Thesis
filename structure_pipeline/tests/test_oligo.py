"""Tests for oligosaccharide registry, parsing, and AF3 builder helpers."""

import json
import subprocess
import tempfile
import textwrap
from pathlib import Path

import pytest

from structure_pipeline.runners.oligo import (
    OligoRegistry,
    OligoSpec,
    build_af3_bonded_atom_pairs,
    build_af3_oligo_ligand_entry,
)
from structure_pipeline.runners.ligand_utils import build_af3_ligand_entries


# ── OligoSpec validation ───────────────────────────────────────────


class TestOligoSpec:
    def test_valid_spec(self):
        spec = OligoSpec(prefix="CEL", monomer="GLC", bond_atom_pair=("C1", "O4"))
        assert spec.prefix == "CEL"
        assert spec.monomer == "GLC"
        assert spec.bond_atom_pair == ("C1", "O4")

    def test_empty_prefix_raises(self):
        with pytest.raises(ValueError, match="prefix must be non-empty"):
            OligoSpec(prefix="", monomer="GLC", bond_atom_pair=("C1", "O4"))

    def test_empty_monomer_raises(self):
        with pytest.raises(ValueError, match="monomer must be non-empty"):
            OligoSpec(prefix="CEL", monomer="", bond_atom_pair=("C1", "O4"))

    def test_empty_atom_name_raises(self):
        with pytest.raises(ValueError, match="atom names must be non-empty"):
            OligoSpec(prefix="CEL", monomer="GLC", bond_atom_pair=("C1", ""))

    def test_wrong_length_bond_atom_pair_raises(self):
        with pytest.raises(ValueError, match="exactly 2 elements"):
            OligoSpec(prefix="CEL", monomer="GLC", bond_atom_pair=("C1", "O4", "X"))


# ── OligoRegistry.parse ───────────────────────────────────────────


class TestOligoRegistryParse:
    @pytest.fixture
    def registry(self):
        return OligoRegistry(
            specs={
                "NAG": OligoSpec(prefix="NAG", monomer="NAG", bond_atom_pair=("C1", "O4")),
                "CEL": OligoSpec(prefix="CEL", monomer="GLC", bond_atom_pair=("C1", "O4")),
                "STA": OligoSpec(prefix="STA", monomer="BGC", bond_atom_pair=("C3", "O6")),
            }
        )

    def test_parse_known_oligo(self, registry):
        result = registry.parse("CEL6")
        assert result is not None
        prefix, n, spec = result
        assert prefix == "CEL"
        assert n == 6
        assert spec.monomer == "GLC"

    def test_parse_case_insensitive_prefix(self, registry):
        """CCD codes from CIF are uppercase, but test lowercase input."""
        result = registry.parse("cel6")
        assert result is not None
        assert result[0] == "CEL"
        assert result[1] == 6

    def test_parse_unknown_prefix(self, registry):
        assert registry.parse("XYZ6") is None

    def test_parse_no_number(self, registry):
        assert registry.parse("CEL") is None

    def test_parse_n_equals_one_returns_none(self, registry):
        """n=1 is a single monomer, not an oligo."""
        assert registry.parse("CEL1") is None

    def test_parse_n_equals_two(self, registry):
        result = registry.parse("NAG2")
        assert result is not None
        assert result[1] == 2

    def test_parse_n_large(self, registry):
        result = registry.parse("STA12")
        assert result is not None
        assert result[1] == 12

    def test_parse_regular_ccd_code(self, registry):
        """Regular CCD codes like CU, ATP should return None."""
        assert registry.parse("CU") is None
        assert registry.parse("ATP") is None

    def test_empty_registry_returns_none(self):
        empty = OligoRegistry.empty()
        assert empty.parse("CEL6") is None


# ── OligoRegistry.from_yaml ───────────────────────────────────────


class TestOligoRegistryFromYaml:
    def test_load_valid_yaml(self):
        yaml_content = textwrap.dedent("""\
            oligo_definitions:
              NAG:
                monomer: NAG
                bond_atom_pair: ["C1", "O4"]
              CEL:
                monomer: GLC
                bond_atom_pair: ["C3", "O6"]
        """)
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False
        ) as f:
            f.write(yaml_content)
            path = Path(f.name)

        try:
            reg = OligoRegistry.from_yaml(path)
            assert len(reg.specs) == 2
            assert reg.specs["NAG"].monomer == "NAG"
            assert reg.specs["NAG"].bond_atom_pair == ("C1", "O4")
            assert reg.specs["CEL"].monomer == "GLC"
            assert reg.specs["CEL"].bond_atom_pair == ("C3", "O6")
        finally:
            path.unlink()

    def test_load_uppercase_conversion(self):
        """Prefix and monomer are uppercased."""
        yaml_content = textwrap.dedent("""\
            oligo_definitions:
              sta:
                monomer: bgc
                bond_atom_pair: ["C1", "O4"]
        """)
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False
        ) as f:
            f.write(yaml_content)
            path = Path(f.name)

        try:
            reg = OligoRegistry.from_yaml(path)
            assert "STA" in reg.specs
            assert reg.specs["STA"].monomer == "BGC"
        finally:
            path.unlink()

    def test_missing_monomer_raises(self):
        yaml_content = textwrap.dedent("""\
            oligo_definitions:
              NAG:
                bond_atom_pair: ["C1", "O4"]
        """)
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False
        ) as f:
            f.write(yaml_content)
            path = Path(f.name)

        try:
            with pytest.raises(ValueError, match="Missing 'monomer'"):
                OligoRegistry.from_yaml(path)
        finally:
            path.unlink()

    def test_missing_bond_atom_pair_raises(self):
        yaml_content = textwrap.dedent("""\
            oligo_definitions:
              NAG:
                monomer: NAG
        """)
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False
        ) as f:
            f.write(yaml_content)
            path = Path(f.name)

        try:
            with pytest.raises(ValueError, match="bond_atom_pair"):
                OligoRegistry.from_yaml(path)
        finally:
            path.unlink()

    def test_empty_yaml_raises(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False
        ) as f:
            f.write("")
            path = Path(f.name)

        try:
            with pytest.raises(ValueError, match="Empty YAML"):
                OligoRegistry.from_yaml(path)
        finally:
            path.unlink()

    def test_no_oligo_definitions_key_raises(self):
        yaml_content = "some_other_key: {}\n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False
        ) as f:
            f.write(yaml_content)
            path = Path(f.name)

        try:
            with pytest.raises(ValueError, match="No 'oligo_definitions'"):
                OligoRegistry.from_yaml(path)
        finally:
            path.unlink()


# ── AF3 oligo ligand entry builder ────────────────────────────────


class TestBuildAf3OligoLigandEntry:
    def test_basic(self):
        spec = OligoSpec(prefix="CEL", monomer="GLC", bond_atom_pair=("C1", "O4"))
        entry = build_af3_oligo_ligand_entry("C", spec, 4)
        assert entry == {"ligand": {"id": "C", "ccdCodes": ["GLC", "GLC", "GLC", "GLC"]}}

    def test_n_equals_2(self):
        spec = OligoSpec(prefix="NAG", monomer="NAG", bond_atom_pair=("C1", "O4"))
        entry = build_af3_oligo_ligand_entry("C", spec, 2)
        assert entry["ligand"]["ccdCodes"] == ["NAG", "NAG"]


# ── AF3 bondedAtomPairs builder ───────────────────────────────────


class TestBuildAf3BondedAtomPairs:
    def test_n2(self):
        spec = OligoSpec(prefix="CEL", monomer="GLC", bond_atom_pair=("C1", "O4"))
        pairs = build_af3_bonded_atom_pairs("C", spec, 2)
        assert pairs == [[["C", 1, "C1"], ["C", 2, "O4"]]]

    def test_n3(self):
        spec = OligoSpec(prefix="CEL", monomer="GLC", bond_atom_pair=("C1", "O4"))
        pairs = build_af3_bonded_atom_pairs("C", spec, 3)
        assert len(pairs) == 2
        assert pairs[0] == [["C", 1, "C1"], ["C", 2, "O4"]]
        assert pairs[1] == [["C", 2, "C1"], ["C", 3, "O4"]]

    def test_n6(self):
        spec = OligoSpec(prefix="STA", monomer="BGC", bond_atom_pair=("C3", "O6"))
        pairs = build_af3_bonded_atom_pairs("C", spec, 6)
        assert len(pairs) == 5
        # Check first and last
        assert pairs[0] == [["C", 1, "C3"], ["C", 2, "O6"]]
        assert pairs[4] == [["C", 5, "C3"], ["C", 6, "O6"]]

    def test_n8(self):
        spec = OligoSpec(prefix="NAG", monomer="NAG", bond_atom_pair=("C1", "O4"))
        pairs = build_af3_bonded_atom_pairs("C", spec, 8)
        assert len(pairs) == 7
        # Verify all residue indices are 1-based and sequential
        for idx, pair in enumerate(pairs):
            assert pair[0][1] == idx + 1  # residue i
            assert pair[1][1] == idx + 2  # residue i+1


# ── Integration: build_af3_ligand_entries with oligo ──────────────


class TestBuildAf3LigandEntriesWithOligo:
    @pytest.fixture
    def registry(self):
        return OligoRegistry(
            specs={
                "NAG": OligoSpec(prefix="NAG", monomer="NAG", bond_atom_pair=("C1", "O4")),
                "CEL": OligoSpec(prefix="CEL", monomer="GLC", bond_atom_pair=("C1", "O4")),
                "STA": OligoSpec(prefix="STA", monomer="BGC", bond_atom_pair=("C3", "O6")),
            }
        )

    def test_oligo_ligand_returns_monomers_and_bap(self, registry):
        entries, bap = build_af3_ligand_entries("CEL6", oligo_registry=registry)
        # Should have CU + oligo ligand
        assert len(entries) == 2
        assert entries[0]["ligand"]["ccdCodes"] == ["CU"]
        assert entries[1]["ligand"]["ccdCodes"] == ["GLC"] * 6
        assert entries[1]["ligand"]["id"] == "C"
        # bondedAtomPairs should have 5 entries
        assert len(bap) == 5

    def test_regular_ligand_returns_single_ccd_no_bap(self, registry):
        entries, bap = build_af3_ligand_entries("ATP", oligo_registry=registry)
        assert len(entries) == 2
        assert entries[1]["ligand"]["ccdCodes"] == ["ATP"]
        assert bap == []

    def test_cu_ligand(self, registry):
        entries, bap = build_af3_ligand_entries("CU", oligo_registry=registry)
        assert len(entries) == 1
        assert entries[0]["ligand"]["ccdCodes"] == ["CU"]
        assert bap == []

    def test_no_registry_backward_compat(self):
        """Without a registry, behaviour is identical to original."""
        entries, bap = build_af3_ligand_entries("CEL6")
        assert len(entries) == 2
        assert entries[1]["ligand"]["ccdCodes"] == ["CEL6"]
        assert bap == []

    def test_none_registry_backward_compat(self):
        entries, bap = build_af3_ligand_entries("CEL6", oligo_registry=None)
        assert len(entries) == 2
        assert entries[1]["ligand"]["ccdCodes"] == ["CEL6"]
        assert bap == []


# ── Golden JSON tests ─────────────────────────────────────────────


class TestGoldenJson:
    """Verify exact JSON structure for known oligo ligands."""

    @pytest.fixture
    def registry(self):
        return OligoRegistry(
            specs={
                "NAG": OligoSpec(prefix="NAG", monomer="NAG", bond_atom_pair=("C1", "O4")),
            }
        )

    def test_nag3_golden(self, registry):
        entries, bap = build_af3_ligand_entries("NAG3", oligo_registry=registry)

        expected_sequences = [
            {"ligand": {"id": "B", "ccdCodes": ["CU"]}},
            {"ligand": {"id": "C", "ccdCodes": ["NAG", "NAG", "NAG"]}},
        ]
        expected_bap = [
            [["C", 1, "C1"], ["C", 2, "O4"]],
            [["C", 2, "C1"], ["C", 3, "O4"]],
        ]

        assert entries == expected_sequences
        assert bap == expected_bap

    def test_nag6_golden_full_json(self, registry):
        """Build a complete AF3-like JSON and verify structure."""
        entries, bap = build_af3_ligand_entries("NAG6", oligo_registry=registry)

        full_json = {
            "name": "P001_NAG6",
            "dialect": "alphafold3",
            "version": 4,
            "modelSeeds": [1, 2, 3],
            "sequences": [
                {"protein": {"id": "A", "sequence": "MVTPE"}},
                *entries,
            ],
        }
        if bap:
            full_json["bondedAtomPairs"] = bap

        # Verify round-trip through JSON serialization
        serialized = json.dumps(full_json, indent=2)
        roundtripped = json.loads(serialized)

        assert roundtripped["sequences"][2]["ligand"]["ccdCodes"] == ["NAG"] * 6
        assert len(roundtripped["bondedAtomPairs"]) == 5
        assert "userCCDPath" not in roundtripped
