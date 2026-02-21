"""Tests for RF3 runner inputs."""

import json
from pathlib import Path
from unittest.mock import MagicMock

from structure_pipeline.cases import Case
from structure_pipeline.runners.rf3 import RF3Runner


def _make_case(ccd_code: str) -> Case:
    return Case(
        case_id="abc123",
        protein_id="P001",
        ligand_id="L001",
        ligand_ccd_code=ccd_code,
        model="rf3",
        msa_source="mmseqs",
        msa_path="",
    )


def _make_runner() -> RF3Runner:
    cfg = MagicMock()
    return RF3Runner(cfg)


class TestRf3MandatoryCu:
    def test_adds_cu_ligand(self, tmp_path: Path):
        runner = _make_runner()
        case = _make_case("CEL6")

        input_file = runner.build_input_file(
            cases=[case],
            protein_sequences={"P001": "AAAA"},
            output_dir=tmp_path,
            msa_paths={},
            ligand_cif_paths={},
        )

        data = json.loads(input_file.read_text())
        components = data[0]["components"]
        ligand_components = [c for c in components if "ccd_code" in c]
        assert [c["ccd_code"] for c in ligand_components] == ["CU", "CEL6"]

    def test_no_duplicate_when_cu_is_main(self, tmp_path: Path):
        runner = _make_runner()
        case = _make_case("CU")

        input_file = runner.build_input_file(
            cases=[case],
            protein_sequences={"P001": "AAAA"},
            output_dir=tmp_path,
            msa_paths={},
            ligand_cif_paths={},
        )

        data = json.loads(input_file.read_text())
        components = data[0]["components"]
        ligand_components = [c for c in components if "ccd_code" in c]
        assert [c["ccd_code"] for c in ligand_components] == ["CU"]
