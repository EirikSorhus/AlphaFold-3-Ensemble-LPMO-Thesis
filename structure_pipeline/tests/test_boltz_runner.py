"""Tests for Boltz runner inputs."""

import yaml
from pathlib import Path
from unittest.mock import MagicMock

from structure_pipeline.cases import Case
from structure_pipeline.runners.boltz import BoltzRunner


def _load_single_yaml(path: Path) -> dict:
	with open(path) as f:
		return yaml.safe_load(f)


def _make_case(ccd_code: str) -> Case:
	return Case(
		case_id="abc123",
		protein_id="P001",
		ligand_id="L001",
		ligand_ccd_code=ccd_code,
		model="boltz",
		msa_source="mmseqs",
		msa_path="",
	)


def _make_runner(use_affinity: bool = False) -> BoltzRunner:
	cfg = MagicMock()
	cfg.boltz.use_affinity = use_affinity
	return BoltzRunner(cfg)


class TestBoltzMandatoryCu:
	def test_adds_cu_ligand(self, tmp_path: Path):
		runner = _make_runner()
		case = _make_case("CEL6")

		out_dir = runner.build_input_file(
			cases=[case],
			protein_sequences={"P001": "AAAA"},
			output_dir=tmp_path,
		)

		data = _load_single_yaml(out_dir / "P001_CEL6.yaml")
		ligands = [entry["ligand"]["ccd"] for entry in data["sequences"] if "ligand" in entry]
		ids = [entry["ligand"]["id"] for entry in data["sequences"] if "ligand" in entry]

		assert ligands == ["CU", "CEL6"]
		assert ids == ["B", "C"]

	def test_no_duplicate_when_cu_is_main(self, tmp_path: Path):
		runner = _make_runner()
		case = _make_case("CU")

		out_dir = runner.build_input_file(
			cases=[case],
			protein_sequences={"P001": "AAAA"},
			output_dir=tmp_path,
		)

		data = _load_single_yaml(out_dir / "P001_CU.yaml")
		ligands = [entry["ligand"]["ccd"] for entry in data["sequences"] if "ligand" in entry]
		ids = [entry["ligand"]["id"] for entry in data["sequences"] if "ligand" in entry]

		assert ligands == ["CU"]
		assert ids == ["B"]

	def test_affinity_binder_tracks_main_ligand(self, tmp_path: Path):
		runner = _make_runner(use_affinity=True)
		case = _make_case("CEL6")

		out_dir = runner.build_input_file(
			cases=[case],
			protein_sequences={"P001": "AAAA"},
			output_dir=tmp_path,
		)

		data = _load_single_yaml(out_dir / "P001_CEL6.yaml")
		assert data["properties"][0]["affinity"]["binder"] == "C"
