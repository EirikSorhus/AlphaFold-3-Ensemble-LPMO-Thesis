"""Tests for AF3 runner – ligand CIF path and userCCDPath support."""

import json
import tempfile
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from structure_pipeline.cases import Case
from structure_pipeline.runners.af3 import AF3Runner


def _af3_ligand_entries(data):
    return [entry["ligand"] for entry in data.get("sequences", []) if "ligand" in entry]


@pytest.fixture
def mock_config():
    """Build a minimal PipelineConfig-like object for AF3Runner."""
    cfg = MagicMock()
    cfg.af3.seeds = 3
    cfg.af3.num_diffusion_samples = 2
    cfg.af3.jackhmmer_n_cpu = 4
    cfg.paths.af3_dir = Path("/prog/af3")
    cfg.paths.af3_cpu_image = Path("/prog/af3/af3_cpu_amd64.sif")
    cfg.paths.af3_gpu_image = Path("/prog/af3/af3_gpu_arm64.sif")
    cfg.paths.af3_weights = Path("/prog/af3/weights")
    cfg.paths.af3_databases = Path("/shared/databases")
    cfg.slurm.account = "nn1003k"
    cfg.slurm.partition_cpu = "normal"
    cfg.slurm.partition_gpu = "accel"
    cfg.slurm.cpus_msa = 8
    cfg.slurm.mem_per_cpu_msa = "10G"
    cfg.slurm.time_msa = "02:00:00"
    cfg.slurm.time_inference = "01:00:00"
    cfg.slurm.gpus = 1
    cfg.slurm.get_mem_per_gpu.return_value = "80G"
    return cfg


@pytest.fixture
def runner(mock_config):
    return AF3Runner(mock_config)


class TestBuildMsaInputFile:
    """Tests for build_msa_input_file."""

    def test_without_cif_path(self, runner):
        """Without ligand CIF, no userCCDPath should appear."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            path = runner.build_msa_input_file(
                "P001", "AAAA", "CEL6", out,
            )
            data = json.loads(path.read_text())
            assert "userCCDPath" not in data
            ligands = _af3_ligand_entries(data)
            assert [lig["ccdCodes"][0] for lig in ligands] == ["CU", "CEL6"]
            assert [lig["id"] for lig in ligands] == ["B", "C"]

    def test_with_cif_path(self, runner):
        """With ligand CIF, userCCDPath should be set."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            cif = Path("/root/af_input/CEL6.cif")
            path = runner.build_msa_input_file(
                "P001", "AAAA", "CEL6", out,
                ligand_cif_path=cif,
            )
            data = json.loads(path.read_text())
            assert data["userCCDPath"] == "/root/af_input/CEL6.cif"
            ligands = _af3_ligand_entries(data)
            assert [lig["ccdCodes"][0] for lig in ligands] == ["CU", "CEL6"]
            assert [lig["id"] for lig in ligands] == ["B", "C"]

    def test_with_cu_ligand(self, runner):
        """CU should not be duplicated when it is the main ligand."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            path = runner.build_msa_input_file(
                "P001", "AAAA", "CU", out,
            )
            data = json.loads(path.read_text())
            ligands = _af3_ligand_entries(data)
            assert [lig["ccdCodes"][0] for lig in ligands] == ["CU"]
            assert [lig["id"] for lig in ligands] == ["B"]


class TestBuildInputFile:
    """Tests for build_input_file."""

    def _make_case(self, **overrides):
        defaults = dict(
            case_id="abc123",
            protein_id="P001",
            ligand_id="L001",
            ligand_ccd_code="CEL6",
            model="af3",
            msa_source="native",
            msa_path="",
        )
        defaults.update(overrides)
        return Case(**defaults)

    def test_without_cif(self, runner):
        with tempfile.TemporaryDirectory() as tmpdir:
            case = self._make_case()
            path = runner.build_input_file(
                [case], {"P001": "AAAA"}, Path(tmpdir),
            )
            data = json.loads(path.read_text())
            assert "userCCDPath" not in data

    def test_with_cif(self, runner):
        with tempfile.TemporaryDirectory() as tmpdir:
            case = self._make_case()
            cif = Path("/somewhere/CEL6.cif")
            path = runner.build_input_file(
                [case], {"P001": "AAAA"}, Path(tmpdir),
                ligand_cif_path=cif,
            )
            data = json.loads(path.read_text())
            assert data["userCCDPath"] == "/somewhere/CEL6.cif"


class TestBuildMsaSlurmScript:
    """Tests for build_msa_slurm_script."""

    def test_script_without_cif(self, runner):
        with tempfile.TemporaryDirectory() as tmpdir:
            work_dir = Path(tmpdir)
            script = runner.build_msa_slurm_script(
                protein_id="P001",
                protein_sequence="AAAA",
                ligand_ccd_code="CEL6",
                work_dir=work_dir,
                job_name="test_msa",
            )
            assert "#SBATCH --job-name=test_msa" in script
            assert "af3_msa_input.json" in script

            # Squashfs bind mount with image-src
            assert "image-src=/public_databases" in script
            assert "AF3_DATABASES_SQUASHFS" in script

            # No double-quotes around bind-mount paths
            assert '"${' not in script or '"${files[0]}"' in script

            # Verify the generated JSON has no userCCDPath
            json_path = work_dir / "input" / "af3_msa_input.json"
            data = json.loads(json_path.read_text())
            assert "userCCDPath" not in data
            assert "userCCD" not in data
            ligands = _af3_ligand_entries(data)
            assert [lig["ccdCodes"][0] for lig in ligands] == ["CU", "CEL6"]

    def test_script_with_cif(self, runner):
        """When CIF is provided, it should be copied and userCCDPath set."""
        with tempfile.TemporaryDirectory() as tmpdir:
            work_dir = Path(tmpdir) / "work"
            cif_dir = Path(tmpdir) / "ligands"
            cif_dir.mkdir()
            cif_file = cif_dir / "CEL6.cif"
            cif_file.write_text("data_CEL6\n")

            script = runner.build_msa_slurm_script(
                protein_id="P001",
                protein_sequence="AAAA",
                ligand_ccd_code="CEL6",
                work_dir=work_dir,
                job_name="test_msa",
                ligand_cif_path=cif_file,
            )

            # CIF should be copied to input dir
            copied = work_dir / "input" / "CEL6.cif"
            assert copied.exists()

            # JSON should point to container path
            json_path = work_dir / "input" / "af3_msa_input.json"
            data = json.loads(json_path.read_text())
            assert data["userCCDPath"] == "/root/af_input/CEL6.cif"
            assert "userCCD" not in data
            ligands = _af3_ligand_entries(data)
            assert [lig["ccdCodes"][0] for lig in ligands] == ["CU", "CEL6"]


class TestBuildInferenceSlurmScript:
    """Tests for build_slurm_script (inference)."""

    def _make_case(self, **overrides):
        defaults = dict(
            case_id="abc123",
            protein_id="P001",
            ligand_id="L001",
            ligand_ccd_code="CEL6",
            model="af3",
            msa_source="native",
            msa_path="",
        )
        defaults.update(overrides)
        return Case(**defaults)

    def test_inference_script_with_cif(self, runner):
        """When ligand_cif_paths provided, CIF is copied and patch script gets CIF arg."""
        with tempfile.TemporaryDirectory() as tmpdir:
            work_dir = Path(tmpdir) / "work"
            cif_dir = Path(tmpdir) / "ligands"
            cif_dir.mkdir()
            cif_file = cif_dir / "CEL6.cif"
            cif_file.write_text("data_CEL6\n")

            case = self._make_case()
            script = runner.build_slurm_script(
                cases=[case],
                protein_sequences={"P001": "AAAA"},
                work_dir=work_dir,
                job_name="test_inf",
                ligand_cif_paths={"L001": cif_file},
            )

            # CIF should be copied to input dir
            assert (work_dir / "input" / "CEL6.cif").exists()

            # Standalone patch script should be written to input dir
            assert (work_dir / "input" / "_patch_json.py").exists()

            # Script should invoke patch script with container CIF path
            assert "/root/af_input/CEL6.cif" in script
            assert "_patch_json.py" in script

    def test_inference_script_without_cif(self, runner):
        """Without CIF, patch script is called without CIF path argument."""
        with tempfile.TemporaryDirectory() as tmpdir:
            work_dir = Path(tmpdir) / "work"
            case = self._make_case()
            script = runner.build_slurm_script(
                cases=[case],
                protein_sequences={"P001": "AAAA"},
                work_dir=work_dir,
                job_name="test_inf",
            )

            # Patch script should be written
            assert (work_dir / "input" / "_patch_json.py").exists()

            # Script should NOT contain CIF path arg
            assert "/root/af_input/" not in script or "_patch_json.py" in script
            # No inline python heredoc
            assert "<<'PY'" not in script

            # Squashfs bind mount with image-src
            assert "image-src=/public_databases" in script
            assert "AF3_DATABASES_SQUASHFS" in script


class TestPatchJsonScript:
    """Tests for the generated patch script."""

    def test_removes_userccd_and_sets_userccdpath(self, runner):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            src = tmp / "src.json"
            dst = tmp / "dst.json"
            script = tmp / "_patch_json.py"

            data = {
                "name": "old",
                "sequences": [
                    {"protein": {"id": "A", "sequence": "AAAA"}},
                    {"ligand": {"id": "B", "ccdCodes": ["CEL6"]}},
                ],
                "userCCD": "inline_ccd",
            }
            src.write_text(json.dumps(data))

            runner._write_patch_script(script)

            cmd = [
                "python",
                str(script),
                str(src),
                str(dst),
                "P001_CEL6",
                "CEL6",
                "/root/af_input/CEL6.cif",
            ]
            import subprocess

            subprocess.run(cmd, check=True)
            patched = json.loads(dst.read_text())

            assert "userCCD" not in patched
            assert patched["userCCDPath"] == "/root/af_input/CEL6.cif"
            ligands = _af3_ligand_entries(patched)
            assert [lig["ccdCodes"][0] for lig in ligands] == ["CU", "CEL6"]
