"""Integration test with fake runners."""

import tempfile
from pathlib import Path

import pytest

from structure_pipeline.config import PipelineConfig
from structure_pipeline.cases import Case
from structure_pipeline.runners.base import RunnerInterface, RunnerResult


class FakeRunner(RunnerInterface):
    """Fake runner for testing."""

    @property
    def model_name(self) -> str:
        return "fake"

    def build_input_file(self, cases, protein_sequences, output_dir):
        output_dir.mkdir(parents=True, exist_ok=True)
        input_file = output_dir / "fake_input.json"
        input_file.write_text('{"test": true}')
        return input_file

    def build_command(self, input_file, output_dir, **kwargs):
        return ["echo", "fake", "command"]

    def build_slurm_script(self, cases, protein_sequences, work_dir, job_name):
        return f"""#!/bin/bash
#SBATCH --job-name={job_name}
echo "Fake job for {len(cases)} cases"
touch {work_dir}/DONE.ok
"""

    def parse_outputs(self, output_dir, cases):
        return [
            RunnerResult(
                success=True,
                case_id=case.case_id,
                output_dir=output_dir,
            )
            for case in cases
        ]


class TestFakeRunner:
    """Tests using the fake runner."""

    def test_build_input_file(self):
        """Test building input file."""
        runner = FakeRunner(config=None)
        cases = [
            Case(
                case_id="test1",
                protein_id="P001",
                ligand_id="L001",
                ligand_ccd_code="AMY",
                model="fake",
                msa_source="none",
                msa_path="",
            )
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            protein_seqs = {"P001": "AAAA"}
            input_file = runner.build_input_file(cases, protein_seqs, output_dir)

            assert input_file.exists()
            assert input_file.suffix == ".json"

    def test_build_slurm_script(self):
        """Test building SLURM script."""
        runner = FakeRunner(config=None)
        cases = [
            Case(
                case_id="test1",
                protein_id="P001",
                ligand_id="L001",
                ligand_ccd_code="AMY",
                model="fake",
                msa_source="none",
                msa_path="",
            )
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            work_dir = Path(tmpdir)
            protein_seqs = {"P001": "AAAA"}
            script = runner.build_slurm_script(cases, protein_seqs, work_dir, "test_job")

            assert "#!/bin/bash" in script
            assert "#SBATCH" in script
            assert "test_job" in script
            assert "DONE.ok" in script
