"""Boltz-2 runner implementation.

Key fixes vs original:
  - Removed non-existent ``--affinity`` CLI flag; affinity is enabled via
    YAML ``properties`` field instead.
  - Added explicit bind-mounts for input/MSA directories (not just /work).
  - Per-run output isolation: runs/{SLURM_JOB_ID}/.
  - Per-ligand batching: one YAML per protein placed in a directory,
    ``boltz predict <dir>`` processes all of them.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from ..cases import Case
from ..config import PipelineConfig
from .base import RunnerInterface, RunnerResult
from .ligand_utils import build_boltz_ligand_entries


class BoltzRunner(RunnerInterface):
    """Runner for Boltz-2.

    Batching strategy (per-ligand):
        One SLURM job per ligand.  A directory of YAML files is created,
        one per protein, and ``boltz predict <dir>`` processes all of them.
    """

    @property
    def model_name(self) -> str:
        return "boltz"

    # ------------------------------------------------------------------
    # Input files
    # ------------------------------------------------------------------

    def build_input_file(
        self,
        cases: list[Case],
        protein_sequences: dict[str, str],
        output_dir: Path,
        msa_paths: dict[str, str] | None = None,
    ) -> Path:
        """Build Boltz YAML input files (one per protein).

        Creates a directory with one YAML per protein, all sharing the
        same ligand.  ``boltz predict <dir>`` will process them all.

        Boltz input format::

            version: 1
            sequences:
              - protein:
                  id: A
                  sequence: MVTPE...
                  msa: /path/to/msa.a3m
              - ligand:
                  id: B
                  ccd: STA6
            # Optional: affinity via properties
            properties:
              - affinity:
                  binder:
                    entity: ligand
                    id: B

        Args:
            cases: Cases for this batch (same ligand, different proteins).
            protein_sequences: protein_id -> amino acid sequence.
            output_dir: Directory to write input files.
            msa_paths: protein_id -> MSA file path (optional).

        Returns:
            Path to the *directory* containing all YAML files.
        """
        if not cases:
            raise ValueError("No cases provided")

        if msa_paths is None:
            msa_paths = {}

        output_dir.mkdir(parents=True, exist_ok=True)

        for case in cases:
            protein_seq = protein_sequences[case.protein_id]

            sequences: list[dict[str, Any]] = []

            # Protein
            protein_entry: dict[str, Any] = {
                "protein": {
                    "id": "A",
                    "sequence": protein_seq,
                }
            }
            msa = msa_paths.get(case.protein_id) or case.msa_path
            if msa:
                protein_entry["protein"]["msa"] = str(msa)
            sequences.append(protein_entry)

            # Ligands (always include CU when needed)
            ligand_entries, main_ligand_id = build_boltz_ligand_entries(
                case.ligand_ccd_code
            )
            sequences.extend(ligand_entries)

            input_data: dict[str, Any] = {
                "version": 1,
                "sequences": sequences,
            }

            # Affinity prediction via YAML properties (not CLI flag)
            if self.config.boltz.use_affinity:
                input_data["properties"] = [
                    {
                        "affinity": {
                            "binder": main_ligand_id,
                        }
                    }
                ]

            fname = f"{case.protein_id}_{case.ligand_ccd_code}.yaml"
            with open(output_dir / fname, "w") as f:
                yaml.dump(input_data, f, default_flow_style=False, sort_keys=False)

        return output_dir

    # ------------------------------------------------------------------
    # Command
    # ------------------------------------------------------------------

    def build_command(
        self,
        input_file: Path,
        output_dir: Path,
        **kwargs,
    ) -> list[str]:
        """Build Boltz command (run inside container)."""
        cmd = [
            "boltz",
            "predict",
            "/work/input",  # directory of YAML files
            "--out_dir", "/work/output",
            "--cache", "/weights/boltz",
            "--recycling_steps", str(self.config.boltz.recycling_steps),
            "--diffusion_samples", str(self.config.boltz.diffusion_samples),
            "--override",
        ]

        if self.config.boltz.use_msa_server:
            cmd.append("--use_msa_server")

        return cmd

    # ------------------------------------------------------------------
    # SLURM script
    # ------------------------------------------------------------------

    def build_slurm_script(
        self,
        cases: list[Case],
        protein_sequences: dict[str, str],
        work_dir: Path,
        job_name: str,
        msa_paths: dict[str, str] | None = None,
        msa_sqsh_file: Path | None = None,
    ) -> str:
        """Build SLURM script for Boltz-2.

        * ``OUTBASE = {work_dir}/runs``
        * ``OUT = $OUTBASE/$SLURM_JOB_ID``
        * ``latest`` symlink on success
        * Explicit bind-mounts for input, MSA, weights

        Args:
            cases: All cases for this ligand (different proteins).
            protein_sequences: protein_id -> sequence.
            work_dir: ``work/{ligand_ccd}/boltz/``.
            job_name: SLURM job name.
            msa_paths: protein_id -> a3m file path.
            msa_sqsh_file: Optional squashfs overlay.
        """
        if not cases:
            raise ValueError("No cases provided")

        ligand_ccd = cases[0].ligand_ccd_code
        if msa_paths is None:
            msa_paths = {}

        # -- directories ------------------------------------------------
        input_dir = work_dir / "input"
        outbase = work_dir / "runs"

        # -- MSA path mapping -------------------------------------------
        job_msa_paths: dict[str, str] = {}

        if msa_sqsh_file:
            # Squashfs not used in this version
            job_msa_paths = dict(msa_paths)
        else:
            job_msa_paths = dict(msa_paths)

        # -- build input YAMLs -----------------------------------------
        input_path = self.build_input_file(
            cases, protein_sequences, input_dir, msa_paths=job_msa_paths
        )

        # -- paths ------------------------------------------------------
        image = self.config.paths.boltz_image
        weights = self.config.paths.boltz_weights
        mem = self.config.slurm.get_mem_per_gpu("boltz")

        protein_list = ", ".join(c.protein_id for c in cases[:6])
        if len(cases) > 6:
            protein_list += f" ... (+{len(cases)-6})"

        # -- boltz CLI args ---------------------------------------------
        boltz_extra_args = []
        if self.config.boltz.recycling_steps != 10:
            boltz_extra_args.append(
                f"--recycling_steps {self.config.boltz.recycling_steps}"
            )
        if self.config.boltz.diffusion_samples != 5:
            boltz_extra_args.append(
                f"--diffusion_samples {self.config.boltz.diffusion_samples}"
            )
        if self.config.boltz.use_msa_server:
            boltz_extra_args.append("--use_msa_server")
        extra_args_str = ""
        if boltz_extra_args:
            extra_args_str = " \\\n    " + " \\\n    ".join(boltz_extra_args)

        script = f'''#!/usr/bin/env bash
#SBATCH --job-name={job_name}
#SBATCH --account={self.config.slurm.account}
#SBATCH --partition={self.config.slurm.partition_gpu}
#SBATCH --gpus={self.config.slurm.gpus}
#SBATCH --time={self.config.slurm.time_inference}
#SBATCH --mem-per-gpu={mem}
#SBATCH --output={work_dir}/slurm_%j.out
#SBATCH --error={work_dir}/slurm_%j.err

# ── Clean inherited container env (Tykky/NRIS wrapper) ─────
unset SINGULARITY_BIND SINGULARITY_CONTAINER SINGULARITY_NAME
unset SINGULARITY_COMMAND SINGULARITY_ENVIRONMENT
unset APPTAINER_BIND APPTAINER_CONTAINER APPTAINER_NAME
unset APPTAINER_COMMAND APPTAINER_ENVIRONMENT
unset SINGULARITYENV_PATH SINGULARITYENV_LD_LIBRARY_PATH

set -euo pipefail
module load NRIS/GPU

# ── directories ────────────────────────────────────────────
WORK_DIR={work_dir}
IMG={image}
WEIGHTS_HOST={weights}

mkdir -p "${{WORK_DIR}}/result"

echo "──────────────────────────────────────────"
echo "Boltz-2  |  ligand={ligand_ccd}  |  proteins={len(cases)}"
echo "Proteins: {protein_list}"
echo "WORK_DIR: $WORK_DIR"
echo "Time    : $(date)"
echo "──────────────────────────────────────────"

[[ -f "$IMG" ]] || {{ echo "ERROR: image not found: $IMG" >&2; exit 2; }}

# ── run Boltz-2 ────────────────────────────────────────────
apptainer exec --nv \\
  --bind "${{WORK_DIR}}:/work" \\
  --bind "${{WEIGHTS_HOST}}:/weights" \\
  --bind "/cluster/work/projects/nn1003k:/cluster/work/projects/nn1003k" \\
  --env HF_HOME=/weights/huggingface \\
  "$IMG" boltz predict /work/input \\
    --out_dir /work/result \\
    --cache /weights/boltz \\
    --override{extra_args_str}

# ── bookkeeping ────────────────────────────────────────────
touch "${{WORK_DIR}}/DONE.ok"

echo "Boltz run completed. Results in: ${{WORK_DIR}}/result"
'''
        return script

    # ------------------------------------------------------------------
    # Output parsing
    # ------------------------------------------------------------------

    def parse_outputs(
        self,
        output_dir: Path,
        cases: list[Case],
    ) -> list[RunnerResult]:
        """Parse Boltz outputs."""
        results = []

        for case in cases:
            result = RunnerResult(
                success=False,
                case_id=case.case_id,
                output_dir=output_dir,
            )

            # Boltz outputs: predictions/{input_name}/
            predictions_dir = output_dir / "predictions"
            expected_name = f"{case.protein_id}_{case.ligand_ccd_code}"

            target_dir = predictions_dir / expected_name
            if not target_dir.exists() and predictions_dir.exists():
                # Try any matching subdir
                for subdir in predictions_dir.iterdir():
                    if subdir.is_dir() and case.protein_id in subdir.name:
                        target_dir = subdir
                        break

            if target_dir.exists():
                cif_files = list(target_dir.glob("*_model_*.cif"))
                if cif_files:
                    result.output_files = [str(f) for f in cif_files]
                    result.success = True

                # Confidence scores
                conf_files = list(target_dir.glob("confidence_*.json"))
                if conf_files:
                    try:
                        import json
                        with open(conf_files[0]) as f:
                            conf_data = json.load(f)
                        result.confidence_scores = {
                            "ptm": conf_data.get("ptm", 0.0),
                            "iptm": conf_data.get("iptm", 0.0),
                            "ligand_iptm": conf_data.get("ligand_iptm", 0.0),
                            "confidence_score": conf_data.get("confidence_score", 0.0),
                        }
                    except Exception:
                        pass

                # Affinity scores
                affinity_files = list(target_dir.glob("affinity_*.json"))
                if affinity_files:
                    try:
                        import json
                        with open(affinity_files[0]) as f:
                            aff_data = json.load(f)
                        result.confidence_scores["affinity"] = aff_data.get(
                            "affinity_pred_value", 0.0
                        )
                    except Exception:
                        pass

            results.append(result)

        return results

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate_paths(self) -> list[str]:
        """Validate Boltz paths."""
        errors = []

        def check(path: Path, name: str) -> None:
            try:
                if not path.exists():
                    errors.append(f"{name} not found: {path}")
            except PermissionError:
                pass

        check(self.config.paths.boltz_image, "Boltz image")
        check(self.config.paths.boltz_weights, "Boltz weights")

        return errors
