"""RoseTTAFold3 runner implementation.

Matches the proven batch-processing pattern from rf3_B6EQJ6_direct_input.sh:
  - OUTBASE / OUT=${OUTBASE}/${SLURM_JOB_ID}  (per-run isolation)
  - Foundry env vars, Triton side-install
  - Explicit bind-mounts computed from all referenced paths
  - One JSON with one example per protein, all sharing the same ligand
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ..cases import Case
from ..config import PipelineConfig
from .base import RunnerInterface, RunnerResult
from .ligand_utils import build_rf3_ligand_components


class RF3Runner(RunnerInterface):
    """Runner for RoseTTAFold3 (RF3).

    Batching strategy (per-ligand):
        One SLURM job per ligand.  All proteins are combined into a single
        JSON input array so that the expensive model-load happens only once.
    """

    @property
    def model_name(self) -> str:
        return "rf3"

    # ------------------------------------------------------------------
    # Input file
    # ------------------------------------------------------------------

    def build_input_file(
        self,
        cases: list[Case],
        protein_sequences: dict[str, str],
        output_dir: Path,
        msa_paths: dict[str, str] | None = None,
        ligand_cif_paths: dict[str, Path] | None = None,
    ) -> Path:
        """Build RF3 JSON input file.

        Creates one example per protein, all sharing the same ligand CCD.

        Args:
            cases: Cases for this batch (same ligand, different proteins).
            protein_sequences: protein_id -> amino acid sequence.
            output_dir: Directory to write input file.
            msa_paths: protein_id -> a3m file path (optional).
            ligand_cif_paths: ligand_id -> CIF file path (optional).

        Returns:
            Path to generated JSON input file.
        """
        if not cases:
            raise ValueError("No cases provided")

        if msa_paths is None:
            msa_paths = {}
        if ligand_cif_paths is None:
            ligand_cif_paths = {}

        examples: list[dict[str, Any]] = []

        for case in cases:
            protein_seq = protein_sequences[case.protein_id]
            components: list[dict[str, Any]] = []

            # Protein component
            protein_component: dict[str, Any] = {
                "seq": protein_seq,
                "chain_id": "A",
            }
            msa = msa_paths.get(case.protein_id) or case.msa_path
            if msa:
                protein_component["msa_path"] = str(msa)
            components.append(protein_component)

            # Ligand components - prefer CIF path, fall back to CCD code
            main_cif_path = None
            if case.ligand_id in ligand_cif_paths:
                main_cif_path = str(ligand_cif_paths[case.ligand_id])
            components.extend(
                build_rf3_ligand_components(case.ligand_ccd_code, main_cif_path)
            )

            examples.append(
                {
                    "name": f"{case.protein_id}_{case.ligand_ccd_code}",
                    "components": components,
                }
            )

        output_dir.mkdir(parents=True, exist_ok=True)
        input_path = output_dir / "rf3_input.json"
        with open(input_path, "w") as f:
            json.dump(examples, f, indent=2)
        return input_path

    # ------------------------------------------------------------------
    # Command
    # ------------------------------------------------------------------

    def build_command(
        self,
        input_file: Path,
        output_dir: Path,
        **kwargs,
    ) -> list[str]:
        """Build Hydra-style RF3 fold command (run inside container)."""
        return [
            "python",
            f'"{self.config.paths.rf3_foundry_root}/models/rf3/src/rf3/cli.py"',
            "fold",
            f'inputs="{input_file}"',
            f'ckpt_path="{self.config.paths.rf3_checkpoint}"',
            f'out_dir="{output_dir}"',
            f"n_recycles={self.config.rf3.n_recycles}",
            f"diffusion_batch_size={self.config.rf3.diffusion_batch_size}",
            f"num_steps={self.config.rf3.num_steps}",
            f"seed={self.config.rf3.seed}",
            f"early_stopping_plddt_threshold={self.config.rf3.early_stopping_plddt_threshold}",
            f"skip_existing={str(self.config.rf3.skip_existing)}",
            "dump_predictions=True",
            "dump_trajectories=False",
        ]

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
        ligand_cif_paths: dict[str, Path] | None = None,
    ) -> str:
        """Build SLURM script for RF3.

        Produces a script that closely matches the proven
        ``rf3_B6EQJ6_direct_input.sh`` pattern:

        * ``OUTBASE = {work_dir}/runs``
        * ``OUT = $OUTBASE/$SLURM_JOB_ID``   (per-run isolation)
        * ``latest`` symlink updated after success
        * Explicit bind-mounts computed from all referenced paths

        Args:
            cases: All cases for this ligand (different proteins).
            protein_sequences: protein_id -> sequence.
            work_dir: ``work/{ligand_ccd}/rf3/``.
            job_name: SLURM job name.
            msa_paths: protein_id -> a3m file path.
            msa_sqsh_file: Optional squashfs overlay for MSAs.
            ligand_cif_paths: ligand_id -> CIF file path.
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

        # -- build input JSON -------------------------------------------
        input_file = self.build_input_file(
            cases,
            protein_sequences,
            input_dir,
            msa_paths=job_msa_paths,
            ligand_cif_paths=ligand_cif_paths,
        )

        # -- paths ------------------------------------------------------
        foundry_root = self.config.paths.rf3_foundry_root
        checkpoint = self.config.paths.rf3_checkpoint
        image = self.config.paths.rf3_image
        tmp_base = self.config.slurm.tmp_base
        mem = self.config.slurm.get_mem_per_gpu("rf3")

        # -- summary for log --------------------------------------------
        protein_list = ", ".join(c.protein_id for c in cases[:6])
        if len(cases) > 6:
            protein_list += f" ... (+{len(cases)-6})"

        triton_version = "3.5.0"
        ckpt_dir = checkpoint.parent

        # Build the script as a plain string (no f-string) for the bash
        # body inside bash -lc "…" to avoid double-escaping issues.
        # All {var} are Python .format() substitutions; all ${var} are
        # literal bash and must be written as ${{var}} in the f-string.
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

require_file() {{
  local path="$1"
  local label="$2"
  if [[ ! -f "${{path}}" ]]; then
    echo "ERROR: Missing ${{label}}: ${{path}}" >&2
    exit 1
  fi
}}

# ── directories ────────────────────────────────────────────
WORK_DIR={work_dir}
INPUT_DIR={input_dir}
INPUT_JSON={input_file}
OUTBASE={outbase}
OUT=${{OUTBASE}}/${{SLURM_JOB_ID}}

FOUNDRY_ROOT={foundry_root}
CKPT_DIR={ckpt_dir}
CKPT_FILE={checkpoint}
SIF={image}
TMPBASE={tmp_base}
TMPDIR=${{TMPBASE}}/foundry_${{USER}}

mkdir -p "${{OUT}}" "${{OUTBASE}}" "${{TMPDIR}}"

echo "──────────────────────────────────────────"
echo "RF3  |  ligand={ligand_ccd}  |  proteins={len(cases)}"
echo "Proteins: {protein_list}"
echo "OUTBASE : $OUTBASE"
echo "OUT     : $OUT"
echo "Time    : $(date)"
echo "──────────────────────────────────────────"

# ── validate files ─────────────────────────────────────────
require_file "${{SIF}}" "Apptainer image"
require_file "${{CKPT_FILE}}" "RF3 checkpoint"
require_file "${{INPUT_JSON}}" "RF3 input JSON"

# ── environment ────────────────────────────────────────────
set -a
[[ -f "${{FOUNDRY_ROOT}}/.env" ]] && source "${{FOUNDRY_ROOT}}/.env"
set +a

export APPTAINERENV_FOUNDRY_CHECKPOINT_DIRS="${{CKPT_DIR}}"
export APPTAINERENV_XDG_CACHE_HOME="${{OUT}}/cache"
export APPTAINERENV_TMPDIR="${{TMPDIR}}"
export APPTAINERENV_PYTHONPATH="${{FOUNDRY_ROOT}}/src:${{FOUNDRY_ROOT}}/models/rf3/src:${{FOUNDRY_ROOT}}/models/mpnn/src:${{FOUNDRY_ROOT}}/models/rfd3/src:${{PYTHONPATH:-}}"
export APPTAINERENV_PROJECT_ROOT="${{FOUNDRY_ROOT}}"
export APPTAINERENV_HYDRA_FULL_ERROR=1
mkdir -p "${{OUT}}/cache"

TRITON_VERSION={triton_version}
TRITON_SITE=${{TMPDIR}}/triton_site
mkdir -p "${{TRITON_SITE}}"

# ── run RF3 ────────────────────────────────────────────────
echo "=== Run RF3 fold ==="
apptainer exec --nv \\
  --bind "${{FOUNDRY_ROOT}}:${{FOUNDRY_ROOT}}" \\
  --bind "${{TMPBASE}}:${{TMPBASE}}" \\
  --bind "/cluster/work/projects/nn1003k:/cluster/work/projects/nn1003k" \\
  "${{SIF}}" bash -lc "
    set -e
    export PYTHONPATH=\\"${{TRITON_SITE}}:\\${{PYTHONPATH:-}}\\"
    if ! python - <<'PY'
from importlib.metadata import version, PackageNotFoundError
import sys
try:
    v = version('triton')
    sys.exit(0 if v.startswith('3.5.') else 1)
except PackageNotFoundError:
    sys.exit(1)
PY
    then
      echo 'Installing triton==${{TRITON_VERSION}} into ${{TRITON_SITE}}...'
      pip install --no-cache-dir --target \\"${{TRITON_SITE}}\\" \\"triton==${{TRITON_VERSION}}\\"
      export PYTHONPATH=\\"${{TRITON_SITE}}:\\${{PYTHONPATH:-}}\\"
    fi
    python \\"${{FOUNDRY_ROOT}}/models/rf3/src/rf3/cli.py\\" fold \\
      inputs=\\"${{INPUT_JSON}}\\" \\
      ckpt_path=\\"${{CKPT_FILE}}\\" \\
      out_dir=\\"${{OUT}}\\" \\
      skip_existing={str(self.config.rf3.skip_existing)} \\
      dump_predictions=True \\
      dump_trajectories=False \\
      n_recycles={self.config.rf3.n_recycles} \\
      diffusion_batch_size={self.config.rf3.diffusion_batch_size} \\
      num_steps={self.config.rf3.num_steps} \\
      seed={self.config.rf3.seed} \\
      early_stopping_plddt_threshold={self.config.rf3.early_stopping_plddt_threshold}
  "

# ── bookkeeping ────────────────────────────────────────────
ln -sfn "${{OUT}}" "${{WORK_DIR}}/latest"
touch "${{WORK_DIR}}/DONE.ok"

echo "Done. Outputs in: ${{OUT}}"
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
        """Parse RF3 outputs for each case."""
        results = []

        for case in cases:
            result = RunnerResult(
                success=False,
                case_id=case.case_id,
                output_dir=output_dir,
            )

            pattern = f"{case.protein_id}_{case.ligand_ccd_code}"

            # Structure files
            cif_files = list(output_dir.glob(f"**/{pattern}*_model_*.cif*"))
            if cif_files:
                result.output_files = [str(f) for f in cif_files]
                result.success = True

            # Metrics CSV
            metrics_files = list(output_dir.glob(f"**/{pattern}*_metrics.csv"))
            if metrics_files:
                try:
                    import csv as csv_mod

                    with open(metrics_files[0], newline="") as fh:
                        reader = csv_mod.DictReader(fh)
                        for row in reader:
                            result.confidence_scores = {
                                "ptm": float(row.get("pTM", 0.0)),
                                "plddt": float(row.get("pLDDT", 0.0)),
                            }
                            break
                except Exception:
                    pass

            # Early-stopping marker
            score_files = list(output_dir.glob(f"**/{pattern}*.score"))
            if score_files:
                try:
                    content = score_files[0].read_text()
                    if "early_stopped" in content.lower():
                        result.metadata["early_stopped"] = True
                except Exception:
                    pass

            results.append(result)

        return results

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate_paths(self) -> list[str]:
        """Validate RF3 paths."""
        errors = []

        def check(path: Path, name: str) -> None:
            try:
                if not path.exists():
                    errors.append(f"{name} not found: {path}")
            except PermissionError:
                pass

        check(self.config.paths.rf3_foundry_root, "RF3 Foundry root")
        check(self.config.paths.rf3_image, "RF3 image")
        check(self.config.paths.rf3_checkpoint, "RF3 checkpoint")

        return errors
