"""AlphaFold3 runner implementation.

Two-stage fan-in architecture:
  1. MSA generation (CPU) -- one job per (protein, ligand) pair
     ``--run_inference=false`` produces ``*_data.json``
  2. Inference (GPU)      -- one job per *ligand* batching all proteins
     Uses ``_data.json`` from stage 1, with ``--dependency=afterok``

Key fixes vs original:
  - Inference now uses the ``_data.json`` produced by MSA stage instead
    of creating a brand-new JSON (which discards computed MSA).
  - Per-run output isolation via runs/{SLURM_JOB_ID}/.
  - Bind-mounts for input/output directories.
"""

from __future__ import annotations

import json
import logging
import textwrap
from pathlib import Path
from typing import Any

from ..cases import Case
from ..config import PipelineConfig
from .base import RunnerInterface, RunnerResult
from .ligand_utils import build_af3_ligand_entries
from .oligo import OligoRegistry

logger = logging.getLogger(__name__)


class AF3Runner(RunnerInterface):
    """Runner for AlphaFold3.

    AF3 has two stages:
      1. MSA generation (CPU) -- ``--run_inference=false``
      2. Inference (GPU)      -- uses pre-computed ``_data.json``

    The ``build_slurm_script`` method generates the *inference* script.
    Use ``build_msa_slurm_script`` for the MSA stage.
    """

    def __init__(
        self,
        config: PipelineConfig,
        oligo_registry: OligoRegistry | None = None,
    ) -> None:
        super().__init__(config)
        self.oligo_registry = oligo_registry or OligoRegistry.empty()

    @property
    def model_name(self) -> str:
        return "af3"

    # ------------------------------------------------------------------
    # Input file builders
    # ------------------------------------------------------------------

    def build_input_file(
        self,
        cases: list[Case],
        protein_sequences: dict[str, str],
        output_dir: Path,
        ligand_cif_path: Path | None = None,
    ) -> Path:
        """Build AF3 JSON input file for a single (protein, ligand) pair.

        AF3 input format (v4)::

            {
                "name": "B6EQJ6_STA6",
                "dialect": "alphafold3",
                "version": 4,
                "modelSeeds": [1, 2, ...],
                "sequences": [
                    {"protein": {"id": "A", "sequence": "..."}},
                    {"ligand":  {"id": "B", "ccdCodes": ["STA6"]}}
                ],
                "userCCDPath": "/path/to/custom.cif"  # optional
            }

        Args:
            cases: Cases for this batch.
            protein_sequences: protein_id -> sequence.
            output_dir: Directory to write input file.
            ligand_cif_path: Optional path to custom CIF file for the
                ligand.  When provided, ``userCCDPath`` is set in the
                output JSON so AF3 uses the custom CCD definition.
        """
        if not cases:
            raise ValueError("No cases provided")

        # For AF3 MSA, we create one JSON per (protein, ligand) pair.
        case = cases[0]
        protein_seq = protein_sequences[case.protein_id]
        seeds = list(range(1, self.config.af3.seeds + 1))

        ligand_entries, bonded_atom_pairs = build_af3_ligand_entries(
            case.ligand_ccd_code, oligo_registry=self.oligo_registry,
        )
        is_oligo = len(bonded_atom_pairs) > 0

        sequences: list[dict[str, Any]] = [
            {"protein": {"id": "A", "sequence": protein_seq}},
            *ligand_entries,
        ]

        input_data: dict[str, Any] = {
            "name": f"{case.protein_id}_{case.ligand_ccd_code}",
            "dialect": "alphafold3",
            "version": 4,
            "modelSeeds": seeds,
            "sequences": sequences,
        }

        if bonded_atom_pairs:
            input_data["bondedAtomPairs"] = bonded_atom_pairs

        # Oligo ligands use CCD monomers, not a custom CIF.
        if ligand_cif_path is not None and not is_oligo:
            input_data["userCCDPath"] = str(ligand_cif_path)
            logger.info("AF3 input: set userCCDPath for %s", case.ligand_ccd_code)
        elif is_oligo and ligand_cif_path is not None:
            logger.info(
                "AF3 input: skipping userCCDPath for oligo %s (using CCD monomers)",
                case.ligand_ccd_code,
            )

        output_dir.mkdir(parents=True, exist_ok=True)
        input_path = output_dir / "af3_input.json"
        with open(input_path, "w") as f:
            json.dump(input_data, f, indent=2)
        return input_path

    def build_msa_input_file(
        self,
        protein_id: str,
        protein_sequence: str,
        ligand_ccd_code: str,
        output_dir: Path,
        ligand_cif_path: Path | None = None,
    ) -> Path:
        """Build AF3 JSON for MSA-only stage.

        Includes the ligand so the resulting ``_data.json`` is complete
        and can be used directly for inference.

        Args:
            protein_id: Protein identifier.
            protein_sequence: Amino acid sequence.
            ligand_ccd_code: CCD code for the ligand.
            output_dir: Directory to write input file.
            ligand_cif_path: Optional path to custom CIF file.  When
                provided, ``userCCDPath`` is added to the JSON.
        """
        ligand_entries, bonded_atom_pairs = build_af3_ligand_entries(
            ligand_ccd_code, oligo_registry=self.oligo_registry,
        )
        is_oligo = len(bonded_atom_pairs) > 0

        input_data: dict[str, Any] = {
            "name": f"{protein_id}_{ligand_ccd_code}",
            "dialect": "alphafold3",
            "version": 4,
            "modelSeeds": list(range(1, self.config.af3.seeds + 1)),
            "sequences": [
                {"protein": {"id": "A", "sequence": protein_sequence}},
                *ligand_entries,
            ],
        }

        if bonded_atom_pairs:
            input_data["bondedAtomPairs"] = bonded_atom_pairs

        # Oligo ligands use CCD monomers, not a custom CIF.
        if ligand_cif_path is not None and not is_oligo:
            input_data["userCCDPath"] = str(ligand_cif_path)
            logger.info("AF3 MSA input: set userCCDPath for %s", ligand_ccd_code)
        elif is_oligo and ligand_cif_path is not None:
            logger.info(
                "AF3 MSA input: skipping userCCDPath for oligo %s (using CCD monomers)",
                ligand_ccd_code,
            )

        output_dir.mkdir(parents=True, exist_ok=True)
        input_path = output_dir / "af3_msa_input.json"
        with open(input_path, "w") as f:
            json.dump(input_data, f, indent=2)
        return input_path

    # ------------------------------------------------------------------
    # Command
    # ------------------------------------------------------------------

    def build_command(
        self,
        input_file: Path,
        output_dir: Path,
        msa_only: bool = False,
        **kwargs,
    ) -> list[str]:
        """Build AF3 command (run inside container)."""
        cmd = [
            "/opt/af3-venv/bin/python",
            "/opt/alphafold3/run_alphafold.py",
            f"--json_path=/root/af_input/{input_file.name}",
            "--model_dir=/root/models",
            "--db_dir=/root/public_databases",
            "--output_dir=/root/af_output",
        ]

        if msa_only:
            cmd.append("--run_inference=false")
            cmd.extend([
                "--jackhmmer_n_cpu=8",
                "--jackhmmer_binary_path=/hmmer/bin/jackhmmer",
                "--nhmmer_binary_path=/hmmer/bin/nhmmer",
                "--hmmalign_binary_path=/hmmer/bin/hmmalign",
                "--hmmsearch_binary_path=/hmmer/bin/hmmsearch",
                "--hmmbuild_binary_path=/hmmer/bin/hmmbuild",
            ])
        else:
            cmd.append(
                f"--num_diffusion_samples={self.config.af3.num_diffusion_samples}"
            )
            if self.config.af3.num_recycles is not None:
                cmd.append(f"--num_recycles={self.config.af3.num_recycles}")

        return cmd

    # ------------------------------------------------------------------
    # MSA SLURM script  (CPU, one per protein+ligand)
    # ------------------------------------------------------------------

    def build_msa_slurm_script(
        self,
        protein_id: str,
        protein_sequence: str,
        ligand_ccd_code: str,
        work_dir: Path,
        job_name: str,
        ligand_cif_path: Path | None = None,
    ) -> str:
        """Build SLURM script for AF3 MSA stage.

        Uses SCRATCH for the MSA run, then copies the ``_data.json``
        result back to ``work_dir/input/`` for the inference stage.

        Directory layout::

            work_dir = work/af3_msa/{protein_id}_{ligand_ccd}/
            work_dir/input/af3_msa_input.json   (original input)
            work_dir/input/produced_msa.json     (copied _data.json)
            work_dir/runs/{SLURM_JOB_ID}/       (scratch copy)

        Args:
            protein_id: Protein identifier.
            protein_sequence: Amino acid sequence.
            ligand_ccd_code: CCD code for the ligand.
            work_dir: Working directory for this MSA job.
            job_name: SLURM job name.
            ligand_cif_path: Optional path to custom CIF file for the
                ligand.  When provided the file is copied to scratch and
                bind-mounted so AF3 can resolve the custom CCD.
        """
        input_dir = work_dir / "input"
        outbase = work_dir / "runs"
        image = self.config.paths.af3_cpu_image
        weights = self.config.paths.af3_weights
        databases = self.config.paths.af3_databases

        # Check if this ligand is an oligosaccharide
        oligo_info = self.oligo_registry.parse(ligand_ccd_code)
        is_oligo = oligo_info is not None

        # Build input file (includes userCCDPath when ligand_cif_path is given)
        # When a custom CIF is provided, copy it into input_dir so it
        # gets rsync'd to scratch.  The JSON's userCCDPath is written
        # to point to the container-local path (/root/af_input/<name>).
        # For oligo ligands we skip CIF handling — the JSON uses CCD monomers.
        container_cif_path: str | None = None
        if ligand_cif_path is not None and not is_oligo:
            import shutil
            input_dir.mkdir(parents=True, exist_ok=True)
            dest = input_dir / ligand_cif_path.name
            shutil.copy2(ligand_cif_path, dest)
            container_cif_path = f"/root/af_input/{ligand_cif_path.name}"

        input_file = self.build_msa_input_file(
            protein_id, protein_sequence, ligand_ccd_code, input_dir,
            ligand_cif_path=Path(container_cif_path) if container_cif_path else None,
        )

        # Build MSA command args (for inside container)
        msa_cmd_args = [
            "/opt/af3-venv/bin/python /opt/alphafold3/run_alphafold.py",
            f"--json_path=/root/af_input/{input_file.name}",
            "--model_dir=/root/models",
            "--db_dir=/root/public_databases",
            "--output_dir=/root/af_output",
            "--run_inference=false",
            f"--jackhmmer_n_cpu={self.config.af3.jackhmmer_n_cpu}",
            "--jackhmmer_binary_path=/hmmer/bin/jackhmmer",
            "--nhmmer_binary_path=/hmmer/bin/nhmmer",
            "--hmmalign_binary_path=/hmmer/bin/hmmalign",
            "--hmmsearch_binary_path=/hmmer/bin/hmmsearch",
            "--hmmbuild_binary_path=/hmmer/bin/hmmbuild",
        ]
        msa_cmd_str = " \\\n    ".join(msa_cmd_args)

        script = f'''#!/usr/bin/env bash
#SBATCH --job-name={job_name}
#SBATCH --account={self.config.slurm.account}
#SBATCH --partition={self.config.slurm.partition_cpu}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={self.config.slurm.cpus_msa}
#SBATCH --mem-per-cpu={self.config.slurm.mem_per_cpu_msa}
#SBATCH --time={self.config.slurm.time_msa}
#SBATCH --output={work_dir}/slurm_%j.out
#SBATCH --error={work_dir}/slurm_%j.err

# ── Clean inherited container env (Tykky/NRIS wrapper) ─────
unset SINGULARITY_BIND SINGULARITY_CONTAINER SINGULARITY_NAME
unset SINGULARITY_COMMAND SINGULARITY_ENVIRONMENT
unset APPTAINER_BIND APPTAINER_CONTAINER APPTAINER_NAME
unset APPTAINER_COMMAND APPTAINER_ENVIRONMENT
unset SINGULARITYENV_PATH SINGULARITYENV_LD_LIBRARY_PATH

module load NRIS/CPU
set -euo pipefail

WORK_DIR={work_dir}
INPUT_DIR={input_dir}
OUTBASE={outbase}
OUT=${{OUTBASE}}/${{SLURM_JOB_ID}}

mkdir -p "${{OUT}}" "${{INPUT_DIR}}"

echo "──────────────────────────────────────────"
echo "AF3 MSA  |  protein={protein_id}  |  ligand={ligand_ccd_code}"
echo "OUT : $OUT"
echo "Time: $(date)"
echo "──────────────────────────────────────────"

# Setup scratch
TMP_ON_SCRATCH="${{SCRATCH}}/af3_tmp_${{SLURM_JOB_ID}}"
mkdir -p "${{SCRATCH}}/input" "${{SCRATCH}}/result" "${{TMP_ON_SCRATCH}}"
export TMPDIR="${{TMP_ON_SCRATCH}}"
export APPTAINERENV_TMPDIR="${{TMPDIR}}"

AF3_DIR={self.config.paths.af3_dir}
AF3_IMAGE=${{AF3_DIR}}/{image.name}
AF3_MODEL_PARAMETERS_DIR=${{AF3_DIR}}/weights
AF3_DATABASES_SQUASHFS={databases}

# Copy input to scratch
rsync -a --delete ${{INPUT_DIR}}/ ${{SCRATCH}}/input/

# Trap to copy _data.json back on exit
trap 'shopt -s nullglob; files=(${{SCRATCH}}/result/*/*_data.json); if (( ${{#files[@]}} )); then cp -f "${{files[0]}}" "${{INPUT_DIR}}/produced_msa.json"; cp -f "${{files[0]}}" "${{OUT}}/"; fi; ln -sfn "${{OUT}}" "${{WORK_DIR}}/latest"; touch "${{WORK_DIR}}/DONE.ok"' EXIT

apptainer exec --cleanenv \\
  --env PATH=/hmmer/bin:/opt/af3-venv/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \\
  --env JAX_PLATFORM_NAME=cpu \\
  --pwd /opt/alphafold3 \\
  --bind ${{SCRATCH}}/input:/root/af_input \\
  --bind ${{SCRATCH}}/result:/root/af_output \\
  --bind ${{TMP_ON_SCRATCH}}:/tmp \\
  --bind ${{AF3_MODEL_PARAMETERS_DIR}}:/root/models \\
  --bind ${{AF3_DATABASES_SQUASHFS}}:/root/public_databases:image-src=/public_databases \\
  ${{AF3_IMAGE}} \\
  {msa_cmd_str}

echo "AF3 MSA completed at $(date)"
'''
        return script

    # ------------------------------------------------------------------
    # JSON patching helper (written to input_dir at generation time)
    # ------------------------------------------------------------------

    @staticmethod
    def _write_patch_script(
        output_path: Path,
        oligo_ccd_codes: list[str] | None = None,
        oligo_bonded_atom_pairs: list[list[list]] | None = None,
    ) -> None:
        """Write a standalone JSON-patching helper to *output_path*.

        The generated script accepts positional args::

            src dst name ligand_ccd [user_ccd_path]

        It reads the MSA ``_data.json`` from *src*, patches
        name / ligand / CCD fields, and writes the result to *dst*.

        When *oligo_ccd_codes* and *oligo_bonded_atom_pairs* are
        provided, the ligand entry uses the multi-monomer ``ccdCodes``
        list and ``bondedAtomPairs`` is written at the top level.

        Using a standalone script (invoked with ``/usr/bin/python3``)
        avoids any PATH-dependent Python resolution — critical on
        mixed-architecture clusters where the submission environment
        may inject an AMD64 Python into ``$PATH``.
        """
        # Embed oligo data as Python literals when applicable.
        if oligo_ccd_codes and oligo_bonded_atom_pairs:
            oligo_codes_repr = repr(oligo_ccd_codes)
            oligo_bap_repr = repr(oligo_bonded_atom_pairs)
            script = textwrap.dedent(f"""\
                #!/usr/bin/env python3
                \"\"\"Patch AF3 MSA _data.json for inference (oligo mode) – auto-generated.\"\"\"
                import json, sys

                src, dst, name, ligand = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

                with open(src) as fh:
                    data = json.load(fh)

                data["name"] = name

                seqs = [entry for entry in data.get("sequences", []) if "ligand" not in entry]

                # CU is always present
                seqs.append({{"ligand": {{"id": "B", "ccdCodes": ["CU"]}}}})

                # Oligo ligand with multi-monomer ccdCodes
                oligo_codes = {oligo_codes_repr}
                seqs.append({{"ligand": {{"id": "C", "ccdCodes": oligo_codes}}}})

                data["sequences"] = seqs

                # bondedAtomPairs for the oligosaccharide chain
                data["bondedAtomPairs"] = {oligo_bap_repr}

                # Remove CIF-related keys – oligo uses CCD monomers
                for key in ("userCCD", "userCCDPath"):
                    if key in data:
                        print(f"Removing {{key}} from AF3 JSON (oligo mode)", file=sys.stderr)
                        del data[key]

                with open(dst, "w") as fh:
                    json.dump(data, fh, indent=2)
            """)
        else:
            script = textwrap.dedent("""\
                #!/usr/bin/env python3
                \"\"\"Patch AF3 MSA _data.json for inference – auto-generated.\"\"\"
                import json, sys

                src, dst, name, ligand = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
                user_ccd_path = sys.argv[5] if len(sys.argv) > 5 else ""

                with open(src) as fh:
                    data = json.load(fh)

                data["name"] = name

                seqs = [entry for entry in data.get("sequences", []) if "ligand" not in entry]

                if ligand == "CU":
                    seqs.append({"ligand": {"id": "B", "ccdCodes": ["CU"]}})
                else:
                    seqs.append({"ligand": {"id": "B", "ccdCodes": ["CU"]}})
                    seqs.append({"ligand": {"id": "C", "ccdCodes": [ligand]}})

                data["sequences"] = seqs

                if "userCCD" in data:
                    print("Removing inline userCCD from AF3 JSON", file=sys.stderr)
                    del data["userCCD"]

                if user_ccd_path:
                    data["userCCDPath"] = user_ccd_path
                    print(f"Set userCCDPath to {user_ccd_path}", file=sys.stderr)
                elif "userCCDPath" in data:
                    print("Removing userCCDPath (no CIF path provided)", file=sys.stderr)
                    del data["userCCDPath"]

                with open(dst, "w") as fh:
                    json.dump(data, fh, indent=2)
            """)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            f.write(script)
        output_path.chmod(0o755)

    # ------------------------------------------------------------------
    # Inference SLURM script  (GPU, one per ligand, all proteins)
    # ------------------------------------------------------------------

    def build_slurm_script(
        self,
        cases: list[Case],
        protein_sequences: dict[str, str],
        work_dir: Path,
        job_name: str,
        msa_work_dirs: dict[str, Path] | None = None,
        ligand_cif_paths: dict[str, Path] | None = None,
    ) -> str:
        """Build SLURM script for AF3 inference stage.

        Processes all proteins for one ligand sequentially in a single
        GPU job.  For each protein, uses the ``_data.json`` produced by
        the corresponding MSA job.

        Args:
            cases: Cases for this ligand (different proteins).
            protein_sequences: protein_id -> sequence.
            work_dir: ``work/{ligand_ccd}/af3/``.
            job_name: SLURM job name.
            msa_work_dirs: protein_id -> MSA work_dir path, where
                ``input/produced_msa.json`` is expected.
            ligand_cif_paths: ligand_id -> CIF file path.  When provided
                the CIF is copied into the input dir and
                ``userCCDPath`` in the inference JSON is rewritten to
                the container-local path.
        """
        if not cases:
            raise ValueError("No cases provided")

        ligand_ccd = cases[0].ligand_ccd_code
        if msa_work_dirs is None:
            msa_work_dirs = {}
        if ligand_cif_paths is None:
            ligand_cif_paths = {}

        input_dir = work_dir / "input"
        outbase = work_dir / "runs"
        image = self.config.paths.af3_gpu_image
        weights = self.config.paths.af3_weights
        databases = self.config.paths.af3_databases
        mem = self.config.slurm.get_mem_per_gpu("af3")

        input_dir.mkdir(parents=True, exist_ok=True)

        # Check if this ligand is an oligosaccharide
        oligo_info = self.oligo_registry.parse(ligand_ccd)
        is_oligo = oligo_info is not None

        # Copy custom CIF files into input_dir for bind-mounting
        # Map ligand_id -> container path (only for ligands in this batch)
        # For oligo ligands we skip CIF — the JSON uses CCD monomers.
        container_cif_map: dict[str, str] = {}
        if not is_oligo:
            batch_ligand_ids = {case.ligand_id for case in cases}
            for lig_id, cif_path in ligand_cif_paths.items():
                if lig_id not in batch_ligand_ids:
                    continue
                import shutil
                dest = input_dir / cif_path.name
                shutil.copy2(cif_path, dest)
                container_cif_map[lig_id] = f"/root/af_input/{cif_path.name}"

        protein_list = ", ".join(c.protein_id for c in cases[:6])
        if len(cases) > 6:
            protein_list += f" ... (+{len(cases)-6})"

        # Build the inference loop: for each protein, copy its _data.json
        # into the input dir and run AF3 inference.
        # Determine container CIF path for this ligand (if any)
        container_cif = ""
        if not is_oligo:
            for case in cases:
                cif = container_cif_map.get(case.ligand_id, "")
                if cif:
                    container_cif = cif
                    break

        # Write standalone JSON-patch helper into input_dir so the bash
        # script can invoke it with /usr/bin/python3 (guaranteed native
        # architecture) instead of a PATH-dependent ``python3``.
        # For oligo ligands, embed the monomer codes and bond pairs.
        patch_script = input_dir / "_patch_json.py"
        if is_oligo:
            _prefix, n, spec = oligo_info  # type: ignore[misc]
            from .ligand_utils import AF3_MAIN_LIGAND_ID
            from .oligo import build_af3_bonded_atom_pairs
            oligo_codes = [spec.monomer] * n
            oligo_bap = build_af3_bonded_atom_pairs(AF3_MAIN_LIGAND_ID, spec, n)
            self._write_patch_script(
                patch_script,
                oligo_ccd_codes=oligo_codes,
                oligo_bonded_atom_pairs=oligo_bap,
            )
        else:
            self._write_patch_script(patch_script)

        inference_blocks = []
        for case in cases:
            pid = case.protein_id
            msa_dir = msa_work_dirs.get(pid)
            if msa_dir:
                data_json = msa_dir / "input" / "produced_msa.json"
            else:
                data_json = Path(f"MISSING_MSA_FOR_{pid}")

            job_input_name = f"{pid}_{ligand_ccd}_data.json"

            # Build patch command with explicit /usr/bin/python3
            patch_args = (
                f'"$DATA_JSON" "${{INPUT_DIR}}/{job_input_name}"'
                f' "{pid}_{ligand_ccd}" "{ligand_ccd}"'
            )
            if container_cif:
                patch_args += f' "{container_cif}"'

            # Build num_recycles arg conditionally
            num_recycles_arg = ""
            if self.config.af3.num_recycles is not None:
                num_recycles_arg = f" \\\n      --num_recycles={self.config.af3.num_recycles}"

            inference_blocks.append(f'''
echo "── Inference: {pid} + {ligand_ccd} ──"
DATA_JSON={data_json}
if [[ ! -f "$DATA_JSON" ]]; then
  echo "WARNING: MSA data not found for {pid}: $DATA_JSON — skipping" >&2
else
  /usr/bin/python3 "${{INPUT_DIR}}/_patch_json.py" {patch_args}

  mkdir -p "${{OUT}}"

  apptainer exec --nv --cleanenv \\
    --pwd /opt/alphafold3 \\
    --bind ${{INPUT_DIR}}:/root/af_input \\
    --bind ${{OUT}}:/root/af_output \\
    --bind ${{AF3_MODEL_PARAMETERS_DIR}}:/root/models \\
    --bind ${{AF3_DATABASES_SQUASHFS}}:/root/public_databases:image-src=/public_databases \\
    ${{AF3_IMAGE}} \\
    /opt/af3-venv/bin/python /opt/alphafold3/run_alphafold.py \\
      --json_path=/root/af_input/{job_input_name} \\
      --model_dir=/root/models \\
      --db_dir=/root/public_databases \\
      --num_diffusion_samples={self.config.af3.num_diffusion_samples}{num_recycles_arg} \\
      --output_dir=/root/af_output

  echo "  done: {pid}"
fi''')

        inference_loop = "\n".join(inference_blocks)

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

WORK_DIR={work_dir}
INPUT_DIR={input_dir}
OUTBASE={outbase}
OUT=${{OUTBASE}}/${{SLURM_JOB_ID}}

AF3_DIR={self.config.paths.af3_dir}
AF3_IMAGE=${{AF3_DIR}}/{image.name}
AF3_MODEL_PARAMETERS_DIR=${{AF3_DIR}}/weights
AF3_DATABASES_SQUASHFS={databases}

mkdir -p "${{OUT}}" "${{INPUT_DIR}}"

echo "──────────────────────────────────────────"
echo "AF3 inference  |  ligand={ligand_ccd}  |  proteins={len(cases)}"
echo "Proteins: {protein_list}"
echo "OUT     : $OUT"
echo "Time    : $(date)"
echo "──────────────────────────────────────────"
{inference_loop}

# ── bookkeeping ────────────────────────────────────────────
ln -sfn "${{OUT}}" "${{WORK_DIR}}/latest"
touch "${{WORK_DIR}}/DONE.ok"

echo "AF3 inference completed at $(date)"
echo "Results in: ${{OUT}}"
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
        """Parse AF3 outputs."""
        results = []

        for case in cases:
            result = RunnerResult(
                success=False,
                case_id=case.case_id,
                output_dir=output_dir,
            )

            # AF3 creates: {output_dir}/{name}/{name}_model.cif
            name = f"{case.protein_id}_{case.ligand_ccd_code}"
            case_dir = output_dir / name

            if case_dir.exists():
                # Confidence
                conf_files = list(case_dir.glob("*_summary_confidences_*.json"))
                if conf_files:
                    try:
                        with open(conf_files[0]) as f:
                            conf_data = json.load(f)
                        result.confidence_scores = {
                            "ptm": conf_data.get("ptm", 0.0),
                            "iptm": conf_data.get("iptm", 0.0),
                            "ranking_score": conf_data.get("ranking_score", 0.0),
                        }
                    except Exception:
                        pass

                # Structure files
                structure_files = list(case_dir.glob("*_model.cif"))
                if structure_files:
                    result.output_files = [str(f) for f in structure_files]
                    result.success = True

            results.append(result)

        return results

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate_paths(self) -> list[str]:
        """Validate AF3 paths."""
        errors = []

        def check(path: Path, name: str) -> None:
            try:
                if not path.exists():
                    errors.append(f"{name} not found: {path}")
            except PermissionError:
                pass

        check(self.config.paths.af3_cpu_image, "AF3 CPU image")
        check(self.config.paths.af3_gpu_image, "AF3 GPU image")
        check(self.config.paths.af3_weights, "AF3 weights")
        check(self.config.paths.af3_databases, "AF3 databases squashfs")

        return errors
