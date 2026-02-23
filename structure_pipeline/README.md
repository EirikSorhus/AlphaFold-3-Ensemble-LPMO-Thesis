# Structure Pipeline

Multi-model structure prediction pipeline for running **AlphaFold3**, **Boltz-2**, and **RoseTTAFold3** on LPMO proteins with ligand binding on the Olivia / Betzy HPC cluster (SLURM).

## Limitations

- **Single-chain proteins only** – multi-chain paired MSAs are not supported.
- **Pre-computed MSAs required** – Boltz-2 and RF3 need `.a3m` files from mmseqs / ColabFold.
- **CCD ligands** – ligands must be provided as CIF files in CCD format.

---

## Installation on Olivia (HPC)

### Option 1: Use Pre-built Container (recommended)

```bash
export PATH="/cluster/work/projects/nn1003k/eirik/conda/structure_pipeline_env/bin:$PATH"
export PATH="/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/bin:$PATH"
structure-pipeline --help
```

### Option 2: Build Your Own Container

```bash
module load NRIS/CPU
module load hpc-container-wrapper

export http_proxy=http://10.63.2.48:3128/
export https_proxy=http://10.63.2.48:3128/

conda-containerize new --mamba --prefix /path/to/your/install env.yml
export PATH="/path/to/your/install/bin:$PATH"
pip install -e .
```

> **Note**: Use `python -m structure_pipeline.cli` or the `bin/structure-pipeline` wrapper rather than the pip-installed entry point (container path issues).

---

## Quick Start

### 1. Setup Configuration

```bash
structure-pipeline init-config --output config/pipeline.yaml
# Edit paths to match your environment
```

### 2. Prepare Input Data

```
structure_pipeline/
├── input/
│   ├── fasta/
│   │   └── proteins.fasta       # FASTA sequences
│   ├── ligand/
│   │   ├── amylose/             # Category subdirectories
│   │   │   └── STA6.cif
│   │   ├── cellulose/
│   │   │   └── CEL6.cif
│   │   └── chitin/
│   │       └── NAG6.cif
│   └── msa/
│       ├── B6EQJ6.a3m           # Pre-computed mmseqs MSAs
│       └── D0EW65.a3m
```

### 3. Generate Manifests

```bash
structure-pipeline manifest -c config/pipeline.yaml
```

Creates `manifest/{proteins,ligands,msa,cases}.csv`.

### 4. Validate

```bash
structure-pipeline validate -c config/pipeline.yaml
```

### 5. Run

```bash
# Submit all jobs
structure-pipeline run -c config/pipeline.yaml

# Dry run
structure-pipeline run -c config/pipeline.yaml --dry-run

# Single model / protein
structure-pipeline run -c config/pipeline.yaml --model rf3
structure-pipeline run -c config/pipeline.yaml --protein B6EQJ6

# Resume (default: --resume, skips completed jobs)
structure-pipeline run -c config/pipeline.yaml --resume
```

### 6. Status

```bash
structure-pipeline status -c config/pipeline.yaml
```

---

## Architecture

### Batching Strategy

**One SLURM job per (ligand, model)** – all proteins are batched together.

| Model | Stage | Grouping |
|-------|-------|----------|
| **RF3** | single GPU job | All proteins for one ligand in one JSON |
| **Boltz-2** | single GPU job | Directory of YAMLs (one per protein) |
| **AF3** | MSA (CPU) + Inference (GPU) | MSA: one job per (protein, ligand); Inference: one job per ligand with `--dependency=afterok` |

### AF3 Fan-in

```
MSA jobs (CPU, per protein+ligand)
  ├── af3_msa_B6EQJ6_STA6  ──→  _data.json
  ├── af3_msa_D0EW65_STA6  ──→  _data.json
  └── af3_msa_Q9RFX5_STA6  ──→  _data.json
              │
              ▼  --dependency=afterok:ID1:ID2:ID3
Inference job (GPU, per ligand)
  └── af3_STA6_inf  (loops through all _data.json files)
```

### Directory Layout

```
work/
├── {ligand_ccd}/                    # e.g. STA6/
│   ├── af3/
│   │   ├── input/
│   │   ├── runs/{SLURM_JOB_ID}/    # per-run isolation
│   │   ├── latest -> runs/12345    # symlink to most recent run
│   │   ├── DONE.ok
│   │   └── status.jsonl
│   ├── boltz/
│   │   ├── input/                   # one YAML per protein
│   │   ├── runs/{SLURM_JOB_ID}/
│   │   ├── latest -> runs/12346
│   │   ├── DONE.ok
│   │   └── status.jsonl
│   └── rf3/
│       ├── input/rf3_input.json
│       ├── runs/{SLURM_JOB_ID}/
│       ├── latest -> runs/12347
│       ├── DONE.ok
│       └── status.jsonl
├── af3_msa/                         # AF3 MSA stage outputs
│   ├── B6EQJ6_STA6/
│   │   ├── input/
│   │   │   ├── af3_msa_input.json
│   │   │   └── produced_msa.json    # _data.json copied here
│   │   ├── runs/{SLURM_JOB_ID}/
│   │   ├── latest -> runs/12340
│   │   └── DONE.ok
│   └── D0EW65_STA6/
│       └── ...
```

### Bind-Mount Optimisation

The pipeline computes **minimal common-prefix bind-mounts** from all referenced paths:

```
Input paths:
  /cluster/projects/nn1003k/prog/foundry/checkpoints/ckpt.pt
  /cluster/projects/nn1003k/prog/foundry/models/rf3/src
  /cluster/work/projects/nn1003k/eirik/tmp

Computed mounts (2 instead of 3):
  --bind /cluster/projects/nn1003k/prog/foundry:/cluster/projects/nn1003k/prog/foundry
  --bind /cluster/work/projects/nn1003k/eirik/tmp:/cluster/work/projects/nn1003k/eirik/tmp
```

Use `--dry-run` to see the mount plan without submitting jobs.

---

## Configuration Reference

### `pipeline.yaml`

```yaml
paths:
  af3_dir: /cluster/projects/nn1003k/prog/af3
  af3_cpu_image: /cluster/projects/nn1003k/prog/af3/af3_cpu_amd64.sif
  af3_gpu_image: /cluster/projects/nn1003k/prog/af3/af3_gpu_arm64.sif
  af3_weights: /cluster/projects/nn1003k/prog/af3/weights
  af3_databases: /cluster/shared/alphafold/public_databases

  boltz_image: /cluster/projects/nn1003k/prog/boltz/boltz2_alt.sif
  boltz_weights: /cluster/projects/nn1003k/prog/boltz/weights

  rf3_foundry_root: /cluster/projects/nn1003k/prog/foundry
  rf3_image: /cluster/projects/nn1003k/prog/foundry/foundry.sif
  rf3_checkpoint: /cluster/projects/nn1003k/prog/foundry/checkpoints/rf3_foundry_01_24_latest_remapped.ckpt

slurm:
  account: nn1003k
  partition_cpu: normal
  partition_gpu: accel
  mem_per_gpu: 20G                # Default for all models
  # mem_per_gpu_af3: 32G          # Per-model override
  # mem_per_gpu_boltz: 20G
  # mem_per_gpu_rf3: 20G
  time_msa: "02:00:00"
  time_inference: "01:00:00"
  gpus: 1
  cpus_msa: 8
  mem_per_cpu_msa: 10G
  tmp_base: /cluster/work/projects/nn1003k/eirik/tmp

af3:
  seeds: 10
  num_diffusion_samples: 5

boltz:
  diffusion_samples: 5
  recycling_steps: 10
  use_affinity: true              # Via YAML 'properties' field (not CLI flag)
  use_msa_server: false

rf3:
  seed: 42
  diffusion_batch_size: 5
  n_recycles: 10
  num_steps: 50
  early_stopping_plddt_threshold: 0.5
  skip_existing: false

inputs:
  proteins_fasta: input/fasta/proteins.fasta
  ligands_dir: input/ligand
  msa_dir: input/msa

outputs:
  work_dir: work
  manifest_dir: manifest
  results_dir: results

models_enabled: [af3, boltz, rf3]
```

### Per-model Memory

```yaml
slurm:
  mem_per_gpu: 20G          # fallback
  mem_per_gpu_af3: 32G      # AF3 needs more VRAM
  mem_per_gpu_boltz: 20G
  mem_per_gpu_rf3: 20G
```

---

## CLI Reference

| Command | Description |
|---------|-------------|
| `init-config` | Generate example config file |
| `manifest` | Parse inputs → generate case matrix |
| `validate` | Check all paths and dependencies |
| `run` | Submit prediction jobs to SLURM |
| `status` | Show job status summary |
| `version` | Show version |

### Run Options

```
  -c, --config PATH               Config YAML
  -n, --dry-run                    Preview jobs without submitting
  -r, --resume       Skip completed jobs (default: resume)
  -m, --model TEXT                 af3 / boltz / rf3
  -p, --protein TEXT               Single protein ID
  --local                          Run as subprocess (no SLURM)
  --af3-seeds INT
  --af3-diffusion-samples INT
  --boltz-diffusion-samples INT
  --boltz-recycling-steps INT
  --rf3-seed INT
  --rf3-diffusion-batch-size INT
  --rf3-n-recycles INT
  --rf3-num-steps INT
```

---

## FASTA Format

```
>UniProtIDs_Q1K4Q1;Q873G1_Neurospora_crassa_AA9_family_LPMO
MVTPETRGRLKGFAAPPPARSGRGQLLLLLLL
```

The first UniProt ID (before `;`) is used as the primary `protein_id`.

## MSA Matching

MSA files are matched to proteins by:

1. Parsing the query header of each `.a3m` file
2. Extracting UniProt-like IDs from headers and filenames
3. Matching against `protein_id`

## Resume Behaviour

- Each SLURM job writes `DONE.ok` on success.
- `--resume` (default) skips groups with existing `DONE.ok`.
- Individual case events are logged in `status.jsonl`.
- `--no-resume` re-submits everything.

## Model-Specific Notes

### AlphaFold3

- **Two-stage**: MSA (CPU) → Inference (GPU) with SLURM dependency.
- MSA stage produces `_data.json`; inference uses it directly (no re-generation).
- Uses native AF3 MSA (jackhmmer/nhmmer), not pre-computed mmseqs.

### Boltz-2

- Single GPU job.  Pre-computed MSA from `msa_dir`.
- **Affinity prediction** via YAML `properties` field (not `--affinity` CLI flag, which does not exist).
- One YAML input file per protein; `boltz predict <dir>` processes all.

### RoseTTAFold3

- Single GPU job via Foundry framework with Hydra CLI.
- Triton is side-installed into a tmpdir if not present.
- Batch JSON with one example per protein amortises model-load cost.
- Uses `foundry.sif` (not `foundry_pyg_24.11.sif`).

---

## Testing

```bash
export PATH="/cluster/work/projects/nn1003k/eirik/conda/structure_pipeline_env/bin:$PATH"
cd /cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline
pytest tests/ -v
```

## Troubleshooting

### "Missing MSA" warnings

Proteins without matching `.a3m` files will skip Boltz/RF3. Check filenames and query headers.

### SLURM job failures

Logs are in:
- `work/{ligand_ccd}/{model}/slurm_*.out`
- `work/{ligand_ccd}/{model}/slurm_*.err`
- AF3 MSA logs: `work/af3_msa/{protein}_{ligand}/slurm_*.out`

### Container path issues

All paths in `pipeline.yaml` must be absolute and accessible from compute nodes.

---

## License

MIT
