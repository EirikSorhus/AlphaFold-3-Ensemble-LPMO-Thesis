"""Configuration management with Pydantic validation."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import yaml
from pydantic import BaseModel, Field, field_validator, model_validator


class PathsConfig(BaseModel):
    """Paths to model weights, containers, and databases."""

    # AF3 paths
    af3_dir: Path = Field(description="Root directory for AF3 installation")
    af3_cpu_image: Path = Field(description="Path to AF3 CPU container image")
    af3_gpu_image: Path = Field(description="Path to AF3 GPU container image")
    af3_weights: Path = Field(description="Path to AF3 model weights")
    af3_databases: Path = Field(
        description="Path to AF3 databases squashfs image"
    )

    # Boltz paths
    boltz_image: Path = Field(description="Path to Boltz container image")
    boltz_weights: Path = Field(description="Path to Boltz weights directory")

    # RF3 paths
    rf3_foundry_root: Path = Field(description="Path to Foundry installation root")
    rf3_image: Path = Field(description="Path to RF3 container image")
    rf3_checkpoint: Path = Field(description="Path to RF3 checkpoint file")

    @field_validator("*", mode="before")
    @classmethod
    def expand_env_vars(cls, v: Any) -> Any:
        """Expand environment variables in paths."""
        if isinstance(v, str):
            return os.path.expandvars(v)
        return v


class SlurmConfig(BaseModel):
    """SLURM job configuration."""

    account: str = Field(default="nn1003k", description="SLURM account")
    partition_cpu: str = Field(default="normal", description="CPU partition name")
    partition_gpu: str = Field(default="accel", description="GPU partition name")
    mem_per_gpu: str = Field(default="80G", description="Default memory per GPU")
    mem_per_gpu_af3: str | None = Field(
        default=None, description="Memory per GPU for AF3 (overrides mem_per_gpu)"
    )
    mem_per_gpu_boltz: str | None = Field(
        default=None, description="Memory per GPU for Boltz (overrides mem_per_gpu)"
    )
    mem_per_gpu_rf3: str | None = Field(
        default=None, description="Memory per GPU for RF3 (overrides mem_per_gpu)"
    )
    time_msa: str = Field(default="02:00:00", description="Time limit for MSA jobs")
    time_inference: str = Field(
        default="01:00:00", description="Time limit for inference jobs"
    )
    gpus: int = Field(default=1, description="Number of GPUs per job")
    cpus_msa: int = Field(default=8, description="CPUs for MSA generation")
    mem_per_cpu_msa: str = Field(default="10G", description="Memory per CPU for MSA")
    tmp_base: Path = Field(
        default=Path("/cluster/work/projects/nn1003k/eirik/tmp"),
        description="Base directory for temporary files on the cluster",
    )

    def get_mem_per_gpu(self, model: str) -> str:
        """Get memory per GPU for a specific model, with fallback."""
        override = getattr(self, f"mem_per_gpu_{model}", None)
        return override if override is not None else self.mem_per_gpu


class AF3Config(BaseModel):
    """AF3-specific configuration."""

    seeds: int = Field(default=10, ge=1, description="Number of model seeds")
    num_diffusion_samples: int = Field(
        default=5, ge=1, description="Diffusion samples per seed"
    )
    num_recycles: int | None = Field(
        default=None, ge=1, description="Number of recycles (None = use AF3 default)"
    )
    # MSA-specific
    jackhmmer_n_cpu: int = Field(default=8, description="CPUs for jackhmmer")


class BoltzConfig(BaseModel):
    """Boltz-2 specific configuration."""

    seeds: int = Field(default=5, ge=1, description="Number of seeds (recycling_steps)")
    diffusion_samples: int = Field(default=5, ge=1, description="Diffusion samples")
    recycling_steps: int = Field(default=10, ge=1, description="Recycling steps")
    sampling_steps: int | None = Field(
        default=None, ge=1, description="Sampling steps (None = use Boltz default)"
    )
    use_potentials: bool | None = Field(
        default=None, description="Use potentials (None = use Boltz default)"
    )
    use_affinity: bool = Field(default=True, description="Enable affinity prediction")
    use_msa_server: bool = Field(
        default=False, description="Use MSA server (vs pre-computed)"
    )


class RF3Config(BaseModel):
    """RF3-specific configuration."""

    seed: int = Field(default=42, description="Model seed")
    diffusion_batch_size: int = Field(
        default=5, ge=1, description="Number of structures per diffusion batch"
    )
    n_recycles: int = Field(default=10, ge=1, description="Number of recycles")
    num_steps: int = Field(default=50, ge=1, description="Diffusion steps")
    early_stopping_plddt_threshold: float = Field(
        default=0.5, ge=0.0, le=1.0, description="Early stopping threshold"
    )
    skip_existing: bool = Field(default=False, description="Skip existing outputs")


class InputsConfig(BaseModel):
    """Input file locations."""

    proteins_fasta: Path = Field(description="Path to proteins FASTA file")
    ligands_dir: Path = Field(description="Path to ligands directory (with subdirs)")
    msa_dir: Path | None = Field(
        default=None, description="Path to MSA directory with a3m files (optional)"
    )
    msa_sqsh_files: list[Path] = Field(
        default_factory=list,
        description="List of squashfs files containing a3m MSAs",
    )
    # For squashfs MSAs, specify mount point for jobs
    msa_sqsh_mount_point: Path = Field(
        default=Path("/msa"),
        description="Mount point for squashfs MSAs inside container/job",
    )
    oligo_definitions: Path | None = Field(
        default=None,
        description=(
            "Path to a YAML file defining oligosaccharide prefixes, "
            "their monomer CCD codes, and bond atom pairs. "
            "When provided, AF3 builds multi-monomer ccdCodes and "
            "bondedAtomPairs for matching ligands instead of using "
            "custom CIF files."
        ),
    )

    @field_validator("msa_sqsh_files", mode="before")
    @classmethod
    def expand_sqsh_paths(cls, v: Any) -> list[Path]:
        """Expand paths in squashfs list."""
        if v is None:
            return []
        if isinstance(v, list):
            return [Path(os.path.expandvars(str(p))) for p in v]
        return [Path(os.path.expandvars(str(v)))]

    @model_validator(mode="after")
    def resolve_relative_paths(self) -> "InputsConfig":
        """Convert relative paths to absolute."""
        if not self.proteins_fasta.is_absolute():
            object.__setattr__(self, 'proteins_fasta', self.proteins_fasta.resolve())
        if not self.ligands_dir.is_absolute():
            object.__setattr__(self, 'ligands_dir', self.ligands_dir.resolve())
        if self.msa_dir and not self.msa_dir.is_absolute():
            object.__setattr__(self, 'msa_dir', self.msa_dir.resolve())
        # Resolve squashfs paths
        resolved_sqsh = []
        for sqsh in self.msa_sqsh_files:
            if not sqsh.is_absolute():
                resolved_sqsh.append(sqsh.resolve())
            else:
                resolved_sqsh.append(sqsh)
        object.__setattr__(self, 'msa_sqsh_files', resolved_sqsh)
        return self


class OutputsConfig(BaseModel):
    """Output configuration."""

    work_dir: Path = Field(description="Working directory for outputs")
    manifest_dir: Path = Field(description="Directory for manifest files")
    results_dir: Path = Field(description="Directory for aggregated results")

    @model_validator(mode="after")
    def resolve_relative_paths(self) -> "OutputsConfig":
        """Convert relative paths to absolute."""
        if not self.work_dir.is_absolute():
            object.__setattr__(self, 'work_dir', self.work_dir.resolve())
        if not self.manifest_dir.is_absolute():
            object.__setattr__(self, 'manifest_dir', self.manifest_dir.resolve())
        if not self.results_dir.is_absolute():
            object.__setattr__(self, 'results_dir', self.results_dir.resolve())
        return self


class PipelineConfig(BaseModel):
    """Root pipeline configuration."""

    paths: PathsConfig
    slurm: SlurmConfig = Field(default_factory=SlurmConfig)
    af3: AF3Config = Field(default_factory=AF3Config)
    boltz: BoltzConfig = Field(default_factory=BoltzConfig)
    rf3: RF3Config = Field(default_factory=RF3Config)
    inputs: InputsConfig
    outputs: OutputsConfig

    # Global settings
    models_enabled: list[str] = Field(
        default=["af3", "boltz", "rf3"],
        description="Models to run (af3, boltz, rf3)",
    )
    dry_run: bool = Field(default=False, description="Dry run mode")
    verbose: bool = Field(default=False, description="Verbose logging")

    @field_validator("models_enabled", mode="after")
    @classmethod
    def validate_models(cls, v: list[str]) -> list[str]:
        """Validate model names."""
        valid = {"af3", "boltz", "rf3"}
        for model in v:
            if model not in valid:
                raise ValueError(f"Invalid model: {model}. Must be one of {valid}")
        return v

    @model_validator(mode="after")
    def create_output_dirs(self) -> "PipelineConfig":
        """Create output directories if they don't exist."""
        for path in [
            self.outputs.work_dir,
            self.outputs.manifest_dir,
            self.outputs.results_dir,
        ]:
            path.mkdir(parents=True, exist_ok=True)
        return self

    @classmethod
    def from_yaml(cls, path: Path) -> "PipelineConfig":
        """Load configuration from YAML file."""
        with open(path) as f:
            data = yaml.safe_load(f)
        return cls(**data)

    def to_yaml(self, path: Path) -> None:
        """Save configuration to YAML file."""
        # Convert to dict with string paths
        data = self.model_dump(mode="json")
        with open(path, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)


def load_config(config_path: Path | None = None) -> PipelineConfig:
    """Load configuration from file or environment.

    Args:
        config_path: Path to config file. If None, tries:
            1. STRUCTURE_PIPELINE_CONFIG env var
            2. ./config/pipeline.yaml
            3. ~/structure_pipeline.yaml

    Returns:
        Validated PipelineConfig
    """
    if config_path is None:
        # Try environment variable
        env_path = os.environ.get("STRUCTURE_PIPELINE_CONFIG")
        if env_path:
            config_path = Path(env_path)
        # Try current directory
        elif Path("config/pipeline.yaml").exists():
            config_path = Path("config/pipeline.yaml")
        # Try home directory
        elif Path.home().joinpath("structure_pipeline.yaml").exists():
            config_path = Path.home() / "structure_pipeline.yaml"
        else:
            raise FileNotFoundError(
                "No config file found. Provide --config or set STRUCTURE_PIPELINE_CONFIG"
            )

    return PipelineConfig.from_yaml(config_path)
