"""Base runner interface for structure prediction models."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..cases import Case
from ..config import PipelineConfig


@dataclass
class RunnerResult:
    """Result from running a prediction job."""

    success: bool
    case_id: str
    output_dir: Path | None = None
    output_files: list[str] = field(default_factory=list)
    confidence_scores: dict[str, float] = field(default_factory=dict)
    error_message: str = ""
    runtime_seconds: float = 0.0
    metadata: dict[str, Any] = field(default_factory=dict)


class RunnerInterface(ABC):
    """Abstract base class for model runners.

    Each runner is responsible for:
    1. Building model-specific input files (JSON/YAML)
    2. Generating the command to run the model
    3. Parsing output and confidence scores
    """

    def __init__(self, config: PipelineConfig):
        """Initialize runner with pipeline configuration.

        Args:
            config: Pipeline configuration
        """
        self.config = config

    @property
    @abstractmethod
    def model_name(self) -> str:
        """Return the model name (af3, boltz, rf3)."""
        pass

    @abstractmethod
    def build_input_file(
        self,
        cases: list[Case],
        protein_sequences: dict[str, str],
        output_dir: Path,
    ) -> Path:
        """Build model-specific input file for a batch of cases.

        All cases should be for the same ligand but different proteins.

        Args:
            cases: List of cases (same ligand, different proteins)
            protein_sequences: Mapping of protein_id -> amino acid sequence
            output_dir: Directory to write input file

        Returns:
            Path to the generated input file
        """
        pass

    @abstractmethod
    def build_command(
        self,
        input_file: Path,
        output_dir: Path,
        **kwargs,
    ) -> list[str]:
        """Build the command to run the model.

        Args:
            input_file: Path to input file
            output_dir: Directory for outputs
            **kwargs: Additional model-specific arguments

        Returns:
            Command as list of strings
        """
        pass

    @abstractmethod
    def build_slurm_script(
        self,
        cases: list[Case],
        protein_sequences: dict[str, str],
        work_dir: Path,
        job_name: str,
    ) -> str:
        """Build complete SLURM job script.

        Args:
            cases: List of cases (same ligand, different proteins)
            protein_sequences: Mapping of protein_id -> amino acid sequence
            work_dir: Working directory for job
            job_name: SLURM job name

        Returns:
            Complete SLURM script as string
        """
        pass

    @abstractmethod
    def parse_outputs(
        self,
        output_dir: Path,
        cases: list[Case],
    ) -> list[RunnerResult]:
        """Parse model outputs and extract results.

        Args:
            output_dir: Directory containing outputs
            cases: List of cases that were run

        Returns:
            List of RunnerResult objects
        """
        pass

    def validate_paths(self) -> list[str]:
        """Validate that required paths exist.

        Returns:
            List of error messages (empty if all valid)
        """
        errors = []
        # Subclasses should override to check specific paths
        return errors
