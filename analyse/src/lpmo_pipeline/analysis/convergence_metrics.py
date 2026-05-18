"""Convergence metrics for ligand reproducibility within one condition."""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from lpmo_pipeline.config import load_defaults_config, load_runtime_paths_config

_RUNTIME_PATHS = load_runtime_paths_config()
_DEFAULTS_CONFIG = load_defaults_config()
_CHAIN_SCHEMA = _DEFAULTS_CONFIG.get("chain_schema") or {}
_DEFAULT_THRESHOLDS_PATH = _RUNTIME_PATHS.pipeline_assets.thresholds_config
DEFAULT_GLYCAN_CHAINS = tuple(str(chain) for chain in (_CHAIN_SCHEMA.get("glycans") or ["B", "C", "D"]))
_DEFAULT_GLYCAN_SELECTION = " ".join(DEFAULT_GLYCAN_CHAINS) or "B"
DEFAULT_LIGAND_SELECTION = f"chainID {_DEFAULT_GLYCAN_SELECTION} and not name H*"

POSE_CONVERGENCE_COLUMNS = [
    "pose_id",
    "condition_id",
    "reference_pose_id",
    "ligand_rmsd_to_reference",
    "convergent_flag",
]

CONDITION_CONVERGENCE_COLUMNS = [
    "condition_id",
    "reference_pose_id",
    "n_qc_pass_poses",
    "convergence_fraction",
    "median_ligand_rmsd",
    "iqr_ligand_rmsd",
    "low_convergence_flag",
]


@dataclass(frozen=True)
class ConvergenceConfig:
    """Thresholds and selections for convergence metrics."""

    alignment_selection: str = "protein and name CA"
    ligand_selection: str = DEFAULT_LIGAND_SELECTION
    convergent_rmsd_max_a: float = 2.0
    low_convergence_flag_threshold: float = 0.3


@dataclass(frozen=True)
class ConvergencePoseInput:
    """Minimal per-pose inputs for pre-clustering convergence metrics."""

    pose_id: str
    condition_id: str
    complex_pdb: Path
    seed: int | None = None
    sample: int | None = None
    ranking_score: float | None = None
    mean_plddt: float | None = None


@dataclass(frozen=True)
class PoseConvergenceMetrics:
    """Per-pose convergence output."""

    pose_id: str
    condition_id: str
    reference_pose_id: str
    ligand_rmsd_to_reference: float
    convergent_flag: bool

    def to_row(self) -> dict[str, Any]:
        return {
            "pose_id": self.pose_id,
            "condition_id": self.condition_id,
            "reference_pose_id": self.reference_pose_id,
            "ligand_rmsd_to_reference": self.ligand_rmsd_to_reference,
            "convergent_flag": self.convergent_flag,
        }


@dataclass(frozen=True)
class ConditionConvergenceSummary:
    """Condition-level aggregate convergence summary."""

    condition_id: str
    reference_pose_id: str
    n_qc_pass_poses: int
    convergence_fraction: float
    median_ligand_rmsd: float
    iqr_ligand_rmsd: float
    low_convergence_flag: bool

    def to_row(self) -> dict[str, Any]:
        return {
            "condition_id": self.condition_id,
            "reference_pose_id": self.reference_pose_id,
            "n_qc_pass_poses": self.n_qc_pass_poses,
            "convergence_fraction": self.convergence_fraction,
            "median_ligand_rmsd": self.median_ligand_rmsd,
            "iqr_ligand_rmsd": self.iqr_ligand_rmsd,
            "low_convergence_flag": self.low_convergence_flag,
        }


@lru_cache(maxsize=4)
def load_convergence_config(config_path: str | None = None) -> ConvergenceConfig:
    """Load convergence settings from thresholds.yaml."""

    path = Path(config_path) if config_path else _DEFAULT_THRESHOLDS_PATH
    with open(path) as handle:
        payload = yaml.safe_load(handle) or {}

    convergence = payload.get("convergence") or {}
    return ConvergenceConfig(
        alignment_selection=str(convergence.get("alignment_selection") or "protein and name CA"),
        ligand_selection=str(convergence.get("ligand_selection") or DEFAULT_LIGAND_SELECTION),
        convergent_rmsd_max_a=float(convergence.get("convergent_rmsd_max_a", 2.0)),
        low_convergence_flag_threshold=float(convergence.get("low_convergence_flag_threshold", 0.3)),
    )


def select_reference_pose(pose_inputs: list[ConvergencePoseInput]) -> ConvergencePoseInput:
    """Select the pre-clustering reference pose.

    Prefer QC-passing seed-1 poses when present. Within that subset, prefer
    higher ranking score or mean pLDDT when available, otherwise the lowest
    sample number for deterministic behavior.
    """

    if not pose_inputs:
        raise ValueError("select_reference_pose requires at least one pose")

    seed_one_inputs = [pose_input for pose_input in pose_inputs if pose_input.seed == 1]
    candidates = seed_one_inputs or list(pose_inputs)
    return min(candidates, key=_reference_sort_key)


def _reference_sort_key(pose_input: ConvergencePoseInput) -> tuple[float, float, int, int, str]:
    ranking_sort = -float(pose_input.ranking_score) if pose_input.ranking_score is not None else math.inf
    plddt_sort = -float(pose_input.mean_plddt) if pose_input.mean_plddt is not None else math.inf
    seed_sort = pose_input.seed if pose_input.seed is not None else math.inf
    sample_sort = pose_input.sample if pose_input.sample is not None else math.inf
    return (ranking_sort, plddt_sort, seed_sort, sample_sort, pose_input.pose_id)


def compute_condition_convergence(
    pose_inputs: list[ConvergencePoseInput],
    *,
    config: ConvergenceConfig | None = None,
) -> tuple[list[PoseConvergenceMetrics], ConditionConvergenceSummary]:
    """Compute pre-clustering convergence metrics for one condition."""

    if not pose_inputs:
        raise ValueError("compute_condition_convergence requires at least one pose")

    resolved_config = config or load_convergence_config()
    condition_ids = {pose_input.condition_id for pose_input in pose_inputs}
    if len(condition_ids) != 1:
        raise ValueError("compute_condition_convergence expects exactly one condition_id")

    reference_input = select_reference_pose(pose_inputs)
    reference_universe = _load_universe(reference_input.complex_pdb)
    reference_ligand = _select_ligand_atoms(reference_universe, resolved_config.ligand_selection)

    results: list[PoseConvergenceMetrics] = []
    for pose_input in pose_inputs:
        mobile_universe = _load_universe(pose_input.complex_pdb)
        _align_protein_to_reference(
            mobile_universe,
            reference_universe,
            resolved_config.alignment_selection,
        )
        mobile_ligand = _select_ligand_atoms(mobile_universe, resolved_config.ligand_selection)
        ligand_rmsd = _compute_atomwise_rmsd(reference_ligand.positions, mobile_ligand.positions)
        results.append(
            PoseConvergenceMetrics(
                pose_id=pose_input.pose_id,
                condition_id=pose_input.condition_id,
                reference_pose_id=reference_input.pose_id,
                ligand_rmsd_to_reference=ligand_rmsd,
                convergent_flag=ligand_rmsd < resolved_config.convergent_rmsd_max_a,
            )
        )

    summary = build_condition_convergence_summary(results, resolved_config)
    return results, summary


def build_condition_convergence_summary(
    pose_metrics: list[PoseConvergenceMetrics],
    config: ConvergenceConfig | None = None,
) -> ConditionConvergenceSummary:
    """Aggregate per-pose convergence metrics to one condition row."""

    if not pose_metrics:
        raise ValueError("build_condition_convergence_summary requires at least one pose metric")

    resolved_config = config or load_convergence_config()
    condition_ids = {metric.condition_id for metric in pose_metrics}
    if len(condition_ids) != 1:
        raise ValueError("build_condition_convergence_summary expects exactly one condition_id")

    reference_pose_ids = {metric.reference_pose_id for metric in pose_metrics}
    if len(reference_pose_ids) != 1:
        raise ValueError("build_condition_convergence_summary expects exactly one reference pose")

    rmsd_values = np.asarray([metric.ligand_rmsd_to_reference for metric in pose_metrics], dtype=float)
    convergence_fraction = float(np.mean([metric.convergent_flag for metric in pose_metrics]))
    iqr = float(np.percentile(rmsd_values, 75) - np.percentile(rmsd_values, 25))

    return ConditionConvergenceSummary(
        condition_id=next(iter(condition_ids)),
        reference_pose_id=next(iter(reference_pose_ids)),
        n_qc_pass_poses=len(pose_metrics),
        convergence_fraction=convergence_fraction,
        median_ligand_rmsd=float(np.median(rmsd_values)),
        iqr_ligand_rmsd=iqr,
        low_convergence_flag=convergence_fraction < resolved_config.low_convergence_flag_threshold,
    )


def write_pose_convergence_tsv(
    pose_metrics: list[PoseConvergenceMetrics],
    output_path: Path,
) -> None:
    """Write per-pose convergence metrics to TSV."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=POSE_CONVERGENCE_COLUMNS, delimiter="\t")
        writer.writeheader()
        for metric in pose_metrics:
            writer.writerow(metric.to_row())


def write_condition_convergence_summary_tsv(
    summary_rows: list[ConditionConvergenceSummary],
    output_path: Path,
) -> None:
    """Write per-condition convergence summaries to TSV."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CONDITION_CONVERGENCE_COLUMNS, delimiter="\t")
        writer.writeheader()
        for summary in summary_rows:
            writer.writerow(summary.to_row())


def _load_universe(complex_pdb: Path) -> Any:
    import MDAnalysis as mda

    if not complex_pdb.exists():
        raise FileNotFoundError(f"Missing convergence input PDB: {complex_pdb}")
    return mda.Universe(str(complex_pdb))


def _align_protein_to_reference(mobile_universe: Any, reference_universe: Any, alignment_selection: str) -> None:
    from MDAnalysis.analysis import align

    mobile_atoms = mobile_universe.select_atoms(alignment_selection)
    reference_atoms = reference_universe.select_atoms(alignment_selection)
    if mobile_atoms.n_atoms == 0 or reference_atoms.n_atoms == 0:
        raise ValueError("Alignment selection produced zero atoms for convergence metrics")
    if mobile_atoms.n_atoms != reference_atoms.n_atoms:
        raise ValueError(
            "Alignment atom count mismatch for convergence metrics: "
            f"{mobile_atoms.n_atoms} != {reference_atoms.n_atoms}"
        )
    align.alignto(mobile_universe, reference_universe, select=alignment_selection)


def _select_ligand_atoms(universe: Any, ligand_selection: str) -> Any:
    ligand_atoms = universe.select_atoms(ligand_selection)
    if ligand_atoms.n_atoms == 0:
        raise ValueError("Ligand selection produced zero atoms for convergence metrics")
    return ligand_atoms


def _compute_atomwise_rmsd(reference_positions: np.ndarray, mobile_positions: np.ndarray) -> float:
    if reference_positions.shape != mobile_positions.shape:
        raise ValueError(
            "Ligand atom count mismatch for convergence metrics: "
            f"{reference_positions.shape} != {mobile_positions.shape}"
        )
    squared = np.sum((mobile_positions - reference_positions) ** 2, axis=1)
    return float(np.sqrt(np.mean(squared)))