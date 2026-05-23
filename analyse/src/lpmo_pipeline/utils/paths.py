# src/lpmo_pipeline/utils/paths.py
"""
Responsibility: Canonical path construction for all pipeline artifacts.
Input:  run_id components (model, protein_id, ligand_id, seed, etc.)
Output: Path objects following the directory convention
"""
from __future__ import annotations

from pathlib import Path


# --------------------------------------------------------------------------- #
# Base layout
# --------------------------------------------------------------------------- #
def results_root(base: Path, mode: str = "production") -> Path:
    """Top-level results directory.

    Args:
        base: Project root (e.g. /cluster/work/projects/.../analyse)
        mode: "tuning" or "production"

    Returns:
        Path to results/{tuning_*|del_a|del_b}
    """
    return base / "results" if mode == "production" else base / "results"


# --------------------------------------------------------------------------- #
# Per-run paths
# --------------------------------------------------------------------------- #
def run_dir(
    base: Path,
    model: str,
    protein_id: str,
    ligand_id: str,
    seed: int,
    *,
    mode: str = "production",
    del_variant: str = "a",
) -> Path:
    """Directory for a single prediction run.

    Convention:
        results/del_{a|b}/{model}/{protein_id}/{ligand_id}/seed_{seed}/
    """
    root = results_root(base, mode)
    return root / f"del_{del_variant}" / model / protein_id / ligand_id / f"seed_{seed}"


def normalized_cif_path(run: Path) -> Path:
    return run / "normalized.cif"


def atom_map_tsv_path(run: Path) -> Path:
    return run / "atom_map.tsv"


def rename_log_path(run: Path) -> Path:
    return run / "rename_log.json"


def privateer_input_path(run: Path) -> Path:
    return run / "privateer_input.cif"


def protonation_dir(run: Path) -> Path:
    return run / "protonated"


def analysis_export_dir(run: Path) -> Path:
    return run / "analysis_export"


def posebusters_pdb_path(run: Path) -> Path:
    return protonation_dir(run) / "for_posebusters.pdb"


def prolif_ligand_pdb_path(run: Path) -> Path:
    return analysis_export_dir(run) / "ligand_only_for_prolif.pdb"


def prolif_complex_pdb_path(run: Path) -> Path:
    return analysis_export_dir(run) / "complex_for_prolif.pdb"


def complex_h_pdb_path(run: Path) -> Path:
    return protonation_dir(run) / "complex_H.pdb"


# --------------------------------------------------------------------------- #
# PLACER paths
# --------------------------------------------------------------------------- #
def placer_ensemble_dir(run: Path) -> Path:
    return run / "placer_ensemble"


def placer_scores_path(run: Path) -> Path:
    return placer_ensemble_dir(run) / "placer_scores.json"


# --------------------------------------------------------------------------- #
# QC paths
# --------------------------------------------------------------------------- #
def qc_report_path(run: Path) -> Path:
    return run / "qc_report.json"


# --------------------------------------------------------------------------- #
# Analysis paths
# --------------------------------------------------------------------------- #
def ifp_matrix_path(run: Path) -> Path:
    return run / "analysis" / "ifp_matrix.csv"


def geometry_metrics_path(run: Path) -> Path:
    return run / "analysis" / "geometry_metrics.json"


# --------------------------------------------------------------------------- #
# Cluster paths (per protein×ligand, not per-seed)
# --------------------------------------------------------------------------- #
def cluster_dir(base: Path, protein_id: str, ligand_id: str, del_variant: str = "a") -> Path:
    return results_root(base) / f"del_{del_variant}" / "clusters" / f"{protein_id}_{ligand_id}"


def cluster_signatures_path(cluster_d: Path) -> Path:
    return cluster_d / "cluster_signatures.json"


# --------------------------------------------------------------------------- #
# Report paths
# --------------------------------------------------------------------------- #
def report_dir(base: Path, del_variant: str = "a") -> Path:
    return results_root(base) / f"del_{del_variant}"


def summary_json_path(base: Path, del_variant: str = "a") -> Path:
    return report_dir(base, del_variant) / "summary.json"


def metrics_csv_path(base: Path, del_variant: str = "a") -> Path:
    return report_dir(base, del_variant) / "metrics.csv"


def report_html_path(base: Path, del_variant: str = "a") -> Path:
    return report_dir(base, del_variant) / "report.html"


# --------------------------------------------------------------------------- #
# Tuning paths
# --------------------------------------------------------------------------- #
def tuning_dir(base: Path, model: str) -> Path:
    return results_root(base) / f"tuning_{model.lower()}_results"


def best_params_path(base: Path, model: str) -> Path:
    return tuning_dir(base, model) / "best_params.yaml"


# --------------------------------------------------------------------------- #
# Manifest
# --------------------------------------------------------------------------- #
def manifest_path(base: Path) -> Path:
    return results_root(base) / "reproducibility" / "run_manifest.json"


# --------------------------------------------------------------------------- #
# Work-root (structure_pipeline prediction artifacts)
# --------------------------------------------------------------------------- #
def work_root(structure_pipeline_base: Path) -> Path:
    """Path to the ``work/`` directory inside a structure_pipeline checkout."""
    return structure_pipeline_base / "work"


def discovery_manifest_path(base: Path) -> Path:
    """Where to write the JSON discovery manifest."""
    return results_root(base) / "discovery_manifest.json"
