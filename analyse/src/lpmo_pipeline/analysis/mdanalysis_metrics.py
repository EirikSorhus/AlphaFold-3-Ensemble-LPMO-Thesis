"""Downstream geometry metrics for pose_geometry.tsv.

This analysis-stage branch reuses the same Cu / histidine-brace identification
rules as hard QC, but computes geometry-only features that are intentionally
kept out of the IFP clustering input.
"""
from __future__ import annotations

import csv
import json
import logging
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from lpmo_pipeline.config import load_defaults_config, load_runtime_paths_config
from lpmo_pipeline.qc.custom_geometry_checks import check_geometry

logger = logging.getLogger(__name__)

_RUNTIME_PATHS = load_runtime_paths_config()
_DEFAULTS_CONFIG = load_defaults_config()
_CHAIN_SCHEMA = _DEFAULTS_CONFIG.get("chain_schema") or {}
DEFAULT_PROTEIN_CHAIN = str(_CHAIN_SCHEMA.get("protein") or "A")
DEFAULT_GLYCAN_CHAINS = tuple(str(chain) for chain in (_CHAIN_SCHEMA.get("glycans") or ["B", "C", "D"]))
DEFAULT_CU_CHAIN = str(_CHAIN_SCHEMA.get("metal") or "E")

GEOMETRY_STATUS_NOT_COMPUTABLE = "geometry_not_computable"
GEOMETRY_STATUS_COMPUTABLE_IMPLAUSIBLE = "geometry_computable_implausible"
GEOMETRY_STATUS_PLAUSIBLE = "geometry_plausible"
GEOMETRY_STATUS_HIGHLY_PLAUSIBLE = "geometry_highly_plausible"

POSE_GEOMETRY_COLUMNS = [
    "pose_id",
    "model",
    "protein_id",
    "ligand_id",
    "proximal_sugar_id",
    "brace_integrity_flag",
    "Cu_C1_distance",
    "Cu_C4_distance",
    "oxyl_H_C1_distance",
    "oxyl_H_C4_distance",
    "Cu_oxyl_H_C1_angle",
    "Cu_oxyl_H_C4_angle",
    "sugar_face_orientation",
    "ring_normal_vs_brace_normal",
    "oxyl_H_score_C1",
    "oxyl_H_score_C4",
    "geometry_status_C1",
    "geometry_status_C4",
    "his_brace_angle_deg",
    "core_rmsd_vs_reference",
    "pocket_rmsd_vs_crystal",
]


@dataclass(frozen=True)
class GeometryThresholds:
    """Locked thresholds for the downstream geometry branch."""

    oxyl_h_min_a: float = 1.5
    oxyl_h_max_a: float = 4.0
    oxyl_h_optimum_a: float = 2.1
    oxyl_h_score_normalization: float = 1.9
    oxyl_h_highly_plausible_min_a: float = 1.8
    oxyl_h_highly_plausible_max_a: float = 2.5
    his_brace_max_search_a: float = 3.0
    cu_oxyl_bond_length_a: float = 1.9
    c_h_bond_length_a: float = 1.1


@dataclass
class PoseGeometryMetrics:
    """All geometry outputs for one pose-level row in pose_geometry.tsv."""

    pose_id: str = ""
    model: str = ""
    protein_id: str = ""
    ligand_id: str = ""

    proximal_sugar_id: str = ""
    brace_integrity_flag: bool = False

    cu_c1_distance: float | None = None
    cu_c4_distance: float | None = None
    oxyl_h_c1_distance: float | None = None
    oxyl_h_c4_distance: float | None = None
    cu_oxyl_h_c1_angle: float | None = None
    cu_oxyl_h_c4_angle: float | None = None
    sugar_face_orientation: str = GEOMETRY_STATUS_NOT_COMPUTABLE
    ring_normal_vs_brace_normal: float | None = None
    oxyl_h_score_c1: float = math.nan
    oxyl_h_score_c4: float = math.nan
    geometry_status_c1: str = GEOMETRY_STATUS_NOT_COMPUTABLE
    geometry_status_c4: str = GEOMETRY_STATUS_NOT_COMPUTABLE

    his_brace_angle_deg: float | None = None
    core_rmsd_vs_reference: float | None = None
    pocket_rmsd_vs_crystal: float | None = None

    original_cu_coordinates: tuple[float, float, float] | None = None
    repositioned_cu_coordinates: tuple[float, float, float] | None = None
    virtual_oxyl_coordinates: tuple[float, float, float] | None = None
    virtual_h_c1_coordinates: tuple[float, float, float] | None = None
    virtual_h_c4_coordinates: tuple[float, float, float] | None = None

    @property
    def min_cu_c1(self) -> float | None:
        return self.cu_c1_distance

    @property
    def min_cu_c4(self) -> float | None:
        return self.cu_c4_distance

    def to_row(self) -> dict[str, Any]:
        """Flatten to the pose_geometry.tsv row contract."""
        return {
            "pose_id": self.pose_id,
            "model": self.model,
            "protein_id": self.protein_id,
            "ligand_id": self.ligand_id,
            "proximal_sugar_id": self.proximal_sugar_id,
            "brace_integrity_flag": self.brace_integrity_flag,
            "Cu_C1_distance": self.cu_c1_distance,
            "Cu_C4_distance": self.cu_c4_distance,
            "oxyl_H_C1_distance": self.oxyl_h_c1_distance,
            "oxyl_H_C4_distance": self.oxyl_h_c4_distance,
            "Cu_oxyl_H_C1_angle": self.cu_oxyl_h_c1_angle,
            "Cu_oxyl_H_C4_angle": self.cu_oxyl_h_c4_angle,
            "sugar_face_orientation": self.sugar_face_orientation,
            "ring_normal_vs_brace_normal": self.ring_normal_vs_brace_normal,
            "oxyl_H_score_C1": self.oxyl_h_score_c1,
            "oxyl_H_score_C4": self.oxyl_h_score_c4,
            "geometry_status_C1": self.geometry_status_c1,
            "geometry_status_C4": self.geometry_status_c4,
            "his_brace_angle_deg": self.his_brace_angle_deg,
            "core_rmsd_vs_reference": self.core_rmsd_vs_reference,
            "pocket_rmsd_vs_crystal": self.pocket_rmsd_vs_crystal,
        }

    def to_legacy_geometry_dict(self) -> dict[str, Any]:
        """Provide the minimal legacy keys still expected by downstream code."""
        return {
            "min_cu_c1": self.cu_c1_distance,
            "min_cu_c4": self.cu_c4_distance,
            "his_brace_angle_deg": self.his_brace_angle_deg,
            "core_rmsd_vs_reference": self.core_rmsd_vs_reference,
            "pocket_rmsd_vs_crystal": self.pocket_rmsd_vs_crystal,
            "geometry_status_C1": self.geometry_status_c1,
            "geometry_status_C4": self.geometry_status_c4,
            "oxyl_H_score_C1": self.oxyl_h_score_c1,
            "oxyl_H_score_C4": self.oxyl_h_score_c4,
        }

    def to_debug_dict(self) -> dict[str, Any]:
        """Extend the row payload with visual-debug coordinates."""
        data = self.to_row()
        data.update(
            {
                "original_cu_coordinates": self.original_cu_coordinates,
                "repositioned_cu_coordinates": self.repositioned_cu_coordinates,
                "virtual_oxyl_coordinates": self.virtual_oxyl_coordinates,
                "virtual_h_c1_coordinates": self.virtual_h_c1_coordinates,
                "virtual_h_c4_coordinates": self.virtual_h_c4_coordinates,
            }
        )
        return data


@dataclass(frozen=True)
class _ProximalSugar:
    chain_name: str
    residue: Any
    min_heavy_atom_distance: float


@dataclass(frozen=True)
class _TargetGeometry:
    cu_c_distance: float | None
    oxyl_h_distance: float | None
    cu_oxyl_h_angle: float | None
    score: float
    status: str
    virtual_h_coordinates: tuple[float, float, float] | None = None


def _default_thresholds_path() -> Path:
    return _RUNTIME_PATHS.pipeline_assets.thresholds_config


@lru_cache(maxsize=4)
def load_geometry_thresholds(config_path: str | None = None) -> GeometryThresholds:
    """Load the locked downstream geometry thresholds from YAML."""
    path = Path(config_path) if config_path else _default_thresholds_path()
    with open(path) as f:
        data = yaml.safe_load(f) or {}

    geometry = data.get("geometry_plausibility") or {}
    rules = data.get("geometry_rules") or {}
    his_brace = rules.get("his_brace") or {}
    virtual_oxyl = rules.get("virtual_oxyl") or {}
    virtual_h = rules.get("virtual_h") or {}

    return GeometryThresholds(
        oxyl_h_min_a=float(geometry.get("oxyl_h_min_a", 1.5)),
        oxyl_h_max_a=float(geometry.get("oxyl_h_max_a", 4.0)),
        oxyl_h_optimum_a=float(geometry.get("oxyl_h_optimum_a", 2.1)),
        oxyl_h_score_normalization=float(geometry.get("oxyl_h_score_normalization", 1.9)),
        oxyl_h_highly_plausible_min_a=float(
            geometry.get("oxyl_h_highly_plausible_min_a", 1.8)
        ),
        oxyl_h_highly_plausible_max_a=float(
            geometry.get("oxyl_h_highly_plausible_max_a", 2.5)
        ),
        his_brace_max_search_a=float(his_brace.get("max_search_dist_a", 3.0)),
        cu_oxyl_bond_length_a=float(virtual_oxyl.get("cu_oxyl_bond_length_a", 1.9)),
        c_h_bond_length_a=float(virtual_h.get("c_h_bond_length_a", 1.1)),
    )


def compute_pose_metrics(
    pdb_path: Path,
    pose_id: str = "",
    protein_chain: str = DEFAULT_PROTEIN_CHAIN,
    glycan_chains: list[str] | None = None,
    cu_chain: str = DEFAULT_CU_CHAIN,
    reference_pdb: Path | None = None,
    crystal_pdb: Path | None = None,
    core_selection: str = "protein and not (resid 1:5 or name H*)",
    pocket_residues: list[int] | None = None,
    model: str = "",
    protein_id: str = "",
    ligand_id: str = "",
    thresholds: GeometryThresholds | None = None,
) -> PoseGeometryMetrics:
    """Compute a downstream geometry row from a structure file.

    This function owns the pose-local virtual oxyl/H geometry fields in
    ``pose_geometry.tsv``. Alignment-derived RMSD values are populated by the
    convergence and crystal-anchoring stages when those comparisons are
    available; they are intentionally left ``None`` here.
    """
    try:
        import gemmi
    except ImportError as exc:
        raise ImportError("gemmi is required to compute downstream geometry metrics") from exc

    structure = gemmi.read_structure(str(pdb_path))
    return compute_pose_metrics_from_structure(
        structure=structure,
        pose_id=pose_id,
        protein_chain=protein_chain,
        glycan_chains=glycan_chains,
        cu_chain=cu_chain,
        reference_pdb=reference_pdb,
        crystal_pdb=crystal_pdb,
        core_selection=core_selection,
        pocket_residues=pocket_residues,
        model=model,
        protein_id=protein_id,
        ligand_id=ligand_id,
        thresholds=thresholds,
    )


def compute_pose_metrics_from_structure(
    structure: Any,
    pose_id: str = "",
    protein_chain: str = DEFAULT_PROTEIN_CHAIN,
    glycan_chains: list[str] | None = None,
    cu_chain: str = DEFAULT_CU_CHAIN,
    reference_pdb: Path | None = None,
    crystal_pdb: Path | None = None,
    core_selection: str = "protein and not (resid 1:5 or name H*)",
    pocket_residues: list[int] | None = None,
    model: str = "",
    protein_id: str = "",
    ligand_id: str = "",
    thresholds: GeometryThresholds | None = None,
) -> PoseGeometryMetrics:
    """Compute one pose_geometry row from an in-memory structure."""
    del protein_chain, reference_pdb, crystal_pdb, core_selection, pocket_residues

    glycan_chains = list(glycan_chains or DEFAULT_GLYCAN_CHAINS)
    thresholds = thresholds or load_geometry_thresholds()

    metrics = PoseGeometryMetrics(
        pose_id=pose_id,
        model=model,
        protein_id=protein_id,
        ligand_id=ligand_id,
    )
    logger.info("Computing downstream geometry metrics for pose %s", pose_id)

    geom_result = check_geometry(
        structure=structure,
        pose_id=pose_id,
        cu_chain=cu_chain,
        glycan_chains=glycan_chains,
        his_brace_max_search_a=thresholds.his_brace_max_search_a,
    )
    metrics.his_brace_angle_deg = geom_result.his_brace_angle

    if not geom_result.cu_found:
        return metrics

    original_cu = np.asarray(geom_result.cu_position, dtype=float)
    metrics.original_cu_coordinates = tuple(original_cu.tolist())
    selected_brace = _selected_brace_positions(geom_result)
    metrics.brace_integrity_flag = (
        len(selected_brace) == 3
        and all(
            m.distance_angstrom <= thresholds.his_brace_max_search_a
            for m in geom_result.cu_his_measurements
        )
    )

    proximal_sugar = _select_proximal_sugar(structure, original_cu, glycan_chains)
    if proximal_sugar is None:
        return metrics

    metrics.proximal_sugar_id = _format_residue_id(
        proximal_sugar.chain_name,
        proximal_sugar.residue,
    )

    brace_normal = _compute_plane_normal(selected_brace) if len(selected_brace) == 3 else None
    repositioned_cu = _compute_centroid(selected_brace) if len(selected_brace) == 3 else None
    if repositioned_cu is not None:
        metrics.repositioned_cu_coordinates = tuple(repositioned_cu.tolist())
    nterminal_n = _find_nterminal_brace_nitrogen(geom_result)
    virtual_oxyl = _place_virtual_oxyl(
        repositioned_cu=repositioned_cu,
        nterminal_n=nterminal_n,
        brace_normal=brace_normal,
        bond_length=thresholds.cu_oxyl_bond_length_a,
    )
    if virtual_oxyl is not None:
        metrics.virtual_oxyl_coordinates = tuple(virtual_oxyl.tolist())

    ring_normal, ring_centroid = _compute_ring_frame(proximal_sugar.residue)
    if ring_normal is not None and brace_normal is not None:
        metrics.ring_normal_vs_brace_normal = _angle_between_vectors_deg(
            ring_normal,
            brace_normal,
        )
    metrics.sugar_face_orientation = _classify_sugar_face_orientation(
        ring_normal=ring_normal,
        ring_centroid=ring_centroid,
        repositioned_cu=repositioned_cu,
    )

    c1_atom = _get_atom_by_name(proximal_sugar.residue, "C1")
    c4_atom = _get_atom_by_name(proximal_sugar.residue, "C4")

    c1_metrics = _compute_target_geometry(
        structure=structure,
        chain_name=proximal_sugar.chain_name,
        residue=proximal_sugar.residue,
        target_atom=c1_atom,
        target_name="C1",
        repositioned_cu=repositioned_cu,
        virtual_oxyl=virtual_oxyl,
        thresholds=thresholds,
    )
    c4_metrics = _compute_target_geometry(
        structure=structure,
        chain_name=proximal_sugar.chain_name,
        residue=proximal_sugar.residue,
        target_atom=c4_atom,
        target_name="C4",
        repositioned_cu=repositioned_cu,
        virtual_oxyl=virtual_oxyl,
        thresholds=thresholds,
    )

    metrics.cu_c1_distance = c1_metrics.cu_c_distance
    metrics.cu_c4_distance = c4_metrics.cu_c_distance
    metrics.oxyl_h_c1_distance = c1_metrics.oxyl_h_distance
    metrics.oxyl_h_c4_distance = c4_metrics.oxyl_h_distance
    metrics.cu_oxyl_h_c1_angle = c1_metrics.cu_oxyl_h_angle
    metrics.cu_oxyl_h_c4_angle = c4_metrics.cu_oxyl_h_angle
    metrics.oxyl_h_score_c1 = c1_metrics.score
    metrics.oxyl_h_score_c4 = c4_metrics.score
    metrics.geometry_status_c1 = c1_metrics.status
    metrics.geometry_status_c4 = c4_metrics.status
    metrics.virtual_h_c1_coordinates = c1_metrics.virtual_h_coordinates
    metrics.virtual_h_c4_coordinates = c4_metrics.virtual_h_coordinates
    return metrics


def write_pose_geometry_tsv(metrics_rows: list[PoseGeometryMetrics], output_path: Path) -> None:
    """Write pose_geometry.tsv from computed pose rows."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=POSE_GEOMETRY_COLUMNS, delimiter="\t")
        writer.writeheader()
        for metrics in metrics_rows:
            writer.writerow(metrics.to_row())

    logger.info("Wrote pose_geometry.tsv (%d rows) to %s", len(metrics_rows), output_path)


def write_geometry_metrics(metrics: PoseGeometryMetrics, output_path: Path) -> None:
    """Write a single-pose debug JSON representation of the geometry row."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(metrics.to_debug_dict(), f, indent=2, allow_nan=True)
    logger.info("Wrote geometry metrics debug JSON to %s", output_path)


def write_geometry_debug_pdb(
    structure: Any,
    metrics: PoseGeometryMetrics,
    output_path: Path,
) -> None:
    """Write a visual debug PDB with virtual geometry atoms appended.

    This is intentionally test-only output for manual inspection.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    structure.write_pdb(str(output_path))
    output_path.write_text(
        _append_virtual_geometry_atoms_to_pdb_text(output_path.read_text(), metrics)
    )
    logger.info("Wrote geometry debug PDB to %s", output_path)


def _selected_brace_positions(geom_result: Any) -> list[np.ndarray]:
    positions: list[np.ndarray] = []
    for measurement in geom_result.cu_his_measurements:
        positions.append(np.asarray(measurement.his_atom_position, dtype=float))
    return positions


def _find_nterminal_brace_nitrogen(geom_result: Any) -> np.ndarray | None:
    for measurement in geom_result.cu_his_measurements:
        if measurement.his_resnum == 1 and measurement.his_atom == "N":
            return np.asarray(measurement.his_atom_position, dtype=float)
    return None


def _select_proximal_sugar(
    structure: Any,
    cu_pos: np.ndarray,
    glycan_chains: list[str],
) -> _ProximalSugar | None:
    best: _ProximalSugar | None = None
    for chain_name in glycan_chains:
        chain = _get_chain(structure, chain_name)
        if chain is None:
            continue
        for residue in chain:
            heavy_atoms = [atom for atom in residue if not _is_hydrogen(atom)]
            if not heavy_atoms:
                continue
            min_distance = min(
                float(np.linalg.norm(_atom_position(atom) - cu_pos))
                for atom in heavy_atoms
            )
            candidate = _ProximalSugar(
                chain_name=chain_name,
                residue=residue,
                min_heavy_atom_distance=min_distance,
            )
            if best is None or candidate.min_heavy_atom_distance < best.min_heavy_atom_distance:
                best = candidate
    return best


def _compute_target_geometry(
    structure: Any,
    chain_name: str,
    residue: Any,
    target_atom: Any | None,
    target_name: str,
    repositioned_cu: np.ndarray | None,
    virtual_oxyl: np.ndarray | None,
    thresholds: GeometryThresholds,
) -> _TargetGeometry:
    if target_atom is None or repositioned_cu is None or virtual_oxyl is None:
        return _TargetGeometry(
            None,
            None,
            None,
            math.nan,
            GEOMETRY_STATUS_NOT_COMPUTABLE,
            None,
        )

    target_pos = _atom_position(target_atom)
    virtual_h = _place_virtual_h(
        structure=structure,
        chain_name=chain_name,
        residue=residue,
        target_atom=target_atom,
        bond_length=thresholds.c_h_bond_length_a,
    )
    if virtual_h is None:
        logger.debug("Virtual H placement failed for %s", target_name)
        return _TargetGeometry(
            None,
            None,
            None,
            math.nan,
            GEOMETRY_STATUS_NOT_COMPUTABLE,
            None,
        )

    cu_c_distance = float(np.linalg.norm(repositioned_cu - target_pos))
    oxyl_h_distance = float(np.linalg.norm(virtual_oxyl - virtual_h))
    cu_oxyl_h_angle = _angle_from_points_deg(repositioned_cu, virtual_oxyl, virtual_h)
    score = _score_oxyl_h_distance(oxyl_h_distance, thresholds)
    status = _assign_geometry_status(oxyl_h_distance, score, thresholds)
    return _TargetGeometry(
        cu_c_distance,
        oxyl_h_distance,
        cu_oxyl_h_angle,
        score,
        status,
        tuple(virtual_h.tolist()),
    )


def _place_virtual_oxyl(
    repositioned_cu: np.ndarray | None,
    nterminal_n: np.ndarray | None,
    brace_normal: np.ndarray | None,
    bond_length: float,
) -> np.ndarray | None:
    if repositioned_cu is None or nterminal_n is None or brace_normal is None:
        return None

    n_to_cu = repositioned_cu - nterminal_n
    planar_direction = n_to_cu - np.dot(n_to_cu, brace_normal) * brace_normal
    unit_direction = _unit_vector(planar_direction)
    if unit_direction is None:
        return None
    return repositioned_cu + unit_direction * bond_length


def _place_virtual_h(
    structure: Any,
    chain_name: str,
    residue: Any,
    target_atom: Any,
    bond_length: float,
) -> np.ndarray | None:
    target_pos = _atom_position(target_atom)
    neighbors = _infer_bonded_heavy_neighbors(
        structure=structure,
        chain_name=chain_name,
        residue=residue,
        target_atom=target_atom,
    )
    if not neighbors:
        return None

    centroid = _compute_centroid([_atom_position(atom) for atom in neighbors])
    direction = _unit_vector(target_pos - centroid)
    if direction is None:
        return None
    return target_pos + direction * bond_length


def _infer_bonded_heavy_neighbors(
    structure: Any,
    chain_name: str,
    residue: Any,
    target_atom: Any,
) -> list[Any]:
    neighbors: dict[int, tuple[float, Any]] = {}
    target_pos = _atom_position(target_atom)

    def _add_neighbor(atom: Any) -> None:
        if atom is target_atom or _is_hydrogen(atom):
            return
        distance = float(np.linalg.norm(_atom_position(atom) - target_pos))
        prev = neighbors.get(id(atom))
        if prev is None or distance < prev[0]:
            neighbors[id(atom)] = (distance, atom)

    for atom in residue:
        distance = float(np.linalg.norm(_atom_position(atom) - target_pos))
        if 0.2 < distance <= 1.85:
            _add_neighbor(atom)

    for atom in _connected_heavy_neighbors_from_struct_conn(
        structure=structure,
        chain_name=chain_name,
        residue=residue,
        target_atom=target_atom,
    ):
        _add_neighbor(atom)

    ordered = sorted(neighbors.values(), key=lambda item: item[0])
    return [atom for _, atom in ordered]


def _connected_heavy_neighbors_from_struct_conn(
    structure: Any,
    chain_name: str,
    residue: Any,
    target_atom: Any,
) -> list[Any]:
    connections = getattr(structure, "connections", None)
    if not connections:
        return []

    linked_atoms: list[Any] = []
    for connection in connections:
        partner1 = getattr(connection, "partner1", None)
        partner2 = getattr(connection, "partner2", None)
        if partner1 is None or partner2 is None:
            continue

        if _atom_address_matches(partner1, chain_name, residue, target_atom):
            linked = _resolve_atom_address(structure, partner2)
        elif _atom_address_matches(partner2, chain_name, residue, target_atom):
            linked = _resolve_atom_address(structure, partner1)
        else:
            continue

        if linked is not None and not _is_hydrogen(linked):
            linked_atoms.append(linked)

    return linked_atoms


def _atom_address_matches(
    atom_address: Any,
    chain_name: str,
    residue: Any,
    atom: Any,
) -> bool:
    seqid = getattr(getattr(atom_address, "res_id", None), "seqid", None)
    resname = getattr(getattr(atom_address, "res_id", None), "name", None)
    return (
        getattr(atom_address, "chain_name", None) == chain_name
        and int(getattr(seqid, "num", -999999)) == int(residue.seqid.num)
        and str(resname) == str(residue.name)
        and str(getattr(atom_address, "atom_name", "")).strip() == str(atom.name).strip()
    )


def _resolve_atom_address(structure: Any, atom_address: Any) -> Any | None:
    chain = _get_chain(structure, getattr(atom_address, "chain_name", ""))
    if chain is None:
        return None

    seqid = getattr(getattr(atom_address, "res_id", None), "seqid", None)
    resnum = int(getattr(seqid, "num", -999999))
    resname = str(getattr(getattr(atom_address, "res_id", None), "name", ""))
    atom_name = str(getattr(atom_address, "atom_name", "")).strip()

    for residue in chain:
        if int(residue.seqid.num) != resnum:
            continue
        if resname and str(residue.name) != resname:
            continue
        return _get_atom_by_name(residue, atom_name)

    return None


def _compute_ring_frame(residue: Any) -> tuple[np.ndarray | None, np.ndarray | None]:
    ring_atom_names = ["O5", "C1", "C2", "C3", "C4", "C5"]
    ring_points = [
        _atom_position(atom)
        for atom_name in ring_atom_names
        if (atom := _get_atom_by_name(residue, atom_name)) is not None
    ]
    if len(ring_points) < 3:
        return None, None

    centroid = _compute_centroid(ring_points)
    normal = _compute_plane_normal(ring_points)
    return normal, centroid


def _compute_plane_normal(points: list[np.ndarray]) -> np.ndarray | None:
    if len(points) < 3:
        return None
    point_matrix = np.vstack(points)
    centered = point_matrix - np.mean(point_matrix, axis=0)
    if np.linalg.matrix_rank(centered) < 2:
        return None
    _, _, vh = np.linalg.svd(centered)
    return _unit_vector(vh[-1])


def _compute_centroid(points: list[np.ndarray]) -> np.ndarray:
    return np.mean(np.vstack(points), axis=0)


def _classify_sugar_face_orientation(
    ring_normal: np.ndarray | None,
    ring_centroid: np.ndarray | None,
    repositioned_cu: np.ndarray | None,
) -> str:
    if ring_normal is None or ring_centroid is None or repositioned_cu is None:
        return GEOMETRY_STATUS_NOT_COMPUTABLE

    facing = float(np.dot(ring_normal, repositioned_cu - ring_centroid))
    if abs(facing) < 1e-6:
        return "ambiguous"
    if facing > 0:
        return "ring_normal_toward_cu"
    return "ring_normal_away_from_cu"


def _score_oxyl_h_distance(distance: float, thresholds: GeometryThresholds) -> float:
    if distance < thresholds.oxyl_h_min_a or distance > thresholds.oxyl_h_max_a:
        return math.nan
    return max(
        1.0 - abs(distance - thresholds.oxyl_h_optimum_a) / thresholds.oxyl_h_score_normalization,
        0.0,
    )


def _assign_geometry_status(
    distance: float | None,
    score: float,
    thresholds: GeometryThresholds,
) -> str:
    if distance is None:
        return GEOMETRY_STATUS_NOT_COMPUTABLE
    if math.isnan(score):
        return GEOMETRY_STATUS_COMPUTABLE_IMPLAUSIBLE
    if (
        thresholds.oxyl_h_highly_plausible_min_a
        <= distance
        <= thresholds.oxyl_h_highly_plausible_max_a
    ):
        return GEOMETRY_STATUS_HIGHLY_PLAUSIBLE
    return GEOMETRY_STATUS_PLAUSIBLE


def _angle_from_points_deg(point_a: np.ndarray, vertex: np.ndarray, point_b: np.ndarray) -> float | None:
    return _angle_between_vectors_deg(point_a - vertex, point_b - vertex)


def _angle_between_vectors_deg(vector_a: np.ndarray, vector_b: np.ndarray) -> float | None:
    unit_a = _unit_vector(vector_a)
    unit_b = _unit_vector(vector_b)
    if unit_a is None or unit_b is None:
        return None

    cosine = float(np.dot(unit_a, unit_b))
    cosine = max(-1.0, min(1.0, cosine))
    return float(np.degrees(np.arccos(cosine)))


def _unit_vector(vector: np.ndarray) -> np.ndarray | None:
    norm = float(np.linalg.norm(vector))
    if norm < 1e-9:
        return None
    return vector / norm


def _format_residue_id(chain_name: str, residue: Any) -> str:
    return f"{chain_name}:{residue.name}{residue.seqid.num}"


def _get_chain(structure: Any, chain_name: str) -> Any | None:
    model = structure[0]
    for chain in model:
        if chain.name == chain_name:
            return chain
    return None


def _get_atom_by_name(residue: Any, atom_name: str) -> Any | None:
    for atom in residue:
        if atom.name.strip() == atom_name:
            return atom
    return None


def _atom_position(atom: Any) -> np.ndarray:
    return np.array([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float)


def _is_hydrogen(atom: Any) -> bool:
    return str(atom.element.name).strip().upper().startswith("H")


def _append_virtual_geometry_atoms_to_pdb_text(
    pdb_text: str,
    metrics: PoseGeometryMetrics,
) -> str:
    lines = pdb_text.splitlines()
    while lines and lines[-1].strip() == "END":
        lines.pop()

    serial = _next_pdb_serial(lines)
    virtual_atoms = _virtual_debug_atoms(metrics)
    if virtual_atoms:
        lines.append("REMARK 700 Virtual geometry atoms below are test-only debug output")
        lines.append("REMARK 700 CUX=repositioned Cu, OX=virtual oxyl, HC1/HC4=virtual H atoms")
    for atom_name, element, coords in virtual_atoms:
        lines.append(
            _format_pdb_atom_line(
                serial=serial,
                atom_name=atom_name,
                resname="DBG",
                chain_id="X",
                resseq=1,
                x=coords[0],
                y=coords[1],
                z=coords[2],
                element=element,
            )
        )
        serial += 1

    lines.append("END")
    return "\n".join(lines) + "\n"


def _virtual_debug_atoms(
    metrics: PoseGeometryMetrics,
) -> list[tuple[str, str, tuple[float, float, float]]]:
    atoms: list[tuple[str, str, tuple[float, float, float]]] = []
    if metrics.repositioned_cu_coordinates is not None:
        atoms.append(("CUX", "Cu", metrics.repositioned_cu_coordinates))
    if metrics.virtual_oxyl_coordinates is not None:
        atoms.append(("OX", "O", metrics.virtual_oxyl_coordinates))
    if metrics.virtual_h_c1_coordinates is not None:
        atoms.append(("HC1", "H", metrics.virtual_h_c1_coordinates))
    if metrics.virtual_h_c4_coordinates is not None:
        atoms.append(("HC4", "H", metrics.virtual_h_c4_coordinates))
    return atoms


def _next_pdb_serial(lines: list[str]) -> int:
    max_serial = 0
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                max_serial = max(max_serial, int(line[6:11].strip()))
            except ValueError:
                continue
    return max_serial + 1


def _format_pdb_atom_line(
    *,
    serial: int,
    atom_name: str,
    resname: str,
    chain_id: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    return (
        f"HETATM{serial:>5} {atom_name:>4} {resname:>3} {chain_id:1}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{0.00:>6.2f}          {element:>2}"
    )
