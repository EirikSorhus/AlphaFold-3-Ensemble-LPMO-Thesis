"""Residue-level contact extraction for pose_residue_contact_table.tsv."""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from lpmo_pipeline.analysis.prolif_ifp import IFPResult


_FEATURE_SEPARATOR = "|"
_PROTEIN_RESIDUE_PATTERN = re.compile(
    r"^(?P<residue_name>[A-Za-z0-9]+?)(?P<residue_number>-?\d+)(?:\.(?P<residue_chain>[A-Za-z0-9]+))?$"
)

POSE_RESIDUE_CONTACT_COLUMNS = [
    "pose_id",
    "protein_id",
    "condition_id",
    "residue_chain",
    "residue_number",
    "residue_name",
    "interaction_type",
    "contact_present",
    "ligand_residue_label",
    "distance_if_available",
    "is_core_region",
    "is_catalytic_surface_region",
    "is_non_core_region",
    "is_cbm_region",
    "is_linker_region",
]


@dataclass(frozen=True)
class ResidueRegionFlags:
    """Optional residue-region annotations for downstream interpretation."""

    is_core_region: bool = False
    is_catalytic_surface_region: bool = False
    is_non_core_region: bool = False
    is_cbm_region: bool = False
    is_linker_region: bool = False


@dataclass(frozen=True)
class ResidueContactRecord:
    """One residue-level interaction row derived directly from a ProLIF feature."""

    pose_id: str
    protein_id: str
    condition_id: str
    residue_chain: str
    residue_number: int
    residue_name: str
    interaction_type: str
    contact_present: int
    ligand_residue_label: str
    distance_if_available: float | None = None
    is_core_region: bool = False
    is_catalytic_surface_region: bool = False
    is_non_core_region: bool = False
    is_cbm_region: bool = False
    is_linker_region: bool = False

    def to_row(self) -> dict[str, Any]:
        return {
            "pose_id": self.pose_id,
            "protein_id": self.protein_id,
            "condition_id": self.condition_id,
            "residue_chain": self.residue_chain,
            "residue_number": self.residue_number,
            "residue_name": self.residue_name,
            "interaction_type": self.interaction_type,
            "contact_present": self.contact_present,
            "ligand_residue_label": self.ligand_residue_label,
            "distance_if_available": "" if self.distance_if_available is None else self.distance_if_available,
            "is_core_region": self.is_core_region,
            "is_catalytic_surface_region": self.is_catalytic_surface_region,
            "is_non_core_region": self.is_non_core_region,
            "is_cbm_region": self.is_cbm_region,
            "is_linker_region": self.is_linker_region,
        }


def _parse_feature_name(feature_name: str) -> tuple[str, str, str]:
    parts = feature_name.split(_FEATURE_SEPARATOR)
    if len(parts) != 3:
        raise ValueError(
            "Unsupported IFP feature format for residue contact extraction: "
            f"{feature_name!r}"
        )
    ligand_residue_label, protein_residue_label, interaction_type = (part.strip() for part in parts)
    if not ligand_residue_label or not protein_residue_label or not interaction_type:
        raise ValueError(
            "IFP feature contains an empty component and cannot be mapped to a residue-contact row: "
            f"{feature_name!r}"
        )
    return ligand_residue_label, protein_residue_label, interaction_type


def _parse_protein_residue_label(protein_residue_label: str) -> tuple[str, int, str]:
    match = _PROTEIN_RESIDUE_PATTERN.match(protein_residue_label)
    if match is None or match.group("residue_chain") is None:
        raise ValueError(
            "Unsupported protein residue label for residue contact extraction: "
            f"{protein_residue_label!r}"
        )

    residue_name = match.group("residue_name")
    residue_number = int(match.group("residue_number"))
    residue_chain = match.group("residue_chain")
    return residue_chain, residue_number, residue_name


def build_pose_residue_contact_rows(
    results: list[IFPResult],
    pose_metadata_by_id: dict[str, dict[str, str]],
    *,
    region_annotations_by_pose_id: dict[str, dict[tuple[str, int, str], ResidueRegionFlags]] | None = None,
) -> list[dict[str, Any]]:
    """Build pose_residue_contact_table.tsv rows from per-pose ProLIF results.

    The extraction is intentionally direct: each emitted row corresponds to one
    existing ProLIF feature bit so the residue-contact table stays consistent
    with the feature labels used by pose_ifp_table.tsv.
    """

    rows: list[dict[str, Any]] = []
    region_annotations_by_pose_id = region_annotations_by_pose_id or {}

    for result in results:
        if result.status != "ok" or not result.feature_names:
            continue
        if len(result.feature_names) != len(result.flat_bitvector):
            raise ValueError(
                "feature_names and flat_bitvector must have the same length for residue contact extraction "
                f"(pose_id={result.pose_id!r})"
            )

        pose_metadata = pose_metadata_by_id.get(result.pose_id)
        if pose_metadata is None:
            raise KeyError(f"Missing pose metadata for residue contact extraction: {result.pose_id!r}")

        pose_regions = region_annotations_by_pose_id.get(result.pose_id, {})
        for feature_name, bit in zip(result.feature_names, result.flat_bitvector, strict=True):
            ligand_residue_label, protein_residue_label, interaction_type = _parse_feature_name(feature_name)
            residue_chain, residue_number, residue_name = _parse_protein_residue_label(protein_residue_label)
            flags = pose_regions.get(
                (residue_chain, residue_number, residue_name),
                ResidueRegionFlags(),
            )
            rows.append(
                ResidueContactRecord(
                    pose_id=result.pose_id,
                    protein_id=pose_metadata["protein_id"],
                    condition_id=pose_metadata["condition_id"],
                    residue_chain=residue_chain,
                    residue_number=residue_number,
                    residue_name=residue_name,
                    interaction_type=interaction_type,
                    contact_present=int(bit),
                    ligand_residue_label=ligand_residue_label,
                    is_core_region=flags.is_core_region,
                    is_catalytic_surface_region=flags.is_catalytic_surface_region,
                    is_non_core_region=flags.is_non_core_region,
                    is_cbm_region=flags.is_cbm_region,
                    is_linker_region=flags.is_linker_region,
                ).to_row()
            )

    return rows


def write_pose_residue_contact_table(rows: list[dict[str, Any]], output_path: Path) -> None:
    """Write the canonical residue-level contact table as TSV."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=POSE_RESIDUE_CONTACT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
