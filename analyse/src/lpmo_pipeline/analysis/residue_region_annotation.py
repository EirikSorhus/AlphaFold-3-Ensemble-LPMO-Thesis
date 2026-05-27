"""Metadata-backed residue region annotation for core/non-core analysis."""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from lpmo_pipeline.analysis.residue_contact_extraction import ResidueRegionFlags
from lpmo_pipeline.analysis.prolif_ifp import IFPResult, parse_ifp_feature_name


_PROTEIN_RESIDUE_PATTERN = re.compile(
    r"^(?P<residue_name>[A-Za-z0-9]+?)(?P<residue_number>-?\d+)(?:\.(?P<residue_chain>[A-Za-z0-9]+))?$"
)


@dataclass(frozen=True)
class ProteinRegionDefinition:
    """Residue ranges resolved to normalized chain-A numbering."""

    protein_id: str
    core_start: int | None = None
    core_end: int | None = None

    def flags_for(self, chain: str, residue_number: int, residue_name: str, construct_type: str) -> ResidueRegionFlags:
        if chain != "A":
            return ResidueRegionFlags()
        normalized_construct = str(construct_type).strip().lower()
        if normalized_construct in {"domain_only", "catalytic_domain"}:
            return ResidueRegionFlags(is_core_region=True)
        is_core = (
            self.core_start is not None
            and self.core_end is not None
            and self.core_start <= residue_number <= self.core_end
        )
        return ResidueRegionFlags(
            is_core_region=is_core,
            is_non_core_region=not is_core,
        )


def _as_int(value: Any) -> int | None:
    if value in (None, ""):
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _first_present(row: dict[str, Any], *keys: str) -> str:
    for key in keys:
        value = row.get(key)
        if value not in (None, ""):
            text = str(value).strip()
            if text:
                return text
    return ""


def load_protein_region_definitions(path: Path | None) -> dict[str, ProteinRegionDefinition]:
    """Load metadata core ranges and convert full-sequence positions to mature-relative numbering."""

    if path is None or not Path(path).exists():
        return {}

    definitions: dict[str, ProteinRegionDefinition] = {}
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            protein_id = _first_present(row, "protein_id", "uniprot_id", "UniProt_ID")
            if not protein_id or protein_id in definitions:
                continue
            signal_end = _as_int(_first_present(row, "Signal_End", "signal_end")) or 0
            core_start_full = _as_int(_first_present(row, "LPMO_Core_Start", "LPMO_CoreStart", "core_start"))
            core_end_full = _as_int(_first_present(row, "LPMO_Core_End", "core_end"))
            core_start = (core_start_full - signal_end) if core_start_full is not None else None
            core_end = (core_end_full - signal_end) if core_end_full is not None else None
            if core_start is not None:
                core_start = max(core_start, 1)
            if core_end is not None:
                core_end = max(core_end, 1)
            definitions[protein_id] = ProteinRegionDefinition(
                protein_id=protein_id,
                core_start=core_start,
                core_end=core_end,
            )
    return definitions


def parse_protein_residue_label(protein_residue_label: str) -> tuple[str, int, str] | None:
    """Parse flattened ProLIF protein residue labels such as ``ASN88.A``."""

    match = _PROTEIN_RESIDUE_PATTERN.match(str(protein_residue_label))
    if match is None or match.group("residue_chain") is None:
        return None
    return (
        str(match.group("residue_chain")),
        int(match.group("residue_number")),
        str(match.group("residue_name")),
    )


def build_region_annotations_for_ifp_results(
    results: list[IFPResult],
    *,
    protein_id: str,
    construct_type: str,
    region_definitions: dict[str, ProteinRegionDefinition] | None = None,
) -> dict[str, dict[tuple[str, int, str], ResidueRegionFlags]]:
    """Build per-pose residue flags for all residues observed in full-chain raw IFP features."""

    region_definition = (region_definitions or {}).get(protein_id)
    annotations: dict[str, dict[tuple[str, int, str], ResidueRegionFlags]] = {}
    for result in results:
        pose_annotations: dict[tuple[str, int, str], ResidueRegionFlags] = {}
        for feature_name in result.feature_names:
            parsed_feature = parse_ifp_feature_name(feature_name)
            parsed_residue = parse_protein_residue_label(parsed_feature[1])
            if parsed_residue is None:
                continue
            chain, residue_number, residue_name = parsed_residue
            if region_definition is None:
                flags = ResidueRegionFlags(is_core_region=True)
            else:
                flags = region_definition.flags_for(chain, residue_number, residue_name, construct_type)
            pose_annotations[(chain, residue_number, residue_name)] = flags
        annotations[result.pose_id] = pose_annotations
    return annotations


def residue_is_core(
    protein_residue_label: str,
    *,
    protein_id: str,
    construct_type: str,
    region_definitions: dict[str, ProteinRegionDefinition] | None = None,
) -> bool:
    """Return whether an IFP protein residue label belongs to the core region."""

    parsed = parse_protein_residue_label(protein_residue_label)
    if parsed is None:
        return False
    chain, residue_number, residue_name = parsed
    definition = (region_definitions or {}).get(protein_id)
    if definition is None:
        return True
    return definition.flags_for(chain, residue_number, residue_name, construct_type).is_core_region
