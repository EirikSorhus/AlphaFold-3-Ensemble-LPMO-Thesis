# src/lpmo_pipeline/mapping/cross_model_atom_mapping.py
"""
Responsibility: Build a cross-model atom mapping for a single structure.
  Atom names are NOT assumed identical between AF3, RF3, Boltz-2.
  Mapping key: [element, residue_ccd, local_bond_graph, 3D_proximity]

Input:  normalized.cif (gemmi.Structure)
        reference bond graphs from CCD cache
Output: atom_map.tsv  — columns: chain, resname, resnum, old_atom, new_atom, match_method
        rename_log.json — detailed log of decisions

HARD RULE: atom_mapping_coverage must be 100%.  Any unmapped atom → GateFailure.
"""
from __future__ import annotations

import csv
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------
@dataclass
class AtomMapping:
    """Single atom mapping entry."""

    chain: str
    resname: str
    resnum: int
    old_atom_name: str
    new_atom_name: str
    element: str
    match_method: str  # "element_bondgraph" | "element_3d_proximity" | "exact"


@dataclass
class MappingResult:
    """Full mapping result for one structure."""

    mappings: list[AtomMapping] = field(default_factory=list)
    unmapped: list[dict[str, Any]] = field(default_factory=list)
    coverage: float = 0.0  # Must be 1.0 to pass gate
    log_entries: list[dict[str, Any]] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Core mapping logic
# ---------------------------------------------------------------------------
def build_atom_mapping(
    structure: Any,  # gemmi.Structure
    ccd_bond_graphs: dict[str, dict[str, Any]],
    reference_atom_names: dict[str, dict[str, str]] | None = None,
) -> MappingResult:
    """Build cross-model atom mapping for all residues in *structure*.

    Strategy (applied per residue, in order of preference):
        1. Exact match: atom name already matches CCD reference → keep
        2. Element + local bond graph: match by element type and bond
           connectivity pattern within the residue
        3. Element + 3D proximity: for remaining atoms, match by element
           and nearest-neighbour in 3D space to reference positions

    Args:
        structure: Parsed gemmi.Structure (normalized chains).
        ccd_bond_graphs: {comp_id: {atom_name: [bonded_atom_names]}} from CCD.
        reference_atom_names: Optional {comp_id: {element_idx: canonical_name}}
            for explicit renaming targets.

    Returns:
        MappingResult with mappings, unmapped list, coverage, and log.
    """
    result = MappingResult()
    total_atoms = 0
    mapped_atoms = 0

    # PSEUDOCODE — iterate over all residues
    for model in structure:  # gemmi: structure[0] for first model
        for chain in model:
            for residue in chain:
                comp_id = residue.name
                ref_graph = ccd_bond_graphs.get(comp_id, {})

                for atom in residue:
                    total_atoms += 1
                    mapping = _map_single_atom(
                        atom=atom,
                        residue=residue,
                        chain_name=chain.name,
                        ref_graph=ref_graph,
                        reference_atom_names=reference_atom_names,
                    )

                    if mapping is not None:
                        result.mappings.append(mapping)
                        mapped_atoms += 1
                        result.log_entries.append({
                            "chain": chain.name,
                            "resname": comp_id,
                            "resnum": residue.seqid.num,
                            "old": atom.name,
                            "new": mapping.new_atom_name,
                            "method": mapping.match_method,
                        })
                    else:
                        result.unmapped.append({
                            "chain": chain.name,
                            "resname": comp_id,
                            "resnum": residue.seqid.num,
                            "atom": atom.name,
                            "element": atom.element.name,
                        })

    result.coverage = mapped_atoms / total_atoms if total_atoms > 0 else 0.0
    logger.info(
        "Atom mapping: %d/%d mapped (%.1f%%). Unmapped: %d",
        mapped_atoms, total_atoms, result.coverage * 100, len(result.unmapped),
    )

    return result


def _map_single_atom(
    atom: Any,  # gemmi.Atom
    residue: Any,  # gemmi.Residue
    chain_name: str,
    ref_graph: dict[str, Any],
    reference_atom_names: dict[str, dict[str, str]] | None,
) -> AtomMapping | None:
    """Attempt to map a single atom using the 3-tier strategy.

    Returns AtomMapping on success, None on failure.
    """
    comp_id = residue.name
    resnum = residue.seqid.num
    element = atom.element.name
    old_name = atom.name

    # --- Tier 1: Exact match ---
    if old_name in ref_graph:
        return AtomMapping(
            chain=chain_name, resname=comp_id, resnum=resnum,
            old_atom_name=old_name, new_atom_name=old_name,
            element=element, match_method="exact",
        )

    # --- Tier 2: Element + local bond graph ---
    # Build local bond graph for this atom from the structure
    local_neighbors = _get_bonded_neighbors(atom, residue)
    local_neighbor_elements = sorted([n.element.name for n in local_neighbors])

    # Find CCD atom with same element and same neighbor-element signature
    for ref_atom_name, ref_neighbors in ref_graph.items():
        ref_element = _extract_element_from_atom_name(ref_atom_name)
        if ref_element != element:
            continue
        ref_neighbor_elements = sorted([
            _extract_element_from_atom_name(n) for n in ref_neighbors
        ])
        if local_neighbor_elements == ref_neighbor_elements:
            return AtomMapping(
                chain=chain_name, resname=comp_id, resnum=resnum,
                old_atom_name=old_name, new_atom_name=ref_atom_name,
                element=element, match_method="element_bondgraph",
            )

    # --- Tier 3: Element + 3D proximity ---
    # Find closest atom of same element in reference (if reference positions available)
    if reference_atom_names and comp_id in reference_atom_names:
        # Match by element + proximity (placeholder for 3D KDTree logic)
        ref_names = reference_atom_names[comp_id]
        for _key, ref_name in ref_names.items():
            if _extract_element_from_atom_name(ref_name) == element:
                return AtomMapping(
                    chain=chain_name, resname=comp_id, resnum=resnum,
                    old_atom_name=old_name, new_atom_name=ref_name,
                    element=element, match_method="element_3d_proximity",
                )

    # Failed to map
    logger.warning(
        "UNMAPPED atom: chain=%s res=%s%d atom=%s elem=%s",
        chain_name, comp_id, resnum, old_name, element,
    )
    return None


def _get_bonded_neighbors(atom: Any, residue: Any) -> list[Any]:
    """Get bonded neighbors of *atom* within *residue*.

    PSEUDOCODE: Use gemmi's bond information or distance-based heuristic.
    """
    neighbors = []
    # Use gemmi topology or distance cutoff (1.0–1.8 Å for covalent bonds)
    BOND_CUTOFF = 1.85  # Å
    for other_atom in residue:
        if other_atom.name == atom.name:
            continue
        dist = atom.pos.dist(other_atom.pos)
        if dist < BOND_CUTOFF:
            neighbors.append(other_atom)
    return neighbors


def _extract_element_from_atom_name(atom_name: str) -> str:
    """Extract element symbol from atom name (e.g. 'CA' → 'C', 'OG1' → 'O')."""
    # Strip digits and whitespace, take first 1-2 chars
    clean = atom_name.strip()
    if len(clean) >= 2 and clean[:2] in {"CU", "FE", "ZN", "MN", "MG", "CA"}:
        # Two-letter elements (careful: CA = C-alpha, not calcium in proteins)
        # Context-dependent; for metals in chain E, treat as 2-letter
        return clean[:2]
    return clean[0] if clean else ""


# ---------------------------------------------------------------------------
# I/O: write atom_map.tsv and rename_log.json
# ---------------------------------------------------------------------------
def write_atom_map_tsv(result: MappingResult, output_path: Path) -> None:
    """Write atom_map.tsv following schema/atom_map_tsv_schema.json."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "chain", "resname", "resnum", "old_atom", "new_atom",
            "element", "match_method",
        ])
        for m in result.mappings:
            writer.writerow([
                m.chain, m.resname, m.resnum, m.old_atom_name,
                m.new_atom_name, m.element, m.match_method,
            ])
    logger.info("Wrote atom_map.tsv (%d entries) to %s", len(result.mappings), output_path)


def write_rename_log(result: MappingResult, output_path: Path) -> None:
    """Write rename_log.json with full decision log."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    log = {
        "coverage": result.coverage,
        "total_mapped": len(result.mappings),
        "total_unmapped": len(result.unmapped),
        "entries": result.log_entries,
        "unmapped": result.unmapped,
    }
    with open(output_path, "w") as f:
        json.dump(log, f, indent=2)
    logger.info("Wrote rename_log.json to %s", output_path)
