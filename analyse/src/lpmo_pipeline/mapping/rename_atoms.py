# src/lpmo_pipeline/mapping/rename_atoms.py
"""
Responsibility: Apply atom renaming from atom_map.tsv to a gemmi.Structure.
Input:  gemmi.Structure + MappingResult (or loaded atom_map.tsv)
Output: Modified gemmi.Structure with canonical atom names
"""
from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from lpmo_pipeline.mapping.cross_model_atom_mapping import AtomMapping

logger = logging.getLogger(__name__)


def load_atom_map_tsv(path: Path) -> list[AtomMapping]:
    """Load atom_map.tsv into a list of AtomMapping objects."""
    mappings: list[AtomMapping] = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            mappings.append(AtomMapping(
                chain=row["chain"],
                resname=row["resname"],
                resnum=int(row["resnum"]),
                old_atom_name=row["old_atom"],
                new_atom_name=row["new_atom"],
                element=row["element"],
                match_method=row["match_method"],
            ))
    logger.info("Loaded %d atom mappings from %s", len(mappings), path)
    return mappings


def apply_rename(
    structure: Any,  # gemmi.Structure
    mappings: list[AtomMapping],
) -> tuple[Any, int]:
    """Apply atom renaming to *structure* in-place.

    Args:
        structure: gemmi.Structure to modify.
        mappings: List of AtomMapping entries (old→new).

    Returns:
        (modified_structure, n_renamed)
    """
    # Build lookup: (chain, resnum, old_atom_name) → new_atom_name
    rename_lookup: dict[tuple[str, int, str], str] = {}
    for m in mappings:
        if m.old_atom_name != m.new_atom_name:
            rename_lookup[(m.chain, m.resnum, m.old_atom_name)] = m.new_atom_name

    n_renamed = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    key = (chain.name, residue.seqid.num, atom.name)
                    if key in rename_lookup:
                        old = atom.name
                        atom.name = rename_lookup[key]
                        n_renamed += 1
                        logger.debug(
                            "Renamed %s:%d %s → %s",
                            chain.name, residue.seqid.num, old, atom.name,
                        )

    logger.info("Applied %d atom renames", n_renamed)
    return structure, n_renamed
