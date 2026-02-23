"""Ligand helpers shared across model runners."""

from __future__ import annotations

import logging
from typing import Any

from .oligo import OligoRegistry, OligoSpec, build_af3_bonded_atom_pairs, build_af3_oligo_ligand_entry

logger = logging.getLogger(__name__)

MANDATORY_CU_CCD = "CU"
AF3_CU_LIGAND_ID = "B"
AF3_MAIN_LIGAND_ID = "C"
BOLTZ_CU_LIGAND_ID = "B"
BOLTZ_MAIN_LIGAND_ID = "C"


def build_af3_ligand_entries(
    main_ccd: str,
    oligo_registry: OligoRegistry | None = None,
) -> tuple[list[dict[str, Any]], list[list[list[Any]]]]:
    """Build AF3 ligand entries, always including CU when needed.

    When *oligo_registry* is provided and *main_ccd* matches a known
    oligosaccharide (e.g. ``"CEL6"``), the main ligand entry is built
    with multiple monomer ``ccdCodes`` and the corresponding
    ``bondedAtomPairs`` are returned.

    Returns:
        A tuple ``(sequences_entries, bonded_atom_pairs)``.
        ``bonded_atom_pairs`` is an empty list for non-oligo ligands.
    """
    bonded_atom_pairs: list[list[list[Any]]] = []

    # Check for oligosaccharide
    oligo_info = None
    if oligo_registry is not None:
        oligo_info = oligo_registry.parse(main_ccd)

    if main_ccd == MANDATORY_CU_CCD:
        return (
            [{"ligand": {"id": AF3_CU_LIGAND_ID, "ccdCodes": [MANDATORY_CU_CCD]}}],
            bonded_atom_pairs,
        )

    cu_entry = {"ligand": {"id": AF3_CU_LIGAND_ID, "ccdCodes": [MANDATORY_CU_CCD]}}

    if oligo_info is not None:
        prefix, n, spec = oligo_info
        logger.info(
            "Building AF3 oligo ligand: %s → %d× %s (bond: %s–%s)",
            main_ccd, n, spec.monomer, spec.bond_atom_pair[0], spec.bond_atom_pair[1],
        )
        oligo_entry = build_af3_oligo_ligand_entry(AF3_MAIN_LIGAND_ID, spec, n)
        bonded_atom_pairs = build_af3_bonded_atom_pairs(AF3_MAIN_LIGAND_ID, spec, n)
        return [cu_entry, oligo_entry], bonded_atom_pairs

    return (
        [cu_entry, {"ligand": {"id": AF3_MAIN_LIGAND_ID, "ccdCodes": [main_ccd]}}],
        bonded_atom_pairs,
    )


def build_boltz_ligand_entries(main_ccd: str) -> tuple[list[dict[str, Any]], str]:
    """Build Boltz ligand entries and return the main ligand id."""
    if main_ccd == MANDATORY_CU_CCD:
        return (
            [
                {"ligand": {"id": BOLTZ_CU_LIGAND_ID, "ccd": MANDATORY_CU_CCD}},
            ],
            BOLTZ_CU_LIGAND_ID,
        )

    return (
        [
            {"ligand": {"id": BOLTZ_CU_LIGAND_ID, "ccd": MANDATORY_CU_CCD}},
            {"ligand": {"id": BOLTZ_MAIN_LIGAND_ID, "ccd": main_ccd}},
        ],
        BOLTZ_MAIN_LIGAND_ID,
    )


def build_rf3_ligand_components(
    main_ccd: str,
    main_cif_path: str | None,
) -> list[dict[str, Any]]:
    """Build RF3 ligand components, always including CU when needed."""
    components: list[dict[str, Any]] = []

    if main_ccd != MANDATORY_CU_CCD:
        components.append({"ccd_code": MANDATORY_CU_CCD})

    if main_cif_path:
        components.append({"path": main_cif_path})
    else:
        components.append({"ccd_code": main_ccd})

    return components
