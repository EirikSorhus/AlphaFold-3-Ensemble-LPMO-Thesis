# src/lpmo_pipeline/io/ccd_lookup.py
"""
Responsibility: Validate and look up CCD (Chemical Component Dictionary) monosaccharide codes.
Input:  comp_id strings from mmCIF _chem_comp
Output: Validation result (valid CCD mono? yes/no) + canonical bond graph

HARD RULE: Privateer input must use valid CCD monosaccharides ONLY.
           Custom oligomer IDs like CEL6 are REJECTED.
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Known CCD monosaccharide codes relevant for LPMO substrates
# This is the canonical allow-list. Extend as needed, but NEVER add
# multi-residue oligomer aliases (CEL6, CBI, etc.).
# ---------------------------------------------------------------------------
VALID_CCD_MONOSACCHARIDES: set[str] = {
    # N-acetylglucosamine (chitin)
    "NAG",
    # Glucose / cellulose monomers
    "BGC",  # beta-D-glucose
    "GLC",  # alpha-D-glucose
    "AGC",  # alpha-D-glucose (alternate)
    # Mannose
    "MAN",  # alpha-D-mannose
    "BMA",  # beta-D-mannose
    # Galactose
    "GAL",  # beta-D-galactose
    "GLA",  # alpha-D-galactose
    # Xylose
    "XYS",  # beta-D-xylose
    "XYP",  # alpha-D-xylopyranose
    # Fucose
    "FUC",  # alpha-L-fucose
    "FCA",  # alpha-L-fucose (alternate)
    # Glucuronic acid
    "GCU",  # alpha-D-glucuronic acid
    "BDP",  # beta-D-glucuronic acid
    # N-acetylmuramic acid
    "MUB",
    # Maltose-type
    "MAL",  # maltose (disaccharide — check if Privateer handles)
    # Amylose-relevant
    "GLC",  # reuse
}

# Known INVALID oligomer aliases that models sometimes output
REJECTED_OLIGOMER_ALIASES: set[str] = {
    "CEL6", "CEL5", "CEL4", "CEL3", "CEL2",
    "CBI", "NAG6", "NAG5", "NAG4", "NAG3",
}


@dataclass
class CCDLookupResult:
    """Result of a CCD monosaccharide lookup."""

    comp_id: str
    is_valid_ccd_mono: bool
    canonical_name: str = ""
    bond_graph: dict[str, Any] = field(default_factory=dict)
    rejection_reason: str = ""


def validate_comp_id(comp_id: str) -> CCDLookupResult:
    """Check if a comp_id is a valid CCD monosaccharide.

    Args:
        comp_id: The 3-letter CCD component ID from mmCIF.

    Returns:
        CCDLookupResult with is_valid_ccd_mono flag.
    """
    comp_id_upper = comp_id.strip().upper()

    if comp_id_upper in REJECTED_OLIGOMER_ALIASES:
        logger.warning(
            "REJECTED oligomer alias '%s' — must decompose into monomers", comp_id
        )
        return CCDLookupResult(
            comp_id=comp_id_upper,
            is_valid_ccd_mono=False,
            rejection_reason=f"Custom oligomer alias '{comp_id_upper}' not allowed; "
            "decompose into constituent CCD monosaccharides.",
        )

    if comp_id_upper in VALID_CCD_MONOSACCHARIDES:
        return CCDLookupResult(
            comp_id=comp_id_upper,
            is_valid_ccd_mono=True,
            canonical_name=comp_id_upper,
        )

    # Unknown code — attempt CCD database query (placeholder)
    logger.info("comp_id '%s' not in built-in list; querying CCD cache...", comp_id)
    result = _query_ccd_cache(comp_id_upper)
    return result


def validate_all_glycan_residues(
    comp_ids: list[str],
) -> tuple[bool, list[CCDLookupResult]]:
    """Validate all glycan comp_ids in a structure.

    Args:
        comp_ids: List of comp_id strings from the glycan chain(s).

    Returns:
        (all_valid, list_of_results).
        Gate: all_valid must be True (100% recognized).
    """
    results = [validate_comp_id(c) for c in comp_ids]
    all_valid = all(r.is_valid_ccd_mono for r in results)

    if not all_valid:
        failed = [r.comp_id for r in results if not r.is_valid_ccd_mono]
        logger.error(
            "GATE FAIL: privateer_recognized_sugars < 100%%. "
            "Unrecognized comp_ids: %s",
            failed,
        )

    return all_valid, results


# ---------------------------------------------------------------------------
# CCD cache (offline / pre-downloaded)
# ---------------------------------------------------------------------------
_CCD_CACHE: dict[str, dict] = {}  # Loaded lazily


def load_ccd_cache(cache_path: Path) -> None:
    """Load a pre-downloaded CCD JSON cache.

    Expected format: {comp_id: {name, type, bond_graph: {atom1-atom2: order}}}
    """
    global _CCD_CACHE
    with open(cache_path) as f:
        _CCD_CACHE = json.load(f)
    logger.info("Loaded CCD cache with %d entries from %s", len(_CCD_CACHE), cache_path)


def _query_ccd_cache(comp_id: str) -> CCDLookupResult:
    """Look up comp_id in the pre-loaded CCD cache."""
    if comp_id in _CCD_CACHE:
        entry = _CCD_CACHE[comp_id]
        # Check if it's a monosaccharide type
        ccd_type = entry.get("type", "")
        is_mono = ccd_type in {"D-saccharide", "L-saccharide", "saccharide",
                                "D-saccharide, alpha linking",
                                "D-saccharide, beta linking",
                                "L-saccharide, alpha linking",
                                "L-saccharide, beta linking"}
        return CCDLookupResult(
            comp_id=comp_id,
            is_valid_ccd_mono=is_mono,
            canonical_name=entry.get("name", comp_id),
            bond_graph=entry.get("bond_graph", {}),
            rejection_reason="" if is_mono else f"CCD type '{ccd_type}' not a monosaccharide",
        )

    return CCDLookupResult(
        comp_id=comp_id,
        is_valid_ccd_mono=False,
        rejection_reason=f"comp_id '{comp_id}' not found in CCD cache",
    )
