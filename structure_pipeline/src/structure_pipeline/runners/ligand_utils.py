"""Ligand helpers shared across model runners."""

from __future__ import annotations

from typing import Any

MANDATORY_CU_CCD = "CU"
AF3_CU_LIGAND_ID = "B"
AF3_MAIN_LIGAND_ID = "C"
BOLTZ_CU_LIGAND_ID = "B"
BOLTZ_MAIN_LIGAND_ID = "C"


def build_af3_ligand_entries(main_ccd: str) -> list[dict[str, Any]]:
    """Build AF3 ligand entries, always including CU when needed."""
    if main_ccd == MANDATORY_CU_CCD:
        return [
            {"ligand": {"id": AF3_CU_LIGAND_ID, "ccdCodes": [MANDATORY_CU_CCD]}},
        ]

    return [
        {"ligand": {"id": AF3_CU_LIGAND_ID, "ccdCodes": [MANDATORY_CU_CCD]}},
        {"ligand": {"id": AF3_MAIN_LIGAND_ID, "ccdCodes": [main_ccd]}},
    ]


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
