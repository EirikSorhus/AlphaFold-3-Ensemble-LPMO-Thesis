"""Family-aligned residue enrichment from Stage 16b outputs.

This module intentionally starts at the aligned-residue layer. Sequence sourcing
and alignment construction are handled elsewhere; callers must provide an
explicit mapping from protein residue labels to alignment columns or other
family-comparable positions.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any


FAMILY_ALIGNED_RESIDUE_COLUMNS = [
    "family_label",
    "family_aggregation_id",
    "alignment_column",
    "protein_id",
    "construct_type",
    "catalytic_core_residue_index",
    "alignment_residue",
    "residue_chain",
    "residue_number",
    "residue_name",
    "residue_label",
    "n_conditions_with_valid_clusters",
    "n_conditions_with_contact",
    "mean_condition_residue_contact_score",
    "max_condition_residue_contact_score",
    "mean_c1_weighted_residue_score",
    "mean_c4_weighted_residue_score",
    "c1_minus_c4_weighted_delta",
    "is_catalytic_surface_region",
    "is_cbm_region",
    "is_linker_region",
]

FAMILY_RESIDUE_ENRICHMENT_COLUMNS = [
    "family_label",
    "alignment_column",
    "n_proteins_observed",
    "n_proteins_with_contact",
    "n_family_aggregation_units_observed",
    "n_family_aggregation_units_with_contact",
    "mean_condition_residue_contact_score",
    "max_condition_residue_contact_score",
    "mean_c1_weighted_residue_score",
    "mean_c4_weighted_residue_score",
    "mean_c1_minus_c4_weighted_delta",
    "n_positive_c1_minus_c4_delta",
    "n_negative_c1_minus_c4_delta",
]


@dataclass(frozen=True)
class FamilyResidueEnrichmentOutputs:
    family_aligned_residue_table: list[dict[str, Any]]
    family_residue_enrichment: list[dict[str, Any]]


def _write_tsv(path: Path, columns: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_family_aligned_residue_table(rows: list[dict[str, Any]], output_path: Path) -> None:
    _write_tsv(output_path, FAMILY_ALIGNED_RESIDUE_COLUMNS, rows)


def write_family_residue_enrichment(rows: list[dict[str, Any]], output_path: Path) -> None:
    _write_tsv(output_path, FAMILY_RESIDUE_ENRICHMENT_COLUMNS, rows)


def _as_float(value: Any) -> float | None:
    if value in {None, "", "None"}:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _as_int(value: Any) -> int:
    numeric = _as_float(value)
    return int(numeric) if numeric is not None else 0


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


def _mean(values: list[float]) -> float | str:
    if not values:
        return ""
    return sum(values) / len(values)


def _family_aggregation_id(row: dict[str, Any]) -> str:
    value = str(row.get("family_aggregation_id", "")).strip()
    if value:
        return value
    return str(row.get("protein_id", "")).strip()


def _residue_label_from_row(row: dict[str, Any]) -> str:
    label = str(row.get("residue_label", "")).strip()
    if label:
        return label
    residue_name = str(row.get("residue_name", "")).strip()
    residue_number = _as_int(row.get("residue_number"))
    residue_chain = str(row.get("residue_chain", "")).strip()
    if not residue_name:
        return ""
    return f"{residue_name}{residue_number}.{residue_chain}"


def _index_protein_metadata(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    indexed: dict[str, dict[str, Any]] = {}
    for row in rows:
        protein_id = str(row.get("protein_id") or row.get("uniprot_id") or "").strip()
        if protein_id and protein_id not in indexed:
            indexed[protein_id] = row
    return indexed


def _resolve_family_label(
    alignment_row: dict[str, Any],
    metadata_row: dict[str, Any] | None,
) -> str:
    for key in ("family_label", "family", "CAZy_family", "cazy_family"):
        value = str(alignment_row.get(key, "")).strip()
        if value:
            return value
    if metadata_row is None:
        return ""
    for key in ("family_label", "family", "CAZy_family", "cazy_family"):
        value = str(metadata_row.get(key, "")).strip()
        if value:
            return value
    return ""


def compute_family_residue_enrichment_outputs(
    *,
    protein_condition_residue_score_rows: list[dict[str, Any]],
    protein_residue_regio_delta_rows: list[dict[str, Any]],
    residue_alignment_rows: list[dict[str, Any]],
    protein_metadata_rows: list[dict[str, Any]] | None = None,
) -> FamilyResidueEnrichmentOutputs:
    protein_metadata_by_id = _index_protein_metadata(protein_metadata_rows or [])

    condition_scores_by_residue: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in protein_condition_residue_score_rows:
        protein_id = str(row.get("protein_id", "")).strip()
        residue_label = _residue_label_from_row(row)
        if protein_id and residue_label:
            condition_scores_by_residue[(protein_id, residue_label)].append(row)

    regio_by_residue: dict[tuple[str, str], dict[str, Any]] = {}
    for row in protein_residue_regio_delta_rows:
        protein_id = str(row.get("protein_id", "")).strip()
        residue_label = _residue_label_from_row(row)
        if protein_id and residue_label:
            regio_by_residue[(protein_id, residue_label)] = row

    family_aligned_rows: list[dict[str, Any]] = []
    for alignment_row in residue_alignment_rows:
        protein_id = str(alignment_row.get("protein_id", "")).strip()
        residue_label = _residue_label_from_row(alignment_row)
        alignment_column = str(alignment_row.get("alignment_column", "")).strip()
        if not protein_id or not residue_label or not alignment_column:
            continue

        condition_rows = condition_scores_by_residue.get((protein_id, residue_label), [])
        regio_row = regio_by_residue.get((protein_id, residue_label), {})
        if not condition_rows and not regio_row:
            continue

        family_label = _resolve_family_label(
            alignment_row,
            protein_metadata_by_id.get(protein_id),
        )
        mean_condition_scores = [
            numeric
            for row in condition_rows
            if (numeric := _as_float(row.get("residue_contact_score"))) is not None
        ]
        max_condition_score = max(mean_condition_scores, default="")

        residue_chain = str(alignment_row.get("residue_chain") or regio_row.get("residue_chain") or "")
        residue_number = _as_int(alignment_row.get("residue_number") or regio_row.get("residue_number"))
        residue_name = str(alignment_row.get("residue_name") or regio_row.get("residue_name") or "")

        family_aligned_rows.append(
            {
                "family_label": family_label,
                "family_aggregation_id": _family_aggregation_id(alignment_row),
                "alignment_column": alignment_column,
                "protein_id": protein_id,
                "construct_type": str(alignment_row.get("construct_type", "")).strip(),
                "catalytic_core_residue_index": _as_int(
                    alignment_row.get("catalytic_core_residue_index") or residue_number
                ),
                "alignment_residue": str(alignment_row.get("alignment_residue", "")).strip(),
                "residue_chain": residue_chain,
                "residue_number": residue_number,
                "residue_name": residue_name,
                "residue_label": residue_label,
                "n_conditions_with_valid_clusters": _as_int(regio_row.get("n_conditions_with_valid_clusters")),
                "n_conditions_with_contact": _as_int(regio_row.get("n_conditions_with_contact")),
                "mean_condition_residue_contact_score": _mean(mean_condition_scores),
                "max_condition_residue_contact_score": max_condition_score,
                "mean_c1_weighted_residue_score": _as_float(regio_row.get("mean_c1_weighted_residue_score"))
                if regio_row
                else "",
                "mean_c4_weighted_residue_score": _as_float(regio_row.get("mean_c4_weighted_residue_score"))
                if regio_row
                else "",
                "c1_minus_c4_weighted_delta": _as_float(regio_row.get("c1_minus_c4_weighted_delta"))
                if regio_row
                else "",
                "is_catalytic_surface_region": _as_bool(
                    alignment_row.get(
                        "is_catalytic_surface_region",
                        regio_row.get("is_catalytic_surface_region", False),
                    )
                ),
                "is_cbm_region": _as_bool(
                    alignment_row.get("is_cbm_region", regio_row.get("is_cbm_region", False))
                ),
                "is_linker_region": _as_bool(
                    alignment_row.get("is_linker_region", regio_row.get("is_linker_region", False))
                ),
            }
        )

    family_enrichment_rows: list[dict[str, Any]] = []
    aligned_rows_by_column: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in family_aligned_rows:
        aligned_rows_by_column[(str(row.get("family_label", "")), str(row.get("alignment_column", "")))].append(row)

    for family_label, alignment_column in sorted(aligned_rows_by_column):
        rows = aligned_rows_by_column[(family_label, alignment_column)]
        aggregation_unit_rows_by_id: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in rows:
            aggregation_unit_rows_by_id[_family_aggregation_id(row)].append(row)

        aggregation_unit_rows: list[dict[str, Any]] = []
        for aggregation_unit_id in sorted(aggregation_unit_rows_by_id):
            member_rows = aggregation_unit_rows_by_id[aggregation_unit_id]
            aggregation_unit_rows.append(
                {
                    "family_aggregation_id": aggregation_unit_id,
                    "mean_condition_residue_contact_score": _mean(
                        [
                            numeric
                            for row in member_rows
                            if (
                                numeric := _as_float(
                                    row.get("mean_condition_residue_contact_score")
                                )
                            ) is not None
                        ]
                    ),
                    "max_condition_residue_contact_score": max(
                        [
                            numeric
                            for row in member_rows
                            if (
                                numeric := _as_float(
                                    row.get("max_condition_residue_contact_score")
                                )
                            ) is not None
                        ],
                        default="",
                    ),
                    "mean_c1_weighted_residue_score": _mean(
                        [
                            numeric
                            for row in member_rows
                            if (
                                numeric := _as_float(
                                    row.get("mean_c1_weighted_residue_score")
                                )
                            ) is not None
                        ]
                    ),
                    "mean_c4_weighted_residue_score": _mean(
                        [
                            numeric
                            for row in member_rows
                            if (
                                numeric := _as_float(
                                    row.get("mean_c4_weighted_residue_score")
                                )
                            ) is not None
                        ]
                    ),
                    "mean_c1_minus_c4_weighted_delta": _mean(
                        [
                            numeric
                            for row in member_rows
                            if (
                                numeric := _as_float(
                                    row.get("c1_minus_c4_weighted_delta")
                                )
                            ) is not None
                        ]
                    ),
                }
            )

        mean_condition_scores = [
            numeric
            for row in aggregation_unit_rows
            if (numeric := _as_float(row.get("mean_condition_residue_contact_score"))) is not None
        ]
        max_condition_scores = [
            numeric
            for row in aggregation_unit_rows
            if (numeric := _as_float(row.get("max_condition_residue_contact_score"))) is not None
        ]
        c1_scores = [
            numeric
            for row in aggregation_unit_rows
            if (numeric := _as_float(row.get("mean_c1_weighted_residue_score"))) is not None
        ]
        c4_scores = [
            numeric
            for row in aggregation_unit_rows
            if (numeric := _as_float(row.get("mean_c4_weighted_residue_score"))) is not None
        ]
        deltas = [
            numeric
            for row in aggregation_unit_rows
            if (numeric := _as_float(row.get("mean_c1_minus_c4_weighted_delta"))) is not None
        ]

        family_enrichment_rows.append(
            {
                "family_label": family_label,
                "alignment_column": alignment_column,
                "n_proteins_observed": len(rows),
                "n_proteins_with_contact": sum(
                    1 for row in rows if (_as_float(row.get("mean_condition_residue_contact_score")) or 0.0) > 0.0
                ),
                "n_family_aggregation_units_observed": len(aggregation_unit_rows),
                "n_family_aggregation_units_with_contact": sum(
                    1
                    for row in aggregation_unit_rows
                    if (_as_float(row.get("mean_condition_residue_contact_score")) or 0.0) > 0.0
                ),
                "mean_condition_residue_contact_score": _mean(mean_condition_scores),
                "max_condition_residue_contact_score": max(max_condition_scores, default=""),
                "mean_c1_weighted_residue_score": _mean(c1_scores),
                "mean_c4_weighted_residue_score": _mean(c4_scores),
                "mean_c1_minus_c4_weighted_delta": _mean(deltas),
                "n_positive_c1_minus_c4_delta": sum(1 for value in deltas if value > 0.0),
                "n_negative_c1_minus_c4_delta": sum(1 for value in deltas if value < 0.0),
            }
        )

    family_aligned_rows.sort(
        key=lambda row: (
            str(row.get("family_label", "")),
            str(row.get("alignment_column", "")),
            str(row.get("family_aggregation_id", "")),
            str(row.get("protein_id", "")),
            str(row.get("residue_label", "")),
        )
    )

    return FamilyResidueEnrichmentOutputs(
        family_aligned_residue_table=family_aligned_rows,
        family_residue_enrichment=family_enrichment_rows,
    )