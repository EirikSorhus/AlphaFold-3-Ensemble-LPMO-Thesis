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

FAMILY_SUBSTRATE_RESIDUE_ENRICHMENT_COLUMNS = [
    "family_label",
    "alignment_column",
    "target_substrate",
    "n_proteins_observed",
    "n_proteins_with_target_substrate_rows",
    "n_proteins_with_other_substrate_rows",
    "n_family_aggregation_units_observed",
    "n_family_aggregation_units_with_target_substrate_rows",
    "n_family_aggregation_units_with_other_substrate_rows",
    "n_target_substrate_conditions",
    "n_other_substrate_conditions",
    "mean_target_substrate_residue_contact_score",
    "mean_other_substrate_residue_contact_score",
    "target_minus_other_substrate_residue_contact_delta",
    "n_positive_target_minus_other_delta",
    "n_negative_target_minus_other_delta",
    "substrate_group_status",
]

FAMILY_WRONG_LIGAND_RESIDUE_ENRICHMENT_COLUMNS = [
    "family_label",
    "alignment_column",
    "active_substrate",
    "wrong_prediction_substrate",
    "n_proteins_observed",
    "n_single_active_proteins_observed",
    "n_dual_active_proteins_excluded",
    "n_unknown_activity_proteins_excluded",
    "n_proteins_with_right_prediction_rows",
    "n_proteins_with_wrong_prediction_rows",
    "n_proteins_with_paired_right_wrong_rows",
    "n_family_aggregation_units_observed",
    "n_family_aggregation_units_with_paired_right_wrong_rows",
    "n_right_prediction_conditions",
    "n_wrong_prediction_conditions",
    "mean_right_ligand_residue_contact_score",
    "mean_wrong_ligand_residue_contact_score",
    "right_minus_wrong_ligand_residue_contact_delta",
    "n_positive_right_minus_wrong_delta",
    "n_negative_right_minus_wrong_delta",
    "wrong_ligand_group_status",
]


@dataclass(frozen=True)
class FamilyResidueEnrichmentOutputs:
    family_aligned_residue_table: list[dict[str, Any]]
    family_residue_enrichment: list[dict[str, Any]]
    family_substrate_residue_enrichment: list[dict[str, Any]]
    family_wrong_ligand_residue_enrichment: list[dict[str, Any]]


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


def write_family_substrate_residue_enrichment(
    rows: list[dict[str, Any]],
    output_path: Path,
) -> None:
    _write_tsv(output_path, FAMILY_SUBSTRATE_RESIDUE_ENRICHMENT_COLUMNS, rows)


def write_family_wrong_ligand_residue_enrichment(
    rows: list[dict[str, Any]],
    output_path: Path,
) -> None:
    _write_tsv(output_path, FAMILY_WRONG_LIGAND_RESIDUE_ENRICHMENT_COLUMNS, rows)


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


def _normalize_substrate(value: Any) -> str:
    return str(value or "").strip().lower()


def _first_present(row: dict[str, Any], keys: tuple[str, ...]) -> str:
    for key in keys:
        value = row.get(key)
        if value not in {None, "", "None"}:
            return str(value)
    return ""


def _ec_tokens(value: Any) -> list[str]:
    text = str(value or "")
    return [token.strip() for token in text.replace(";", " ").replace(",", " ").split() if token.strip()]


def _active_substrate_set(metadata_row: dict[str, Any] | None) -> set[str]:
    if metadata_row is None:
        return set()

    direct = _normalize_substrate(
        _first_present(
            metadata_row,
            (
                "mapped_substrate_class",
                "substrate_class",
                "substrate_label",
                "experimental_substrate_label",
            ),
        )
    )
    active: set[str] = set()
    if direct:
        if "cellulose" in direct:
            active.add("cellulose")
        if "chitin" in direct:
            active.add("chitin")
        if active:
            return active

    for token in _ec_tokens(_first_present(metadata_row, ("EC_Number", "ec_number", "ec_numbers"))):
        if token in {"1.14.99.54", "1.14.99.56"}:
            active.add("cellulose")
        elif token == "1.14.99.53":
            active.add("chitin")
    return active


def _substrate_from_condition_score_row(row: dict[str, Any]) -> str:
    explicit = _normalize_substrate(row.get("substrate_class"))
    if explicit:
        return explicit
    condition_id = str(row.get("condition_id", "")).strip()
    if "__" not in condition_id or "_DP" not in condition_id:
        return ""
    substrate_dp = condition_id.rsplit("__", 1)[-1]
    substrate = substrate_dp.rsplit("_DP", 1)[0]
    return _normalize_substrate(substrate)


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
        for key in ("protein_id", "uniprot_id", "UniProt_ID", "Uniprot_ID", "UniProtID"):
            protein_id = str(row.get(key) or "").strip()
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
    substrate_classes: tuple[str, ...] = ("cellulose", "chitin"),
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

    family_substrate_enrichment_rows = _build_family_substrate_enrichment_rows(
        family_aligned_rows=family_aligned_rows,
        condition_scores_by_residue=condition_scores_by_residue,
        substrate_classes=tuple(_normalize_substrate(value) for value in substrate_classes),
    )
    family_wrong_ligand_enrichment_rows = _build_family_wrong_ligand_enrichment_rows(
        family_aligned_rows=family_aligned_rows,
        condition_scores_by_residue=condition_scores_by_residue,
        protein_metadata_by_id=protein_metadata_by_id,
        substrate_classes=tuple(_normalize_substrate(value) for value in substrate_classes),
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
        family_substrate_residue_enrichment=family_substrate_enrichment_rows,
        family_wrong_ligand_residue_enrichment=family_wrong_ligand_enrichment_rows,
    )


def _substrate_group_status(
    *,
    n_target_units: int,
    n_other_units: int,
    n_target_conditions: int,
    n_other_conditions: int,
) -> str:
    statuses: list[str] = []
    if n_target_units == 0 or n_target_conditions == 0:
        statuses.append("missing_target_substrate")
    if n_other_units == 0 or n_other_conditions == 0:
        statuses.append("target_only_no_contrast")
    if 0 < min(n_target_units, n_other_units) < 3:
        statuses.append("small_group")
    if min(n_target_units, n_other_units) > 0 and max(n_target_units, n_other_units) / min(
        n_target_units,
        n_other_units,
    ) >= 3:
        statuses.append("unbalanced_groups")
    return ";".join(statuses) if statuses else "ok"


def _build_family_substrate_enrichment_rows(
    *,
    family_aligned_rows: list[dict[str, Any]],
    condition_scores_by_residue: dict[tuple[str, str], list[dict[str, Any]]],
    substrate_classes: tuple[str, ...],
) -> list[dict[str, Any]]:
    rows_by_column: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in family_aligned_rows:
        rows_by_column[(str(row.get("family_label", "")), str(row.get("alignment_column", "")))].append(row)

    out: list[dict[str, Any]] = []
    for family_label, alignment_column in sorted(rows_by_column):
        aligned_rows = rows_by_column[(family_label, alignment_column)]
        for target_substrate in substrate_classes:
            if not target_substrate:
                continue

            proteins_observed: set[str] = set()
            proteins_with_target: set[str] = set()
            proteins_with_other: set[str] = set()
            target_condition_ids: set[str] = set()
            other_condition_ids: set[str] = set()
            unit_rows: dict[str, dict[str, Any]] = {}

            for aligned_row in aligned_rows:
                protein_id = str(aligned_row.get("protein_id", "")).strip()
                residue_label = str(aligned_row.get("residue_label", "")).strip()
                unit_id = _family_aggregation_id(aligned_row)
                if protein_id:
                    proteins_observed.add(protein_id)
                score_rows = condition_scores_by_residue.get((protein_id, residue_label), [])
                target_scores: list[float] = []
                other_scores: list[float] = []
                for score_row in score_rows:
                    score = _as_float(score_row.get("residue_contact_score"))
                    if score is None:
                        continue
                    condition_id = str(score_row.get("condition_id", "")).strip()
                    row_substrate = _substrate_from_condition_score_row(score_row)
                    if row_substrate == target_substrate:
                        target_scores.append(score)
                        proteins_with_target.add(protein_id)
                        if condition_id:
                            target_condition_ids.add(condition_id)
                    elif row_substrate:
                        other_scores.append(score)
                        proteins_with_other.add(protein_id)
                        if condition_id:
                            other_condition_ids.add(condition_id)

                entry = unit_rows.setdefault(
                    unit_id,
                    {
                        "target_scores": [],
                        "other_scores": [],
                    },
                )
                if target_scores:
                    entry["target_scores"].append(_mean(target_scores))
                if other_scores:
                    entry["other_scores"].append(_mean(other_scores))

            unit_summaries: list[dict[str, Any]] = []
            for unit_id in sorted(unit_rows):
                target_unit_scores = [
                    value
                    for value in unit_rows[unit_id]["target_scores"]
                    if isinstance(value, (int, float))
                ]
                other_unit_scores = [
                    value
                    for value in unit_rows[unit_id]["other_scores"]
                    if isinstance(value, (int, float))
                ]
                target_mean = _mean(target_unit_scores)
                other_mean = _mean(other_unit_scores)
                delta = (
                    target_mean - other_mean
                    if isinstance(target_mean, (int, float)) and isinstance(other_mean, (int, float))
                    else ""
                )
                unit_summaries.append(
                    {
                        "family_aggregation_id": unit_id,
                        "target_mean": target_mean,
                        "other_mean": other_mean,
                        "delta": delta,
                    }
                )

            target_means = [
                value
                for row in unit_summaries
                if isinstance(value := row.get("target_mean"), (int, float))
            ]
            other_means = [
                value
                for row in unit_summaries
                if isinstance(value := row.get("other_mean"), (int, float))
            ]
            deltas = [
                value
                for row in unit_summaries
                if isinstance(value := row.get("delta"), (int, float))
            ]
            units_with_target = {
                str(row["family_aggregation_id"])
                for row in unit_summaries
                if isinstance(row.get("target_mean"), (int, float))
            }
            units_with_other = {
                str(row["family_aggregation_id"])
                for row in unit_summaries
                if isinstance(row.get("other_mean"), (int, float))
            }
            mean_target = _mean(target_means)
            mean_other = _mean(other_means)
            group_delta = (
                mean_target - mean_other
                if isinstance(mean_target, (int, float)) and isinstance(mean_other, (int, float))
                else ""
            )

            out.append(
                {
                    "family_label": family_label,
                    "alignment_column": alignment_column,
                    "target_substrate": target_substrate,
                    "n_proteins_observed": len(proteins_observed),
                    "n_proteins_with_target_substrate_rows": len(proteins_with_target),
                    "n_proteins_with_other_substrate_rows": len(proteins_with_other),
                    "n_family_aggregation_units_observed": len(unit_rows),
                    "n_family_aggregation_units_with_target_substrate_rows": len(units_with_target),
                    "n_family_aggregation_units_with_other_substrate_rows": len(units_with_other),
                    "n_target_substrate_conditions": len(target_condition_ids),
                    "n_other_substrate_conditions": len(other_condition_ids),
                    "mean_target_substrate_residue_contact_score": mean_target,
                    "mean_other_substrate_residue_contact_score": mean_other,
                    "target_minus_other_substrate_residue_contact_delta": group_delta,
                    "n_positive_target_minus_other_delta": sum(1 for value in deltas if value > 0.0),
                    "n_negative_target_minus_other_delta": sum(1 for value in deltas if value < 0.0),
                    "substrate_group_status": _substrate_group_status(
                        n_target_units=len(units_with_target),
                        n_other_units=len(units_with_other),
                        n_target_conditions=len(target_condition_ids),
                        n_other_conditions=len(other_condition_ids),
                    ),
                }
            )
    return out


def _wrong_ligand_status(
    *,
    n_paired_units: int,
    n_right_conditions: int,
    n_wrong_conditions: int,
    n_dual_active_excluded: int,
    n_unknown_excluded: int,
) -> str:
    statuses: list[str] = []
    if n_paired_units == 0:
        statuses.append("no_paired_right_wrong_predictions")
    if n_right_conditions == 0:
        statuses.append("missing_right_ligand_predictions")
    if n_wrong_conditions == 0:
        statuses.append("missing_wrong_ligand_predictions")
    if 0 < n_paired_units < 3:
        statuses.append("small_paired_group")
    if n_dual_active_excluded:
        statuses.append("dual_active_excluded")
    if n_unknown_excluded:
        statuses.append("unknown_activity_excluded")
    return ";".join(statuses) if statuses else "ok"


def _build_family_wrong_ligand_enrichment_rows(
    *,
    family_aligned_rows: list[dict[str, Any]],
    condition_scores_by_residue: dict[tuple[str, str], list[dict[str, Any]]],
    protein_metadata_by_id: dict[str, dict[str, Any]],
    substrate_classes: tuple[str, ...],
) -> list[dict[str, Any]]:
    requested = set(substrate_classes)
    if not {"cellulose", "chitin"}.issubset(requested):
        return []

    rows_by_column: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in family_aligned_rows:
        rows_by_column[(str(row.get("family_label", "")), str(row.get("alignment_column", "")))].append(row)

    out: list[dict[str, Any]] = []
    for family_label, alignment_column in sorted(rows_by_column):
        aligned_rows = rows_by_column[(family_label, alignment_column)]
        for active_substrate, wrong_substrate in (
            ("cellulose", "chitin"),
            ("chitin", "cellulose"),
        ):
            proteins_observed: set[str] = set()
            single_active_proteins: set[str] = set()
            dual_active_excluded: set[str] = set()
            unknown_activity_excluded: set[str] = set()
            proteins_with_right: set[str] = set()
            proteins_with_wrong: set[str] = set()
            paired_proteins: set[str] = set()
            right_condition_ids: set[str] = set()
            wrong_condition_ids: set[str] = set()
            single_active_unit_ids: set[str] = set()
            unit_rows: dict[str, dict[str, Any]] = {}

            for aligned_row in aligned_rows:
                protein_id = str(aligned_row.get("protein_id", "")).strip()
                residue_label = str(aligned_row.get("residue_label", "")).strip()
                unit_id = _family_aggregation_id(aligned_row)
                if protein_id:
                    proteins_observed.add(protein_id)
                active = _active_substrate_set(protein_metadata_by_id.get(protein_id))
                if active == {"cellulose", "chitin"}:
                    dual_active_excluded.add(protein_id)
                    continue
                if active != {active_substrate}:
                    if not active:
                        unknown_activity_excluded.add(protein_id)
                    continue
                single_active_proteins.add(protein_id)
                single_active_unit_ids.add(unit_id)

                right_scores: list[float] = []
                wrong_scores: list[float] = []
                for score_row in condition_scores_by_residue.get((protein_id, residue_label), []):
                    score = _as_float(score_row.get("residue_contact_score"))
                    if score is None:
                        continue
                    condition_id = str(score_row.get("condition_id", "")).strip()
                    row_substrate = _substrate_from_condition_score_row(score_row)
                    if row_substrate == active_substrate:
                        right_scores.append(score)
                        proteins_with_right.add(protein_id)
                        if condition_id:
                            right_condition_ids.add(condition_id)
                    elif row_substrate == wrong_substrate:
                        wrong_scores.append(score)
                        proteins_with_wrong.add(protein_id)
                        if condition_id:
                            wrong_condition_ids.add(condition_id)

                if right_scores and wrong_scores:
                    paired_proteins.add(protein_id)
                    entry = unit_rows.setdefault(
                        unit_id,
                        {
                            "right_scores": [],
                            "wrong_scores": [],
                        },
                    )
                    entry["right_scores"].append(_mean(right_scores))
                    entry["wrong_scores"].append(_mean(wrong_scores))

            unit_summaries: list[dict[str, Any]] = []
            for unit_id in sorted(unit_rows):
                right_mean = _mean(
                    [
                        value
                        for value in unit_rows[unit_id]["right_scores"]
                        if isinstance(value, (int, float))
                    ]
                )
                wrong_mean = _mean(
                    [
                        value
                        for value in unit_rows[unit_id]["wrong_scores"]
                        if isinstance(value, (int, float))
                    ]
                )
                delta = (
                    right_mean - wrong_mean
                    if isinstance(right_mean, (int, float)) and isinstance(wrong_mean, (int, float))
                    else ""
                )
                unit_summaries.append(
                    {
                        "family_aggregation_id": unit_id,
                        "right_mean": right_mean,
                        "wrong_mean": wrong_mean,
                        "delta": delta,
                    }
                )

            right_means = [
                value
                for row in unit_summaries
                if isinstance(value := row.get("right_mean"), (int, float))
            ]
            wrong_means = [
                value
                for row in unit_summaries
                if isinstance(value := row.get("wrong_mean"), (int, float))
            ]
            mean_right = _mean(right_means)
            mean_wrong = _mean(wrong_means)
            group_delta = (
                mean_right - mean_wrong
                if isinstance(mean_right, (int, float)) and isinstance(mean_wrong, (int, float))
                else ""
            )
            deltas = [
                value
                for row in unit_summaries
                if isinstance(value := row.get("delta"), (int, float))
            ]

            out.append(
                {
                    "family_label": family_label,
                    "alignment_column": alignment_column,
                    "active_substrate": active_substrate,
                    "wrong_prediction_substrate": wrong_substrate,
                    "n_proteins_observed": len(proteins_observed),
                    "n_single_active_proteins_observed": len(single_active_proteins),
                    "n_dual_active_proteins_excluded": len(dual_active_excluded),
                    "n_unknown_activity_proteins_excluded": len(unknown_activity_excluded),
                    "n_proteins_with_right_prediction_rows": len(proteins_with_right),
                    "n_proteins_with_wrong_prediction_rows": len(proteins_with_wrong),
                    "n_proteins_with_paired_right_wrong_rows": len(paired_proteins),
                    "n_family_aggregation_units_observed": len(single_active_unit_ids),
                    "n_family_aggregation_units_with_paired_right_wrong_rows": len(unit_rows),
                    "n_right_prediction_conditions": len(right_condition_ids),
                    "n_wrong_prediction_conditions": len(wrong_condition_ids),
                    "mean_right_ligand_residue_contact_score": mean_right,
                    "mean_wrong_ligand_residue_contact_score": mean_wrong,
                    "right_minus_wrong_ligand_residue_contact_delta": group_delta,
                    "n_positive_right_minus_wrong_delta": sum(1 for value in deltas if value > 0.0),
                    "n_negative_right_minus_wrong_delta": sum(1 for value in deltas if value < 0.0),
                    "wrong_ligand_group_status": _wrong_ligand_status(
                        n_paired_units=len(unit_rows),
                        n_right_conditions=len(right_condition_ids),
                        n_wrong_conditions=len(wrong_condition_ids),
                        n_dual_active_excluded=len(dual_active_excluded),
                        n_unknown_excluded=len(unknown_activity_excluded),
                    ),
                }
            )
    return out
