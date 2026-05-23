"""Stage 16b residue-importance aggregation from Stage 16 cluster signatures."""
from __future__ import annotations

import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any


PROTEIN_CONDITION_RESIDUE_SCORE_COLUMNS = [
    "protein_id",
    "condition_id",
    "construct_type",
    "substrate_class",
    "dp",
    "residue_chain",
    "residue_number",
    "residue_name",
    "residue_label",
    "residue_contact_score",
    "c1_weighted_residue_score",
    "c4_weighted_residue_score",
    "c1_minus_c4_weighted_delta",
    "cluster_support_count",
    "total_cluster_occupancy_with_contact",
    "max_cluster_residue_frequency",
    "is_catalytic_surface_region",
    "is_non_core_region",
    "is_cbm_region",
    "is_linker_region",
]

PROTEIN_RESIDUE_REGIO_DELTA_COLUMNS = [
    "protein_id",
    "residue_chain",
    "residue_number",
    "residue_name",
    "residue_label",
    "n_conditions_with_valid_clusters",
    "n_conditions_with_contact",
    "mean_c1_weighted_residue_score",
    "mean_c4_weighted_residue_score",
    "c1_minus_c4_weighted_delta",
    "is_catalytic_surface_region",
    "is_non_core_region",
    "is_cbm_region",
    "is_linker_region",
]

CONDITION_PATCH_SUMMARY_COLUMNS = [
    "protein_id",
    "condition_id",
    "construct_type",
    "substrate_class",
    "dp",
    "any_valid_cluster",
    "n_clusters_considered",
    "total_nonnoise_cluster_occupancy",
    "weighted_contact_feature_mass",
    "aromatic_contact_fraction",
    "polar_contact_fraction",
    "charged_contact_fraction",
    "hbond_contact_fraction",
    "catalytic_surface_contact_fraction",
    "non_core_contact_fraction",
    "cbm_contact_fraction",
    "linker_contact_fraction",
]

PROTEIN_PATCH_SUMMARY_COLUMNS = [
    "protein_id",
    "n_conditions_total",
    "n_conditions_with_valid_clusters",
    "mean_aromatic_contact_fraction",
    "mean_polar_contact_fraction",
    "mean_charged_contact_fraction",
    "mean_hbond_contact_fraction",
    "mean_catalytic_surface_contact_fraction",
    "mean_non_core_contact_fraction",
    "mean_cbm_contact_fraction",
    "mean_linker_contact_fraction",
]

_HBOND_INTERACTIONS = {"ImplicitHBAcceptor", "ImplicitHBDonor"}
_AROMATIC_RESIDUES = {"PHE", "TRP", "TYR", "HIS"}
_POLAR_RESIDUES = {"SER", "THR", "ASN", "GLN", "CYS"}
_CHARGED_RESIDUES = {"ASP", "GLU", "LYS", "ARG"}
_CONDITION_PATTERN = re.compile(
    r"^(?P<protein_id>.+?)__(?P<construct_type>.+?)__(?P<substrate_class>.+?)_DP(?P<dp>-?\d+)$"
)
_PROTEIN_RESIDUE_PATTERN = re.compile(
    r"^(?P<residue_name>[A-Za-z0-9]+?)(?P<residue_number>-?\d+)(?:\.(?P<residue_chain>[A-Za-z0-9]+))?$"
)


@dataclass(frozen=True)
class ResidueImportanceOutputs:
    protein_condition_residue_scores: list[dict[str, Any]]
    protein_residue_regio_delta: list[dict[str, Any]]
    condition_patch_summary: list[dict[str, Any]]
    protein_patch_summary: list[dict[str, Any]]


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes"}
    if isinstance(value, (int, float)):
        return bool(value)
    return False


def _as_int(value: Any, default: int = 0) -> int:
    if value in {None, "", "None"}:
        return default
    return int(value)


def _as_float(value: Any, default: float = 0.0) -> float:
    if value in {None, "", "None"}:
        return default
    return float(value)


def _parse_condition_id(
    condition_id: str,
    explicit_metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    match = _CONDITION_PATTERN.match(condition_id)
    parsed = {
        "protein_id": match.group("protein_id") if match is not None else "",
        "construct_type": match.group("construct_type") if match is not None else "",
        "substrate_class": match.group("substrate_class") if match is not None else "",
        "dp": int(match.group("dp")) if match is not None else 0,
    }

    if explicit_metadata is None:
        return parsed

    return {
        "protein_id": str(explicit_metadata.get("protein_id") or parsed["protein_id"]),
        "construct_type": str(explicit_metadata.get("construct_type") or parsed["construct_type"]),
        "substrate_class": str(explicit_metadata.get("substrate_class") or parsed["substrate_class"]),
        "dp": _as_int(explicit_metadata.get("dp"), default=parsed["dp"]),
    }


def _residue_key(row: dict[str, Any]) -> tuple[str, int, str]:
    return (
        str(row.get("residue_chain", "")),
        _as_int(row.get("residue_number"), default=0),
        str(row.get("residue_name", "")),
    )


def _residue_label(residue_key: tuple[str, int, str]) -> str:
    residue_chain, residue_number, residue_name = residue_key
    return f"{residue_name}{residue_number}.{residue_chain}"


def _merge_region_flags(existing: dict[str, bool], row: dict[str, Any]) -> dict[str, bool]:
    return {
        "is_catalytic_surface_region": existing["is_catalytic_surface_region"]
        or _as_bool(row.get("is_catalytic_surface_region", False)),
        "is_non_core_region": existing["is_non_core_region"]
        or _as_bool(row.get("is_non_core_region", False))
        or _as_bool(row.get("is_cbm_region", False))
        or _as_bool(row.get("is_linker_region", False)),
        "is_cbm_region": existing["is_cbm_region"] or _as_bool(row.get("is_cbm_region", False)),
        "is_linker_region": existing["is_linker_region"]
        or _as_bool(row.get("is_linker_region", False)),
    }


def _empty_region_flags() -> dict[str, bool]:
    return {
        "is_catalytic_surface_region": False,
        "is_non_core_region": False,
        "is_cbm_region": False,
        "is_linker_region": False,
    }


def _empty_patch_totals() -> dict[str, float]:
    return {
        "total": 0.0,
        "aromatic": 0.0,
        "polar": 0.0,
        "charged": 0.0,
        "hbond": 0.0,
        "catalytic_surface": 0.0,
        "non_core": 0.0,
        "cbm": 0.0,
        "linker": 0.0,
    }


def _empty_region_totals() -> dict[str, float]:
    return {
        "total": 0.0,
        "catalytic_surface": 0.0,
        "non_core": 0.0,
        "cbm": 0.0,
        "linker": 0.0,
    }


def _condition_protein_id(
    condition_id: str,
    parsed_metadata: dict[str, Any],
    condition_protein_id_by_id: dict[str, str],
) -> str:
    return parsed_metadata["protein_id"] or str(condition_protein_id_by_id.get(condition_id, ""))


def _new_residue_score_row(
    *,
    protein_id: str,
    condition_id: str,
    parsed_metadata: dict[str, Any],
    residue_key: tuple[str, int, str],
    flags: dict[str, bool],
) -> dict[str, Any]:
    return {
        "protein_id": protein_id,
        "condition_id": condition_id,
        "construct_type": parsed_metadata["construct_type"],
        "substrate_class": parsed_metadata["substrate_class"],
        "dp": parsed_metadata["dp"],
        "residue_chain": residue_key[0],
        "residue_number": residue_key[1],
        "residue_name": residue_key[2],
        "residue_label": _residue_label(residue_key),
        "residue_contact_score": 0.0,
        "c1_weighted_residue_score": 0.0,
        "c4_weighted_residue_score": 0.0,
        "cluster_support_count": 0,
        "total_cluster_occupancy_with_contact": 0.0,
        "max_cluster_residue_frequency": 0.0,
        **flags,
    }


def _coerce_cluster_summary_rows(cluster_summary_rows: list[dict[str, Any]] | dict[str, Any]) -> list[dict[str, Any]]:
    if isinstance(cluster_summary_rows, dict):
        return list(cluster_summary_rows.get("clusters", []))
    return list(cluster_summary_rows)


def _cluster_key(row: dict[str, Any]) -> tuple[str, int]:
    return str(row.get("condition_id", "")), _as_int(row.get("cluster_id"), default=-1)


def _parse_protein_residue_label(protein_residue_label: str) -> tuple[str, int, str]:
    match = _PROTEIN_RESIDUE_PATTERN.match(protein_residue_label)
    if match is None or match.group("residue_chain") is None:
        return "", 0, protein_residue_label
    return (
        match.group("residue_chain"),
        int(match.group("residue_number")),
        match.group("residue_name"),
    )


def _cluster_summaries_by_key(
    cluster_summary_rows: list[dict[str, Any]] | dict[str, Any],
) -> dict[tuple[str, int], dict[str, Any]]:
    summaries: dict[tuple[str, int], dict[str, Any]] = {}
    for row in _coerce_cluster_summary_rows(cluster_summary_rows):
        key = _cluster_key(row)
        if not key[0] or key[1] < 0:
            continue
        summaries[key] = row
    return summaries


def compute_residue_importance_outputs(
    cluster_residue_signature_rows: list[dict[str, Any]],
    cluster_summary_rows: list[dict[str, Any]] | dict[str, Any],
    *,
    cluster_ifp_signature_rows: list[dict[str, Any]] | None = None,
    observed_residue_contact_rows: list[dict[str, Any]] | None = None,
    condition_metadata_by_id: dict[str, dict[str, Any]] | None = None,
) -> ResidueImportanceOutputs:
    """Build Stage 16b outputs from Stage 16 cluster annotation surfaces.

    `cluster_residue_signature_rows` and `cluster_ifp_signature_rows` should be
    the rows written by Stage 16. `observed_residue_contact_rows` is only used to
    keep explicit zero-score residue rows for conditions with no valid clusters.
    """

    condition_metadata_by_id = condition_metadata_by_id or {}
    cluster_ifp_signature_rows = cluster_ifp_signature_rows or []
    observed_residue_contact_rows = observed_residue_contact_rows or []
    cluster_summary_by_key = _cluster_summaries_by_key(cluster_summary_rows)
    condition_ids: set[str] = set(condition_metadata_by_id)
    condition_protein_id_by_id: dict[str, str] = {}
    observed_residue_keys_by_condition: dict[str, set[tuple[str, int, str]]] = defaultdict(set)
    observed_residue_flags: dict[tuple[str, tuple[str, int, str]], dict[str, bool]] = {}

    for row in cluster_summary_by_key.values():
        condition_id = str(row.get("condition_id", ""))
        protein_id = str(row.get("protein_id", ""))
        if condition_id:
            condition_ids.add(condition_id)
        if condition_id and protein_id and condition_id not in condition_protein_id_by_id:
            condition_protein_id_by_id[condition_id] = protein_id

    for row in observed_residue_contact_rows:
        condition_id = str(row.get("condition_id", ""))
        protein_id = str(row.get("protein_id", ""))
        if condition_id:
            condition_ids.add(condition_id)
        if condition_id and protein_id and condition_id not in condition_protein_id_by_id:
            condition_protein_id_by_id[condition_id] = protein_id
        if condition_id and _as_bool(row.get("contact_present", False)):
            residue_key = _residue_key(row)
            observed_residue_keys_by_condition[condition_id].add(residue_key)
            observed_key = (condition_id, residue_key)
            observed_residue_flags[observed_key] = _merge_region_flags(
                observed_residue_flags.get(observed_key, _empty_region_flags()),
                row,
            )

    residue_scores_by_condition: dict[tuple[str, tuple[str, int, str]], dict[str, Any]] = {}
    condition_region_weighted_totals: dict[str, dict[str, float]] = defaultdict(_empty_region_totals)

    for row in cluster_residue_signature_rows:
        condition_id = str(row.get("condition_id", ""))
        cluster_id = _as_int(row.get("cluster_id"), default=-1)
        condition_ids.add(condition_id)
        if not condition_id or cluster_id < 0:
            continue
        cluster_summary = cluster_summary_by_key.get((condition_id, cluster_id), {})
        occupancy = _as_float(row.get("occupancy"), default=_as_float(cluster_summary.get("occupancy")))
        residue_frequency = _as_float(row.get("contact_frequency"), default=0.0)
        c1_fraction = _as_float(cluster_summary.get("c1_plausible_fraction"), default=0.0)
        c4_fraction = _as_float(cluster_summary.get("c4_plausible_fraction"), default=0.0)
        protein_id = str(row.get("protein_id", ""))
        if condition_id and protein_id and condition_id not in condition_protein_id_by_id:
            condition_protein_id_by_id[condition_id] = protein_id

        residue_key = _residue_key(row)
        aggregate_key = (condition_id, residue_key)

        aggregate = residue_scores_by_condition.get(aggregate_key)
        if aggregate is None:
            explicit_metadata = condition_metadata_by_id.get(condition_id)
            parsed_metadata = _parse_condition_id(condition_id, explicit_metadata)
            flags = _merge_region_flags(_empty_region_flags(), row)
            aggregate = _new_residue_score_row(
                protein_id=_condition_protein_id(
                    condition_id,
                    parsed_metadata,
                    condition_protein_id_by_id,
                ),
                condition_id=condition_id,
                parsed_metadata=parsed_metadata,
                residue_key=residue_key,
                flags=flags,
            )
            residue_scores_by_condition[aggregate_key] = aggregate

        aggregate["residue_contact_score"] += occupancy * residue_frequency
        aggregate["c1_weighted_residue_score"] += occupancy * residue_frequency * c1_fraction
        aggregate["c4_weighted_residue_score"] += occupancy * residue_frequency * c4_fraction
        aggregate["cluster_support_count"] += 1
        aggregate["total_cluster_occupancy_with_contact"] += occupancy
        aggregate["max_cluster_residue_frequency"] = max(
            aggregate["max_cluster_residue_frequency"],
            residue_frequency,
        )

        weighted_residue_contact = occupancy * residue_frequency
        if weighted_residue_contact > 0.0:
            region_totals = condition_region_weighted_totals[condition_id]
            region_totals["total"] += weighted_residue_contact
            if _as_bool(row.get("is_catalytic_surface_region", False)):
                region_totals["catalytic_surface"] += weighted_residue_contact
            if (
                _as_bool(row.get("is_non_core_region", False))
                or _as_bool(row.get("is_cbm_region", False))
                or _as_bool(row.get("is_linker_region", False))
            ):
                region_totals["non_core"] += weighted_residue_contact
            if _as_bool(row.get("is_cbm_region", False)):
                region_totals["cbm"] += weighted_residue_contact
            if _as_bool(row.get("is_linker_region", False)):
                region_totals["linker"] += weighted_residue_contact

    condition_patch_weighted_totals: dict[str, dict[str, float]] = defaultdict(_empty_patch_totals)
    for row in cluster_ifp_signature_rows:
        condition_id = str(row.get("condition_id", ""))
        cluster_id = _as_int(row.get("cluster_id"), default=-1)
        if not condition_id or cluster_id < 0:
            continue
        cluster_summary = cluster_summary_by_key.get((condition_id, cluster_id), {})
        occupancy = _as_float(row.get("occupancy"), default=_as_float(cluster_summary.get("occupancy")))
        contact_frequency = _as_float(row.get("contact_frequency"), default=0.0)
        weighted_contact = occupancy * contact_frequency
        if weighted_contact <= 0.0:
            continue
        protein_residue_label = str(row.get("protein_residue_label", ""))
        _, _, residue_name = _parse_protein_residue_label(protein_residue_label)
        residue_name = residue_name.upper()
        interaction_type = str(row.get("interaction_type", ""))
        weighted_totals = condition_patch_weighted_totals[condition_id]
        weighted_totals["total"] += weighted_contact
        if residue_name in _AROMATIC_RESIDUES:
            weighted_totals["aromatic"] += weighted_contact
        if residue_name in _POLAR_RESIDUES:
            weighted_totals["polar"] += weighted_contact
        if residue_name in _CHARGED_RESIDUES:
            weighted_totals["charged"] += weighted_contact
        if interaction_type in _HBOND_INTERACTIONS:
            weighted_totals["hbond"] += weighted_contact

    for condition_id, residue_keys in observed_residue_keys_by_condition.items():
        if any(key[0] == condition_id for key in cluster_summary_by_key):
            continue
        explicit_metadata = condition_metadata_by_id.get(condition_id)
        parsed_metadata = _parse_condition_id(condition_id, explicit_metadata)
        protein_id = _condition_protein_id(condition_id, parsed_metadata, condition_protein_id_by_id)
        for residue_key in sorted(
            residue_keys,
            key=lambda item: (str(item[0]), int(item[1]), str(item[2])),
        ):
            aggregate_key = (condition_id, residue_key)
            if aggregate_key in residue_scores_by_condition:
                continue
            residue_scores_by_condition[aggregate_key] = _new_residue_score_row(
                protein_id=protein_id,
                condition_id=condition_id,
                parsed_metadata=parsed_metadata,
                residue_key=residue_key,
                flags=observed_residue_flags.get((condition_id, residue_key), _empty_region_flags()),
            )

    protein_condition_residue_scores = []
    for row in sorted(
        residue_scores_by_condition.values(),
        key=lambda item: (
            str(item["protein_id"]),
            str(item["condition_id"]),
            str(item["residue_chain"]),
            int(item["residue_number"]),
            str(item["residue_name"]),
        ),
    ):
        protein_condition_residue_scores.append(
            {
                **row,
                "c1_minus_c4_weighted_delta": row["c1_weighted_residue_score"] - row["c4_weighted_residue_score"],
            }
        )

    condition_patch_summary = []
    for condition_id in sorted(condition_ids):
        explicit_metadata = condition_metadata_by_id.get(condition_id)
        parsed_metadata = _parse_condition_id(condition_id, explicit_metadata)
        weighted_totals = condition_patch_weighted_totals.get(condition_id, _empty_patch_totals())
        total_mass = weighted_totals["total"]
        condition_cluster_summaries = [
            row for key, row in cluster_summary_by_key.items() if key[0] == condition_id
        ]
        any_valid_cluster = bool(condition_cluster_summaries)

        def _fraction(key: str) -> float:
            return (weighted_totals[key] / total_mass) if total_mass > 0 else 0.0

        region_totals = condition_region_weighted_totals.get(condition_id, _empty_region_totals())
        region_total_mass = region_totals["total"]

        def _region_fraction(key: str) -> float:
            return (region_totals[key] / region_total_mass) if region_total_mass > 0 else 0.0

        condition_patch_summary.append(
            {
                "protein_id": _condition_protein_id(
                    condition_id,
                    parsed_metadata,
                    condition_protein_id_by_id,
                ),
                "condition_id": condition_id,
                "construct_type": parsed_metadata["construct_type"],
                "substrate_class": parsed_metadata["substrate_class"],
                "dp": parsed_metadata["dp"],
                "any_valid_cluster": any_valid_cluster,
                "n_clusters_considered": len(condition_cluster_summaries),
                "total_nonnoise_cluster_occupancy": sum(
                    _as_float(row.get("occupancy")) for row in condition_cluster_summaries
                ),
                "weighted_contact_feature_mass": total_mass,
                "aromatic_contact_fraction": _fraction("aromatic"),
                "polar_contact_fraction": _fraction("polar"),
                "charged_contact_fraction": _fraction("charged"),
                "hbond_contact_fraction": _fraction("hbond"),
                "catalytic_surface_contact_fraction": _region_fraction("catalytic_surface"),
                "non_core_contact_fraction": _region_fraction("non_core"),
                "cbm_contact_fraction": _region_fraction("cbm"),
                "linker_contact_fraction": _region_fraction("linker"),
            }
        )

    valid_condition_counts_by_protein: dict[str, int] = defaultdict(int)
    all_condition_counts_by_protein: dict[str, int] = defaultdict(int)
    protein_patch_rows: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in condition_patch_summary:
        protein_id = str(row["protein_id"])
        all_condition_counts_by_protein[protein_id] += 1
        if _as_bool(row["any_valid_cluster"]):
            valid_condition_counts_by_protein[protein_id] += 1
        protein_patch_rows[protein_id].append(row)

    protein_residue_accumulator: dict[tuple[str, tuple[str, int, str]], dict[str, Any]] = {}
    for row in protein_condition_residue_scores:
        protein_id = str(row["protein_id"])
        residue_key = (
            str(row["residue_chain"]),
            _as_int(row["residue_number"]),
            str(row["residue_name"]),
        )
        aggregate_key = (protein_id, residue_key)
        aggregate = protein_residue_accumulator.get(aggregate_key)
        if aggregate is None:
            aggregate = {
                "protein_id": protein_id,
                "residue_chain": residue_key[0],
                "residue_number": residue_key[1],
                "residue_name": residue_key[2],
                "residue_label": _residue_label(residue_key),
                "n_conditions_with_contact": 0,
                "sum_c1_weighted_residue_score": 0.0,
                "sum_c4_weighted_residue_score": 0.0,
                "is_catalytic_surface_region": _as_bool(row["is_catalytic_surface_region"]),
                "is_non_core_region": _as_bool(row.get("is_non_core_region", False)),
                "is_cbm_region": _as_bool(row["is_cbm_region"]),
                "is_linker_region": _as_bool(row["is_linker_region"]),
            }
            protein_residue_accumulator[aggregate_key] = aggregate

        if _as_float(row.get("residue_contact_score")) > 0.0:
            aggregate["n_conditions_with_contact"] += 1
        aggregate["sum_c1_weighted_residue_score"] += _as_float(row["c1_weighted_residue_score"])
        aggregate["sum_c4_weighted_residue_score"] += _as_float(row["c4_weighted_residue_score"])

    protein_residue_regio_delta = []
    for aggregate in sorted(
        protein_residue_accumulator.values(),
        key=lambda item: (
            str(item["protein_id"]),
            str(item["residue_chain"]),
            int(item["residue_number"]),
            str(item["residue_name"]),
        ),
    ):
        protein_id = str(aggregate["protein_id"])
        n_valid_conditions = valid_condition_counts_by_protein.get(protein_id, 0)
        mean_c1 = (
            aggregate["sum_c1_weighted_residue_score"] / n_valid_conditions if n_valid_conditions else 0.0
        )
        mean_c4 = (
            aggregate["sum_c4_weighted_residue_score"] / n_valid_conditions if n_valid_conditions else 0.0
        )
        protein_residue_regio_delta.append(
            {
                "protein_id": protein_id,
                "residue_chain": aggregate["residue_chain"],
                "residue_number": aggregate["residue_number"],
                "residue_name": aggregate["residue_name"],
                "residue_label": aggregate["residue_label"],
                "n_conditions_with_valid_clusters": n_valid_conditions,
                "n_conditions_with_contact": aggregate["n_conditions_with_contact"],
                "mean_c1_weighted_residue_score": mean_c1,
                "mean_c4_weighted_residue_score": mean_c4,
                "c1_minus_c4_weighted_delta": mean_c1 - mean_c4,
                "is_catalytic_surface_region": aggregate["is_catalytic_surface_region"],
                "is_non_core_region": aggregate["is_non_core_region"],
                "is_cbm_region": aggregate["is_cbm_region"],
                "is_linker_region": aggregate["is_linker_region"],
            }
        )

    protein_patch_summary = []
    for protein_id in sorted(protein_patch_rows):
        rows = protein_patch_rows[protein_id]
        n_rows = len(rows)

        def _mean(column: str) -> float:
            return sum(_as_float(row[column]) for row in rows) / n_rows if n_rows else 0.0

        protein_patch_summary.append(
            {
                "protein_id": protein_id,
                "n_conditions_total": all_condition_counts_by_protein.get(protein_id, 0),
                "n_conditions_with_valid_clusters": valid_condition_counts_by_protein.get(protein_id, 0),
                "mean_aromatic_contact_fraction": _mean("aromatic_contact_fraction"),
                "mean_polar_contact_fraction": _mean("polar_contact_fraction"),
                "mean_charged_contact_fraction": _mean("charged_contact_fraction"),
                "mean_hbond_contact_fraction": _mean("hbond_contact_fraction"),
                "mean_catalytic_surface_contact_fraction": _mean("catalytic_surface_contact_fraction"),
                "mean_non_core_contact_fraction": _mean("non_core_contact_fraction"),
                "mean_cbm_contact_fraction": _mean("cbm_contact_fraction"),
                "mean_linker_contact_fraction": _mean("linker_contact_fraction"),
            }
        )

    return ResidueImportanceOutputs(
        protein_condition_residue_scores=protein_condition_residue_scores,
        protein_residue_regio_delta=protein_residue_regio_delta,
        condition_patch_summary=condition_patch_summary,
        protein_patch_summary=protein_patch_summary,
    )


def _write_rows(output_path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_protein_condition_residue_scores(rows: list[dict[str, Any]], output_path: Path) -> None:
    _write_rows(output_path, PROTEIN_CONDITION_RESIDUE_SCORE_COLUMNS, rows)


def write_protein_residue_regio_delta(rows: list[dict[str, Any]], output_path: Path) -> None:
    _write_rows(output_path, PROTEIN_RESIDUE_REGIO_DELTA_COLUMNS, rows)


def write_condition_patch_summary(rows: list[dict[str, Any]], output_path: Path) -> None:
    _write_rows(output_path, CONDITION_PATCH_SUMMARY_COLUMNS, rows)


def write_protein_patch_summary(rows: list[dict[str, Any]], output_path: Path) -> None:
    _write_rows(output_path, PROTEIN_PATCH_SUMMARY_COLUMNS, rows)
