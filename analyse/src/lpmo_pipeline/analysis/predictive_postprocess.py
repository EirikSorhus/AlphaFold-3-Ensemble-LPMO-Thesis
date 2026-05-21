"""Predictive analysis postprocess over condition-level summary tables."""

from __future__ import annotations

import csv
import json
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from lpmo_pipeline.analysis.predictive_models import (
    BinaryClassificationCVResult,
    run_grouped_binary_logistic_cv,
)

logger = logging.getLogger(__name__)

C1_C4_FEATURE_COLUMNS = [
    "contact_eligible_fraction",
    "noise_fraction",
    "top_cluster_occupancy",
    "hbond_contact_fraction",
    "aromatic_contact_fraction",
    "median_n_non_vdw_interactions",
    "occupancy_weighted_any_plausibility",
    "occupancy_weighted_C1_minus_C4_plausibility",
]

SUBSTRATE_FEATURE_COLUMNS = [
    "contact_eligible_fraction",
    "noise_fraction",
    "top_cluster_occupancy",
    "hbond_contact_fraction",
    "aromatic_contact_fraction",
    "median_n_non_vdw_interactions",
    "occupancy_weighted_any_plausibility",
    "occupancy_weighted_C1_minus_C4_plausibility",
]

_C1_C4_DIRECT_COLUMNS = [
    "contact_eligible_fraction",
    "noise_fraction",
    "top_cluster_occupancy",
    "hbond_contact_fraction",
    "aromatic_contact_fraction",
    "median_n_non_vdw_interactions",
]

_SUBSTRATE_DIRECT_COLUMNS = [
    "contact_eligible_fraction",
    "noise_fraction",
    "top_cluster_occupancy",
    "hbond_contact_fraction",
    "aromatic_contact_fraction",
    "median_n_non_vdw_interactions",
]


@dataclass
class PredictivePostprocessResult:
    """Predictive postprocess output paths."""

    output_dir: Path
    summary_path: Path
    modeling_table_paths: dict[str, Path] = field(default_factory=dict)
    metrics_paths: dict[str, Path] = field(default_factory=dict)
    predictions_paths: dict[str, Path] = field(default_factory=dict)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in fieldnames})


def _as_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _mean_numeric(values: Iterable[float | None]) -> float | str:
    numeric_values = [value for value in values if value is not None]
    if not numeric_values:
        return ""
    return sum(numeric_values) / len(numeric_values)


def _normalize_construct_type(value: Any) -> str:
    return str(value or "").strip().lower()


def _normalize_label(value: Any) -> str:
    return str(value or "").strip().lower()


def _first_present(row: dict[str, Any], keys: tuple[str, ...], default: str = "") -> str:
    for key in keys:
        value = row.get(key)
        if value not in (None, ""):
            return str(value)
    return default


def _protein_id_from_row(row: dict[str, Any]) -> str:
    return _first_present(
        row,
        ("protein_id", "uniprot_id", "UniProt_ID", "Uniprot_ID", "UniProtID"),
    ).strip()


def _choose_construct_type(rows: list[dict[str, Any]]) -> str:
    available = sorted(
        {
            _normalize_construct_type(row.get("construct_type"))
            for row in rows
            if str(row.get("protein_id", "")).strip()
        }
    )
    if not available:
        return ""
    for preferred in ("domain_only", "catalytic_domain"):
        if preferred in available:
            return preferred
    return available[0]


def filter_condition_rows_by_construct(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Keep one construct representation per protein."""

    by_protein: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        protein_id = _protein_id_from_row(row)
        if protein_id:
            by_protein.setdefault(protein_id, []).append(row)

    selected: list[dict[str, Any]] = []
    for protein_id in sorted(by_protein):
        protein_rows = by_protein[protein_id]
        chosen_construct = _choose_construct_type(protein_rows)
        if not chosen_construct:
            continue
        selected.extend(
            row
            for row in protein_rows
            if _normalize_construct_type(row.get("construct_type")) == chosen_construct
        )
    return selected


def _index_by_protein(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    indexed: dict[str, dict[str, Any]] = {}
    for row in rows:
        protein_id = _protein_id_from_row(row)
        if protein_id and protein_id not in indexed:
            indexed[protein_id] = row
    return indexed


def _normalize_ec_tokens(ec_raw: Any) -> list[str]:
    text = str(ec_raw or "").strip()
    if not text:
        return []
    return re.findall(r"1\.14\.99\.[0-9-]+", text)


def _map_ec_token(token: str, family_raw: Any) -> dict[str, str] | None:
    family = str(family_raw or "").strip().upper()
    if token == "1.14.99.54":
        return {
            "substrate_class": "cellulose",
            "regio_class": "C1",
            "activity_label": "cellulose_C1_hydroxylating",
            "mapping_rule": "exact_1.14.99.54",
        }
    if token == "1.14.99.56":
        return {
            "substrate_class": "cellulose",
            "regio_class": "C4",
            "activity_label": "cellulose_C4_dehydrogenating",
            "mapping_rule": "exact_1.14.99.56",
        }
    if token == "1.14.99.53":
        return {
            "substrate_class": "chitin",
            "regio_class": "mixed",
            "activity_label": "chitin_C1_C4_mixed",
            "mapping_rule": "exact_1.14.99.53",
        }
    if token == "1.14.99.55":
        return {
            "substrate_class": "starch",
            "regio_class": "C1",
            "activity_label": "starch_C1_hydroxylating",
            "mapping_rule": "exact_1.14.99.55",
        }
    if token == "1.14.99.-" and "AA17" in family:
        return {
            "substrate_class": "homogalacturonan",
            "regio_class": "C4",
            "activity_label": "homogalacturonan_C4_oxidation",
            "mapping_rule": "special_1.14.99.-_AA17",
        }
    if token == "1.14.99.-":
        return {
            "substrate_class": "xylan_or_other",
            "regio_class": "unknown",
            "activity_label": "xylan_like_oxidative",
            "mapping_rule": "special_1.14.99.-_non_AA17",
        }
    return None


def _resolve_activity_annotation(metadata_row: dict[str, Any]) -> dict[str, str]:
    direct_regio = _first_present(
        metadata_row,
        ("mapped_regio_class", "regio_class", "regio_label", "experimental_regio_label"),
    )
    direct_substrate = _first_present(
        metadata_row,
        ("mapped_substrate_class", "substrate_class", "substrate_label", "experimental_substrate_label"),
    )
    direct_activity = _first_present(
        metadata_row,
        ("mapped_activity_label", "activity_label", "experimental_activity_label"),
    )
    direct_rule = _first_present(metadata_row, ("mapping_rule", "activity_mapping_rule"))

    ec_tokens = _normalize_ec_tokens(
        _first_present(metadata_row, ("EC_Number", "ec_number", "ec_numbers"))
    )
    family = _first_present(metadata_row, ("CAZy_family", "family", "family_label", "aa_family"))
    derived = {
        "regio_label": "unknown",
        "substrate_label": "unknown",
        "activity_label": "unknown",
        "mapping_rule": "missing_ec",
    }
    if not ec_tokens:
        derived = {
            "regio_label": "unknown",
            "substrate_label": "unknown",
            "activity_label": "unknown",
            "mapping_rule": "missing_ec",
        }
    else:
        recognized: list[dict[str, str]] = []
        for token in ec_tokens:
            mapped = _map_ec_token(token, family)
            if mapped is not None:
                recognized.append(mapped)
        if recognized:
            substrate_set = {item["substrate_class"] for item in recognized}
            regio_set: set[str] = set()
            for item in recognized:
                regio = _normalize_label(item["regio_class"])
                if regio == "mixed":
                    regio_set.update({"c1", "c4"})
                elif regio in {"c1", "c4"}:
                    regio_set.add(regio)
            if regio_set == {"c1"}:
                regio_label = "C1"
            elif regio_set == {"c4"}:
                regio_label = "C4"
            elif regio_set == {"c1", "c4"}:
                regio_label = "C1/C4"
            else:
                regio_label = "unknown"

            prioritized_substrates = [
                substrate
                for substrate in ("chitin", "cellulose", "starch", "homogalacturonan", "xylan_or_other")
                if substrate in substrate_set
            ]
            derived = {
                "regio_label": regio_label,
                "substrate_label": "+".join(prioritized_substrates) if prioritized_substrates else "unknown",
                "activity_label": "+".join(sorted({item["activity_label"] for item in recognized})) or "unknown",
                "mapping_rule": "+".join(sorted({item["mapping_rule"] for item in recognized})) or "derived_from_ec",
            }
        else:
            derived = {
                "regio_label": "unknown",
                "substrate_label": "unknown",
                "activity_label": "unknown",
                "mapping_rule": "unmapped_ec",
            }

    return {
        "regio_label": direct_regio or derived["regio_label"],
        "substrate_label": direct_substrate or derived["substrate_label"],
        "activity_label": direct_activity or derived["activity_label"],
        "mapping_rule": direct_rule or ("preannotated_metadata" if (direct_regio or direct_substrate or direct_activity) else derived["mapping_rule"]),
    }


def _parse_regio_targets(label: Any) -> tuple[int, int] | None:
    normalized = _normalize_label(label)
    if normalized == "c1":
        return (1, 0)
    if normalized == "c4":
        return (0, 1)
    if normalized in {"c1/c4", "c1+c4", "mixed"}:
        return (1, 1)
    return None


def _parse_substrate_set(label: Any) -> set[str]:
    normalized = _normalize_label(label)
    active: set[str] = set()
    if "chitin" in normalized:
        active.add("chitin")
    if "cellulose" in normalized:
        active.add("cellulose")
    if "starch" in normalized or "amylose" in normalized:
        active.add("starch")
    return active


def _occupancy_weighted_c1_minus_c4_plausibility(row: dict[str, Any]) -> float | None:
    c1 = _as_float(row.get("occupancy_weighted_c1_plausible_fraction"))
    c4 = _as_float(row.get("occupancy_weighted_c4_plausible_fraction"))
    if c1 is None or c4 is None:
        return None
    return c1 - c4


def _occupancy_weighted_any_plausibility(row: dict[str, Any]) -> float | None:
    c1 = _as_float(row.get("occupancy_weighted_c1_plausible_fraction"))
    c4 = _as_float(row.get("occupancy_weighted_c4_plausible_fraction"))
    values = [value for value in (c1, c4) if value is not None]
    if not values:
        return None
    return max(values)


def _occupancy_weighted_delta_attack_angle(row: dict[str, Any]) -> float | None:
    c1 = _as_float(row.get("occupancy_weighted_attack_angle_C1_median"))
    c4 = _as_float(row.get("occupancy_weighted_attack_angle_C4_median"))
    if c1 is None or c4 is None:
        return None
    return c1 - c4


def _occupancy_weighted_delta_cu_distance(row: dict[str, Any]) -> float | None:
    c1 = _as_float(row.get("occupancy_weighted_Cu_C1_distance_median"))
    c4 = _as_float(row.get("occupancy_weighted_Cu_C4_distance_median"))
    if c1 is None or c4 is None:
        return None
    return c4 - c1


def _aggregate_feature(rows: list[dict[str, Any]], column_name: str) -> float | str:
    return _mean_numeric(_as_float(row.get(column_name)) for row in rows)


def build_c1_c4_modeling_rows(
    condition_rows: list[dict[str, Any]],
    protein_metadata_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Build protein-substrate modeling rows for C1/C4 prediction."""

    selected_rows = filter_condition_rows_by_construct(condition_rows)
    metadata_by_protein = _index_by_protein(protein_metadata_rows)
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for row in selected_rows:
        protein_id = str(row.get("protein_id", "")).strip()
        substrate_class = str(row.get("substrate_class", "")).strip()
        if protein_id and substrate_class:
            grouped.setdefault((protein_id, substrate_class), []).append(row)

    out: list[dict[str, Any]] = []
    for (protein_id, substrate_class), rows in sorted(grouped.items()):
        metadata = metadata_by_protein.get(protein_id, {})
        annotation = _resolve_activity_annotation(metadata)
        targets = _parse_regio_targets(annotation["regio_label"])
        if targets is None:
            continue
        construct_type = rows[0].get("construct_type", "")
        row_out = {
            "model_row_id": f"{protein_id}__{substrate_class}",
            "protein_id": protein_id,
            "substrate_class": substrate_class,
            "construct_type": construct_type,
            "n_source_rows": len(rows),
            "has_C1_activity": targets[0],
            "has_C4_activity": targets[1],
            "occupancy_weighted_any_plausibility": _mean_numeric(
                _occupancy_weighted_any_plausibility(row) for row in rows
            ),
            "occupancy_weighted_C1_minus_C4_plausibility": _mean_numeric(
                _occupancy_weighted_c1_minus_c4_plausibility(row) for row in rows
            ),
        }
        for column_name in _C1_C4_DIRECT_COLUMNS:
            row_out[column_name] = _aggregate_feature(rows, column_name)
        out.append(row_out)
    return out


def build_substrate_modeling_rows(
    condition_rows: list[dict[str, Any]],
    protein_metadata_rows: list[dict[str, Any]],
    *,
    substrate_class: str,
) -> list[dict[str, Any]]:
    """Build protein-level modeling rows for one substrate-specific model."""

    selected_rows = filter_condition_rows_by_construct(condition_rows)
    metadata_by_protein = _index_by_protein(protein_metadata_rows)
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in selected_rows:
        protein_id = str(row.get("protein_id", "")).strip()
        row_substrate = str(row.get("substrate_class", "")).strip().lower()
        if protein_id and row_substrate == substrate_class:
            grouped.setdefault(protein_id, []).append(row)

    out: list[dict[str, Any]] = []
    for protein_id, rows in sorted(grouped.items()):
        metadata = metadata_by_protein.get(protein_id, {})
        annotation = _resolve_activity_annotation(metadata)
        active_substrates = _parse_substrate_set(annotation["substrate_label"])
        if not active_substrates:
            continue
        row_out = {
            "model_row_id": protein_id,
            "protein_id": protein_id,
            "prediction_substrate": substrate_class,
            "construct_type": rows[0].get("construct_type", ""),
            "n_source_rows": len(rows),
            "is_active_on_substrate": 1 if substrate_class in active_substrates else 0,
            "occupancy_weighted_any_plausibility": _mean_numeric(
                _occupancy_weighted_any_plausibility(row) for row in rows
            ),
            "occupancy_weighted_C1_minus_C4_plausibility": _mean_numeric(
                _occupancy_weighted_c1_minus_c4_plausibility(row) for row in rows
            ),
        }
        for column_name in _SUBSTRATE_DIRECT_COLUMNS:
            row_out[column_name] = _aggregate_feature(rows, column_name)
        out.append(row_out)
    return out


def _write_metrics_files(
    *,
    result: BinaryClassificationCVResult,
    metrics_path: Path,
    fold_metrics_path: Path,
    predictions_path: Path,
) -> None:
    summary_row = {
        "model_name": result.model_name,
        "target_column": result.target_column,
        "n_rows": result.n_rows,
        "n_positive": result.n_positive,
        "n_negative": result.n_negative,
        "n_folds": result.n_folds,
        "mean_pr_auc": result.mean_metrics.get("pr_auc", ""),
        "std_pr_auc": result.std_metrics.get("pr_auc", ""),
        "mean_balanced_accuracy": result.mean_metrics.get("balanced_accuracy", ""),
        "std_balanced_accuracy": result.std_metrics.get("balanced_accuracy", ""),
        "mean_recall_at_k": result.mean_metrics.get("recall_at_k", ""),
        "std_recall_at_k": result.std_metrics.get("recall_at_k", ""),
        "feature_importances_json": json.dumps(result.feature_importances, sort_keys=True),
    }
    _write_tsv(metrics_path, [summary_row], list(summary_row.keys()))
    if result.fold_metrics:
        _write_tsv(fold_metrics_path, result.fold_metrics, list(result.fold_metrics[0].keys()))
    else:
        _write_tsv(fold_metrics_path, [], [
            "fold_id",
            "n_test_rows",
            "n_test_positive",
            "n_test_negative",
            "pr_auc",
            "balanced_accuracy",
            "recall_at_k",
            "k",
        ])
    if result.predictions:
        _write_tsv(predictions_path, result.predictions, list(result.predictions[0].keys()))
    else:
        _write_tsv(predictions_path, [], [
            "fold_id",
            "row_id",
            "protein_id",
            "y_true",
            "y_score",
            "y_pred",
        ])


def run_predictive_postprocess(
    *,
    condition_table_path: Path,
    protein_metadata_path: Path,
    output_dir: Path,
    task: str = "all",
    n_folds: int = 5,
    random_state: int = 42,
) -> PredictivePostprocessResult:
    """Run predictive postprocess over condition-level summary tables."""

    condition_rows = read_tsv(condition_table_path)
    protein_metadata_rows = read_tsv(protein_metadata_path)

    predictive_root = output_dir / "10_predictive"
    modeling_dir = predictive_root / "modeling_tables"
    cv_results_dir = predictive_root / "cv_results"
    summary_path = predictive_root / "predictive_summary.json"
    result = PredictivePostprocessResult(output_dir=predictive_root, summary_path=summary_path)

    summary_data: dict[str, Any] = {
        "task": task,
        "condition_table": str(condition_table_path),
        "protein_metadata": str(protein_metadata_path),
        "modeling_tables": {},
        "metrics": {},
        "predictions": {},
    }

    if task in {"all", "c1_c4"}:
        c1_c4_rows = build_c1_c4_modeling_rows(condition_rows, protein_metadata_rows)
        c1_c4_modeling_table = modeling_dir / "c1_c4_modeling_table.tsv"
        c1_c4_columns = [
            "model_row_id",
            "protein_id",
            "substrate_class",
            "construct_type",
            "n_source_rows",
            "has_C1_activity",
            "has_C4_activity",
            *C1_C4_FEATURE_COLUMNS,
        ]
        _write_tsv(c1_c4_modeling_table, c1_c4_rows, c1_c4_columns)
        result.modeling_table_paths["c1_c4"] = c1_c4_modeling_table
        summary_data["modeling_tables"]["c1_c4"] = str(c1_c4_modeling_table)

        for model_name, target_column in (
            ("c1_activity", "has_C1_activity"),
            ("c4_activity", "has_C4_activity"),
        ):
            cv_result = run_grouped_binary_logistic_cv(
                c1_c4_rows,
                feature_columns=C1_C4_FEATURE_COLUMNS,
                target_column=target_column,
                group_column="protein_id",
                row_id_column="model_row_id",
                model_name=model_name,
                n_folds=n_folds,
                random_state=random_state,
            )
            metrics_path = cv_results_dir / f"{model_name}_metrics.tsv"
            fold_metrics_path = cv_results_dir / f"{model_name}_fold_metrics.tsv"
            predictions_path = cv_results_dir / f"{model_name}_predictions.tsv"
            _write_metrics_files(
                result=cv_result,
                metrics_path=metrics_path,
                fold_metrics_path=fold_metrics_path,
                predictions_path=predictions_path,
            )
            result.metrics_paths[model_name] = metrics_path
            result.predictions_paths[model_name] = predictions_path
            summary_data["metrics"][model_name] = str(metrics_path)
            summary_data["predictions"][model_name] = str(predictions_path)

    if task in {"all", "substrate"}:
        for substrate_class in ("cellulose", "chitin", "starch"):
            substrate_rows = build_substrate_modeling_rows(
                condition_rows,
                protein_metadata_rows,
                substrate_class=substrate_class,
            )
            model_name = f"{substrate_class}_activity"
            modeling_table_path = modeling_dir / f"{model_name}_modeling_table.tsv"
            modeling_columns = [
                "model_row_id",
                "protein_id",
                "prediction_substrate",
                "construct_type",
                "n_source_rows",
                "is_active_on_substrate",
                *SUBSTRATE_FEATURE_COLUMNS,
            ]
            _write_tsv(modeling_table_path, substrate_rows, modeling_columns)
            result.modeling_table_paths[model_name] = modeling_table_path
            summary_data["modeling_tables"][model_name] = str(modeling_table_path)

            cv_result = run_grouped_binary_logistic_cv(
                substrate_rows,
                feature_columns=SUBSTRATE_FEATURE_COLUMNS,
                target_column="is_active_on_substrate",
                group_column="protein_id",
                row_id_column="model_row_id",
                model_name=model_name,
                n_folds=n_folds,
                random_state=random_state,
            )
            metrics_path = cv_results_dir / f"{model_name}_metrics.tsv"
            fold_metrics_path = cv_results_dir / f"{model_name}_fold_metrics.tsv"
            predictions_path = cv_results_dir / f"{model_name}_predictions.tsv"
            _write_metrics_files(
                result=cv_result,
                metrics_path=metrics_path,
                fold_metrics_path=fold_metrics_path,
                predictions_path=predictions_path,
            )
            result.metrics_paths[model_name] = metrics_path
            result.predictions_paths[model_name] = predictions_path
            summary_data["metrics"][model_name] = str(metrics_path)
            summary_data["predictions"][model_name] = str(predictions_path)

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w") as handle:
        json.dump(summary_data, handle, indent=2)

    return result