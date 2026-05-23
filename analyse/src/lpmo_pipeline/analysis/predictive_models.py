# src/lpmo_pipeline/analysis/predictive_models.py
"""
Responsibility: Final compact sklearn backend for exploratory predictive models.
Input:  condition-derived modeling tables or optional predictive_cluster_table rows.
Output: Grouped-CV metrics, predictions, and feature importances.

Status:
    The primary backend is sklearn LogisticRegression with L2 regularization,
    fold-local numeric imputation/standardization, balanced class weighting, and
    protein-level grouped CV. Decision-tree/random-forest helpers remain
    sensitivity or compatibility paths, not the primary reported backend.

Anti p-hack rules:
  - HDBSCAN params LOCKED (from tuning, not re-tuned here)
    - Grouped CV: fold by protein_id (all clusters from a protein remain together)
    - Stratified CV: maintain class balance (C1/C4/mixed, polymer type) at protein level
  - CV hierarchy defined in configs/cv_hierarchy.yaml
  - Default CV target is 5 folds when enough protein groups exist.

Main rule:
    - Keep all rows from the same protein in the same CV fold.
    - Treat outputs as exploratory structure-derived associations, not causal
      activity evidence.
"""
from __future__ import annotations

import json
import logging
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier

logger = logging.getLogger(__name__)


@dataclass
class CVFold:
    """A single cross-validation fold."""

    fold_id: int
    train_protein_ids: list[str] = field(default_factory=list)
    test_protein_ids: list[str] = field(default_factory=list)


@dataclass
class CVResult:
    """Cross-validation result."""

    n_folds: int = 0
    metric_name: str = ""  # e.g. "balanced_accuracy", "f1_macro"
    model_type: str = ""
    per_fold_scores: list[float] = field(default_factory=list)
    mean_score: float = 0.0
    std_score: float = 0.0
    feature_importances: dict[str, float] = field(default_factory=dict)
    fold_test_protein_ids: list[list[str]] = field(default_factory=list)


@dataclass
class BinaryClassificationCVResult:
    """Grouped CV result for one binary classification model."""

    model_name: str
    target_column: str
    feature_columns: list[str]
    n_rows: int = 0
    n_positive: int = 0
    n_negative: int = 0
    n_folds: int = 0
    fold_metrics: list[dict[str, Any]] = field(default_factory=list)
    mean_metrics: dict[str, float] = field(default_factory=dict)
    std_metrics: dict[str, float] = field(default_factory=dict)
    feature_importances: dict[str, float] = field(default_factory=dict)
    predictions: list[dict[str, Any]] = field(default_factory=list)


def _as_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(numeric):
        return None
    return numeric


def _as_binary_int(value: Any) -> int | None:
    if isinstance(value, bool):
        return int(value)
    if value in (0, 1):
        return int(value)
    normalized = str(value).strip().lower()
    if normalized in {"0", "false", "no"}:
        return 0
    if normalized in {"1", "true", "yes"}:
        return 1
    return None


def _prepare_numeric_matrix(rows: list[dict[str, Any]], feature_columns: list[str]) -> np.ndarray:
    return np.asarray(
        [
            [_as_float(row.get(column)) if _as_float(row.get(column)) is not None else np.nan for column in feature_columns]
            for row in rows
        ],
        dtype=float,
    )


def _fit_numeric_preprocessor(x_train: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    medians = np.nanmedian(x_train, axis=0)
    medians = np.where(np.isnan(medians), 0.0, medians)
    imputed = np.where(np.isnan(x_train), medians, x_train)
    means = np.mean(imputed, axis=0)
    stds = np.std(imputed, axis=0)
    stds = np.where(stds == 0.0, 1.0, stds)
    return medians, means, stds


def _transform_numeric_matrix(
    x: np.ndarray,
    *,
    medians: np.ndarray,
    means: np.ndarray,
    stds: np.ndarray,
) -> np.ndarray:
    imputed = np.where(np.isnan(x), medians, x)
    return (imputed - means) / stds


def _recall_at_k(y_true: np.ndarray, y_score: np.ndarray, k: int) -> float:
    if k <= 0:
        return 0.0
    positive_total = int(np.sum(y_true))
    if positive_total <= 0:
        return 0.0
    top_k = min(k, len(y_true))
    ranking = np.argsort(-y_score, kind="stable")[:top_k]
    return float(np.sum(y_true[ranking]) / positive_total)


def _binary_balanced_accuracy(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    recalls = []
    for label in (0, 1):
        label_mask = y_true == label
        label_total = int(np.sum(label_mask))
        if label_total > 0:
            recalls.append(float(np.sum(y_pred[label_mask] == label) / label_total))
    return float(np.mean(recalls)) if recalls else 0.0


def _label_balanced_accuracy(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    recalls = []
    for label in sorted(set(y_true.tolist())):
        label_mask = y_true == label
        label_total = int(np.sum(label_mask))
        if label_total > 0:
            recalls.append(float(np.sum(y_pred[label_mask] == label) / label_total))
    return float(np.mean(recalls)) if recalls else 0.0


def build_grouped_stratified_cv(
    protein_ids: list[str],
    activity_labels: dict[str, str],     # {protein_id: "C1"|"C4"|"C1+C4"}
    hierarchy: dict[str, str],           # {protein_id: structure_group}
    n_folds: int = 5,
    random_state: int = 42,
) -> list[CVFold]:
    """Build grouped + stratified CV folds.

    Groups: proteins from the same structure group never split across folds.
    Stratification: class balance (C1/C4/mixed) maintained per fold.

    Args:
        protein_ids: All protein IDs in the dataset.
        activity_labels: Activity class per protein.
        hierarchy: Grouping structure (e.g. by crystal family).
        n_folds: Number of CV folds.

    Returns:
        List of CVFold objects.
    """
    # --- Step 1: Group enzymes by hierarchy ---
    groups: dict[str, list[str]] = {}
    for protein_id in protein_ids:
        group = hierarchy.get(protein_id, protein_id)
        groups.setdefault(group, []).append(protein_id)

    group_keys = sorted(groups.keys())
    if not group_keys:
        return []
    n_actual_folds = max(1, min(n_folds, len(group_keys)))

    group_labels: dict[str, str] = {}
    for gkey, ids in groups.items():
        labels = [activity_labels.get(protein_id, "unknown") for protein_id in ids]
        most_common = Counter(labels).most_common(1)[0][0]
        group_labels[gkey] = most_common

    rng = np.random.default_rng(random_state)
    groups_by_label: dict[str, list[str]] = defaultdict(list)
    for group_key, label in group_labels.items():
        groups_by_label[label].append(group_key)
    test_groups_by_fold: list[list[str]] = [[] for _ in range(n_actual_folds)]
    fold_sizes = [0 for _ in range(n_actual_folds)]
    for label in sorted(groups_by_label):
        label_groups = groups_by_label[label]
        rng.shuffle(label_groups)
        for group_key in label_groups:
            fold_index = min(range(n_actual_folds), key=lambda idx: fold_sizes[idx])
            test_groups_by_fold[fold_index].append(group_key)
            fold_sizes[fold_index] += len(groups[group_key])

    folds: list[CVFold] = []
    all_groups = set(group_keys)
    for fold_id, test_groups in enumerate(test_groups_by_fold):
        test_group_set = set(test_groups)
        train_groups = sorted(all_groups - test_group_set)
        folds.append(
            CVFold(
                fold_id=fold_id,
                train_protein_ids=sorted(
                    protein_id for group in train_groups for protein_id in groups[group]
                ),
                test_protein_ids=sorted(
                    protein_id for group in sorted(test_group_set) for protein_id in groups[group]
                ),
            )
        )

    logger.info("Built %d grouped+stratified CV folds from %d groups", len(folds), len(groups))
    return folds


def run_binary_predictive_baselines(
    rows: list[dict[str, Any]],
    *,
    feature_columns: list[str],
    label_column: str = "experimental_regio_label",
    group_column: str = "protein_id",
    positive_label: str = "C1",
    negative_label: str = "C4",
    n_folds: int = 5,
    random_state: int = 42,
) -> list[CVResult]:
    """Run a narrow leakage-safe binary baseline suite on predictive rows."""

    filtered_rows = [
        row
        for row in rows
        if str(row.get(label_column, "")) in {positive_label, negative_label}
        and str(row.get(group_column, ""))
        and all(_as_float(row.get(column)) is not None for column in feature_columns)
    ]
    if not filtered_rows:
        return []

    protein_ids = sorted({str(row[group_column]) for row in filtered_rows})
    activity_labels = {
        str(row[group_column]): str(row[label_column])
        for row in filtered_rows
    }
    folds = build_grouped_stratified_cv(
        protein_ids,
        activity_labels,
        hierarchy={protein_id: protein_id for protein_id in protein_ids},
        n_folds=n_folds,
        random_state=random_state,
    )

    model_specs = [
        (
            "logistic_regression",
            LogisticRegression(max_iter=1000, class_weight="balanced", random_state=random_state),
        ),
        (
            "decision_tree",
            DecisionTreeClassifier(max_depth=3, class_weight="balanced", random_state=random_state),
        ),
    ]
    results: list[CVResult] = []
    for model_type, estimator in model_specs:
        result = CVResult(
            n_folds=len(folds),
            metric_name="balanced_accuracy",
            model_type=model_type,
        )
        importance_accumulator = np.zeros(len(feature_columns), dtype=float)
        n_importance_folds = 0
        for fold in folds:
            train_proteins = set(fold.train_protein_ids)
            test_proteins = set(fold.test_protein_ids)
            train_rows = [row for row in filtered_rows if str(row[group_column]) in train_proteins]
            test_rows = [row for row in filtered_rows if str(row[group_column]) in test_proteins]
            if not train_rows or not test_rows:
                continue
            y_train = np.asarray(
                [1 if str(row[label_column]) == positive_label else 0 for row in train_rows],
                dtype=int,
            )
            y_test = np.asarray(
                [1 if str(row[label_column]) == positive_label else 0 for row in test_rows],
                dtype=int,
            )
            if len(set(y_train.tolist())) < 2 or len(set(y_test.tolist())) < 2:
                continue
            x_train = np.asarray(
                [[_as_float(row[column]) for column in feature_columns] for row in train_rows],
                dtype=float,
            )
            x_test = np.asarray(
                [[_as_float(row[column]) for column in feature_columns] for row in test_rows],
                dtype=float,
            )
            estimator.fit(x_train, y_train)
            y_pred = estimator.predict(x_test)
            result.per_fold_scores.append(_binary_balanced_accuracy(y_test, y_pred))
            result.fold_test_protein_ids.append(fold.test_protein_ids)
            if hasattr(estimator, "coef_"):
                importance_accumulator += np.abs(estimator.coef_[0])
                n_importance_folds += 1
            elif hasattr(estimator, "feature_importances_"):
                importance_accumulator += estimator.feature_importances_
                n_importance_folds += 1

        result.n_folds = len(result.per_fold_scores)
        result.mean_score = float(np.mean(result.per_fold_scores)) if result.per_fold_scores else 0.0
        result.std_score = float(np.std(result.per_fold_scores)) if result.per_fold_scores else 0.0
        if n_importance_folds:
            result.feature_importances = {
                feature: float(value)
                for feature, value in zip(
                    feature_columns,
                    importance_accumulator / n_importance_folds,
                    strict=True,
                )
            }
        results.append(result)
    return results


def run_grouped_binary_logistic_cv(
    rows: list[dict[str, Any]],
    *,
    feature_columns: list[str],
    target_column: str,
    group_column: str = "protein_id",
    row_id_column: str = "model_row_id",
    model_name: str = "logistic_regression",
    n_folds: int = 5,
    random_state: int = 42,
    c_value: float = 1.0,
) -> BinaryClassificationCVResult:
    """Run grouped 5-fold binary logistic regression with fold-local preprocessing."""

    filtered_rows = []
    for row in rows:
        group_value = str(row.get(group_column, "")).strip()
        target_value = _as_binary_int(row.get(target_column))
        if not group_value or target_value is None:
            continue
        filtered_rows.append(row)

    result = BinaryClassificationCVResult(
        model_name=model_name,
        target_column=target_column,
        feature_columns=list(feature_columns),
        n_rows=len(filtered_rows),
    )
    if not filtered_rows:
        return result

    protein_ids = sorted({str(row[group_column]) for row in filtered_rows})
    targets_by_protein: dict[str, list[int]] = defaultdict(list)
    for row in filtered_rows:
        target_value = _as_binary_int(row[target_column])
        if target_value is not None:
            targets_by_protein[str(row[group_column])].append(target_value)

    activity_labels = {}
    for protein_id in protein_ids:
        protein_targets = set(targets_by_protein[protein_id])
        if protein_targets == {1}:
            activity_labels[protein_id] = "positive"
        elif protein_targets == {0}:
            activity_labels[protein_id] = "negative"
        else:
            activity_labels[protein_id] = "mixed"

    row_targets = [target for targets in targets_by_protein.values() for target in targets]
    result.n_positive = sum(1 for target in row_targets if target == 1)
    result.n_negative = sum(1 for target in row_targets if target == 0)
    if result.n_positive == 0 or result.n_negative == 0:
        return result

    folds = build_grouped_stratified_cv(
        protein_ids,
        activity_labels,
        hierarchy={protein_id: protein_id for protein_id in protein_ids},
        n_folds=n_folds,
        random_state=random_state,
    )
    if not folds:
        return result

    importance_accumulator = np.zeros(len(feature_columns), dtype=float)
    n_importance_folds = 0

    for fold in folds:
        train_proteins = set(fold.train_protein_ids)
        test_proteins = set(fold.test_protein_ids)
        train_rows = [row for row in filtered_rows if str(row[group_column]) in train_proteins]
        test_rows = [row for row in filtered_rows if str(row[group_column]) in test_proteins]
        if not train_rows or not test_rows:
            continue

        y_train = np.asarray([_as_binary_int(row[target_column]) for row in train_rows], dtype=int)
        y_test = np.asarray([_as_binary_int(row[target_column]) for row in test_rows], dtype=int)
        if len(set(y_train.tolist())) < 2:
            continue

        x_train_raw = _prepare_numeric_matrix(train_rows, feature_columns)
        x_test_raw = _prepare_numeric_matrix(test_rows, feature_columns)
        medians, means, stds = _fit_numeric_preprocessor(x_train_raw)
        x_train = _transform_numeric_matrix(x_train_raw, medians=medians, means=means, stds=stds)
        x_test = _transform_numeric_matrix(x_test_raw, medians=medians, means=means, stds=stds)

        estimator = LogisticRegression(
            C=c_value,
            class_weight="balanced",
            solver="liblinear",
            max_iter=1000,
            random_state=random_state,
        )
        estimator.fit(x_train, y_train)
        y_score = estimator.predict_proba(x_test)[:, 1]
        y_pred = estimator.predict(x_test)

        k = int(np.sum(y_test))
        pr_auc = float(average_precision_score(y_test, y_score)) if k > 0 else 0.0
        balanced_accuracy = _binary_balanced_accuracy(y_test, y_pred)
        recall_at_k = _recall_at_k(y_test, y_score, k)
        result.fold_metrics.append(
            {
                "fold_id": fold.fold_id,
                "n_test_rows": len(test_rows),
                "n_test_positive": int(np.sum(y_test)),
                "n_test_negative": int(len(y_test) - np.sum(y_test)),
                "pr_auc": pr_auc,
                "balanced_accuracy": balanced_accuracy,
                "recall_at_k": recall_at_k,
                "k": k,
            }
        )
        for row, true_value, score_value, pred_value in zip(test_rows, y_test, y_score, y_pred, strict=True):
            result.predictions.append(
                {
                    "fold_id": fold.fold_id,
                    "row_id": row.get(row_id_column, ""),
                    "protein_id": row.get(group_column, ""),
                    "y_true": int(true_value),
                    "y_score": float(score_value),
                    "y_pred": int(pred_value),
                }
            )
        if hasattr(estimator, "coef_"):
            importance_accumulator += np.abs(estimator.coef_[0])
            n_importance_folds += 1

    result.n_folds = len(result.fold_metrics)
    for metric_name in ("pr_auc", "balanced_accuracy", "recall_at_k"):
        metric_values = [float(fold_metric[metric_name]) for fold_metric in result.fold_metrics]
        if metric_values:
            result.mean_metrics[metric_name] = float(np.mean(metric_values))
            result.std_metrics[metric_name] = float(np.std(metric_values))
        else:
            result.mean_metrics[metric_name] = 0.0
            result.std_metrics[metric_name] = 0.0
    if n_importance_folds:
        result.feature_importances = {
            feature: float(value)
            for feature, value in zip(
                feature_columns,
                importance_accumulator / n_importance_folds,
                strict=True,
            )
        }
    return result


def run_predictive_model(
    features: dict[str, list[float]],   # {cluster_id: feature_vector}
    labels: dict[str, str],             # {cluster_id: activity_class}
    cluster_to_protein: dict[str, str],  # {cluster_id: protein_id}
    folds: list[CVFold],
    model_type: str = "logistic_regression",
) -> CVResult:
    """Run predictive modeling with grouped CV.

    Args:
        features: Feature vectors per cluster row.
        labels: Activity labels per cluster row.
        cluster_to_protein: Mapping to enforce grouped CV at protein level.
        folds: Pre-built protein-level CV folds.
        model_type: "logistic_regression" | "decision_tree" | "random_forest"

    Returns:
        CVResult with per-fold scores and feature importances.
    """
    result = CVResult(metric_name="balanced_accuracy", model_type=model_type)
    valid_cluster_ids = sorted(
        cluster_id
        for cluster_id in cluster_to_protein
        if cluster_id in features and cluster_id in labels
    )
    if not valid_cluster_ids or not folds:
        return result

    n_features = len(features[valid_cluster_ids[0]])
    if any(len(features[cluster_id]) != n_features for cluster_id in valid_cluster_ids):
        raise ValueError("All feature vectors must have the same length")

    importance_accumulator = np.zeros(n_features, dtype=float)
    n_importance_folds = 0
    for fold in folds:
        train_proteins = set(fold.train_protein_ids)
        test_proteins = set(fold.test_protein_ids)
        train_clusters = [
            cluster_id
            for cluster_id in valid_cluster_ids
            if cluster_to_protein[cluster_id] in train_proteins
        ]
        test_clusters = [
            cluster_id
            for cluster_id in valid_cluster_ids
            if cluster_to_protein[cluster_id] in test_proteins
        ]
        if not train_clusters or not test_clusters:
            continue

        y_train = np.asarray([str(labels[cluster_id]) for cluster_id in train_clusters])
        y_test = np.asarray([str(labels[cluster_id]) for cluster_id in test_clusters])
        if len(set(y_train.tolist())) < 2:
            continue

        x_train_raw = np.asarray([features[cluster_id] for cluster_id in train_clusters], dtype=float)
        x_test_raw = np.asarray([features[cluster_id] for cluster_id in test_clusters], dtype=float)
        medians, means, stds = _fit_numeric_preprocessor(x_train_raw)
        x_train = _transform_numeric_matrix(x_train_raw, medians=medians, means=means, stds=stds)
        x_test = _transform_numeric_matrix(x_test_raw, medians=medians, means=means, stds=stds)

        if model_type == "logistic_regression":
            estimator = LogisticRegression(
                C=1.0,
                class_weight="balanced",
                solver="liblinear",
                max_iter=1000,
                random_state=42,
            )
        elif model_type == "decision_tree":
            estimator = DecisionTreeClassifier(
                max_depth=3,
                class_weight="balanced",
                random_state=42,
            )
        elif model_type == "random_forest":
            estimator = RandomForestClassifier(
                n_estimators=100,
                max_depth=3,
                class_weight="balanced",
                random_state=42,
            )
        else:
            raise ValueError(f"Unsupported model_type: {model_type}")

        estimator.fit(x_train, y_train)
        y_pred = estimator.predict(x_test)
        result.per_fold_scores.append(_label_balanced_accuracy(y_test, y_pred))
        result.fold_test_protein_ids.append(fold.test_protein_ids)

        if hasattr(estimator, "coef_"):
            coefficients = np.asarray(estimator.coef_, dtype=float)
            if coefficients.ndim == 1:
                importance_accumulator += np.abs(coefficients)
            else:
                importance_accumulator += np.mean(np.abs(coefficients), axis=0)
            n_importance_folds += 1
        elif hasattr(estimator, "feature_importances_"):
            importance_accumulator += estimator.feature_importances_
            n_importance_folds += 1

    result.n_folds = len(result.per_fold_scores)
    result.mean_score = float(np.mean(result.per_fold_scores)) if result.per_fold_scores else 0.0
    result.std_score = float(np.std(result.per_fold_scores)) if result.per_fold_scores else 0.0
    if n_importance_folds:
        result.feature_importances = {
            f"feature_{index}": float(value)
            for index, value in enumerate(importance_accumulator / n_importance_folds)
        }

    logger.info(
        "Predictive model (cluster rows): %s CV score = %.3f ± %.3f",
        model_type, result.mean_score, result.std_score,
    )
    return result


def write_cv_result(result: CVResult, output_path: Path) -> None:
    """Write CV result to JSON."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "n_folds": result.n_folds,
        "metric": result.metric_name,
        "model_type": result.model_type,
        "mean_score": result.mean_score,
        "std_score": result.std_score,
        "per_fold": result.per_fold_scores,
        "feature_importances": result.feature_importances,
        "fold_test_protein_ids": result.fold_test_protein_ids,
    }
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
