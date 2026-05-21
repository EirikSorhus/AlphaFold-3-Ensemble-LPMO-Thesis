# src/lpmo_pipeline/analysis/predictive_models.py
"""
Responsibility: Provisional predictive-modeling scaffold from cluster features.
Input:  predictive_cluster_table (cluster-level rows + activity labels)
Output: Model performance metrics, feature importances

Status:
    This module is not the final predictive-analysis implementation. It
    currently provides tested leakage-safety utilities and narrow baseline
    runners while the project plan is revised to select final model families
    and predictor variables.

Anti p-hack rules:
  - HDBSCAN params LOCKED (from tuning, not re-tuned here)
    - Grouped CV: fold by protein_id (all clusters from a protein remain together)
    - Stratified CV: maintain class balance (C1/C4/mixed, polymer type) at protein level
  - CV hierarchy defined in configs/cv_hierarchy.yaml
  - Default CV target is 5 folds when enough protein groups exist.

Main rule:
    - Keep cluster rows as primary modeling rows.
    - Enzyme-level aggregation is sensitivity analysis only.
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
from sklearn.metrics import balanced_accuracy_score
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
            result.per_fold_scores.append(float(balanced_accuracy_score(y_test, y_pred)))
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


def run_predictive_model(
    features: dict[str, list[float]],   # {cluster_id: feature_vector}
    labels: dict[str, str],             # {cluster_id: activity_class}
    cluster_to_protein: dict[str, str],  # {cluster_id: protein_id}
    folds: list[CVFold],
    model_type: str = "random_forest",
) -> CVResult:
    """Run predictive modeling with grouped CV.

    PSEUDOCODE — wraps sklearn classifiers.

    Args:
        features: Feature vectors per cluster row.
        labels: Activity labels per cluster row.
        cluster_to_protein: Mapping to enforce grouped CV at protein level.
        folds: Pre-built protein-level CV folds.
        model_type: "random_forest" | "logistic_regression"

    Returns:
        CVResult with per-fold scores and feature importances.
    """
    result = CVResult(n_folds=len(folds), metric_name="balanced_accuracy")

    # PSEUDOCODE:
    # from sklearn.ensemble import RandomForestClassifier
    # from sklearn.metrics import balanced_accuracy_score
    #
    # importances_accumulator = np.zeros(n_features)
    # for fold in folds:
    #     train_clusters = [
    #         cid for cid, protein_id in cluster_to_protein.items()
    #         if protein_id in fold.train_protein_ids
    #     ]
    #     test_clusters = [
    #         cid for cid, protein_id in cluster_to_protein.items()
    #         if protein_id in fold.test_protein_ids
    #     ]
    #     X_train = np.array([features[cid] for cid in train_clusters])
    #     y_train = np.array([labels[cid] for cid in train_clusters])
    #     X_test = np.array([features[cid] for cid in test_clusters])
    #     y_test = np.array([labels[cid] for cid in test_clusters])
    #
    #     clf = RandomForestClassifier(n_estimators=100, random_state=42)
    #     clf.fit(X_train, y_train)
    #     y_pred = clf.predict(X_test)
    #     score = balanced_accuracy_score(y_test, y_pred)
    #     result.per_fold_scores.append(score)
    #     importances_accumulator += clf.feature_importances_
    #
    # result.mean_score = np.mean(result.per_fold_scores)
    # result.std_score = np.std(result.per_fold_scores)
    # result.feature_importances = {f"feature_{i}": v for i, v in enumerate(importances_accumulator / len(folds))}

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
