# src/lpmo_pipeline/analysis/predictive_models.py
"""
Responsibility: Optional predictive modeling from cluster features.
Input:  activity_features.json (protein-level features + activity labels)
Output: Model performance metrics, feature importances

Anti p-hack rules:
  - HDBSCAN params LOCKED (from tuning, not re-tuned here)
  - Grouped CV: fold by structure, not by pose
  - Stratified CV: maintain class balance (C1/C4/mixed, polymer type)
  - CV hierarchy defined in configs/cv_hierarchy.yaml
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

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
    per_fold_scores: list[float] = field(default_factory=list)
    mean_score: float = 0.0
    std_score: float = 0.0
    feature_importances: dict[str, float] = field(default_factory=dict)


def build_grouped_stratified_cv(
    protein_ids: list[str],
    activity_labels: dict[str, str],     # {protein_id: "C1"|"C4"|"C1+C4"}
    hierarchy: dict[str, str],           # {protein_id: structure_group}
    n_folds: int = 5,
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
    # --- Step 1: Group proteins by hierarchy ---
    groups: dict[str, list[str]] = {}
    for pid in protein_ids:
        group = hierarchy.get(pid, pid)  # Default: each protein is its own group
        groups.setdefault(group, []).append(pid)

    group_keys = sorted(groups.keys())

    # --- Step 2: Assign group labels for stratification ---
    # Label each group by majority class of its proteins
    group_labels: dict[str, str] = {}
    for gkey, pids in groups.items():
        labels = [activity_labels.get(pid, "unknown") for pid in pids]
        # Majority vote
        from collections import Counter
        most_common = Counter(labels).most_common(1)[0][0]
        group_labels[gkey] = most_common

    # --- Step 3: Stratified group split ---
    # PSEUDOCODE: Use sklearn StratifiedKFold on group-level
    # from sklearn.model_selection import StratifiedKFold
    # skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)
    # group_array = np.array(group_keys)
    # label_array = np.array([group_labels[g] for g in group_keys])
    # folds = []
    # for fold_id, (train_idx, test_idx) in enumerate(skf.split(group_array, label_array)):
    #     train_groups = group_array[train_idx]
    #     test_groups = group_array[test_idx]
    #     train_pids = [pid for g in train_groups for pid in groups[g]]
    #     test_pids = [pid for g in test_groups for pid in groups[g]]
    #     folds.append(CVFold(fold_id=fold_id, train_protein_ids=train_pids, test_protein_ids=test_pids))

    folds: list[CVFold] = []  # PSEUDOCODE placeholder
    logger.info("Built %d grouped+stratified CV folds from %d groups", n_folds, len(groups))
    return folds


def run_predictive_model(
    features: dict[str, list[float]],   # {protein_id: feature_vector}
    labels: dict[str, str],             # {protein_id: activity_class}
    folds: list[CVFold],
    model_type: str = "random_forest",
) -> CVResult:
    """Run predictive modeling with grouped CV.

    PSEUDOCODE — wraps sklearn classifiers.

    Args:
        features: Feature vectors per protein.
        labels: Activity labels per protein.
        folds: Pre-built CV folds.
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
    #     X_train = np.array([features[pid] for pid in fold.train_protein_ids])
    #     y_train = np.array([labels[pid] for pid in fold.train_protein_ids])
    #     X_test = np.array([features[pid] for pid in fold.test_protein_ids])
    #     y_test = np.array([labels[pid] for pid in fold.test_protein_ids])
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
        "Predictive model: %s CV score = %.3f ± %.3f",
        model_type, result.mean_score, result.std_score,
    )
    return result


def write_cv_result(result: CVResult, output_path: Path) -> None:
    """Write CV result to JSON."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "n_folds": result.n_folds,
        "metric": result.metric_name,
        "mean_score": result.mean_score,
        "std_score": result.std_score,
        "per_fold": result.per_fold_scores,
        "feature_importances": result.feature_importances,
    }
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
