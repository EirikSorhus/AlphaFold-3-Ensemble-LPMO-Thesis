from __future__ import annotations

from lpmo_pipeline.analysis.predictive_models import (
    build_grouped_stratified_cv,
    run_grouped_binary_logistic_cv,
    run_binary_predictive_baselines,
)


def test_grouped_cv_does_not_split_protein_groups() -> None:
    folds = build_grouped_stratified_cv(
        ["P1", "P2", "P3", "P4"],
        {"P1": "C1", "P2": "C1", "P3": "C4", "P4": "C4"},
        {"P1": "G1", "P2": "G1", "P3": "G2", "P4": "G3"},
        n_folds=2,
        random_state=1,
    )

    assert len(folds) == 2
    for fold in folds:
        assert set(fold.train_protein_ids).isdisjoint(fold.test_protein_ids)
        assert not ({"P1", "P2"} & set(fold.train_protein_ids) and {"P1", "P2"} & set(fold.test_protein_ids))


def test_grouped_cv_defaults_to_five_folds_when_enough_groups() -> None:
    protein_ids = [f"P{index}" for index in range(10)]
    labels = {
        protein_id: "C1" if index % 2 == 0 else "C4"
        for index, protein_id in enumerate(protein_ids)
    }

    folds = build_grouped_stratified_cv(
        protein_ids,
        labels,
        {protein_id: protein_id for protein_id in protein_ids},
    )

    assert len(folds) == 5
    for fold in folds:
        assert set(fold.train_protein_ids).isdisjoint(fold.test_protein_ids)


def test_binary_predictive_baselines_use_grouped_splits_without_leakage() -> None:
    rows = [
        {
            "protein_id": "P1",
            "experimental_regio_label": "C1",
            "occupancy": 0.9,
            "median_Cu_C1": 2.5,
            "median_Cu_C4": 6.0,
        },
        {
            "protein_id": "P2",
            "experimental_regio_label": "C1",
            "occupancy": 0.8,
            "median_Cu_C1": 2.8,
            "median_Cu_C4": 5.8,
        },
        {
            "protein_id": "P3",
            "experimental_regio_label": "C4",
            "occupancy": 0.7,
            "median_Cu_C1": 6.1,
            "median_Cu_C4": 2.7,
        },
        {
            "protein_id": "P4",
            "experimental_regio_label": "C4",
            "occupancy": 0.6,
            "median_Cu_C1": 5.7,
            "median_Cu_C4": 2.4,
        },
    ]

    results = run_binary_predictive_baselines(
        rows,
        feature_columns=["occupancy", "median_Cu_C1", "median_Cu_C4"],
        n_folds=2,
        random_state=4,
    )

    assert {result.model_type for result in results} == {"logistic_regression", "decision_tree"}
    assert all(result.n_folds == 2 for result in results)
    labels_by_protein = {row["protein_id"]: row["experimental_regio_label"] for row in rows}
    for result in results:
        assert all(0.0 <= score <= 1.0 for score in result.per_fold_scores)
        assert set(result.feature_importances) == {"occupancy", "median_Cu_C1", "median_Cu_C4"}
        for test_proteins in result.fold_test_protein_ids:
            assert len(test_proteins) == 2
            assert {labels_by_protein[protein_id] for protein_id in test_proteins} == {"C1", "C4"}


def test_grouped_binary_logistic_cv_emits_predictions_and_metrics() -> None:
    rows = [
        {"model_row_id": "P1", "protein_id": "P1", "target": 1, "f1": 2.0, "f2": 1.0},
        {"model_row_id": "P2", "protein_id": "P2", "target": 1, "f1": 1.8, "f2": 0.9},
        {"model_row_id": "P3", "protein_id": "P3", "target": 0, "f1": -1.7, "f2": -1.0},
        {"model_row_id": "P4", "protein_id": "P4", "target": 0, "f1": -2.1, "f2": -0.8},
    ]

    result = run_grouped_binary_logistic_cv(
        rows,
        feature_columns=["f1", "f2"],
        target_column="target",
        model_name="toy_model",
        n_folds=2,
        random_state=7,
    )

    assert result.model_name == "toy_model"
    assert result.target_column == "target"
    assert result.n_rows == 4
    assert result.n_positive == 2
    assert result.n_negative == 2
    assert result.n_folds == 2
    assert len(result.fold_metrics) == 2
    assert len(result.predictions) == 4
    assert set(result.feature_importances) == {"f1", "f2"}
    assert 0.0 <= result.mean_metrics["balanced_accuracy"] <= 1.0
    assert 0.0 <= result.mean_metrics["pr_auc"] <= 1.0
    assert 0.0 <= result.mean_metrics["recall_at_k"] <= 1.0
