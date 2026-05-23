from __future__ import annotations

import pytest

from lpmo_pipeline.analysis.clustering_pilot import (
    PilotConditionIFP,
    build_feature_prevalence,
    build_filtered_ifp_matrix,
    build_interaction_type_prevalence,
    select_main_clustering_features,
    select_main_clustering_features_for_condition,
    summarize_condition_matrices,
)
from lpmo_pipeline.analysis.prolif_ifp import IFPBatch, IFPResult


def _make_condition(
    condition_id: str,
    *,
    pose_ids: list[str],
    feature_names: list[str],
    matrix: list[list[int]],
    interaction_counts: list[dict[str, int]],
    selected_pose_ids: tuple[str, ...] | None = None,
) -> PilotConditionIFP:
    results = [
        IFPResult(
            pose_id=pose_id,
            status="ok",
            feature_names=feature_names,
            flat_bitvector=row,
            n_total_contacts=sum(row),
            interaction_counts=counts,
        )
        for pose_id, row, counts in zip(pose_ids, matrix, interaction_counts, strict=True)
    ]
    batch = IFPBatch(
        protein_id="P",
        ligand_id="L",
        model="af3",
        results=results,
        matrix=matrix,
        feature_names=feature_names,
    )
    return PilotConditionIFP(
        condition_id=condition_id,
        batch=batch,
        selected_pose_ids=selected_pose_ids,
    )


def test_build_interaction_type_prevalence_uses_selected_pose_subsets() -> None:
    condition_a = _make_condition(
        "cond-a",
        pose_ids=["a0", "a1"],
        feature_names=[
            "NAG1.B|ASN10.A|ImplicitHBDonor",
            "NAG1.B|ASN10.A|VdWContact",
            "NAG2.B|SER11.A|ImplicitHBAcceptor",
        ],
        matrix=[[1, 1, 0], [0, 0, 1]],
        interaction_counts=[
            {"ImplicitHBDonor": 1, "VdWContact": 1},
            {"ImplicitHBAcceptor": 1},
        ],
    )
    condition_b = _make_condition(
        "cond-b",
        pose_ids=["b0", "b1"],
        feature_names=[
            "NAG1.B|ASN10.A|ImplicitHBDonor",
            "NAG2.B|TYR12.A|Hydrophobic",
        ],
        matrix=[[1, 1], [0, 0]],
        interaction_counts=[
            {"ImplicitHBDonor": 2, "Hydrophobic": 1},
            {},
        ],
        selected_pose_ids=("b0",),
    )

    rows = {
        row.interaction_type: row
        for row in build_interaction_type_prevalence([condition_a, condition_b])
    }

    assert rows["ImplicitHBDonor"].n_selected_poses_with_type == 2
    assert rows["ImplicitHBDonor"].pose_prevalence == pytest.approx(2 / 3)
    assert rows["ImplicitHBDonor"].n_conditions_with_type == 2
    assert rows["ImplicitHBDonor"].condition_prevalence == pytest.approx(1.0)
    assert rows["ImplicitHBDonor"].median_count_per_pose_when_present == pytest.approx(1.5)
    assert rows["Hydrophobic"].n_selected_poses_with_type == 1
    assert rows["Hydrophobic"].condition_prevalence == pytest.approx(0.5)


def test_build_feature_prevalence_tracks_pose_and_condition_prevalence() -> None:
    condition_a = _make_condition(
        "cond-a",
        pose_ids=["a0", "a1"],
        feature_names=[
            "NAG1.B|ASN10.A|ImplicitHBDonor",
            "NAG1.B|ASN10.A|VdWContact",
            "NAG2.B|SER11.A|ImplicitHBAcceptor",
        ],
        matrix=[[1, 1, 0], [0, 0, 1]],
        interaction_counts=[
            {"ImplicitHBDonor": 1, "VdWContact": 1},
            {"ImplicitHBAcceptor": 1},
        ],
    )
    condition_b = _make_condition(
        "cond-b",
        pose_ids=["b0", "b1"],
        feature_names=[
            "NAG1.B|ASN10.A|ImplicitHBDonor",
            "NAG2.B|TYR12.A|Hydrophobic",
        ],
        matrix=[[1, 1], [0, 0]],
        interaction_counts=[
            {"ImplicitHBDonor": 1, "Hydrophobic": 1},
            {},
        ],
        selected_pose_ids=("b0",),
    )

    rows = {
        row.feature_name: row
        for row in build_feature_prevalence([condition_a, condition_b])
    }

    donor = rows["NAG1.B|ASN10.A|ImplicitHBDonor"]
    assert donor.interaction_type == "ImplicitHBDonor"
    assert donor.n_selected_poses_with_feature == 2
    assert donor.pose_prevalence == pytest.approx(2 / 3)
    assert donor.n_conditions_with_feature == 2
    assert donor.condition_prevalence == pytest.approx(1.0)
    assert donor.within_condition_max_prevalence == pytest.approx(1.0)
    assert donor.within_condition_variance == pytest.approx(0.0625)

    hydrophobic = rows["NAG2.B|TYR12.A|Hydrophobic"]
    assert hydrophobic.pose_prevalence == pytest.approx(1 / 3)
    assert hydrophobic.condition_prevalence == pytest.approx(0.5)


def test_select_main_clustering_features_applies_interaction_and_rarity_rules() -> None:
    condition = _make_condition(
        "cond-a",
        pose_ids=["a0", "a1", "a2"],
        feature_names=[
            "NAG1.B|ASN10.A|ImplicitHBDonor",
            "NAG1.B|ASN10.A|VdWContact",
            "NAG2.B|TYR12.A|Hydrophobic",
            "NAG3.B|GLU15.A|ImplicitHBAcceptor",
        ],
        matrix=[
            [1, 1, 1, 0],
            [1, 1, 1, 0],
            [0, 0, 0, 1],
        ],
        interaction_counts=[
            {"ImplicitHBDonor": 1, "VdWContact": 1, "Hydrophobic": 1},
            {"ImplicitHBDonor": 1, "VdWContact": 1, "Hydrophobic": 1},
            {"ImplicitHBAcceptor": 1},
        ],
    )

    rows = build_feature_prevalence([condition])
    default_selected = select_main_clustering_features(
        rows,
        rare_feature_pose_prevalence_lt=0.34,
        rare_feature_condition_prevalence_lt_n_conditions=2,
    )
    hydrophobic_selected = select_main_clustering_features(
        rows,
        extra_include_interaction_types={"Hydrophobic"},
        rare_feature_pose_prevalence_lt=0.34,
        rare_feature_condition_prevalence_lt_n_conditions=1,
    )

    assert "NAG1.B|ASN10.A|ImplicitHBDonor" in default_selected
    assert "NAG1.B|ASN10.A|VdWContact" not in default_selected
    assert "NAG2.B|TYR12.A|Hydrophobic" not in default_selected
    assert "NAG3.B|GLU15.A|ImplicitHBAcceptor" not in default_selected
    assert "NAG2.B|TYR12.A|Hydrophobic" in hydrophobic_selected


def test_build_filtered_ifp_matrix_subsets_poses_and_preserves_feature_order() -> None:
    condition = _make_condition(
        "cond-a",
        pose_ids=["a0", "a1"],
        feature_names=[
            "NAG1.B|ASN10.A|ImplicitHBDonor",
            "NAG1.B|ASN10.A|VdWContact",
            "NAG2.B|SER11.A|ImplicitHBAcceptor",
        ],
        matrix=[[1, 1, 0], [0, 0, 1]],
        interaction_counts=[
            {"ImplicitHBDonor": 1, "VdWContact": 1},
            {"ImplicitHBAcceptor": 1},
        ],
    )

    filtered = build_filtered_ifp_matrix(
        condition.batch,
        selected_pose_ids=("a1",),
        allowed_feature_names={
            "NAG2.B|SER11.A|ImplicitHBAcceptor",
            "NAG1.B|ASN10.A|ImplicitHBDonor",
        },
    )

    assert filtered.pose_ids == ["a1"]
    assert filtered.feature_names == [
        "NAG1.B|ASN10.A|ImplicitHBDonor",
        "NAG2.B|SER11.A|ImplicitHBAcceptor",
    ]
    assert filtered.matrix == [[0, 1]]


def test_select_main_clustering_features_for_condition_uses_selected_pose_signal() -> None:
    condition = _make_condition(
        "cond-a",
        pose_ids=["a0", "a1"],
        feature_names=[
            "NAG1.B|ASN10.A|ImplicitHBDonor",
            "NAG1.B|ASN10.A|VdWContact",
            "NAG2.B|TYR12.A|Hydrophobic",
            "NAG3.B|GLU15.A|ImplicitHBAcceptor",
        ],
        matrix=[[1, 1, 1, 0], [0, 0, 0, 1]],
        interaction_counts=[
            {"ImplicitHBDonor": 1, "VdWContact": 1, "Hydrophobic": 1},
            {"ImplicitHBAcceptor": 1},
        ],
    )

    selected = select_main_clustering_features_for_condition(
        condition.batch,
        selected_pose_ids=("a0",),
    )
    selected_with_hydrophobic = select_main_clustering_features_for_condition(
        condition.batch,
        selected_pose_ids=("a0",),
        extra_include_interaction_types={"Hydrophobic"},
    )

    assert selected == ["NAG1.B|ASN10.A|ImplicitHBDonor"]
    assert selected_with_hydrophobic == [
        "NAG1.B|ASN10.A|ImplicitHBDonor",
        "NAG2.B|TYR12.A|Hydrophobic",
    ]


def test_summarize_condition_matrices_marks_insufficient_clusterable_signal() -> None:
    condition = _make_condition(
        "cond-a",
        pose_ids=["a0", "a1"],
        feature_names=[
            "NAG1.B|ASN10.A|HBDonor",
            "NAG2.B|SER11.A|HBAcceptor",
        ],
        matrix=[[1, 0], [0, 1]],
        interaction_counts=[
            {"HBDonor": 1},
            {"HBAcceptor": 1},
        ],
    )

    raw_matrix = build_filtered_ifp_matrix(condition.batch)
    main_matrix = build_filtered_ifp_matrix(
        condition.batch,
        allowed_feature_names={"NAG1.B|ASN10.A|HBDonor"},
    )
    summary = summarize_condition_matrices(
        condition.condition_id,
        raw_matrix,
        main_matrix,
        minimum_clusterable_n=3,
    )

    assert summary.condition_id == "cond-a"
    assert summary.n_selected_poses == 2
    assert summary.raw_feature_count == 2
    assert summary.main_feature_count == 1
    assert summary.clustering_status == "insufficient_clusterable_signal"
    assert summary.formal_clustering_allowed is False


def test_summarize_condition_matrices_marks_empty_main_matrix() -> None:
    condition = _make_condition(
        "cond-b",
        pose_ids=["b0", "b1", "b2"],
        feature_names=["NAG1.B|ASN10.A|VdWContact"],
        matrix=[[1], [1], [0]],
        interaction_counts=[
            {"VdWContact": 1},
            {"VdWContact": 1},
            {},
        ],
    )

    raw_matrix = build_filtered_ifp_matrix(condition.batch)
    main_matrix = build_filtered_ifp_matrix(condition.batch, allowed_feature_names=set())
    summary = summarize_condition_matrices(
        condition.condition_id,
        raw_matrix,
        main_matrix,
        minimum_clusterable_n=3,
    )

    assert summary.clustering_status == "empty_main_matrix"
    assert summary.formal_clustering_allowed is False
