"""Utilities for clustering-pilot feature audits and filtered IFP matrices.

This module is intentionally additive: it builds pilot-specific summaries and
clustering inputs from the existing ``IFPBatch`` contract without changing the
standard production path.
"""
from __future__ import annotations

from dataclasses import dataclass
from statistics import median, pvariance
from typing import Collection, Iterable, Sequence

from lpmo_pipeline.analysis.prolif_ifp import ContactEligibility, IFPBatch, parse_ifp_feature_name

DEFAULT_MAIN_CLUSTERING_INTERACTION_TYPES = (
    "HBDonor",
    "HBAcceptor",
    "PiStacking",
)
DEFAULT_EXCLUDED_CLUSTERING_INTERACTION_TYPES = ("VdWContact",)
DEFAULT_MINIMUM_CLUSTERABLE_N = 10
DEFAULT_INSUFFICIENT_CLUSTERABLE_SIGNAL_LABEL = "insufficient_clusterable_signal"


@dataclass(frozen=True)
class PilotConditionIFP:
    """One condition's aligned IFP batch plus the selected pilot pose subset."""

    condition_id: str
    batch: IFPBatch
    selected_pose_ids: tuple[str, ...] | None = None
    contact_eligibilities: tuple[ContactEligibility, ...] = ()
    n_qc_pass_poses: int | None = None


@dataclass(frozen=True)
class InteractionTypePrevalence:
    """Pilot-wide prevalence summary for one interaction type."""

    interaction_type: str
    n_selected_poses_with_type: int
    pose_prevalence: float
    n_conditions_with_type: int
    condition_prevalence: float
    median_count_per_pose_when_present: float

    def to_row(self) -> dict[str, float | int | str]:
        return {
            "interaction_type": self.interaction_type,
            "n_selected_poses_with_type": self.n_selected_poses_with_type,
            "pose_prevalence": self.pose_prevalence,
            "n_conditions_with_type": self.n_conditions_with_type,
            "condition_prevalence": self.condition_prevalence,
            "median_count_per_pose_when_present": self.median_count_per_pose_when_present,
        }


@dataclass(frozen=True)
class FeaturePrevalence:
    """Pilot-wide prevalence summary for one flattened IFP feature."""

    feature_name: str
    interaction_type: str
    n_selected_poses_with_feature: int
    pose_prevalence: float
    n_conditions_with_feature: int
    condition_prevalence: float
    within_condition_max_prevalence: float
    within_condition_variance: float

    def to_row(self) -> dict[str, float | int | str]:
        return {
            "feature_name": self.feature_name,
            "interaction_type": self.interaction_type,
            "n_selected_poses_with_feature": self.n_selected_poses_with_feature,
            "pose_prevalence": self.pose_prevalence,
            "n_conditions_with_feature": self.n_conditions_with_feature,
            "condition_prevalence": self.condition_prevalence,
            "within_condition_max_prevalence": self.within_condition_max_prevalence,
            "within_condition_variance": self.within_condition_variance,
        }


@dataclass(frozen=True)
class FilteredIFPMatrix:
    """Pilot-specific clustering matrix built from a filtered feature set."""

    pose_ids: list[str]
    feature_names: list[str]
    matrix: list[list[int]]


@dataclass(frozen=True)
class PilotConditionMatrixSummary:
    """Pilot-level summary for whether a condition can enter formal clustering."""

    condition_id: str
    n_selected_poses: int
    raw_feature_count: int
    main_feature_count: int
    minimum_clusterable_n: int
    clustering_status: str
    formal_clustering_allowed: bool

    def to_row(self) -> dict[str, int | bool | str]:
        return {
            "condition_id": self.condition_id,
            "n_selected_poses": self.n_selected_poses,
            "raw_feature_count": self.raw_feature_count,
            "main_feature_count": self.main_feature_count,
            "minimum_clusterable_n": self.minimum_clusterable_n,
            "clustering_status": self.clustering_status,
            "formal_clustering_allowed": self.formal_clustering_allowed,
        }


def _resolve_selected_indices(batch: IFPBatch, selected_pose_ids: Sequence[str] | None) -> list[int]:
    if selected_pose_ids is None:
        return list(range(len(batch.results)))

    selected_set = set(selected_pose_ids)
    available_pose_ids = {result.pose_id for result in batch.results}
    missing_pose_ids = sorted(selected_set - available_pose_ids)
    if missing_pose_ids:
        raise ValueError(f"Selected pose_ids missing from batch: {missing_pose_ids}")

    return [index for index, result in enumerate(batch.results) if result.pose_id in selected_set]


def _prepare_conditions(
    conditions: Iterable[PilotConditionIFP],
) -> list[tuple[PilotConditionIFP, list[int]]]:
    prepared: list[tuple[PilotConditionIFP, list[int]]] = []
    for condition in conditions:
        selected_indices = _resolve_selected_indices(condition.batch, condition.selected_pose_ids)
        if selected_indices:
            prepared.append((condition, selected_indices))
    return prepared


def build_interaction_type_prevalence(
    conditions: Iterable[PilotConditionIFP],
) -> list[InteractionTypePrevalence]:
    """Aggregate pilot-wide prevalence for each interaction type.

    The summaries are computed over the explicitly selected pose subsets for
    each condition.
    """

    prepared = _prepare_conditions(conditions)
    total_selected_poses = sum(len(indices) for _, indices in prepared)
    total_conditions = len(prepared)
    if total_selected_poses == 0 or total_conditions == 0:
        return []

    interaction_types: set[str] = set()
    for condition, selected_indices in prepared:
        for index in selected_indices:
            result = condition.batch.results[index]
            interaction_types.update(result.interaction_counts.keys())
            interaction_types.update(result.interaction_types)

    summaries: list[InteractionTypePrevalence] = []
    for interaction_type in sorted(interaction_types):
        n_selected_poses_with_type = 0
        n_conditions_with_type = 0
        counts_when_present: list[int] = []

        for condition, selected_indices in prepared:
            condition_has_type = False
            for index in selected_indices:
                count = int(condition.batch.results[index].interaction_counts.get(interaction_type, 0))
                if count <= 0:
                    continue
                n_selected_poses_with_type += 1
                counts_when_present.append(count)
                condition_has_type = True
            if condition_has_type:
                n_conditions_with_type += 1

        if n_selected_poses_with_type == 0:
            continue

        summaries.append(
            InteractionTypePrevalence(
                interaction_type=interaction_type,
                n_selected_poses_with_type=n_selected_poses_with_type,
                pose_prevalence=n_selected_poses_with_type / total_selected_poses,
                n_conditions_with_type=n_conditions_with_type,
                condition_prevalence=n_conditions_with_type / total_conditions,
                median_count_per_pose_when_present=float(median(counts_when_present)),
            )
        )

    return summaries


def build_feature_prevalence(
    conditions: Iterable[PilotConditionIFP],
) -> list[FeaturePrevalence]:
    """Aggregate pilot-wide prevalence for each flattened IFP feature."""

    prepared = _prepare_conditions(conditions)
    total_selected_poses = sum(len(indices) for _, indices in prepared)
    total_conditions = len(prepared)
    if total_selected_poses == 0 or total_conditions == 0:
        return []

    active_feature_names: set[str] = set()
    for condition, selected_indices in prepared:
        for index in selected_indices:
            row = condition.batch.matrix[index]
            for feature_name, value in zip(condition.batch.feature_names, row, strict=True):
                if int(value):
                    active_feature_names.add(feature_name)

    summaries: list[FeaturePrevalence] = []
    for feature_name in sorted(active_feature_names):
        _, _, interaction_type = parse_ifp_feature_name(feature_name)
        n_selected_poses_with_feature = 0
        n_conditions_with_feature = 0
        within_condition_prevalences: list[float] = []

        for condition, selected_indices in prepared:
            try:
                feature_index = condition.batch.feature_names.index(feature_name)
            except ValueError:
                feature_index = None

            feature_count = 0
            if feature_index is not None:
                feature_count = sum(int(condition.batch.matrix[index][feature_index]) for index in selected_indices)

            if feature_count > 0:
                n_selected_poses_with_feature += feature_count
                n_conditions_with_feature += 1

            within_condition_prevalences.append(feature_count / len(selected_indices))

        summaries.append(
            FeaturePrevalence(
                feature_name=feature_name,
                interaction_type=interaction_type,
                n_selected_poses_with_feature=n_selected_poses_with_feature,
                pose_prevalence=n_selected_poses_with_feature / total_selected_poses,
                n_conditions_with_feature=n_conditions_with_feature,
                condition_prevalence=n_conditions_with_feature / total_conditions,
                within_condition_max_prevalence=max(within_condition_prevalences, default=0.0),
                within_condition_variance=float(
                    pvariance(within_condition_prevalences) if within_condition_prevalences else 0.0
                ),
            )
        )

    return summaries


def select_main_clustering_features(
    feature_prevalence: Iterable[FeaturePrevalence],
    *,
    default_include_interaction_types: Collection[str] = DEFAULT_MAIN_CLUSTERING_INTERACTION_TYPES,
    extra_include_interaction_types: Collection[str] = (),
    excluded_interaction_types: Collection[str] = DEFAULT_EXCLUDED_CLUSTERING_INTERACTION_TYPES,
    rare_feature_pose_prevalence_lt: float = 0.01,
    rare_feature_condition_prevalence_lt_n_conditions: int = 2,
) -> list[str]:
    """Select the pilot's main clustering features from prevalence summaries."""

    included_interaction_types = set(default_include_interaction_types) | set(extra_include_interaction_types)
    excluded_interaction_types = set(excluded_interaction_types)

    selected_feature_names: list[str] = []
    for feature in feature_prevalence:
        if feature.interaction_type in excluded_interaction_types:
            continue
        if feature.interaction_type not in included_interaction_types:
            continue
        if (
            feature.pose_prevalence < rare_feature_pose_prevalence_lt
            and feature.n_conditions_with_feature < rare_feature_condition_prevalence_lt_n_conditions
        ):
            continue
        selected_feature_names.append(feature.feature_name)

    return sorted(selected_feature_names)


def build_filtered_ifp_matrix(
    batch: IFPBatch,
    *,
    selected_pose_ids: Sequence[str] | None = None,
    allowed_feature_names: Collection[str] | None = None,
) -> FilteredIFPMatrix:
    """Build a pose-subsetted, feature-filtered matrix for pilot clustering."""

    selected_indices = _resolve_selected_indices(batch, selected_pose_ids)
    allowed_feature_set = set(allowed_feature_names) if allowed_feature_names is not None else None
    feature_indices = [
        index
        for index, feature_name in enumerate(batch.feature_names)
        if allowed_feature_set is None or feature_name in allowed_feature_set
    ]

    return FilteredIFPMatrix(
        pose_ids=[batch.results[index].pose_id for index in selected_indices],
        feature_names=[batch.feature_names[index] for index in feature_indices],
        matrix=[
            [int(batch.matrix[row_index][feature_index]) for feature_index in feature_indices]
            for row_index in selected_indices
        ],
    )


def summarize_condition_matrices(
    condition_id: str,
    raw_matrix: FilteredIFPMatrix,
    main_matrix: FilteredIFPMatrix,
    *,
    minimum_clusterable_n: int = DEFAULT_MINIMUM_CLUSTERABLE_N,
    insufficient_clusterable_signal_label: str = DEFAULT_INSUFFICIENT_CLUSTERABLE_SIGNAL_LABEL,
) -> PilotConditionMatrixSummary:
    """Summarize whether a pilot condition can enter formal clustering."""

    n_selected_poses = len(raw_matrix.pose_ids)
    if n_selected_poses < minimum_clusterable_n:
        clustering_status = insufficient_clusterable_signal_label
        formal_clustering_allowed = False
    elif len(main_matrix.feature_names) == 0:
        clustering_status = "empty_main_matrix"
        formal_clustering_allowed = False
    else:
        clustering_status = "ok"
        formal_clustering_allowed = True

    return PilotConditionMatrixSummary(
        condition_id=condition_id,
        n_selected_poses=n_selected_poses,
        raw_feature_count=len(raw_matrix.feature_names),
        main_feature_count=len(main_matrix.feature_names),
        minimum_clusterable_n=minimum_clusterable_n,
        clustering_status=clustering_status,
        formal_clustering_allowed=formal_clustering_allowed,
    )