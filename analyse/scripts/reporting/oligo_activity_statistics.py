#!/usr/bin/env python3
"""Summarize prediction outcomes by experimentally observed oligo activity."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

try:
    from scipy import stats
except ImportError:  # pragma: no cover - scipy is available in the analysis env, but keep outputs robust.
    stats = None

try:
    from sklearn.metrics import average_precision_score, roc_auc_score
except ImportError:  # pragma: no cover - sklearn is available in the analysis env, but keep outputs robust.
    average_precision_score = None
    roc_auc_score = None


CONSTRUCT_LABELS = {
    "domain_only": "Catalytic domain",
    "full_length": "Full-length",
    "all_constructs": "All constructs",
}


NUMERIC_INPUT_COLUMNS = [
    "n_generated",
    "n_stage1_pass",
    "n_stage1_hard_fail",
    "n_ifp_success",
    "n_contact_eligible",
    "n_ifp_clustered",
    "n_noise",
    "n_clusters",
    "top_cluster_occupancy",
    "cluster_entropy",
    "convergence_fraction",
    "median_ligand_rmsd",
    "contact_eligible_fraction",
    "noise_fraction",
    "occupancy_weighted_c1_plausible_fraction",
    "occupancy_weighted_c4_plausible_fraction",
    "occupancy_weighted_Cu_C1_distance_median",
    "occupancy_weighted_Cu_C4_distance_median",
    "occupancy_weighted_oxyl_H_C1_distance_median",
    "occupancy_weighted_oxyl_H_C4_distance_median",
    "occupancy_weighted_Cu_oxyl_H_C1_angle_median",
    "occupancy_weighted_Cu_oxyl_H_C4_angle_median",
    "occupancy_weighted_oxyl_H_score_C1_median",
    "occupancy_weighted_oxyl_H_score_C4_median",
    "weighted_contact_feature_mass",
    "aromatic_contact_fraction",
    "polar_contact_fraction",
    "hbond_contact_fraction",
    "catalytic_surface_contact_fraction",
]


METRIC_INFO = {
    "n_conditions": ("Number of predicted substrate conditions per protein/construct", "count", "descriptive"),
    "hard_qc_pass_fraction": ("Hard-QC pass fraction", "fraction", "higher is better"),
    "hard_qc_fail_fraction": ("Hard-QC hard-fail fraction", "fraction", "lower is better"),
    "ifp_success_fraction": ("IFP success fraction after hard QC", "fraction", "higher is better"),
    "contact_eligible_fraction_weighted": ("Contact-eligible pose fraction", "fraction", "higher is better"),
    "clusterable_pose_fraction": ("Clusterable pose fraction among IFP-success poses", "fraction", "higher is better"),
    "formal_clusterable_condition_fraction": ("Fraction of conditions allowed into formal clustering", "fraction", "higher is better"),
    "valid_cluster_condition_fraction": ("Fraction of conditions with at least one valid cluster", "fraction", "higher is better"),
    "noise_fraction_weighted": ("Noise fraction among clusterable poses", "fraction", "lower is better"),
    "n_clusters_mean": ("Mean number of clusters per condition", "clusters", "descriptive"),
    "n_clusters_median": ("Median number of clusters per condition", "clusters", "descriptive"),
    "top_cluster_occupancy_median": ("Median top-cluster occupancy", "fraction", "higher means more dominant cluster"),
    "cluster_entropy_median": ("Median cluster entropy", "entropy", "higher means more diffuse clustering"),
    "convergence_fraction_median": ("Median convergence fraction", "fraction", "higher is better"),
    "median_substrate_rmsd_median": ("Median within-condition substrate RMSD", "Å", "lower is better"),
    "c1_plausible_fraction_median": ("Median C1 plausible fraction", "fraction", "higher is better"),
    "c4_plausible_fraction_median": ("Median C4 plausible fraction", "fraction", "higher is better"),
    "cu_c1_distance_median": ("Median Cu-C1 distance", "Å", "lower is better"),
    "cu_c4_distance_median": ("Median Cu-C4 distance", "Å", "lower is better"),
    "oxyl_h_c1_distance_median": ("Median Oxyl-H(C1) distance", "Å", "lower is better"),
    "oxyl_h_c4_distance_median": ("Median Oxyl-H(C4) distance", "Å", "lower is better"),
    "cu_oxyl_h_c1_angle_median": ("Median Cu-Oxyl-H(C1) angle", "degrees", "geometry dependent"),
    "cu_oxyl_h_c4_angle_median": ("Median Cu-Oxyl-H(C4) angle", "degrees", "geometry dependent"),
    "oxyl_h_c1_score_median": ("Median Oxyl-H(C1) geometry score", "score", "higher is better"),
    "oxyl_h_c4_score_median": ("Median Oxyl-H(C4) geometry score", "score", "higher is better"),
    "weighted_contact_feature_mass_median": ("Median weighted contact feature mass", "score", "higher means more contact signal"),
    "aromatic_contact_fraction_median": ("Median aromatic contact fraction", "fraction", "descriptive"),
    "polar_contact_fraction_median": ("Median polar contact fraction", "fraction", "descriptive"),
    "hbond_contact_fraction_median": ("Median hydrogen-bond contact fraction", "fraction", "descriptive"),
    "catalytic_surface_contact_fraction_median": ("Median catalytic-surface contact fraction", "fraction", "descriptive"),
    "best_plausible_fraction_median": ("Median best C1/C4 plausible fraction", "fraction", "higher is better"),
    "best_cu_distance_median": ("Median best Cu-C1/C4 distance", "Å", "lower is better"),
    "best_oxyl_h_distance_median": ("Median best Oxyl-H(C1/C4) distance", "Å", "lower is better"),
    "best_oxyl_h_score_median": ("Median best Oxyl-H(C1/C4) geometry score", "score", "higher is better"),
}


SHORT_MEDIAN_COLUMNS = {
    "top_cluster_occupancy": "top_cluster_occupancy_median",
    "cluster_entropy": "cluster_entropy_median",
    "convergence_fraction": "convergence_fraction_median",
    "median_ligand_rmsd": "median_substrate_rmsd_median",
    "occupancy_weighted_c1_plausible_fraction": "c1_plausible_fraction_median",
    "occupancy_weighted_c4_plausible_fraction": "c4_plausible_fraction_median",
    "occupancy_weighted_Cu_C1_distance_median": "cu_c1_distance_median",
    "occupancy_weighted_Cu_C4_distance_median": "cu_c4_distance_median",
    "occupancy_weighted_oxyl_H_C1_distance_median": "oxyl_h_c1_distance_median",
    "occupancy_weighted_oxyl_H_C4_distance_median": "oxyl_h_c4_distance_median",
    "occupancy_weighted_Cu_oxyl_H_C1_angle_median": "cu_oxyl_h_c1_angle_median",
    "occupancy_weighted_Cu_oxyl_H_C4_angle_median": "cu_oxyl_h_c4_angle_median",
    "occupancy_weighted_oxyl_H_score_C1_median": "oxyl_h_c1_score_median",
    "occupancy_weighted_oxyl_H_score_C4_median": "oxyl_h_c4_score_median",
    "weighted_contact_feature_mass": "weighted_contact_feature_mass_median",
    "aromatic_contact_fraction": "aromatic_contact_fraction_median",
    "polar_contact_fraction": "polar_contact_fraction_median",
    "hbond_contact_fraction": "hbond_contact_fraction_median",
    "catalytic_surface_contact_fraction": "catalytic_surface_contact_fraction_median",
    "best_plausible_fraction_condition": "best_plausible_fraction_median",
    "best_cu_distance_condition": "best_cu_distance_median",
    "best_oxyl_h_distance_condition": "best_oxyl_h_distance_median",
    "best_oxyl_h_score_condition": "best_oxyl_h_score_median",
}

PRESENTABLE_METRICS = [
    "hard_qc_pass_fraction",
    "clusterable_pose_fraction",
    "noise_fraction_weighted",
    "valid_cluster_condition_fraction",
    "best_plausible_fraction_median",
    "best_cu_distance_median",
    "best_oxyl_h_distance_median",
    "best_oxyl_h_score_median",
]

PRESENTABLE_LABELS = {
    "Oligo activity true": "Oligo active",
    "Oligo activity false": "Not oligo active",
    "No updated oligo activity record": "No updated record",
}

N_PERMUTATIONS = 10000
PERMUTATION_SEED = 20260603

OLIGO_2X2_ORDER = [
    "Oligo active demonstrated activity",
    "Oligo active not demonstrated activity",
    "All other demonstrated activity",
    "All other no demonstrated activity",
]

OLIGO_2X2_SHORT_LABELS = {
    "Oligo active demonstrated activity": "Oligo active\nDemonstrated",
    "Oligo active not demonstrated activity": "Oligo active\nNot\ndemonstrated",
    "All other demonstrated activity": "All other\nDemonstrated",
    "All other no demonstrated activity": "All other\nNot\ndemonstrated",
}

OLIGO_2X2_PALETTE = {
    "Oligo active demonstrated activity": "#1f77b4",
    "Oligo active not demonstrated activity": "#86b6d8",
    "All other demonstrated activity": "#9b5b2e",
    "All other no demonstrated activity": "#c8a36a",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--analyse-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def semicolon_join(values: Iterable[object]) -> str:
    cleaned = sorted({str(v).strip() for v in values if str(v).strip() and str(v).strip().lower() != "nan"})
    return "; ".join(cleaned)


def to_bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(["true", "1", "yes", "y"])


def to_ja_bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(["ja", "true", "1", "yes", "y"])


def safe_ratio(numerator: float, denominator: float) -> float:
    if not np.isfinite(numerator) or not np.isfinite(denominator) or denominator == 0:
        return float("nan")
    return float(numerator / denominator)


def describe(values: pd.Series) -> dict[str, float | int]:
    arr = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size == 0:
        return {
            "n": 0,
            "mean": np.nan,
            "sd": np.nan,
            "median": np.nan,
            "q1": np.nan,
            "q3": np.nan,
            "min": np.nan,
            "max": np.nan,
        }
    return {
        "n": int(arr.size),
        "mean": float(np.mean(arr)),
        "sd": float(np.std(arr, ddof=1)) if arr.size > 1 else np.nan,
        "median": float(np.median(arr)),
        "q1": float(np.quantile(arr, 0.25)),
        "q3": float(np.quantile(arr, 0.75)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
    }


def cliffs_delta(x: pd.Series, y: pd.Series) -> float:
    x_arr = pd.to_numeric(x, errors="coerce").dropna().to_numpy(dtype=float)
    y_arr = pd.to_numeric(y, errors="coerce").dropna().to_numpy(dtype=float)
    if x_arr.size == 0 or y_arr.size == 0:
        return float("nan")
    greater = 0
    less = 0
    for value in x_arr:
        greater += int(np.sum(value > y_arr))
        less += int(np.sum(value < y_arr))
    return float((greater - less) / (x_arr.size * y_arr.size))


def make_protein_labels(
    condition: pd.DataFrame,
    pdb_data: pd.DataFrame,
    metadata: pd.DataFrame,
    oligo_activity: pd.DataFrame,
) -> pd.DataFrame:
    run_proteins = pd.DataFrame({"protein_id": sorted(condition["protein_id"].dropna().astype(str).unique())})

    pdb = pdb_data.copy()
    pdb = pdb.rename(columns={"Uniprot_ID": "protein_id", "Family": "pdb_family", "EC#": "pdb_ec_numbers"})
    pdb["protein_id"] = pdb["protein_id"].astype(str)
    pdb["oligo_activity_norm"] = pdb["Oligo_Activity"].fillna("Unknown").astype(str).str.strip()
    pdb["oligo_activity_true_row"] = pdb["oligo_activity_norm"].str.lower().eq("true")
    pdb["resolution_numeric"] = pd.to_numeric(pdb["Resolution"], errors="coerce")

    pdb_by_protein = (
        pdb.groupby("protein_id", dropna=False)
        .agg(
            pdb_family=("pdb_family", semicolon_join),
            pdb_ec_numbers=("pdb_ec_numbers", semicolon_join),
            pdb_ids=("PDB", semicolon_join),
            carbohydrate_ligands=("Carbohydrate_Ligands", semicolon_join),
            pdb_dp_values=("DP", semicolon_join),
            oligo_activity_values=("oligo_activity_norm", semicolon_join),
            pdb_oligo_activity_true=("oligo_activity_true_row", "max"),
            n_pdb_rows=("PDB", "size"),
            best_resolution=("resolution_numeric", "min"),
            median_resolution=("resolution_numeric", "median"),
        )
        .reset_index()
    )

    meta = metadata.rename(columns={"UniProt_ID": "protein_id", "CAZy_family": "metadata_family", "EC_Number": "metadata_ec_numbers"}).copy()
    meta["protein_id"] = meta["protein_id"].astype(str)
    meta_by_protein = (
        meta.groupby("protein_id", dropna=False)
        .agg(
            metadata_family=("metadata_family", semicolon_join),
            metadata_ec_numbers=("metadata_ec_numbers", semicolon_join),
            protein_name=("Protein_Name", semicolon_join),
            organism=("Organism", semicolon_join),
        )
        .reset_index()
    )

    oligo = oligo_activity.rename(
        columns={
            "UniProt_ID": "protein_id",
            "Family": "oligo_family",
            "EC#": "oligo_ec_numbers",
            "Ligandkategori_EC": "oligo_ligand_category_ec",
            "Oligo_active": "oligo_activity_value_updated",
            "Oligo_active_classes": "oligo_active_classes",
        }
    ).copy()
    oligo["protein_id"] = oligo["protein_id"].astype(str)
    oligo["oligo_activity_true_row"] = to_ja_bool_series(oligo["oligo_activity_value_updated"])
    oligo_by_protein = (
        oligo.groupby("protein_id", dropna=False)
        .agg(
            oligo_family=("oligo_family", semicolon_join),
            oligo_ec_numbers=("oligo_ec_numbers", semicolon_join),
            oligo_ligand_category_ec=("oligo_ligand_category_ec", semicolon_join),
            oligo_activity_values_updated=("oligo_activity_value_updated", semicolon_join),
            oligo_active_classes=("oligo_active_classes", semicolon_join),
            oligo_activity_true=("oligo_activity_true_row", "max"),
            n_oligo_activity_rows=("protein_id", "size"),
        )
        .reset_index()
    )

    labels = (
        run_proteins.merge(pdb_by_protein, on="protein_id", how="left")
        .merge(meta_by_protein, on="protein_id", how="left")
        .merge(oligo_by_protein, on="protein_id", how="left")
    )
    labels["has_pdb_oligo_record"] = labels["n_pdb_rows"].notna()
    labels["has_updated_oligo_activity_record"] = labels["n_oligo_activity_rows"].notna()
    labels["oligo_activity_true"] = labels["oligo_activity_true"].fillna(False).astype(bool)
    labels["oligo_activity_group"] = np.select(
        [labels["oligo_activity_true"], labels["has_updated_oligo_activity_record"]],
        ["Oligo activity true", "Oligo activity false"],
        default="No updated oligo activity record",
    )
    labels["family"] = (
        labels["oligo_family"].replace("", np.nan).fillna(labels["pdb_family"].replace("", np.nan)).fillna(labels["metadata_family"])
    )
    return labels.sort_values(["oligo_activity_group", "protein_id"]).reset_index(drop=True)


def active_substrates_from_text(text: object) -> set[str]:
    raw = str(text).lower()
    out: set[str] = set()
    if "chitin" in raw or "nag" in raw:
        out.add("chitin")
    if "cellulose" in raw or "cell" in raw or "xylo" in raw:
        out.add("cellulose")
    if "starch" in raw or "amylose" in raw or "malto" in raw:
        out.add("starch")
    return out


def ec_active_substrates(text: object) -> set[str]:
    raw = str(text)
    out: set[str] = set()
    if "1.14.99.53" in raw:
        out.add("chitin")
    if "1.14.99.54" in raw or "1.14.99.56" in raw:
        out.add("cellulose")
    if "1.14.99.55" in raw:
        out.add("starch")
    return out


def make_activity_table(oligo_activity: pd.DataFrame) -> pd.DataFrame:
    rows = []
    oligo = oligo_activity.rename(
        columns={
            "UniProt_ID": "protein_id",
            "Family": "activity_family",
            "EC#": "ec_numbers",
            "Ligandkategori_EC": "ligand_category_ec",
            "Oligo_active": "oligo_activity_value",
            "Oligo_active_classes": "oligo_active_classes",
        }
    ).copy()
    for _, row in oligo.iterrows():
        active = active_substrates_from_text(row.get("ligand_category_ec")) | ec_active_substrates(row.get("ec_numbers"))
        if str(row.get("oligo_activity_value")).strip().lower() == "ja":
            active |= active_substrates_from_text(row.get("oligo_active_classes"))
        rows.append(
            {
                "protein_id": row.get("protein_id"),
                "activity_family": row.get("activity_family"),
                "active_substrates": ",".join(sorted(active)),
                "activity_source": "oligo_active_updated.tsv",
            }
        )
    out = pd.DataFrame(rows).dropna(subset=["protein_id"])
    out["protein_id"] = out["protein_id"].astype(str)
    return out.drop_duplicates("protein_id", keep="first")


def normalize_substrate_class(value: object) -> str:
    substrate = str(value)
    return "starch" if substrate == "amylose" else substrate


def add_condition_metrics(condition: pd.DataFrame, labels: pd.DataFrame, activity: pd.DataFrame) -> pd.DataFrame:
    df = condition.copy()
    for col in NUMERIC_INPUT_COLUMNS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in ["formal_clustering_allowed", "any_valid_cluster"]:
        if col in df.columns:
            df[col] = to_bool_series(df[col]).astype(int)

    df["hard_qc_pass_fraction_condition"] = df.apply(lambda row: safe_ratio(row["n_stage1_pass"], row["n_generated"]), axis=1)
    df["hard_qc_fail_fraction_condition"] = df.apply(lambda row: safe_ratio(row["n_stage1_hard_fail"], row["n_generated"]), axis=1)
    df["ifp_success_fraction_condition"] = df.apply(lambda row: safe_ratio(row["n_ifp_success"], row["n_stage1_pass"]), axis=1)
    df["clusterable_pose_fraction_condition"] = df.apply(lambda row: safe_ratio(row["n_ifp_clustered"], row["n_ifp_success"]), axis=1)
    df["best_plausible_fraction_condition"] = df[
        ["occupancy_weighted_c1_plausible_fraction", "occupancy_weighted_c4_plausible_fraction"]
    ].max(axis=1, skipna=True)
    df["best_cu_distance_condition"] = df[
        ["occupancy_weighted_Cu_C1_distance_median", "occupancy_weighted_Cu_C4_distance_median"]
    ].min(axis=1, skipna=True)
    df["best_oxyl_h_distance_condition"] = df[
        ["occupancy_weighted_oxyl_H_C1_distance_median", "occupancy_weighted_oxyl_H_C4_distance_median"]
    ].min(axis=1, skipna=True)
    df["best_oxyl_h_score_condition"] = df[
        ["occupancy_weighted_oxyl_H_score_C1_median", "occupancy_weighted_oxyl_H_score_C4_median"]
    ].max(axis=1, skipna=True)

    keep = [
        "protein_id",
        "family",
        "oligo_activity_group",
        "oligo_activity_true",
        "has_updated_oligo_activity_record",
        "oligo_activity_values_updated",
        "oligo_active_classes",
        "oligo_ligand_category_ec",
        "has_pdb_oligo_record",
        "oligo_activity_values",
        "n_pdb_rows",
        "best_resolution",
        "pdb_ids",
        "carbohydrate_ligands",
    ]
    out = df.merge(labels[keep], on="protein_id", how="left").merge(
        activity[["protein_id", "active_substrates", "activity_source"]],
        on="protein_id",
        how="left",
    )
    active_sets = out["active_substrates"].fillna("").map(lambda value: set(filter(None, str(value).split(","))))
    substrates = out["substrate_class"].map(normalize_substrate_class)
    out["demonstrated_activity_condition"] = [
        substrate in active_set if active_set else np.nan
        for substrate, active_set in zip(substrates, active_sets, strict=False)
    ]
    out["activity_condition_group"] = [
        "Demonstrated activity"
        if value is True
        else ("No demonstrated activity" if value is False else "No activity metadata")
        for value in out["demonstrated_activity_condition"]
    ]
    return out


def aggregate_condition_block(sub: pd.DataFrame, construct_group: str) -> dict[str, object]:
    row: dict[str, object] = {
        "protein_id": sub["protein_id"].iloc[0],
        "construct_group": construct_group,
        "oligo_activity_group": sub["oligo_activity_group"].iloc[0],
        "oligo_activity_true": bool(sub["oligo_activity_true"].iloc[0]),
        "has_pdb_oligo_record": bool(sub["has_pdb_oligo_record"].iloc[0]),
        "family": semicolon_join(sub["family"]),
        "pdb_ids": semicolon_join(sub["pdb_ids"]),
        "oligo_activity_values": semicolon_join(sub["oligo_activity_values"]),
        "n_conditions": int(len(sub)),
        "n_substrate_classes": int(sub["substrate_class"].nunique(dropna=True)),
        "substrate_classes": semicolon_join(sub["substrate_class"]),
        "dp_values": semicolon_join(sub["dp"]),
    }

    sums = {
        "n_generated_total": "n_generated",
        "n_stage1_pass_total": "n_stage1_pass",
        "n_stage1_hard_fail_total": "n_stage1_hard_fail",
        "n_ifp_success_total": "n_ifp_success",
        "n_contact_eligible_total": "n_contact_eligible",
        "n_ifp_clustered_total": "n_ifp_clustered",
        "n_noise_total": "n_noise",
    }
    for out_col, in_col in sums.items():
        row[out_col] = float(pd.to_numeric(sub[in_col], errors="coerce").sum())

    row["hard_qc_pass_fraction"] = safe_ratio(row["n_stage1_pass_total"], row["n_generated_total"])
    row["hard_qc_fail_fraction"] = safe_ratio(row["n_stage1_hard_fail_total"], row["n_generated_total"])
    row["ifp_success_fraction"] = safe_ratio(row["n_ifp_success_total"], row["n_stage1_pass_total"])
    row["contact_eligible_fraction_weighted"] = safe_ratio(row["n_contact_eligible_total"], row["n_ifp_success_total"])
    row["clusterable_pose_fraction"] = safe_ratio(row["n_ifp_clustered_total"], row["n_ifp_success_total"])
    row["noise_fraction_weighted"] = safe_ratio(row["n_noise_total"], row["n_ifp_clustered_total"])
    row["formal_clusterable_condition_fraction"] = float(sub["formal_clustering_allowed"].mean())
    row["valid_cluster_condition_fraction"] = float(sub["any_valid_cluster"].mean())
    row["n_clusters_mean"] = float(pd.to_numeric(sub["n_clusters"], errors="coerce").mean())
    row["n_clusters_median"] = float(pd.to_numeric(sub["n_clusters"], errors="coerce").median())

    for in_col, out_col in SHORT_MEDIAN_COLUMNS.items():
        row[out_col] = float(pd.to_numeric(sub[in_col], errors="coerce").median())
    return row


def aggregate_protein_construct(condition_labeled: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (protein_id, construct_type), sub in condition_labeled.groupby(["protein_id", "construct_type"], observed=True):
        rows.append(aggregate_condition_block(sub, CONSTRUCT_LABELS.get(str(construct_type), str(construct_type))))
    for protein_id, sub in condition_labeled.groupby("protein_id", observed=True):
        rows.append(aggregate_condition_block(sub, "All constructs"))
    return pd.DataFrame(rows).sort_values(["construct_group", "oligo_activity_group", "protein_id"]).reset_index(drop=True)


def make_group_counts(labels: pd.DataFrame, condition_labeled: pd.DataFrame, protein_construct: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for construct_group, sub in protein_construct.groupby("construct_group", observed=True):
        for group, group_sub in sub.groupby("oligo_activity_group", observed=True):
            rows.append(
                {
                    "construct_group": construct_group,
                    "oligo_activity_group": group,
                    "n_proteins": group_sub["protein_id"].nunique(),
                    "n_protein_construct_rows": len(group_sub),
                    "n_condition_rows": int(
                        condition_labeled[
                            condition_labeled["protein_id"].isin(group_sub["protein_id"])
                            & (
                                condition_labeled["construct_type"].map(CONSTRUCT_LABELS).fillna(condition_labeled["construct_type"]) == construct_group
                                if construct_group != "All constructs"
                                else True
                            )
                        ].shape[0]
                    ),
                }
            )
    label_counts = labels.groupby("oligo_activity_group", observed=True)["protein_id"].nunique().reset_index(name="n_run_proteins_with_label")
    counts = pd.DataFrame(rows).merge(label_counts, on="oligo_activity_group", how="left")
    return counts.sort_values(["construct_group", "oligo_activity_group"]).reset_index(drop=True)


def make_group_descriptive(protein_construct: pd.DataFrame) -> pd.DataFrame:
    rows = []
    metrics = list(METRIC_INFO)
    for construct_group, construct_sub in protein_construct.groupby("construct_group", observed=True):
        for group, group_sub in construct_sub.groupby("oligo_activity_group", observed=True):
            for metric in metrics:
                metric_label, unit, direction = METRIC_INFO[metric]
                desc = describe(group_sub[metric])
                rows.append(
                    {
                        "construct_group": construct_group,
                        "oligo_activity_group": group,
                        "metric": metric,
                        "metric_label": metric_label,
                        "unit": unit,
                        "direction": direction,
                        **desc,
                    }
                )
    return pd.DataFrame(rows)


def make_comparisons(protein_construct: pd.DataFrame) -> pd.DataFrame:
    rows = []
    comparator_specs = [
        ("All non-true run proteins", lambda df: ~df["oligo_activity_true"]),
        ("Updated oligo-active Nei proteins", lambda df: df["oligo_activity_group"].eq("Oligo activity false")),
        (
            "No updated oligo activity record",
            lambda df: df["oligo_activity_group"].eq("No updated oligo activity record"),
        ),
    ]
    for construct_group, construct_sub in protein_construct.groupby("construct_group", observed=True):
        true_sub = construct_sub[construct_sub["oligo_activity_true"]].copy()
        for comparator_name, comparator_mask in comparator_specs:
            comp_sub = construct_sub[comparator_mask(construct_sub)].copy()
            for metric, (metric_label, unit, direction) in METRIC_INFO.items():
                true_desc = describe(true_sub[metric])
                comp_desc = describe(comp_sub[metric])
                rows.append(
                    {
                        "construct_group": construct_group,
                        "comparison": f"Oligo activity true vs {comparator_name}",
                        "metric": metric,
                        "metric_label": metric_label,
                        "unit": unit,
                        "direction": direction,
                        "n_true": true_desc["n"],
                        "true_mean": true_desc["mean"],
                        "true_sd": true_desc["sd"],
                        "true_median": true_desc["median"],
                        "true_q1": true_desc["q1"],
                        "true_q3": true_desc["q3"],
                        "n_comparator": comp_desc["n"],
                        "comparator_mean": comp_desc["mean"],
                        "comparator_sd": comp_desc["sd"],
                        "comparator_median": comp_desc["median"],
                        "comparator_q1": comp_desc["q1"],
                        "comparator_q3": comp_desc["q3"],
                        "mean_difference_true_minus_comparator": (
                            true_desc["mean"] - comp_desc["mean"] if true_desc["n"] and comp_desc["n"] else np.nan
                        ),
                        "median_difference_true_minus_comparator": (
                            true_desc["median"] - comp_desc["median"] if true_desc["n"] and comp_desc["n"] else np.nan
                        ),
                        "cliffs_delta_true_vs_comparator": cliffs_delta(true_sub[metric], comp_sub[metric]),
                    }
                )
    return pd.DataFrame(rows)


def make_presentable_table(comparisons: pd.DataFrame) -> pd.DataFrame:
    sub = comparisons[
        (comparisons["construct_group"].eq("All constructs"))
        & (comparisons["comparison"].eq("Oligo activity true vs All non-true run proteins"))
        & (comparisons["metric"].isin(PRESENTABLE_METRICS))
    ].copy()
    rows = []
    for _, row in sub.iterrows():
        delta = row["median_difference_true_minus_comparator"]
        direction = row["direction"]
        if pd.isna(delta):
            reading = ""
        elif "higher is better" in direction:
            reading = "Better for oligo-active" if delta > 0 else "Lower for oligo-active"
        elif "lower is better" in direction:
            reading = "Better for oligo-active" if delta < 0 else "Higher/worse for oligo-active"
        else:
            reading = "Higher for oligo-active" if delta > 0 else "Lower for oligo-active"
        rows.append(
            {
                "metric": row["metric"],
                "metric_label": row["metric_label"],
                "unit": row["unit"],
                "direction": direction,
                "n_oligo_active": row["n_true"],
                "median_oligo_active": row["true_median"],
                "n_not_oligo_active": row["n_comparator"],
                "median_not_oligo_active": row["comparator_median"],
                "median_delta_oligo_minus_not": delta,
                "cliffs_delta": row["cliffs_delta_true_vs_comparator"],
                "short_reading": reading,
            }
        )
    out = pd.DataFrame(rows)
    order = {metric: idx for idx, metric in enumerate(PRESENTABLE_METRICS)}
    out["_order"] = out["metric"].map(order)
    return out.sort_values("_order").drop(columns="_order").reset_index(drop=True)


def write_markdown_table(table: pd.DataFrame, path: Path) -> None:
    display = table.copy()
    numeric = [
        "median_oligo_active",
        "median_not_oligo_active",
        "median_delta_oligo_minus_not",
        "cliffs_delta",
    ]
    for col in numeric:
        display[col] = pd.to_numeric(display[col], errors="coerce").map(lambda value: "" if pd.isna(value) else f"{value:.3f}")
    display = display[
        [
            "metric_label",
            "direction",
            "n_oligo_active",
            "median_oligo_active",
            "n_not_oligo_active",
            "median_not_oligo_active",
            "median_delta_oligo_minus_not",
            "cliffs_delta",
            "short_reading",
        ]
    ]
    headers = list(display.columns)
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for _, row in display.iterrows():
        lines.append("| " + " | ".join(str(row[col]) for col in headers) + " |")
    path.write_text("\n".join(lines) + "\n")


def write_generic_markdown_table(table: pd.DataFrame, path: Path) -> None:
    display = table.copy()
    for col in display.columns:
        if pd.api.types.is_numeric_dtype(display[col]):
            display[col] = display[col].map(lambda value: "" if pd.isna(value) else (str(int(value)) if float(value).is_integer() else f"{value:.3f}"))
    headers = list(display.columns)
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for _, row in display.iterrows():
        lines.append("| " + " | ".join(str(row[col]) for col in headers) + " |")
    path.write_text("\n".join(lines) + "\n")


def make_presentable_group_medians(protein_construct: pd.DataFrame) -> pd.DataFrame:
    data = protein_construct[protein_construct["construct_group"].eq("All constructs")].copy()
    rows = []
    for metric in PRESENTABLE_METRICS:
        metric_label, unit, direction = METRIC_INFO[metric]
        row = {"metric": metric, "metric_label": metric_label, "unit": unit, "direction": direction}
        for group, label in PRESENTABLE_LABELS.items():
            values = pd.to_numeric(data.loc[data["oligo_activity_group"].eq(group), metric], errors="coerce").dropna()
            key = label.lower().replace(" ", "_")
            row[f"n_{key}"] = int(values.shape[0])
            row[f"median_{key}"] = float(values.median()) if len(values) else np.nan
            row[f"q1_{key}"] = float(values.quantile(0.25)) if len(values) else np.nan
            row[f"q3_{key}"] = float(values.quantile(0.75)) if len(values) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def plot_demonstrated_activity_comparison(protein_construct: pd.DataFrame, output_dir: Path) -> None:
    data = protein_construct[protein_construct["construct_group"].eq("All constructs")].copy()
    data = data[data["oligo_activity_group"].isin(PRESENTABLE_LABELS)]
    long = data.melt(
        id_vars=["protein_id", "oligo_activity_group"],
        value_vars=PRESENTABLE_METRICS,
        var_name="metric",
        value_name="value",
    ).dropna(subset=["value"])
    long["metric_label"] = long["metric"].map(lambda metric: METRIC_INFO[metric][0])
    long["group_label"] = long["oligo_activity_group"].map(PRESENTABLE_LABELS)

    order = ["Oligo active", "Not oligo active", "No updated record"]
    palette = {"Oligo active": "#1f77b4", "Not oligo active": "#9b3a3a", "No updated record": "#9a7b35"}
    fig, axes = plt.subplots(2, 4, figsize=(13.5, 6.4), constrained_layout=True)
    for ax, metric in zip(axes.flat, PRESENTABLE_METRICS, strict=False):
        sub = long[long["metric"].eq(metric)]
        sns.boxplot(
            data=sub,
            x="group_label",
            y="value",
            hue="group_label",
            order=order,
            hue_order=order,
            palette=palette,
            showfliers=False,
            width=0.5,
            ax=ax,
            legend=False,
        )
        sns.stripplot(
            data=sub,
            x="group_label",
            y="value",
            order=order,
            color="#1b1b1b",
            alpha=0.35,
            size=3,
            jitter=0.18,
            ax=ax,
        )
        ax.set_title(METRIC_INFO[metric][0], fontsize=9, color="black")
        ax.set_xlabel("")
        ax.set_ylabel(METRIC_INFO[metric][1].capitalize(), color="black")
        ax.grid(False)
        ax.tick_params(axis="both", colors="black", length=3, width=0.8)
        ax.tick_params(axis="x", rotation=25)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_color("black")
            label.set_fontsize(8)
        for spine in ["left", "bottom"]:
            ax.spines[spine].set_color("black")
            ax.spines[spine].set_linewidth(0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    for path in [
        output_dir / "oligo_activity_demonstrated_condition_comparison.png",
        output_dir / "oligo_activity_demonstrated_condition_comparison.pdf",
    ]:
        fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_demonstrated_activity_outputs(condition_labeled: pd.DataFrame, output_dir: Path) -> None:
    demonstrated = condition_labeled[condition_labeled["activity_condition_group"].eq("Demonstrated activity")].copy()
    demonstrated_metrics = aggregate_protein_construct(demonstrated)
    demonstrated_descriptive = make_group_descriptive(demonstrated_metrics)
    demonstrated_comparisons = make_comparisons(demonstrated_metrics)
    presentable = make_presentable_table(demonstrated_comparisons)
    group_medians = make_presentable_group_medians(demonstrated_metrics)

    demonstrated.to_csv(output_dir / "oligo_activity_demonstrated_conditions.tsv", sep="\t", index=False)
    demonstrated_metrics.to_csv(output_dir / "oligo_activity_demonstrated_condition_protein_construct_metrics.tsv", sep="\t", index=False)
    demonstrated_descriptive.to_csv(output_dir / "oligo_activity_demonstrated_condition_group_descriptive_statistics.tsv", sep="\t", index=False)
    demonstrated_comparisons.to_csv(output_dir / "oligo_activity_demonstrated_condition_comparative_statistics.tsv", sep="\t", index=False)
    presentable.to_csv(output_dir / "oligo_activity_demonstrated_condition_presentable_comparison.tsv", sep="\t", index=False)
    group_medians.to_csv(output_dir / "oligo_activity_demonstrated_condition_presentable_group_medians.tsv", sep="\t", index=False)
    write_markdown_table(presentable, output_dir / "oligo_activity_demonstrated_condition_presentable_comparison.md")
    write_generic_markdown_table(group_medians, output_dir / "oligo_activity_demonstrated_condition_presentable_group_medians.md")
    plot_demonstrated_activity_comparison(demonstrated_metrics, output_dir)


def aggregate_activity_state_block(sub: pd.DataFrame) -> dict[str, object]:
    row: dict[str, object] = {
        "protein_id": sub["protein_id"].iloc[0],
        "oligo_activity_true": bool(sub["oligo_activity_true"].iloc[0]),
        "activity_condition_group": sub["activity_condition_group"].iloc[0],
        "n_conditions": int(len(sub)),
        "n_construct_types": int(sub["construct_type"].nunique(dropna=True)),
        "construct_types": semicolon_join(sub["construct_type"]),
        "substrate_classes": semicolon_join(sub["substrate_class"]),
        "dp_values": semicolon_join(sub["dp"]),
    }
    totals = {
        "n_generated_total": "n_generated",
        "n_stage1_pass_total": "n_stage1_pass",
        "n_stage1_hard_fail_total": "n_stage1_hard_fail",
        "n_ifp_success_total": "n_ifp_success",
        "n_contact_eligible_total": "n_contact_eligible",
        "n_ifp_clustered_total": "n_ifp_clustered",
        "n_noise_total": "n_noise",
    }
    for out_col, in_col in totals.items():
        row[out_col] = float(pd.to_numeric(sub[in_col], errors="coerce").sum())

    row["hard_qc_pass_fraction"] = safe_ratio(row["n_stage1_pass_total"], row["n_generated_total"])
    row["hard_qc_fail_fraction"] = safe_ratio(row["n_stage1_hard_fail_total"], row["n_generated_total"])
    row["ifp_success_fraction"] = safe_ratio(row["n_ifp_success_total"], row["n_stage1_pass_total"])
    row["contact_eligible_fraction_weighted"] = safe_ratio(row["n_contact_eligible_total"], row["n_ifp_success_total"])
    row["clusterable_pose_fraction"] = safe_ratio(row["n_ifp_clustered_total"], row["n_ifp_success_total"])
    row["noise_fraction_weighted"] = safe_ratio(row["n_noise_total"], row["n_ifp_clustered_total"])
    row["formal_clusterable_condition_fraction"] = float(pd.to_numeric(sub["formal_clustering_allowed"], errors="coerce").mean())
    row["valid_cluster_condition_fraction"] = float(pd.to_numeric(sub["any_valid_cluster"], errors="coerce").mean())
    row["n_clusters_mean"] = float(pd.to_numeric(sub["n_clusters"], errors="coerce").mean())
    row["n_clusters_median"] = float(pd.to_numeric(sub["n_clusters"], errors="coerce").median())
    for in_col, out_col in SHORT_MEDIAN_COLUMNS.items():
        row[out_col] = float(pd.to_numeric(sub[in_col], errors="coerce").median())
    return row


def make_oligo_2x2_metrics(condition_labeled: pd.DataFrame) -> pd.DataFrame:
    eligible = condition_labeled[
        condition_labeled["activity_condition_group"].isin(["Demonstrated activity", "No demonstrated activity"])
    ].copy()
    rows = []
    for _, sub in eligible.groupby(["protein_id", "oligo_activity_true", "activity_condition_group"], observed=True):
        row = aggregate_activity_state_block(sub)
        if row["oligo_activity_true"] and row["activity_condition_group"] == "Demonstrated activity":
            group = "Oligo active demonstrated activity"
        elif row["oligo_activity_true"] and row["activity_condition_group"] == "No demonstrated activity":
            group = "Oligo active not demonstrated activity"
        elif (not row["oligo_activity_true"]) and row["activity_condition_group"] == "Demonstrated activity":
            group = "All other demonstrated activity"
        else:
            group = "All other no demonstrated activity"
        row["oligo_2x2_group"] = group
        row["oligo_factor"] = int(row["oligo_activity_true"])
        row["demonstrated_factor"] = int(row["activity_condition_group"] == "Demonstrated activity")
        rows.append(row)
    out = pd.DataFrame(rows)
    out["oligo_2x2_group"] = pd.Categorical(out["oligo_2x2_group"], categories=OLIGO_2X2_ORDER, ordered=True)
    return out.sort_values(["oligo_2x2_group", "protein_id"]).reset_index(drop=True)


def ols_sse(y: np.ndarray, x: np.ndarray) -> tuple[float, int, int]:
    mask = np.isfinite(y) & np.isfinite(x).all(axis=1)
    y = y[mask]
    x = x[mask]
    if len(y) == 0:
        return np.nan, 0, 0
    beta, _, rank, _ = np.linalg.lstsq(x, y, rcond=None)
    residuals = y - x @ beta
    return float(np.sum(residuals**2)), int(len(y)), int(rank)


def partial_f_pvalue(sse_full: float, rank_full: int, sse_reduced: float, rank_reduced: int, n: int) -> tuple[float, float, int, int]:
    df_num = rank_full - rank_reduced
    df_den = n - rank_full
    if (
        stats is None
        or not np.isfinite(sse_full)
        or not np.isfinite(sse_reduced)
        or df_num <= 0
        or df_den <= 0
        or sse_full < 0
    ):
        return np.nan, np.nan, df_num, df_den
    ms_num = max((sse_reduced - sse_full) / df_num, 0.0)
    ms_den = sse_full / df_den
    if ms_den <= 0:
        return np.nan, np.nan, df_num, df_den
    f_value = ms_num / ms_den
    return float(f_value), float(stats.f.sf(f_value, df_num, df_den)), df_num, df_den


def format_p(value: float) -> str:
    if pd.isna(value):
        return "NA"
    if value < 0.001:
        return "<0.001"
    return f"{value:.3f}"


def two_way_anova(metric_df: pd.DataFrame, metric: str) -> list[dict[str, object]]:
    sub = metric_df.dropna(subset=[metric, "oligo_factor", "demonstrated_factor"]).copy()
    y = pd.to_numeric(sub[metric], errors="coerce").to_numpy(dtype=float)
    oligo = pd.to_numeric(sub["oligo_factor"], errors="coerce").to_numpy(dtype=float)
    demo = pd.to_numeric(sub["demonstrated_factor"], errors="coerce").to_numpy(dtype=float)
    interaction = oligo * demo
    full_x = np.column_stack([np.ones(len(sub)), oligo, demo, interaction])
    sse_full, n, rank_full = ols_sse(y, full_x)
    reduced = {
        "Oligo active": np.column_stack([np.ones(len(sub)), demo, interaction]),
        "Demonstrated activity": np.column_stack([np.ones(len(sub)), oligo, interaction]),
        "Interaction": np.column_stack([np.ones(len(sub)), oligo, demo]),
    }
    rows = []
    for term, x_reduced in reduced.items():
        sse_reduced, _, rank_reduced = ols_sse(y, x_reduced)
        f_value, p_value, df_num, df_den = partial_f_pvalue(sse_full, rank_full, sse_reduced, rank_reduced, n)
        rows.append(
            {
                "metric": metric,
                "metric_label": METRIC_INFO[metric][0],
                "term": term,
                "n": n,
                "df_num": df_num,
                "df_den": df_den,
                "f_value": f_value,
                "p_value": p_value,
                "model": "OLS two-way ANOVA: value ~ Oligo_active + Demonstrated_activity + interaction",
            }
        )
    return rows


def pairwise_oligo_vs_rest(metric_df: pd.DataFrame, metric: str) -> list[dict[str, object]]:
    rows = []
    for activity_group in ["Demonstrated activity", "No demonstrated activity"]:
        sub = metric_df[metric_df["activity_condition_group"].eq(activity_group)].copy()
        oligo_values = pd.to_numeric(sub.loc[sub["oligo_activity_true"], metric], errors="coerce").dropna()
        rest_values = pd.to_numeric(sub.loc[~sub["oligo_activity_true"], metric], errors="coerce").dropna()
        p_value = np.nan
        statistic = np.nan
        if stats is not None and len(oligo_values) >= 2 and len(rest_values) >= 2:
            test = stats.ttest_ind(oligo_values, rest_values, equal_var=False, nan_policy="omit")
            statistic = float(test.statistic)
            p_value = float(test.pvalue)
        rows.append(
            {
                "metric": metric,
                "metric_label": METRIC_INFO[metric][0],
                "activity_condition_group": activity_group,
                "n_oligo_active": len(oligo_values),
                "median_oligo_active": float(oligo_values.median()) if len(oligo_values) else np.nan,
                "n_all_other": len(rest_values),
                "median_all_other": float(rest_values.median()) if len(rest_values) else np.nan,
                "median_delta_oligo_minus_all_other": (
                    float(oligo_values.median() - rest_values.median()) if len(oligo_values) and len(rest_values) else np.nan
                ),
                "welch_t": statistic,
                "welch_p_value": p_value,
                "cliffs_delta": cliffs_delta(oligo_values, rest_values),
            }
        )
    return rows


def make_oligo_2x2_summary(metric_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric in PRESENTABLE_METRICS:
        metric_label, unit, direction = METRIC_INFO[metric]
        for group in OLIGO_2X2_ORDER:
            values = pd.to_numeric(metric_df.loc[metric_df["oligo_2x2_group"].eq(group), metric], errors="coerce").dropna()
            desc = describe(values)
            rows.append(
                {
                    "metric": metric,
                    "metric_label": metric_label,
                    "unit": unit,
                    "direction": direction,
                    "oligo_2x2_group": group,
                    **desc,
                }
            )
    return pd.DataFrame(rows)


def plot_oligo_2x2(metric_df: pd.DataFrame, anova: pd.DataFrame, output_dir: Path) -> None:
    long = metric_df.melt(
        id_vars=["protein_id", "oligo_2x2_group", "oligo_activity_true", "activity_condition_group"],
        value_vars=PRESENTABLE_METRICS,
        var_name="metric",
        value_name="value",
    ).dropna(subset=["value"])
    long["oligo_2x2_group"] = pd.Categorical(long["oligo_2x2_group"], categories=OLIGO_2X2_ORDER, ordered=True)
    fig, axes = plt.subplots(2, 4, figsize=(15.2, 7.0), constrained_layout=True)
    for ax, metric in zip(axes.flat, PRESENTABLE_METRICS, strict=False):
        sub = long[long["metric"].eq(metric)].copy()
        counts = sub.groupby("oligo_2x2_group", observed=False)["protein_id"].nunique().to_dict()
        labels = [f"{OLIGO_2X2_SHORT_LABELS[group]}\nn={int(counts.get(group, 0))}" for group in OLIGO_2X2_ORDER]
        sns.boxplot(
            data=sub,
            x="oligo_2x2_group",
            y="value",
            hue="oligo_2x2_group",
            order=OLIGO_2X2_ORDER,
            hue_order=OLIGO_2X2_ORDER,
            palette=OLIGO_2X2_PALETTE,
            showfliers=False,
            width=0.48,
            ax=ax,
            legend=False,
        )
        sns.stripplot(
            data=sub,
            x="oligo_2x2_group",
            y="value",
            order=OLIGO_2X2_ORDER,
            color="#1f1f1f",
            alpha=0.42,
            size=3,
            jitter=0.16,
            ax=ax,
        )
        ax.set_xticks(range(len(OLIGO_2X2_ORDER)))
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_title(METRIC_INFO[metric][0], fontsize=9, color="black")
        ax.set_xlabel("")
        ax.set_ylabel(METRIC_INFO[metric][1].capitalize(), color="black")
        metric_anova = anova[anova["metric"].eq(metric)]
        p_oligo = metric_anova.loc[metric_anova["term"].eq("Oligo active"), "p_value"]
        p_demo = metric_anova.loc[metric_anova["term"].eq("Demonstrated activity"), "p_value"]
        p_int = metric_anova.loc[metric_anova["term"].eq("Interaction"), "p_value"]
        text = (
            f"ANOVA p\n"
            f"Oligo={format_p(p_oligo.iloc[0] if len(p_oligo) else np.nan)}\n"
            f"Activity={format_p(p_demo.iloc[0] if len(p_demo) else np.nan)}\n"
            f"Int={format_p(p_int.iloc[0] if len(p_int) else np.nan)}"
        )
        ax.text(
            0.98,
            0.98,
            text,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=7,
            color="black",
            bbox={"facecolor": "white", "edgecolor": "#777777", "linewidth": 0.5, "alpha": 0.88},
        )
        ax.grid(False)
        ax.tick_params(axis="both", colors="black", length=3, width=0.8)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_color("black")
        for spine in ["left", "bottom"]:
            ax.spines[spine].set_color("black")
            ax.spines[spine].set_linewidth(0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    for path in [
        output_dir / "oligo_active_vs_all_other_activity_2x2_boxplots.png",
        output_dir / "oligo_active_vs_all_other_activity_2x2_boxplots.pdf",
    ]:
        fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_oligo_2x2_outputs(condition_labeled: pd.DataFrame, output_dir: Path) -> None:
    metric_df = make_oligo_2x2_metrics(condition_labeled)
    anova_rows: list[dict[str, object]] = []
    pairwise_rows: list[dict[str, object]] = []
    for metric in PRESENTABLE_METRICS:
        anova_rows.extend(two_way_anova(metric_df, metric))
        pairwise_rows.extend(pairwise_oligo_vs_rest(metric_df, metric))
    anova = pd.DataFrame(anova_rows)
    pairwise = pd.DataFrame(pairwise_rows)
    summary = make_oligo_2x2_summary(metric_df)

    metric_df.to_csv(output_dir / "oligo_active_vs_all_other_activity_2x2_protein_metrics.tsv", sep="\t", index=False)
    summary.to_csv(output_dir / "oligo_active_vs_all_other_activity_2x2_group_summary.tsv", sep="\t", index=False)
    anova.to_csv(output_dir / "oligo_active_vs_all_other_activity_2x2_anova.tsv", sep="\t", index=False)
    pairwise.to_csv(output_dir / "oligo_active_vs_all_other_activity_2x2_pairwise.tsv", sep="\t", index=False)
    write_generic_markdown_table(summary, output_dir / "oligo_active_vs_all_other_activity_2x2_group_summary.md")
    write_generic_markdown_table(anova, output_dir / "oligo_active_vs_all_other_activity_2x2_anova.md")
    plot_oligo_2x2(metric_df, anova, output_dir)


def direction_adjusted_values(metric: str, values: pd.Series) -> pd.Series:
    direction = METRIC_INFO[metric][2]
    numeric = pd.to_numeric(values, errors="coerce")
    if "lower is better" in direction:
        return -numeric
    return numeric


def direction_adjusted_effect(metric: str, active_median: float, other_median: float) -> float:
    if pd.isna(active_median) or pd.isna(other_median):
        return np.nan
    direction = METRIC_INFO[metric][2]
    if "lower is better" in direction:
        return float(other_median - active_median)
    return float(active_median - other_median)


def permutation_p_value(
    values: np.ndarray,
    labels: np.ndarray,
    observed_effect: float,
    rng: np.random.Generator,
    n_permutations: int = N_PERMUTATIONS,
) -> float:
    if (
        values.size == 0
        or labels.size == 0
        or labels.sum() == 0
        or labels.sum() == labels.size
        or not np.isfinite(observed_effect)
    ):
        return np.nan
    count = 0
    active_n = int(labels.sum())
    for _ in range(n_permutations):
        permuted = rng.permutation(labels)
        active = values[permuted.astype(bool)]
        other = values[~permuted.astype(bool)]
        if active.size != active_n or other.size == 0:
            continue
        perm_effect = float(np.median(active) - np.median(other))
        if abs(perm_effect) >= abs(observed_effect):
            count += 1
    return float((count + 1) / (n_permutations + 1))


def pr_auc_score(labels: np.ndarray, scores: np.ndarray) -> float:
    if average_precision_score is None or labels.sum() == 0 or labels.sum() == labels.size:
        return np.nan
    try:
        return float(average_precision_score(labels, scores))
    except ValueError:
        return np.nan


def roc_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    if roc_auc_score is None or labels.sum() == 0 or labels.sum() == labels.size:
        return np.nan
    try:
        return float(roc_auc_score(labels, scores))
    except ValueError:
        return np.nan


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    p = pd.to_numeric(p_values, errors="coerce")
    q = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.dropna().sort_values()
    m = len(valid)
    if m == 0:
        return q
    adjusted = valid.to_numpy(dtype=float) * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    q.loc[valid.index] = adjusted
    return q


def make_permutation_table(
    protein_construct: pd.DataFrame,
    scope_label: str,
    aa17_policy: str,
    seed_offset: int,
) -> pd.DataFrame:
    data = protein_construct[protein_construct["construct_group"].eq("All constructs")].copy()
    if aa17_policy == "aa17_excluded":
        data = data[~data["family"].fillna("").astype(str).str.contains(r"\bAA17\b")].copy()

    rng = np.random.default_rng(PERMUTATION_SEED + seed_offset)
    rows = []
    for metric in METRIC_INFO:
        metric_label, unit, direction = METRIC_INFO[metric]
        sub = data[["protein_id", "oligo_activity_true", metric]].copy()
        sub[metric] = pd.to_numeric(sub[metric], errors="coerce")
        sub = sub.dropna(subset=[metric, "oligo_activity_true"])
        active_values = sub.loc[sub["oligo_activity_true"], metric]
        other_values = sub.loc[~sub["oligo_activity_true"], metric]
        active_desc = describe(active_values)
        other_desc = describe(other_values)
        effect = direction_adjusted_effect(metric, active_desc["median"], other_desc["median"])

        score_series = direction_adjusted_values(metric, sub[metric])
        score_mask = score_series.notna()
        labels = sub.loc[score_mask, "oligo_activity_true"].astype(int).to_numpy()
        scores = score_series.loc[score_mask].to_numpy(dtype=float)
        baseline = float(labels.mean()) if labels.size else np.nan
        raw_values = pd.to_numeric(sub.loc[score_mask, metric], errors="coerce").to_numpy(dtype=float)
        raw_observed = float(np.median(raw_values[labels.astype(bool)]) - np.median(raw_values[~labels.astype(bool)])) if (
            labels.size and labels.sum() and labels.sum() < labels.size
        ) else np.nan

        rows.append(
            {
                "Ligand scope": scope_label,
                "AA17 policy": "AA17 counted" if aa17_policy == "aa17_counted" else "AA17 excluded",
                "Metric": metric,
                "Metric label": metric_label,
                "Unit": unit,
                "Direction": direction,
                "n active": active_desc["n"],
                "n other": other_desc["n"],
                "Median active": active_desc["median"],
                "Median other": other_desc["median"],
                "Effect": effect,
                "Pr-AUC (baseline)": (
                    f"{pr_auc_score(labels, scores):.6g} ({baseline:.6g})"
                    if labels.size and np.isfinite(baseline)
                    else ""
                ),
                "ROC-AUC": roc_auc(labels, scores),
                "p-value": permutation_p_value(raw_values, labels, raw_observed, rng),
            }
        )
    out = pd.DataFrame(rows)
    out["FDR q-value"] = benjamini_hochberg(out["p-value"])
    return out.sort_values(["FDR q-value", "p-value", "Metric"], na_position="last").reset_index(drop=True)


def make_permutation_outputs(condition_labeled: pd.DataFrame, output_dir: Path) -> None:
    scopes = [
        ("all_ligands", "All predicted ligands", condition_labeled.copy()),
        (
            "demonstrated_ligands",
            "Predicted ligands with demonstrated activity",
            condition_labeled[condition_labeled["demonstrated_activity_condition"].eq(True)].copy(),
        ),
    ]
    policies = ["aa17_counted", "aa17_excluded"]
    seed_offset = 0
    for scope_slug, scope_label, scope_data in scopes:
        metrics = aggregate_protein_construct(scope_data) if len(scope_data) else pd.DataFrame(columns=["construct_group", "family"])
        for policy in policies:
            table = make_permutation_table(metrics, scope_label, policy, seed_offset)
            seed_offset += 1
            stem = f"oligo_activity_permutation_{scope_slug}_{policy}"
            table.to_csv(output_dir / f"{stem}.tsv", sep="\t", index=False)
            write_generic_markdown_table(table, output_dir / f"{stem}.md")


def write_notes(path: Path) -> None:
    path.write_text(
        """Oligo activity descriptive and comparative statistics

Grouping:
- Oligo activity true: oligo_active_updated.tsv has Oligo_active=Ja for the UniProt ID.
- Oligo activity false: oligo_active_updated.tsv has Oligo_active=Nei for the UniProt ID.
- No updated oligo activity record: the UniProt ID is present in the prediction run but absent from oligo_active_updated.tsv.
- Demonstrated activity condition: the predicted substrate class matches an active substrate for that protein. This is derived from oligo_active_updated.tsv, using Ligandkategori_EC, EC numbers, and Oligo_active_classes. EC mapping: 1.14.99.53 = chitin, 1.14.99.54/1.14.99.56 = cellulose, 1.14.99.55 = starch. Predicted amylose conditions are treated as starch.
- Oligo active vs all other 2x2 outputs ignore updated-record availability as a grouping variable. "All other" means all proteins in the prediction run without Oligo_active=Ja.
- The permutation outputs are written for all predicted ligands and for predicted ligands with demonstrated activity only. Each ligand scope has one table with AA17 counted and one where AA17 proteins are excluded before aggregation/statistics.

Aggregation:
- The main comparison table uses protein/construct-level rows, not raw condition rows.
- Construct-specific rows are calculated for Catalytic domain and Full-length separately.
- All constructs aggregates all predicted conditions for a protein across construct types.
- Fractions based on counts are calculated after summing counts across that protein/construct.
- Geometry metrics are medians across predicted substrate conditions for that protein/construct.

Key metrics:
- Hard-QC pass fraction = total n_stage1_pass / total n_generated.
- Hard-QC hard-fail fraction = total n_stage1_hard_fail / total n_generated.
- IFP success fraction = total n_ifp_success / total n_stage1_pass.
- Contact-eligible pose fraction = total n_contact_eligible / total n_ifp_success.
- Clusterable pose fraction = total n_ifp_clustered / total n_ifp_success.
- Noise fraction = total n_noise / total n_ifp_clustered, matching the interpretation of noise among clusterable poses.
- Formal clusterable condition fraction = fraction of substrate conditions with formal_clustering_allowed=True.
- Valid cluster condition fraction = fraction of substrate conditions with any_valid_cluster=True.
- N clusters is the number of pose clusters reported per predicted substrate condition.
- C1/C4 plausible fractions are occupancy-weighted fractions of clusters/poses compatible with C1 or C4 geometry.
- Cu-C1, Cu-C4, Oxyl-H(C1), and Oxyl-H(C4) distances are occupancy-weighted median distances in angstrom.
- Cu-Oxyl-H angles are occupancy-weighted median angles in degrees.
- Oxyl-H scores are geometry scores derived from Oxyl-H distance/geometry; higher values indicate stronger compatibility in the pipeline scoring.
- Best C1/C4 plausible fraction is max(C1 plausible fraction, C4 plausible fraction) per condition before aggregation.
- Best Cu distance is min(Cu-C1 distance, Cu-C4 distance) per condition before aggregation.
- Best Oxyl-H distance is min(Oxyl-H(C1) distance, Oxyl-H(C4) distance) per condition before aggregation.
- Best Oxyl-H score is max(Oxyl-H(C1) score, Oxyl-H(C4) score) per condition before aggregation.

Comparison columns:
- mean_difference_true_minus_comparator and median_difference_true_minus_comparator are positive when the Oligo_Activity=True group has higher values.
- cliffs_delta_true_vs_comparator ranges from -1 to 1. Positive means values tend to be higher in the Oligo_Activity=True group; negative means they tend to be lower.
- The 2x2 ANOVA table uses protein-level rows and the model value ~ Oligo_active + Demonstrated_activity + Oligo_active:Demonstrated_activity. P-values are descriptive because the Oligo_Activity=True group is small.
- The permutation tables use All constructs protein-level rows. Effect is median active minus median other, except metrics marked "lower is better" are direction-adjusted as median other minus median active. P-values are two-sided label-permutation p-values on the median difference. FDR q-values use Benjamini-Hochberg correction within each table.

Caution:
- Oligo_Activity=True is a small group in this run, so differences are descriptive and should not be treated as robust statistical inference without follow-up.
"""
    )


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    condition = pd.read_csv(args.results_dir / "combined" / "condition_table.tsv", sep="\t", low_memory=False)
    pdb_data = pd.read_csv(args.analyse_dir / "input_data" / "pdb_structure_data.csv", sep=";")
    metadata = pd.read_csv(args.analyse_dir / "input_data" / "metadata_final_ec_fixed.tsv", sep="\t")
    oligo_activity = pd.read_csv(args.analyse_dir / "input_data" / "oligo_active_updated.tsv", sep="\t")

    labels = make_protein_labels(condition, pdb_data, metadata, oligo_activity)
    activity = make_activity_table(oligo_activity)
    condition_labeled = add_condition_metrics(condition, labels, activity)
    protein_construct = aggregate_protein_construct(condition_labeled)
    group_counts = make_group_counts(labels, condition_labeled, protein_construct)
    group_descriptive = make_group_descriptive(protein_construct)
    comparisons = make_comparisons(protein_construct)
    make_demonstrated_activity_outputs(condition_labeled, args.output_dir)
    make_oligo_2x2_outputs(condition_labeled, args.output_dir)
    make_permutation_outputs(condition_labeled, args.output_dir)

    labels.to_csv(args.output_dir / "oligo_activity_protein_labels.tsv", sep="\t", index=False)
    condition_labeled.to_csv(args.output_dir / "oligo_activity_condition_metrics.tsv", sep="\t", index=False)
    protein_construct.to_csv(args.output_dir / "oligo_activity_protein_construct_metrics.tsv", sep="\t", index=False)
    group_counts.to_csv(args.output_dir / "oligo_activity_group_counts.tsv", sep="\t", index=False)
    group_descriptive.to_csv(args.output_dir / "oligo_activity_group_descriptive_statistics.tsv", sep="\t", index=False)
    comparisons.to_csv(args.output_dir / "oligo_activity_comparative_statistics.tsv", sep="\t", index=False)
    write_notes(args.output_dir / "oligo_activity_statistics_notes.txt")

    print(f"Wrote oligo activity statistics to {args.output_dir}")
    print(f"Oligo_Activity=True proteins in run: {int(labels['oligo_activity_true'].sum())}")
    print(f"Run proteins with PDB oligo records: {int(labels['has_pdb_oligo_record'].sum())}")


if __name__ == "__main__":
    main()
