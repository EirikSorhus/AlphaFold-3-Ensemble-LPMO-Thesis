#!/usr/bin/env python3
"""Correlate cellulose regio labels with C1/C4 geometric proximity bias."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pandas as pd


GEOMETRY_PAIRS = {
    "Cu_C1_minus_C4_distance": (
        "occupancy_weighted_Cu_C1_distance_median",
        "occupancy_weighted_Cu_C4_distance_median",
        "Cu-C1 minus Cu-C4 distance (Å)",
    ),
    "OxylH_C1_minus_C4_distance": (
        "occupancy_weighted_oxyl_H_C1_distance_median",
        "occupancy_weighted_oxyl_H_C4_distance_median",
        "Oxyl-H(C1) minus Oxyl-H(C4) distance (Å)",
    ),
}

CONSTRUCT_LABELS = {"domain_only": "Catalytic domain", "full_length": "Full-length", "all_constructs": "All constructs"}


def has_ec(value: object, ec_number: str) -> bool:
    return ec_number in str(value or "")


def infer_cellulose_regio_from_ec(value: object) -> str:
    has_c1 = has_ec(value, "1.14.99.54")
    has_c4 = has_ec(value, "1.14.99.56")
    if has_c1 and has_c4:
        return "C1+C4"
    if has_c1:
        return "C1"
    if has_c4:
        return "C4"
    return ""


def add_label_columns(active: pd.DataFrame) -> pd.DataFrame:
    active = active.copy()
    active["experimental_C1"] = active["regio_norm"].str.contains("C1").astype(int)
    active["experimental_C4"] = active["regio_norm"].str.contains("C4").astype(int)
    active["exclusive_regio_label"] = active["regio_norm"].map({"C1": 1, "C4": 0})
    active["regio_group"] = active["regio_norm"].map({"C1": "C1-only", "C4": "C4-only"}).fillna("Mixed C1/C4")
    active["chitin_and_cellulose_active"] = (
        active["ec_numbers"].fillna("").astype(str).map(lambda value: has_ec(value, "1.14.99.53"))
        & active["regio_norm"].isin(["C1", "C4", "C1+C4", "C4+C1", "C1/C4"])
    )
    return active


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--analyse-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def pearsonr(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    x_arr = pd.to_numeric(x, errors="coerce").to_numpy(dtype=float)
    y_arr = pd.to_numeric(y, errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(x_arr) & np.isfinite(y_arr)
    x_arr = x_arr[mask]
    y_arr = y_arr[mask]
    n = len(x_arr)
    if n < 3 or np.nanstd(x_arr) == 0 or np.nanstd(y_arr) == 0:
        return float("nan"), float("nan")
    r = float(np.corrcoef(x_arr, y_arr)[0, 1])
    r = max(-0.999999999, min(0.999999999, r))
    z = math.atanh(r) * math.sqrt(n - 3)
    p = 2.0 * (1.0 - NormalDist().cdf(abs(z)))
    return r, p


def rank_series(values: pd.Series) -> pd.Series:
    return values.rank(method="average")


def corr_rows(df: pd.DataFrame, *, construct_group: str, metric: str, label_col: str, analysis_set: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    sub = df.dropna(subset=[metric, label_col]).copy()
    n = len(sub)
    n_proteins = sub["protein_id"].nunique()
    if n:
        pearson_r, pearson_p = pearsonr(sub[metric], sub[label_col])
        spearman_r, spearman_p = pearsonr(rank_series(sub[metric]), rank_series(sub[label_col]))
        rows.append(
            {
                "analysis_set": analysis_set,
                "construct_group": construct_group,
                "metric": metric,
                "experimental_label": label_col,
                "comparison": "continuous_delta_vs_label",
                "n_rows": n,
                "n_proteins": n_proteins,
                "correlation_type": "Pearson/point-biserial",
                "correlation": pearson_r,
                "approx_p_value": pearson_p,
                "spearman_correlation": spearman_r,
                "spearman_approx_p_value": spearman_p,
                "interpretation": "Negative correlation means label=1 is associated with smaller C1-C4 distance, i.e. C1 is closer.",
            }
        )
    closer_col = f"{metric}_C1_closer"
    sub = df.dropna(subset=[closer_col, label_col]).copy()
    if len(sub):
        phi, phi_p = pearsonr(sub[closer_col], sub[label_col])
        tp = int(((sub[closer_col] == 1) & (sub[label_col] == 1)).sum())
        tn = int(((sub[closer_col] == 0) & (sub[label_col] == 0)).sum())
        fp = int(((sub[closer_col] == 1) & (sub[label_col] == 0)).sum())
        fn = int(((sub[closer_col] == 0) & (sub[label_col] == 1)).sum())
        rows.append(
            {
                "analysis_set": analysis_set,
                "construct_group": construct_group,
                "metric": metric,
                "experimental_label": label_col,
                "comparison": "C1_closer_binary_vs_label",
                "n_rows": len(sub),
                "n_proteins": sub["protein_id"].nunique(),
                "correlation_type": "Phi/Pearson binary-binary",
                "correlation": phi,
                "approx_p_value": phi_p,
                "spearman_correlation": phi,
                "spearman_approx_p_value": phi_p,
                "true_positive_C1_closer_and_label1": tp,
                "true_negative_C4_closer_and_label0": tn,
                "false_positive_C1_closer_label0": fp,
                "false_negative_C4_closer_label1": fn,
                "accuracy": (tp + tn) / len(sub) if len(sub) else float("nan"),
                "interpretation": "Positive phi means C1-closer geometry agrees with label=1.",
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    table_dir = args.output_dir / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)

    condition = pd.read_csv(args.results_dir / "combined" / "condition_table.tsv", sep="\t", low_memory=False)
    run_proteins = set(condition["protein_id"].dropna().astype(str))

    full_meta = pd.read_csv(args.analyse_dir / "input_data" / "metadata_final_ec_fixed.tsv", sep="\t")
    full_meta = full_meta.rename(columns={"UniProt_ID": "protein_id", "CAZy_family": "cazy_family", "EC_Number": "ec_numbers"})
    full_meta = full_meta[full_meta["protein_id"].astype(str).isin(run_proteins)].copy()
    full_meta["regio_norm"] = full_meta["ec_numbers"].map(infer_cellulose_regio_from_ec)
    ec_text = full_meta["ec_numbers"].fillna("").astype(str)
    cellulose_active = full_meta[
        ec_text.map(lambda value: has_ec(value, "1.14.99.54") or has_ec(value, "1.14.99.56"))
        & full_meta["regio_norm"].isin(["C1", "C4", "C1+C4"])
    ].copy()
    cellulose_active["activity_source"] = "metadata_ec_cellulose"

    active_cols = ["protein_id", "cazy_family", "ec_numbers", "regio_norm", "activity_source"]
    cellulose_active = (
        cellulose_active[active_cols]
        .sort_values(["protein_id"])
        .drop_duplicates(subset=["protein_id"], keep="first")
        .reset_index(drop=True)
    )
    cellulose_active = add_label_columns(cellulose_active)

    cols = ["protein_id", "construct_type", "condition_id", "substrate_class", "dp", *[col for pair in GEOMETRY_PAIRS.values() for col in pair[:2]]]
    df = condition[cols].copy()
    df = df[df["substrate_class"].astype(str).str.lower().eq("cellulose")]
    for col in [col for pair in GEOMETRY_PAIRS.values() for col in pair[:2]]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.merge(
        cellulose_active[
            [
                "protein_id",
                "cazy_family",
                "ec_numbers",
                "activity_source",
                "chitin_and_cellulose_active",
                "regio_norm",
                "regio_group",
                "experimental_C1",
                "experimental_C4",
                "exclusive_regio_label",
            ]
        ],
        on="protein_id",
        how="inner",
    )

    for metric, (c1_col, c4_col, _) in GEOMETRY_PAIRS.items():
        df[metric] = df[c1_col] - df[c4_col]
        df[f"{metric}_C1_closer"] = np.where(df[metric].notna(), (df[metric] < 0).astype(int), np.nan)
        df[f"{metric}_closer_site"] = np.where(df[metric] < 0, "C1", np.where(df[metric] > 0, "C4", "Tie"))

    agg_rows = []
    construct_groups = {
        "Catalytic domain": df[df["construct_type"].eq("domain_only")].copy(),
        "Full-length": df[df["construct_type"].eq("full_length")].copy(),
        "All constructs": df.copy(),
    }
    agg_cols = list(GEOMETRY_PAIRS.keys())
    for construct_group, sub in construct_groups.items():
        if sub.empty:
            continue
        grouped = (
            sub.groupby(
                [
                    "protein_id",
                    "cazy_family",
                    "ec_numbers",
                    "activity_source",
                    "chitin_and_cellulose_active",
                    "regio_norm",
                    "regio_group",
                    "experimental_C1",
                    "experimental_C4",
                    "exclusive_regio_label",
                ],
                dropna=False,
            )[agg_cols]
            .median()
            .reset_index()
        )
        grouped["construct_group"] = construct_group
        for metric in GEOMETRY_PAIRS:
            grouped[f"{metric}_C1_closer"] = np.where(grouped[metric].notna(), (grouped[metric] < 0).astype(int), np.nan)
            grouped[f"{metric}_closer_site"] = np.select(
                [grouped[metric].isna(), grouped[metric] < 0, grouped[metric] > 0],
                ["", "C1", "C4"],
                default="Tie",
            )
        agg_rows.append(grouped)
    per_protein = pd.concat(agg_rows, ignore_index=True)

    summary_rows = []
    for construct_group, sub in per_protein.groupby("construct_group", observed=True):
        for metric in GEOMETRY_PAIRS:
            summary_rows.extend(corr_rows(sub, construct_group=construct_group, metric=metric, label_col="exclusive_regio_label", analysis_set="exclusive_C1_vs_C4_only"))
            summary_rows.extend(corr_rows(sub, construct_group=construct_group, metric=metric, label_col="experimental_C1", analysis_set="all_cellulose_active_C1_label"))
            summary_rows.extend(corr_rows(sub, construct_group=construct_group, metric=metric, label_col="experimental_C4", analysis_set="all_cellulose_active_C4_label"))
    summary = pd.DataFrame(summary_rows)

    group_rows = []
    for construct_group, sub in per_protein.groupby("construct_group", observed=True):
        for metric, (_, _, label) in GEOMETRY_PAIRS.items():
            stats = (
                sub.dropna(subset=[metric])
                .groupby("regio_group", observed=True)
                .agg(
                    n_proteins=("protein_id", "nunique"),
                    median_delta=(metric, "median"),
                    mean_delta=(metric, "mean"),
                    n_C1_closer=(f"{metric}_C1_closer", "sum"),
                )
                .reset_index()
            )
            for row in stats.to_dict("records"):
                row["construct_group"] = construct_group
                row["metric"] = metric
                row["metric_label"] = label
                row["fraction_C1_closer"] = row["n_C1_closer"] / row["n_proteins"] if row["n_proteins"] else np.nan
                group_rows.append(row)
    group_summary = pd.DataFrame(group_rows)

    per_protein.to_csv(table_dir / "cellulose_regio_geometry_per_protein.tsv", sep="\t", index=False)
    summary.to_csv(table_dir / "cellulose_regio_geometry_correlation_summary.tsv", sep="\t", index=False)
    group_summary.to_csv(table_dir / "cellulose_regio_geometry_group_summary.tsv", sep="\t", index=False)

    notes = """Cellulose regio geometry correlation analysis

Input:
- combined/condition_table.tsv from the allow_single_cluster rerun.
- input_data/metadata_final_ec_fixed.tsv for EC-derived cellulose activity and regio labels.

Filtering:
- Only predicted cellulose conditions are used.
- Proteins are retained when metadata_final_ec_fixed.tsv has cellulose EC 1.14.99.54 and/or 1.14.99.56.
- The cellulose regio label is inferred from cellulose EC numbers only: 1.14.99.54 = C1, 1.14.99.56 = C4, both = C1+C4. Chitin EC 1.14.99.53 is C1 for chitin, but is not used to label cellulose regio in this cellulose-only analysis.
- Primary rows are protein-level summaries, not condition-level rows. For each protein and construct group, cellulose DP conditions are aggregated by median.

Geometry deltas:
- Cu_C1_minus_C4_distance = occupancy_weighted_Cu_C1_distance_median - occupancy_weighted_Cu_C4_distance_median.
- OxylH_C1_minus_C4_distance = occupancy_weighted_oxyl_H_C1_distance_median - occupancy_weighted_oxyl_H_C4_distance_median.
- Negative delta means C1 is closer. Positive delta means C4 is closer.

Correlation:
- exclusive_C1_vs_C4_only uses only C1-only and C4-only proteins; C1-only is encoded as 1 and C4-only as 0.
- all_cellulose_active_C1_label and all_cellulose_active_C4_label include mixed C1+C4 proteins as positive for both labels.
- Continuous correlations use Pearson/point-biserial and Spearman approximations.
- Binary correlations compare whether C1 is closer against the experimental binary label and report phi plus confusion counts.

Caution:
- The explicit cellulose regio-labelled set is small. Treat p-values as descriptive approximations, not strong inferential evidence.
"""
    (args.output_dir / "cellulose_regio_geometry_method_notes.txt").write_text(notes)


if __name__ == "__main__":
    main()
