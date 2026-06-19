#!/usr/bin/env python3
"""Render a curated holo-reference crystal anchoring Tanimoto plot.

This figure is intentionally narrower than the main Figure 7 renderer:

* keep only selected holo references;
* compare only matching substrate classes, while allowing all predicted DP
  variants for that substrate;
* use existing valid cluster medoids only;
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import MaxNLocator


CURATED_REFERENCES = {
    ("A0A0S2GKZ1", "5ACI"): {"substrate": "cellulose", "crystal_ligand": "CEL", "crystal_dp": 6},
    ("A0A0S2GKZ1", "7PXW"): {"substrate": "cellulose", "crystal_ligand": "CEL", "crystal_dp": 4},
    ("A0A223GEC9", "6YDC"): {"substrate": "cellulose", "crystal_ligand": "CEL", "crystal_dp": 6},
}

SUBSTRATE_FROM_LIGAND = {
    "CEL": "cellulose",
    "STA": "amylose",
    "GLC": "amylose",
}

REFERENCE_DP_PALETTE = {
    4: "#2f8f46",
    6: "#3b6ea8",
}

DP_MARKERS = {
    4: "o",
    6: "s",
    8: "^",
}

CONSTRUCT_PALETTE = {"domain_only": "#3b6ea8", "full_length": "#d17a22"}
CONSTRUCT_DISPLAY = {"domain_only": "Catalytic Domain", "full_length": "Full-Length"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--figure-name", default="figure7_crystal_anchoring_curated_holo.png")
    return parser.parse_args()


def parse_ligand_id(ligand_id: Any) -> tuple[str, int | None]:
    match = re.match(r"^([A-Za-z]+)(\d+)$", str(ligand_id or "").strip())
    if match is None:
        return "", None
    return match.group(1).upper(), int(match.group(2))


def build_curated_table(results_dir: Path) -> pd.DataFrame:
    anchor = pd.read_csv(results_dir / "combined" / "crystal_anchor_table.tsv", sep="\t", low_memory=False)
    condition = pd.read_csv(results_dir / "combined" / "condition_table.tsv", sep="\t", low_memory=False)

    condition_keep = condition[["condition_id", "any_valid_cluster", "substrate_class", "dp"]].copy()
    condition_keep["any_valid_cluster_bool"] = condition_keep["any_valid_cluster"].astype(str).str.lower() == "true"
    df = anchor.merge(
        condition_keep.drop_duplicates("condition_id"),
        on="condition_id",
        how="left",
        suffixes=("", "_condition"),
    )

    parsed = df["ligand_id"].map(parse_ligand_id)
    df["ligand_code"] = [item[0] for item in parsed]
    df["predicted_dp"] = [item[1] for item in parsed]
    df["predicted_substrate"] = df["ligand_code"].map(SUBSTRATE_FROM_LIGAND)
    df["curated_key"] = list(zip(df["protein_id"].astype(str), df["crystal_reference_id"].astype(str)))
    df["curated_reference"] = df["curated_key"].map(CURATED_REFERENCES)
    df = df[df["curated_reference"].notna()].copy()

    df["crystal_substrate"] = df["curated_reference"].map(lambda item: item["substrate"])
    df["crystal_ligand"] = df["curated_reference"].map(lambda item: item["crystal_ligand"])
    df["crystal_dp"] = df["curated_reference"].map(lambda item: int(item["crystal_dp"]))
    df = df[
        (df["representative_role"].astype(str) == "cluster_medoid")
        & (df["any_valid_cluster_bool"])
        & (df["predicted_substrate"] == df["crystal_substrate"])
        & (df["comparison_status"].astype(str) == "ok")
    ].copy()

    df["ifp_tanimoto"] = pd.to_numeric(df["ifp_tanimoto"], errors="coerce")
    df["local_pocket_rmsd"] = pd.to_numeric(df["local_pocket_rmsd"], errors="coerce")
    df = df[df["ifp_tanimoto"].notna() & df["local_pocket_rmsd"].notna()].copy()

    df["exact_holo_dp_match"] = df["predicted_dp"] == df["crystal_dp"]
    cluster_numeric = pd.to_numeric(df["cluster_id"], errors="coerce")
    df["cluster_label"] = df["cluster_id"].astype(str)
    df.loc[cluster_numeric.notna(), "cluster_label"] = (
        (cluster_numeric[cluster_numeric.notna()].astype(int) + 1).astype(str)
    )
    df["point_label"] = (
        df["protein_id"].astype(str)
        + " C"
        + df["cluster_label"].astype(str)
        + " vs "
        + df["crystal_reference_id"].astype(str)
    )
    return df.sort_values(["protein_id", "ligand_id", "crystal_reference_id", "representative_pose_id"])


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin(["true", "1", "yes", "ja"])


def finite_numeric(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = df.copy()
    for column in columns:
        out[column] = pd.to_numeric(out[column], errors="coerce")
        out = out[out[column].notna() & np.isfinite(out[column])].copy()
    return out


def render_plot(curated_df: pd.DataFrame, results_dir: Path, output_dir: Path, figure_name: str) -> None:
    sns.set_theme(
        context="paper",
        style="ticks",
        font_scale=1.15,
        rc={
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#333333",
            "axes.labelcolor": "#222222",
            "xtick.color": "#222222",
            "ytick.color": "#222222",
            "font.family": "DejaVu Sans",
        },
    )
    anchor = pd.read_csv(results_dir / "combined" / "crystal_anchor_table.tsv", sep="\t", low_memory=False)
    condition = pd.read_csv(results_dir / "combined" / "condition_table.tsv", sep="\t", low_memory=False)
    merged = anchor.merge(
        condition[["condition_id", "substrate_class", "dp", "construct_type"]].drop_duplicates("condition_id"),
        on="condition_id",
        how="left",
    )

    fig = plt.figure(figsize=(10.2, 7.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 0.68], width_ratios=[0.35, 1.0, 1.0, 0.35])
    ax_rmsd = fig.add_subplot(gs[0, 0:2])
    ax_ifp = fig.add_subplot(gs[0, 2:4])
    ax_cluster = fig.add_subplot(gs[1, 1:3])

    rmsd = finite_numeric(merged, ["local_pocket_rmsd"])
    rmsd["substrate_display"] = rmsd["substrate_class"].astype(str).str.capitalize()
    rmsd["construct_display"] = rmsd["construct_type"].astype(str).map(CONSTRUCT_DISPLAY).fillna(
        rmsd["construct_type"].astype(str).str.replace("_", " ").str.title()
    )
    construct_palette = {CONSTRUCT_DISPLAY[key]: value for key, value in CONSTRUCT_PALETTE.items()}
    sns.boxplot(
        data=rmsd,
        x="substrate_display",
        y="local_pocket_rmsd",
        hue="construct_display",
        showfliers=False,
        palette=construct_palette,
        ax=ax_rmsd,
    )
    sns.stripplot(
        data=rmsd,
        x="substrate_display",
        y="local_pocket_rmsd",
        hue="construct_display",
        dodge=True,
        color="#111111",
        alpha=0.45,
        size=2.2,
        ax=ax_rmsd,
        legend=False,
    )
    ax_rmsd.set_xlabel("Substrate")
    ax_rmsd.set_ylabel("RMSD (A)")
    ax_rmsd.legend(
        title="Construct",
        loc="upper center",
        bbox_to_anchor=(0.5, 1.13),
        ncol=2,
        frameon=False,
    )
    ax_rmsd.margins(x=0.02)
    ax_rmsd.text(-0.10, 1.02, "A", transform=ax_rmsd.transAxes, fontsize=16, fontweight="bold")

    for reference_dp, color in REFERENCE_DP_PALETTE.items():
        for dp, marker in DP_MARKERS.items():
            subset = curated_df[
                (curated_df["crystal_dp"] == reference_dp)
                & (curated_df["predicted_dp"] == dp)
            ]
            if subset.empty:
                continue
            ax_ifp.scatter(
                subset["local_pocket_rmsd"],
                subset["ifp_tanimoto"],
                s=70,
                marker=marker,
                color=color,
                edgecolor="#202020",
                linewidth=0.7,
                alpha=0.88,
            )

    ax_ifp.set_xlabel("Local Pocket RMSD (A)")
    ax_ifp.set_ylabel("IFP Tanimoto Similarity")
    ax_ifp.text(-0.10, 1.02, "B", transform=ax_ifp.transAxes, fontsize=16, fontweight="bold")
    y_max = float(pd.to_numeric(curated_df["ifp_tanimoto"], errors="coerce").max())
    ax_ifp.set_ylim(-0.035, max(0.5, y_max + 0.08))
    sns.despine(ax=ax_rmsd)
    sns.despine(ax=ax_ifp)

    label_rows = curated_df.sort_values(
        ["ifp_tanimoto", "local_pocket_rmsd", "point_label"],
        ascending=[False, True, True],
    ).reset_index(drop=True)
    x_values = pd.to_numeric(curated_df["local_pocket_rmsd"], errors="coerce")
    ax_ifp.set_xlim(float(x_values.min()) - 0.006, float(x_values.max()) + 0.006)

    top5 = label_rows.head(5).copy()
    remaining = label_rows.drop(top5.index)
    left4 = remaining.sort_values(["local_pocket_rmsd", "ifp_tanimoto"]).head(4).copy()
    remaining = remaining.drop(left4.index)

    label_specs: list[tuple[pd.Series, float, float, str]] = []
    for y_text, (_, row) in zip(np.linspace(0.91, 0.67, len(top5)), top5.iterrows()):
        label_specs.append((row, 0.78, float(y_text), "right"))
    for y_text, (_, row) in zip(np.linspace(0.47, 0.29, len(left4)), left4.iterrows()):
        label_specs.append((row, 0.06, float(y_text), "left"))
    for y_text, (_, row) in zip(np.linspace(0.58, 0.08, len(remaining)), remaining.iterrows()):
        label_specs.append((row, 1.03, float(y_text), "left"))

    for row, x_text, y_text, ha in label_specs:
        ax_ifp.annotate(
            str(row["point_label"]),
            (row["local_pocket_rmsd"], row["ifp_tanimoto"]),
            xytext=(x_text, y_text),
            xycoords="data",
            textcoords=ax_ifp.transAxes,
            ha=ha,
            va="center",
            fontsize=8.2,
            color="#222222",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.92, "pad": 0.6},
            arrowprops={"arrowstyle": "-", "color": "#777777", "linewidth": 0.5},
            annotation_clip=False,
            zorder=5,
        )

    color_handles = [
        plt.Line2D([0], [0], marker="o", color="w", label=f"Reference DP{dp}", markerfacecolor=color,
                   markeredgecolor="#202020", markersize=8)
        for dp, color in REFERENCE_DP_PALETTE.items()
        if dp in set(curated_df["crystal_dp"])
    ]
    marker_handles = [
        plt.Line2D([0], [0], marker=marker, color="#444444", label=f"Predicted DP{dp}",
                   markerfacecolor="white", markeredgecolor="#202020", linestyle="", markersize=8)
        for dp, marker in DP_MARKERS.items()
        if dp in set(curated_df["predicted_dp"])
    ]
    ax_ifp.legend(
        handles=color_handles + marker_handles,
        title="Reference DP / Predicted DP",
        loc="upper left",
        bbox_to_anchor=(0.95, 1.0),
        frameon=False,
        fontsize=9.5,
        title_fontsize=9.5,
    )

    pdb = pd.read_csv(results_dir.parent.parent / "input_data" / "metadata_pdb_count.tsv", sep="\t", low_memory=False)
    pdb = pdb.rename(columns={"UniProt_ID": "protein_id"})
    pdb["PDB_Before_AF3_Cutoff_Count"] = pd.to_numeric(pdb["PDB_Before_AF3_Cutoff_Count"], errors="coerce")
    rates = (
        condition.assign(valid=as_bool(condition["any_valid_cluster"]))
        .groupby("protein_id", observed=True)["valid"]
        .mean()
        .reset_index()
    )
    rates = rates.merge(
        pdb[["protein_id", "PDB_Before_AF3_Cutoff_Count"]],
        on="protein_id",
        how="left",
    ).fillna({"PDB_Before_AF3_Cutoff_Count": 0})
    sns.scatterplot(
        data=rates,
        x="PDB_Before_AF3_Cutoff_Count",
        y="valid",
        color="#3b6ea8",
        ax=ax_cluster,
    )
    x = rates["PDB_Before_AF3_Cutoff_Count"].astype(float).to_numpy()
    y = rates["valid"].astype(float).to_numpy()
    r2 = np.nan
    if len(rates) >= 2 and np.nanstd(x) > 0:
        slope, intercept = np.polyfit(x, y, 1)
        xs = np.linspace(np.nanmin(x), np.nanmax(x), 100)
        ax_cluster.plot(xs, slope * xs + intercept, color="#555555", linestyle="--", linewidth=1.2)
        pred = slope * x + intercept
        ss_res = np.nansum((y - pred) ** 2)
        ss_tot = np.nansum((y - np.nanmean(y)) ** 2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot else np.nan
    ax_cluster.set_ylim(-0.03, 1.03)
    if len(rates):
        xmin = float(np.nanmin(x))
        xmax = float(np.nanmax(x))
        ax_cluster.set_xlim(xmin - 0.35, xmax + 0.35)
    ax_cluster.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax_cluster.set_xlabel("X-Ray Crystal References Before AF3 Cutoff")
    ax_cluster.set_ylabel("Valid-Cluster Fraction")
    ax_cluster.text(
        0.98,
        0.04,
        f"R$^2$ = {r2:.2f}" if np.isfinite(r2) else "R$^2$ = NA",
        transform=ax_cluster.transAxes,
        ha="right",
        va="bottom",
        fontsize=10,
    )
    ax_cluster.text(-0.10, 1.02, "C", transform=ax_cluster.transAxes, fontsize=16, fontweight="bold")
    sns.despine(ax=ax_cluster)

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / figure_name
    fig.savefig(output_path, dpi=600, bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    results_dir = args.results_dir.resolve()
    output_dir = (args.output_dir or results_dir / "postprocess" / "final_figures").resolve()
    table_dir = output_dir / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)

    df = build_curated_table(results_dir)
    table_path = table_dir / "figure7_crystal_anchoring_curated_holo.tsv"
    df.to_csv(table_path, sep="\t", index=False)
    render_plot(df, results_dir, output_dir, args.figure_name)

    print(f"Rows plotted: {len(df)}")
    print(f"Conditions plotted: {df['condition_id'].nunique()}")
    print(f"Table: {table_path}")
    print(f"Figure: {output_dir / args.figure_name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
