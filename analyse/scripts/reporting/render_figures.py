#!/usr/bin/env python3
"""Render downstream figures from frozen LPMO analysis result tables."""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


SUBSTRATE_ORDER = ["chitin", "cellulose", "starch"]
CONSTRUCT_ORDER = ["domain_only", "full_length"]
CONSTRUCT_PALETTE = {
    "domain_only": "#3b6ea8",
    "full_length": "#d17a22",
}
SUBSTRATE_PALETTE = {
    "chitin": "#2b8cbe",
    "cellulose": "#41ab5d",
    "starch": "#e34a33",
    "amylose": "#e34a33",
}
REGIO_PALETTE = {
    "C1": "#3b6ea8",
    "C4": "#d17a22",
}


@dataclass(frozen=True)
class FigureRecord:
    figure_id: str
    title: str
    source_files: str
    filters: str
    n_rows_after_filter: int
    output_path: str
    status: str = "ok"
    message: str = ""


class FigureRenderer:
    def __init__(self, results_dir: Path, output_dir: Path | None = None) -> None:
        self.results_dir = results_dir.resolve()
        self.output_dir = (output_dir or self.results_dir / "postprocess" / "figures").resolve()
        self.records: list[FigureRecord] = []
        self.tables: dict[Path, pd.DataFrame] = {}
        self.output_dir.mkdir(parents=True, exist_ok=True)
        for subdir in [
            "descriptive",
            "crystal",
            "family_residue",
            "predictive",
            "cbm_paired",
            "structures",
        ]:
            (self.output_dir / subdir).mkdir(parents=True, exist_ok=True)

    def rel(self, path: Path) -> str:
        try:
            return str(path.resolve().relative_to(self.results_dir))
        except ValueError:
            return str(path)

    def read_tsv(self, rel_path: str) -> pd.DataFrame:
        path = self.results_dir / rel_path
        if path in self.tables:
            return self.tables[path].copy()
        if not path.is_file():
            raise FileNotFoundError(path)
        df = pd.read_csv(path, sep="\t", low_memory=False)
        self.tables[path] = df
        return df.copy()

    def read_json(self, rel_path: str) -> dict:
        path = self.results_dir / rel_path
        if not path.is_file():
            raise FileNotFoundError(path)
        return json.loads(path.read_text())

    def save(
        self,
        fig: plt.Figure,
        figure_id: str,
        title: str,
        source_files: list[str],
        filters: str,
        n_rows_after_filter: int,
        dpi: int = 220,
    ) -> None:
        output_path = self.output_dir / figure_id
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        self.records.append(
            FigureRecord(
                figure_id=figure_id,
                title=title,
                source_files=";".join(source_files),
                filters=filters,
                n_rows_after_filter=int(n_rows_after_filter),
                output_path=self.rel(output_path),
            )
        )

    def skip(
        self,
        figure_id: str,
        title: str,
        source_files: list[str],
        message: str,
    ) -> None:
        output_path = self.output_dir / figure_id
        if output_path.exists():
            output_path.unlink()
        self.records.append(
            FigureRecord(
                figure_id=figure_id,
                title=title,
                source_files=";".join(source_files),
                filters="not_rendered",
                n_rows_after_filter=0,
                output_path=self.rel(output_path),
                status="skipped",
                message=message,
            )
        )

    def manifest(self) -> None:
        manifest_path = self.output_dir / "figure_manifest.tsv"
        df = pd.DataFrame([record.__dict__ for record in self.records])
        df.to_csv(manifest_path, sep="\t", index=False)


def configure_style() -> None:
    sns.set_theme(
        context="paper",
        style="whitegrid",
        font_scale=0.9,
        rc={
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#222222",
            "axes.labelcolor": "#222222",
            "xtick.color": "#222222",
            "ytick.color": "#222222",
            "text.color": "#222222",
            "savefig.facecolor": "white",
            "axes.titleweight": "bold",
        },
    )


def require_columns(df: pd.DataFrame, columns: list[str]) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise KeyError("missing columns: " + ", ".join(missing))


def coerce_numeric(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    for column in columns:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")
    return df


def add_plot_categories(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "construct_type" not in df.columns and "_construct_type" in df.columns:
        df["construct_type"] = df["_construct_type"]
    if "construct_type" in df.columns:
        df["construct_type"] = df["construct_type"].fillna(df.get("_construct_type"))
        known = [value for value in CONSTRUCT_ORDER if value in set(df["construct_type"].dropna())]
        extra = sorted(set(df["construct_type"].dropna()) - set(known))
        df["construct_type"] = pd.Categorical(df["construct_type"], known + extra, ordered=True)
    if "substrate_class" in df.columns:
        known = [value for value in SUBSTRATE_ORDER if value in set(df["substrate_class"].dropna())]
        extra = sorted(set(df["substrate_class"].dropna()) - set(known))
        df["substrate_class"] = pd.Categorical(df["substrate_class"], known + extra, ordered=True)
    if "dp" in df.columns:
        df["dp"] = pd.to_numeric(df["dp"], errors="coerce")
    if {"substrate_class", "dp"}.issubset(df.columns):
        dp_label = df["dp"].map(lambda value: f"DP{int(value)}" if pd.notna(value) else "DP?")
        df["substrate_dp"] = df["substrate_class"].astype(str) + " " + dp_label
    return df


def substrate_order(df: pd.DataFrame) -> list[str]:
    if "substrate_class" not in df.columns:
        return []
    values = [str(value) for value in df["substrate_class"].dropna().unique()]
    known = [value for value in SUBSTRATE_ORDER if value in values]
    extra = sorted(set(values) - set(known))
    return known + extra


def substrate_palette(df: pd.DataFrame) -> dict[str, str]:
    order = substrate_order(df)
    fallback = sns.color_palette("Set2", n_colors=max(1, len(order))).as_hex()
    return {
        value: SUBSTRATE_PALETTE.get(value, fallback[idx % len(fallback)])
        for idx, value in enumerate(order)
    }


def finite_df(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = df.copy()
    for column in columns:
        out = out[pd.to_numeric(out[column], errors="coerce").notna()]
        out = out[np.isfinite(pd.to_numeric(out[column], errors="coerce"))]
    return out


def nice_metric_name(name: str) -> str:
    name = re.sub(r"^(delta_|mean_|median_)", "", name)
    name = name.replace("occupancy_weighted_", "")
    name = name.replace("_", " ")
    return name


def heatmap_size(n_rows: int, n_cols: int, min_w: float = 6.0, min_h: float = 4.0) -> tuple[float, float]:
    width = max(min_w, min(16.0, 0.45 * n_cols + 3.0))
    height = max(min_h, min(18.0, 0.22 * n_rows + 2.8))
    return width, height


def plot_qc_attrition(renderer: FigureRenderer) -> None:
    figure_id = "descriptive/qc_attrition_by_construct_substrate_dp.png"
    title = "QC attrition by construct, substrate and DP"
    source = "combined/condition_table.tsv"
    try:
        df = add_plot_categories(renderer.read_tsv(source))
        require_columns(
            df,
            [
                "construct_type",
                "substrate_class",
                "dp",
                "n_generated",
                "n_prep_error",
                "n_stage1_hard_fail",
                "n_stage1_pass",
            ],
        )
        count_cols = ["n_prep_error", "n_stage1_hard_fail", "n_stage1_pass"]
        df = coerce_numeric(df, ["n_generated", *count_cols])
        grouped = (
            df.groupby(["construct_type", "substrate_class", "dp"], observed=True)[count_cols]
            .sum()
            .reset_index()
        )
        grouped = add_plot_categories(grouped)
        grouped = grouped.sort_values(["construct_type", "substrate_class", "dp"])
        fig, axes = plt.subplots(
            1,
            max(1, grouped["construct_type"].nunique()),
            figsize=(13, 4.8),
            sharey=True,
            squeeze=False,
        )
        colors = {
            "n_stage1_pass": "#3a923a",
            "n_stage1_hard_fail": "#c44e52",
            "n_prep_error": "#8172b3",
        }
        labels = {
            "n_stage1_pass": "Stage 1 pass",
            "n_stage1_hard_fail": "Hard fail",
            "n_prep_error": "Prep error",
        }
        for ax, construct in zip(axes.flat, grouped["construct_type"].dropna().unique(), strict=False):
            sub = grouped[grouped["construct_type"] == construct].copy()
            x = np.arange(len(sub))
            bottom = np.zeros(len(sub))
            for col in ["n_stage1_pass", "n_stage1_hard_fail", "n_prep_error"]:
                ax.bar(x, sub[col].fillna(0), bottom=bottom, color=colors[col], label=labels[col])
                bottom += sub[col].fillna(0).to_numpy()
            ax.set_title(str(construct))
            ax.set_xticks(x)
            ax.set_xticklabels(sub["substrate_dp"], rotation=45, ha="right")
            ax.set_xlabel("")
            ax.set_ylabel("Number of generated poses")
        axes.flat[0].legend(frameon=False, loc="upper left")
        fig.suptitle(title)
        renderer.save(fig, figure_id, title, [source], "summed per construct/substrate/DP", len(grouped))
    except Exception as exc:
        renderer.skip(figure_id, title, [source], str(exc))


def plot_cluster_count_distribution(renderer: FigureRenderer) -> None:
    figure_id = "descriptive/cluster_count_distribution.png"
    title = "Cluster count distribution"
    source = "combined/condition_table.tsv"
    try:
        df = add_plot_categories(renderer.read_tsv(source))
        require_columns(df, ["construct_type", "substrate_class", "dp", "n_clusters"])
        df = finite_df(coerce_numeric(df, ["n_clusters"]), ["n_clusters"])
        fig, ax = plt.subplots(figsize=(10, 5))
        sns.violinplot(
            data=df,
            x="substrate_class",
            y="n_clusters",
            hue="construct_type",
            order=substrate_order(df),
            hue_order=CONSTRUCT_ORDER,
            palette=CONSTRUCT_PALETTE,
            cut=0,
            inner="quartile",
            ax=ax,
        )
        sns.stripplot(
            data=df,
            x="substrate_class",
            y="n_clusters",
            hue="construct_type",
            order=substrate_order(df),
            hue_order=CONSTRUCT_ORDER,
            dodge=True,
            color="#222222",
            alpha=0.25,
            size=2,
            ax=ax,
            legend=False,
        )
        ax.set_title(title)
        ax.set_xlabel("Substrate")
        ax.set_ylabel("Number of clusters")
        ax.legend(title="Construct", frameon=False)
        renderer.save(fig, figure_id, title, [source], "non-missing n_clusters", len(df))
    except Exception as exc:
        renderer.skip(figure_id, title, [source], str(exc))


def plot_top_cluster_occupancy(renderer: FigureRenderer) -> None:
    figure_id = "descriptive/top_cluster_occupancy_distribution.png"
    title = "Top-cluster occupancy distribution"
    source = "combined/condition_table.tsv"
    try:
        df = add_plot_categories(renderer.read_tsv(source))
        require_columns(df, ["construct_type", "substrate_class", "top_cluster_occupancy"])
        df = finite_df(coerce_numeric(df, ["top_cluster_occupancy"]), ["top_cluster_occupancy"])
        fig, ax = plt.subplots(figsize=(10, 5))
        sns.boxplot(
            data=df,
            x="substrate_class",
            y="top_cluster_occupancy",
            hue="construct_type",
            order=substrate_order(df),
            hue_order=CONSTRUCT_ORDER,
            palette=CONSTRUCT_PALETTE,
            showfliers=False,
            ax=ax,
        )
        sns.stripplot(
            data=df,
            x="substrate_class",
            y="top_cluster_occupancy",
            hue="construct_type",
            order=substrate_order(df),
            hue_order=CONSTRUCT_ORDER,
            dodge=True,
            color="#222222",
            alpha=0.25,
            size=2,
            ax=ax,
            legend=False,
        )
        ax.set_ylim(-0.03, 1.03)
        ax.set_title(title)
        ax.set_xlabel("Substrate")
        ax.set_ylabel("Top-cluster occupancy")
        ax.legend(title="Construct", frameon=False)
        renderer.save(fig, figure_id, title, [source], "non-missing top_cluster_occupancy", len(df))
    except Exception as exc:
        renderer.skip(figure_id, title, [source], str(exc))


def plot_c1_c4_plausibility(renderer: FigureRenderer) -> None:
    figure_id = "descriptive/c1_c4_plausibility_by_substrate_dp.png"
    title = "C1/C4 plausibility by substrate and DP"
    source = "combined/condition_table.tsv"
    try:
        df = add_plot_categories(renderer.read_tsv(source))
        c1 = "occupancy_weighted_c1_plausible_fraction"
        c4 = "occupancy_weighted_c4_plausible_fraction"
        require_columns(df, ["construct_type", "substrate_class", "dp", c1, c4])
        df = coerce_numeric(df, [c1, c4])
        long = df.melt(
            id_vars=["construct_type", "substrate_class", "dp", "substrate_dp"],
            value_vars=[c1, c4],
            var_name="site",
            value_name="plausible_fraction",
        )
        long["site"] = long["site"].map({c1: "C1", c4: "C4"})
        long = finite_df(long, ["plausible_fraction"])
        g = sns.catplot(
            data=long,
            x="substrate_dp",
            y="plausible_fraction",
            hue="site",
            col="construct_type",
            kind="point",
            errorbar=("pi", 50),
            palette=REGIO_PALETTE,
            height=4.5,
            aspect=1.35,
            sharey=True,
        )
        g.set_axis_labels("", "Occupancy-weighted plausible fraction")
        g.set_titles("{col_name}")
        for ax in g.axes.flat:
            ax.tick_params(axis="x", rotation=45)
            ax.set_ylim(-0.03, 1.03)
        g.fig.suptitle(title, y=1.04)
        renderer.save(g.fig, figure_id, title, [source], "non-missing C1/C4 plausible fractions", len(long))
    except Exception as exc:
        renderer.skip(figure_id, title, [source], str(exc))


def plot_cu_c1_vs_cu_c4(renderer: FigureRenderer) -> None:
    figure_id = "descriptive/cu_c1_vs_cu_c4_cluster_scatter.png"
    title = "Cluster Cu-C1 versus Cu-C4 distances"
    source = "combined/cluster_table.tsv"
    try:
        df = add_plot_categories(renderer.read_tsv(source))
        require_columns(
            df,
            [
                "construct_type",
                "substrate_class",
                "Cu_C1_distance_median",
                "Cu_C4_distance_median",
                "occupancy",
            ],
        )
        df = finite_df(coerce_numeric(df, ["Cu_C1_distance_median", "Cu_C4_distance_median", "occupancy"]), ["Cu_C1_distance_median", "Cu_C4_distance_median"])
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        sns.scatterplot(
            data=df,
            x="Cu_C1_distance_median",
            y="Cu_C4_distance_median",
            hue="substrate_class",
            style="construct_type",
            size="occupancy",
            sizes=(20, 220),
            palette=substrate_palette(df),
            ax=ax,
        )
        xy = [
            min(df["Cu_C1_distance_median"].min(), df["Cu_C4_distance_median"].min()),
            max(df["Cu_C1_distance_median"].max(), df["Cu_C4_distance_median"].max()),
        ]
        ax.plot(xy, xy, color="#555555", linewidth=1, linestyle="--")
        ax.set_title(title)
        ax.set_xlabel("Median Cu-C1 distance")
        ax.set_ylabel("Median Cu-C4 distance")
        ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
        renderer.save(fig, figure_id, title, [source], "non-missing C1/C4 distance medians", len(df))
    except Exception as exc:
        renderer.skip(figure_id, title, [source], str(exc))


def plot_top_residue_contact_heatmap(renderer: FigureRenderer) -> None:
    figure_id = "descriptive/top_residue_contact_heatmap.png"
    title = "Top residue contact scores"
    source = "combined/protein_condition_residue_scores.tsv"
    try:
        df = add_plot_categories(renderer.read_tsv(source))
        require_columns(df, ["residue_label", "substrate_class", "residue_contact_score"])
        df = finite_df(coerce_numeric(df, ["residue_contact_score"]), ["residue_contact_score"])
        top_residues = (
            df.groupby("residue_label", observed=True)["residue_contact_score"]
            .mean()
            .sort_values(ascending=False)
            .head(45)
            .index
        )
        sub = df[df["residue_label"].isin(top_residues)].copy()
        pivot = sub.pivot_table(
            index="residue_label",
            columns="substrate_class",
            values="residue_contact_score",
            aggfunc="mean",
            observed=True,
        )
        pivot = pivot.loc[top_residues]
        width, height = heatmap_size(len(pivot), len(pivot.columns))
        fig, ax = plt.subplots(figsize=(width, height))
        sns.heatmap(pivot, cmap="viridis", ax=ax, cbar_kws={"label": "Mean contact score"})
        ax.set_title(title)
        ax.set_xlabel("Substrate")
        ax.set_ylabel("Residue")
        renderer.save(fig, figure_id, title, [source], "top 45 residues by mean contact score", len(sub))
    except Exception as exc:
        renderer.skip(figure_id, title, [source], str(exc))


def plot_contact_chemistry(renderer: FigureRenderer) -> None:
    figure_id = "descriptive/contact_chemistry_by_substrate_dp.png"
    title = "Contact chemistry by substrate and DP"
    source = "combined/condition_patch_summary.tsv"
    try:
        df = add_plot_categories(renderer.read_tsv(source))
        metrics = [
            "aromatic_contact_fraction",
            "polar_contact_fraction",
            "charged_contact_fraction",
            "hbond_contact_fraction",
            "catalytic_surface_contact_fraction",
            "non_core_contact_fraction",
            "cbm_contact_fraction",
            "linker_contact_fraction",
        ]
        available = [metric for metric in metrics if metric in df.columns]
        require_columns(df, ["construct_type", "substrate_class", "dp", *available])
        df = coerce_numeric(df, available)
        long = df.melt(
            id_vars=["construct_type", "substrate_class", "dp", "substrate_dp"],
            value_vars=available,
            var_name="metric",
            value_name="fraction",
        )
        long = finite_df(long, ["fraction"])
        long["metric"] = long["metric"].map(nice_metric_name)
        g = sns.catplot(
            data=long,
            x="substrate_dp",
            y="fraction",
            hue="construct_type",
            col="metric",
            col_wrap=4,
            kind="point",
            errorbar=("pi", 50),
            palette=CONSTRUCT_PALETTE,
            height=3.1,
            aspect=1.15,
            sharey=True,
        )
        g.set_axis_labels("", "Fraction")
        g.set_titles("{col_name}")
        for ax in g.axes.flat:
            ax.tick_params(axis="x", rotation=45)
            ax.set_ylim(-0.03, 1.03)
        g.fig.suptitle(title, y=1.03)
        renderer.save(g.fig, figure_id, title, [source], "non-missing contact chemistry fractions", len(long))
    except Exception as exc:
        renderer.skip(figure_id, title, [source], str(exc))


def plot_pocket_rmsd_vs_ifp(renderer: FigureRenderer) -> None:
    figure_id = "crystal/pocket_rmsd_vs_ifp_similarity.png"
    title = "Pocket RMSD versus IFP similarity"
    source = "combined/crystal_anchor_table.tsv"
    try:
        anchor = renderer.read_tsv(source)
        cond = add_plot_categories(renderer.read_tsv("combined/condition_table.tsv"))
        require_columns(anchor, ["condition_id", "local_pocket_rmsd", "ifp_tanimoto"])
        meta_cols = [col for col in ["condition_id", "construct_type", "substrate_class", "dp"] if col in cond.columns]
        df = anchor.merge(cond[meta_cols].drop_duplicates("condition_id"), on="condition_id", how="left")
        df = add_plot_categories(df)
        df = finite_df(coerce_numeric(df, ["local_pocket_rmsd", "ifp_tanimoto"]), ["local_pocket_rmsd"])
        fig, ax = plt.subplots(figsize=(7, 5.5))
        sns.scatterplot(
            data=df,
            x="local_pocket_rmsd",
            y="ifp_tanimoto",
            hue="substrate_class" if "substrate_class" in df.columns else None,
            style="construct_type" if "construct_type" in df.columns else None,
            palette=substrate_palette(df),
            s=65,
            ax=ax,
        )
        ax.set_title(title)
        ax.set_xlabel("Local pocket RMSD")
        ax.set_ylabel("IFP Tanimoto")
        ax.set_ylim(-0.03, 1.03)
        ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
        renderer.save(fig, figure_id, title, [source, "combined/condition_table.tsv"], "non-missing pocket RMSD; IFP missing kept as NaN", len(df))
    except Exception as exc:
        renderer.skip(figure_id, title, [source], str(exc))


def plot_clustering_rate_vs_crystal_coverage(renderer: FigureRenderer) -> None:
    figure_id = "crystal/clustering_rate_vs_crystal_coverage.png"
    title = "Clustering rate versus crystal-reference coverage"
    sources = ["combined/condition_table.tsv", "combined/crystal_anchor_table.tsv"]
    try:
        cond = add_plot_categories(renderer.read_tsv(sources[0]))
        anchor = renderer.read_tsv(sources[1])
        require_columns(cond, ["protein_id", "any_valid_cluster"])
        require_columns(anchor, ["protein_id", "crystal_reference_id"])
        cond["any_valid_cluster"] = cond["any_valid_cluster"].astype(str).str.lower().isin(["true", "1", "yes"])
        coverage = (
            anchor.dropna(subset=["crystal_reference_id"])
            .groupby("protein_id", observed=True)["crystal_reference_id"]
            .nunique()
            .rename("n_crystal_refs")
            .reset_index()
        )
        rates = cond.groupby("protein_id", observed=True)["any_valid_cluster"].mean().rename("valid_cluster_rate").reset_index()
        df = rates.merge(coverage, on="protein_id", how="left").fillna({"n_crystal_refs": 0})
        df["has_crystal_reference"] = np.where(df["n_crystal_refs"] > 0, "has crystal ref", "no crystal ref")
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), gridspec_kw={"width_ratios": [1.3, 1]})
        sns.scatterplot(
            data=df,
            x="n_crystal_refs",
            y="valid_cluster_rate",
            hue="has_crystal_reference",
            palette={"has crystal ref": "#3b6ea8", "no crystal ref": "#888888"},
            s=65,
            ax=axes[0],
        )
        axes[0].set_xlabel("Number of crystal references")
        axes[0].set_ylabel("Fraction of conditions with valid cluster")
        axes[0].set_ylim(-0.03, 1.03)
        axes[0].legend(frameon=False)
        sns.boxplot(data=df, x="has_crystal_reference", y="valid_cluster_rate", color="#d9d9d9", ax=axes[1])
        sns.stripplot(data=df, x="has_crystal_reference", y="valid_cluster_rate", color="#222222", alpha=0.45, ax=axes[1])
        axes[1].set_xlabel("")
        axes[1].set_ylabel("")
        axes[1].set_ylim(-0.03, 1.03)
        fig.suptitle(title)
        renderer.save(fig, figure_id, title, sources, "protein-level valid-cluster rate and unique crystal refs", len(df))
    except Exception as exc:
        renderer.skip(figure_id, title, sources, str(exc))


def plot_family_residue_enrichment(renderer: FigureRenderer) -> None:
    figure_id = "family_residue/family_residue_enrichment_heatmap.png"
    title = "Family-aligned residue enrichment"
    source = "postprocess/08_family_residue_enrichment/family_residue_enrichment.tsv"
    try:
        df = renderer.read_tsv(source)
        require_columns(df, ["family_label", "alignment_column", "n_proteins_observed", "mean_condition_residue_contact_score"])
        df = coerce_numeric(df, ["alignment_column", "n_proteins_observed", "mean_condition_residue_contact_score"])
        df = df[df["n_proteins_observed"] >= 2].copy()
        df["row_label"] = df["family_label"].astype(str) + ":" + df["alignment_column"].astype("Int64").astype(str)
        top = df.sort_values("mean_condition_residue_contact_score", ascending=False).head(60)
        pivot = top.pivot_table(
            index="row_label",
            columns="family_label",
            values="mean_condition_residue_contact_score",
            aggfunc="mean",
            observed=True,
        )
        width, height = heatmap_size(len(pivot), len(pivot.columns))
        fig, ax = plt.subplots(figsize=(width, height))
        sns.heatmap(pivot, cmap="mako", ax=ax, cbar_kws={"label": "Mean contact score"})
        ax.set_title(title)
        ax.set_xlabel("Family")
        ax.set_ylabel("Family:alignment column")
        renderer.save(fig, figure_id, title, [source], "n_proteins_observed >= 2; top 60 rows", len(top))
    except Exception as exc:
        renderer.skip(figure_id, title, [source], str(exc))


def plot_cv_performance(renderer: FigureRenderer) -> None:
    figure_id = "predictive/cv_performance_summary.png"
    title = "Exploratory grouped-CV performance summary"
    source_glob = "postprocess/10_predictive/cv_results/*_metrics.tsv"
    try:
        paths = sorted((renderer.results_dir / "postprocess/10_predictive/cv_results").glob("*_metrics.tsv"))
        if not paths:
            raise FileNotFoundError(source_glob)
        frames = []
        for path in paths:
            frame = pd.read_csv(path, sep="\t")
            frame["metric_file"] = path.name
            frames.append(frame)
        df = pd.concat(frames, ignore_index=True)
        require_columns(df, ["model_name", "mean_pr_auc", "std_pr_auc", "mean_balanced_accuracy", "std_balanced_accuracy"])
        df = coerce_numeric(df, ["mean_pr_auc", "std_pr_auc", "mean_balanced_accuracy", "std_balanced_accuracy"])
        long = pd.concat(
            [
                df[["model_name", "mean_pr_auc", "std_pr_auc"]].rename(
                    columns={"mean_pr_auc": "mean", "std_pr_auc": "std"}
                ).assign(metric="PR AUC"),
                df[["model_name", "mean_balanced_accuracy", "std_balanced_accuracy"]].rename(
                    columns={"mean_balanced_accuracy": "mean", "std_balanced_accuracy": "std"}
                ).assign(metric="Balanced accuracy"),
            ],
            ignore_index=True,
        )
        long = finite_df(long, ["mean"])
        fig, ax = plt.subplots(figsize=(9.5, 5.2))
        x_labels = list(long["model_name"].drop_duplicates())
        x_base = np.arange(len(x_labels))
        offsets = {"PR AUC": -0.17, "Balanced accuracy": 0.17}
        colors = {"PR AUC": "#3b6ea8", "Balanced accuracy": "#d17a22"}
        for metric, sub in long.groupby("metric", observed=True):
            xs = [x_base[x_labels.index(name)] + offsets[metric] for name in sub["model_name"]]
            ax.errorbar(
                xs,
                sub["mean"],
                yerr=sub["std"].fillna(0),
                fmt="o",
                color=colors[metric],
                label=metric,
                capsize=3,
            )
        ax.set_xticks(x_base)
        ax.set_xticklabels(x_labels, rotation=35, ha="right")
        ax.set_ylim(-0.03, 1.03)
        ax.set_ylabel("Mean score +/- fold SD")
        ax.set_xlabel("")
        ax.set_title(title)
        ax.legend(frameon=False)
        renderer.save(fig, figure_id, title, [source_glob], "all *_metrics.tsv files", len(df))
    except Exception as exc:
        renderer.skip(figure_id, title, [source_glob], str(exc))


def plot_cbm_paired_delta(renderer: FigureRenderer) -> None:
    figure_id = "cbm_paired/paired_delta_plot_primary_metrics.png"
    title = "CBM paired primary metrics"
    sources = [
        "postprocess/15_cbm_paired_analysis/cbm_paired_comparison_table.tsv",
        "postprocess/15_cbm_paired_analysis/cbm_primary_metric_summary.tsv",
    ]
    try:
        paired = renderer.read_tsv(sources[0])
        summary = renderer.read_tsv(sources[1])
        require_columns(summary, ["source_column", "median", "bootstrap_ci_low", "bootstrap_ci_high"])
        metric_cols = [col for col in summary["source_column"].dropna().astype(str) if col in paired.columns]
        if not metric_cols:
            raise KeyError("no primary metric source columns found in paired table")
        paired = coerce_numeric(paired, metric_cols)
        long = paired.melt(value_vars=metric_cols, var_name="metric", value_name="value")
        long = finite_df(long, ["value"])
        long["metric_label"] = long["metric"].map(nice_metric_name)
        fig, ax = plt.subplots(figsize=(10, 5.4))
        sns.stripplot(data=long, x="metric_label", y="value", color="#222222", alpha=0.22, size=2, ax=ax)
        summ = summary[summary["source_column"].isin(metric_cols)].copy()
        summ = coerce_numeric(summ, ["median", "bootstrap_ci_low", "bootstrap_ci_high"])
        labels = [nice_metric_name(col) for col in metric_cols]
        for idx, col in enumerate(metric_cols):
            row = summ[summ["source_column"] == col].head(1)
            if row.empty or pd.isna(row["median"].iloc[0]):
                continue
            median = row["median"].iloc[0]
            low = row["bootstrap_ci_low"].iloc[0]
            high = row["bootstrap_ci_high"].iloc[0]
            if pd.notna(low) and pd.notna(high):
                ax.errorbar(idx, median, yerr=[[median - low], [high - median]], fmt="D", color="#c44e52", capsize=4)
            else:
                ax.scatter([idx], [median], marker="D", color="#c44e52", zorder=4)
        ax.axhline(0, color="#555555", linewidth=1, linestyle="--")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.set_ylabel("Delta / full-length metric value")
        ax.set_xlabel("")
        ax.set_title(title)
        renderer.save(fig, figure_id, title, sources, "source_column metrics present in paired table", len(long))
    except Exception as exc:
        renderer.skip(figure_id, title, sources, str(exc))


def plot_cbm_domain_full_paired_lines(renderer: FigureRenderer) -> None:
    figure_id = "cbm_paired/domain_only_vs_full_length_paired_lines.png"
    title = "Domain-only versus full-length paired metrics"
    sources = [
        "postprocess/15_cbm_paired_analysis/cbm_construct_condition_summary.tsv",
        "postprocess/15_cbm_paired_analysis/cbm_paired_comparison_table.tsv",
    ]
    try:
        summary = add_plot_categories(renderer.read_tsv(sources[0]))
        paired = renderer.read_tsv(sources[1])
        require_columns(paired, ["domain_only_condition_id", "full_length_condition_id", "protein_id"])
        metrics = [
            "top_cluster_occupancy",
            "cluster_entropy",
            "geometry_plausible_fraction",
            "C4_minus_C1_geometry_bias",
            "catalytic_surface_contact_fraction",
            "aromatic_contact_fraction",
            "polar_contact_fraction",
        ]
        metrics = [metric for metric in metrics if metric in summary.columns]
        require_columns(summary, ["condition_id", "construct_type", *metrics])
        pair_map = paired[["protein_id", "domain_only_condition_id", "full_length_condition_id"]].copy()
        pair_map["pair_id"] = np.arange(len(pair_map))
        long_ids = pd.concat(
            [
                pair_map[["pair_id", "protein_id", "domain_only_condition_id"]].rename(columns={"domain_only_condition_id": "condition_id"}).assign(expected_construct="domain_only"),
                pair_map[["pair_id", "protein_id", "full_length_condition_id"]].rename(columns={"full_length_condition_id": "condition_id"}).assign(expected_construct="full_length"),
            ],
            ignore_index=True,
        )
        matched = long_ids.merge(summary[["condition_id", "construct_type", *metrics]], on="condition_id", how="left")
        value_long = matched.melt(
            id_vars=["pair_id", "protein_id", "condition_id", "expected_construct", "construct_type"],
            value_vars=metrics,
            var_name="metric",
            value_name="value",
        )
        value_long = finite_df(coerce_numeric(value_long, ["value"]), ["value"])
        value_long["metric_label"] = value_long["metric"].map(nice_metric_name)
        n_metrics = min(4, value_long["metric_label"].nunique())
        selected = list(value_long["metric_label"].drop_duplicates())[:n_metrics]
        sub = value_long[value_long["metric_label"].isin(selected)].copy()
        fig, axes = plt.subplots(1, n_metrics, figsize=(3.2 * n_metrics, 4.7), sharex=True)
        if n_metrics == 1:
            axes = [axes]
        xmap = {"domain_only": 0, "full_length": 1}
        for ax, metric in zip(axes, selected, strict=False):
            metric_df = sub[sub["metric_label"] == metric]
            for _, pair in metric_df.groupby("pair_id", observed=True):
                if pair["expected_construct"].nunique() < 2:
                    continue
                pair = pair.sort_values("expected_construct", key=lambda s: s.map(xmap))
                ax.plot(pair["expected_construct"].map(xmap), pair["value"], color="#999999", alpha=0.18, linewidth=0.8)
            sns.pointplot(
                data=metric_df,
                x="expected_construct",
                y="value",
                order=CONSTRUCT_ORDER,
                errorbar=("pi", 50),
                color="#222222",
                ax=ax,
            )
            ax.set_title(metric)
            ax.set_xlabel("")
            ax.set_ylabel("Value")
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["domain", "full"], rotation=0)
        fig.suptitle(title)
        renderer.save(fig, figure_id, title, sources, "matched condition IDs; first four available metrics", len(sub))
    except Exception as exc:
        renderer.skip(figure_id, title, sources, str(exc))


def plot_cbm_bridge_fraction(renderer: FigureRenderer) -> None:
    figure_id = "cbm_paired/cbm_bridge_fraction_by_substrate_dp.png"
    title = "Full-length CBM bridge fraction by substrate and DP"
    source = "postprocess/15_cbm_paired_analysis/cbm_construct_condition_summary.tsv"
    try:
        df = add_plot_categories(renderer.read_tsv(source))
        require_columns(df, ["construct_type", "substrate_class", "dp", "bridge_fraction"])
        df = df[df["construct_type"].astype(str) == "full_length"].copy()
        df = finite_df(coerce_numeric(df, ["bridge_fraction"]), ["bridge_fraction"])
        fig, ax = plt.subplots(figsize=(9.5, 4.8))
        sns.pointplot(
            data=df,
            x="substrate_dp",
            y="bridge_fraction",
            hue="substrate_class",
            palette=substrate_palette(df),
            errorbar=("pi", 50),
            ax=ax,
        )
        sns.stripplot(data=df, x="substrate_dp", y="bridge_fraction", color="#222222", alpha=0.22, size=2, ax=ax)
        ax.set_ylim(-0.03, 1.03)
        ax.set_xlabel("")
        ax.set_ylabel("Bridge fraction")
        ax.set_title(title)
        ax.tick_params(axis="x", rotation=45)
        ax.legend(frameon=False, title="Substrate")
        renderer.save(fig, figure_id, title, [source], "full_length rows with non-missing bridge_fraction", len(df))
    except Exception as exc:
        renderer.skip(figure_id, title, [source], str(exc))


def plot_cbm_active_site_ifp_change(renderer: FigureRenderer) -> None:
    figure_id = "cbm_paired/active_site_ifp_change_heatmap.png"
    title = "Active-site IFP change in CBM pairs"
    source = "postprocess/15_cbm_paired_analysis/cbm_paired_comparison_table.tsv"
    try:
        df = add_plot_categories(renderer.read_tsv(source))
        require_columns(df, ["protein_id", "substrate_class", "dp", "catalytic_domain_ifp_jaccard_distance"])
        df = coerce_numeric(df, ["catalytic_domain_ifp_jaccard_distance"])
        df["column"] = df["substrate_class"].astype(str) + " DP" + df["dp"].astype(str)
        pivot = df.pivot_table(
            index="protein_id",
            columns="column",
            values="catalytic_domain_ifp_jaccard_distance",
            aggfunc="mean",
            observed=True,
        )
        pivot = pivot.dropna(how="all")
        if pivot.empty:
            raise ValueError("no non-missing catalytic_domain_ifp_jaccard_distance values")
        width, height = heatmap_size(len(pivot), len(pivot.columns))
        fig, ax = plt.subplots(figsize=(width, height))
        sns.heatmap(pivot, cmap="rocket_r", vmin=0, vmax=1, ax=ax, cbar_kws={"label": "Jaccard distance"})
        ax.set_title(title)
        ax.set_xlabel("Substrate/DP")
        ax.set_ylabel("Protein")
        renderer.save(fig, figure_id, title, [source], "mean per protein/substrate/DP; missing kept blank", int(df["catalytic_domain_ifp_jaccard_distance"].notna().sum()))
    except Exception as exc:
        renderer.skip(figure_id, title, [source], str(exc))


def plot_cbm_geometry_bias_delta(renderer: FigureRenderer) -> None:
    figure_id = "cbm_paired/geometry_bias_delta_heatmap.png"
    title = "CBM delta in C4 minus C1 geometry bias"
    source = "postprocess/15_cbm_paired_analysis/cbm_paired_comparison_table.tsv"
    try:
        df = add_plot_categories(renderer.read_tsv(source))
        require_columns(df, ["protein_id", "substrate_class", "dp", "delta_C4_minus_C1_geometry_bias"])
        df = coerce_numeric(df, ["delta_C4_minus_C1_geometry_bias"])
        df["column"] = df["substrate_class"].astype(str) + " DP" + df["dp"].astype(str)
        pivot = df.pivot_table(
            index="protein_id",
            columns="column",
            values="delta_C4_minus_C1_geometry_bias",
            aggfunc="mean",
            observed=True,
        )
        order = pivot.mean(axis=1).sort_values(ascending=False).index
        pivot = pivot.loc[order]
        max_abs = float(np.nanmax(np.abs(pivot.to_numpy()))) if pivot.size else 1.0
        if not math.isfinite(max_abs) or max_abs == 0:
            max_abs = 1.0
        width, height = heatmap_size(len(pivot), len(pivot.columns))
        fig, ax = plt.subplots(figsize=(width, height))
        sns.heatmap(
            pivot,
            cmap="vlag",
            center=0,
            vmin=-max_abs,
            vmax=max_abs,
            ax=ax,
            cbar_kws={"label": "Delta C4-C1 geometry bias"},
        )
        ax.set_title(title)
        ax.set_xlabel("Substrate/DP")
        ax.set_ylabel("Protein")
        renderer.save(fig, figure_id, title, [source], "mean per protein/substrate/DP; sorted by mean delta", int(df["delta_C4_minus_C1_geometry_bias"].notna().sum()))
    except Exception as exc:
        renderer.skip(figure_id, title, [source], str(exc))


FIGURE_FUNCTIONS: list[Callable[[FigureRenderer], None]] = [
    plot_qc_attrition,
    plot_cluster_count_distribution,
    plot_top_cluster_occupancy,
    plot_c1_c4_plausibility,
    plot_cu_c1_vs_cu_c4,
    plot_top_residue_contact_heatmap,
    plot_contact_chemistry,
    plot_pocket_rmsd_vs_ifp,
    plot_clustering_rate_vs_crystal_coverage,
    plot_family_residue_enrichment,
    plot_cv_performance,
    plot_cbm_paired_delta,
    plot_cbm_domain_full_paired_lines,
    plot_cbm_bridge_fraction,
    plot_cbm_active_site_ifp_change,
    plot_cbm_geometry_bias_delta,
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results/2full_pipeline_array_rerun_downstream_20260528_164450"),
        help="Root directory containing combined/ and postprocess/ result tables.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Figure output directory. Defaults to <results-dir>/postprocess/figures.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Return a non-zero exit code if any figure is skipped.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_style()
    renderer = FigureRenderer(args.results_dir, args.output_dir)
    for plot_func in FIGURE_FUNCTIONS:
        plot_func(renderer)
    renderer.manifest()
    ok = sum(record.status == "ok" for record in renderer.records)
    skipped = sum(record.status != "ok" for record in renderer.records)
    print(f"Rendered {ok} figures; skipped {skipped}.")
    print(f"Manifest: {renderer.output_dir / 'figure_manifest.tsv'}")
    return 2 if args.strict and skipped else 0


if __name__ == "__main__":
    raise SystemExit(main())
