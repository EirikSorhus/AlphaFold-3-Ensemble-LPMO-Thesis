#!/usr/bin/env python3
"""Render the revised figure2 package from frozen LPMO result tables."""

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
from scipy import stats


SUBSTRATE_ORDER = ["chitin", "cellulose", "starch"]
CONSTRUCT_ORDER = ["domain_only", "full_length"]
CONSTRUCT_PALETTE = {"domain_only": "#3b6ea8", "full_length": "#d17a22"}
SUBSTRATE_PALETTE = {"chitin": "#2b8cbe", "cellulose": "#41ab5d", "starch": "#e34a33", "amylose": "#e34a33"}
RELATION_PALETTE = {"right": "#2f7f4f", "wrong": "#9b3a3a"}
REGIO_PALETTE = {"C1": "#3b6ea8", "C4": "#d17a22"}


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


class Renderer:
    def __init__(self, results_dir: Path, output_dir: Path, analyse_dir: Path) -> None:
        self.results_dir = results_dir.resolve()
        self.output_dir = output_dir.resolve()
        self.analyse_dir = analyse_dir.resolve()
        self.records: list[FigureRecord] = []
        self.tables: dict[Path, pd.DataFrame] = {}
        self.output_dir.mkdir(parents=True, exist_ok=True)
        for sub in ["cbm_paired", "crystal", "descriptive", "predictive", "family_residue", "tables"]:
            (self.output_dir / sub).mkdir(parents=True, exist_ok=True)

    def rel(self, path: Path) -> str:
        try:
            return str(path.resolve().relative_to(self.results_dir))
        except ValueError:
            return str(path)

    def read_tsv(self, rel_path: str, *, base: str = "results") -> pd.DataFrame:
        root = self.results_dir if base == "results" else self.analyse_dir
        path = root / rel_path
        if path in self.tables:
            return self.tables[path].copy()
        if not path.is_file():
            raise FileNotFoundError(path)
        df = pd.read_csv(path, sep="\t", low_memory=False)
        self.tables[path] = df
        return df.copy()

    def save(self, fig: plt.Figure, figure_id: str, title: str, sources: list[str], filters: str, n: int) -> None:
        out = self.output_dir / figure_id
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=220, bbox_inches="tight")
        plt.close(fig)
        self.records.append(FigureRecord(figure_id, title, ";".join(sources), filters, int(n), self.rel(out)))

    def skip(self, figure_id: str, title: str, sources: list[str], message: str) -> None:
        out = self.output_dir / figure_id
        if out.exists():
            out.unlink()
        self.records.append(
            FigureRecord(figure_id, title, ";".join(sources), "not_rendered", 0, self.rel(out), "skipped", message)
        )

    def write_table(self, df: pd.DataFrame, rel_path: str) -> None:
        out = self.output_dir / rel_path
        out.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, sep="\t", index=False)

    def manifest(self) -> None:
        pd.DataFrame([r.__dict__ for r in self.records]).to_csv(
            self.output_dir / "figure_manifest.tsv", sep="\t", index=False
        )


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


def require(df: pd.DataFrame, cols: list[str]) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError("missing columns: " + ", ".join(missing))


def num(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    df = df.copy()
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def bool_col(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin(["true", "1", "yes", "ja"])


def add_categories(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "construct_type" not in df.columns and "_construct_type" in df.columns:
        df["construct_type"] = df["_construct_type"]
    if "construct_type" in df.columns:
        known = [x for x in CONSTRUCT_ORDER if x in set(df["construct_type"].dropna().astype(str))]
        extra = sorted(set(df["construct_type"].dropna().astype(str)) - set(known))
        df["construct_type"] = pd.Categorical(df["construct_type"].astype(str), known + extra, ordered=True)
    if "substrate_class" in df.columns:
        known = [x for x in SUBSTRATE_ORDER if x in set(df["substrate_class"].dropna().astype(str))]
        extra = sorted(set(df["substrate_class"].dropna().astype(str)) - set(known))
        df["substrate_class"] = pd.Categorical(df["substrate_class"].astype(str), known + extra, ordered=True)
    if "dp" in df.columns:
        df["dp"] = pd.to_numeric(df["dp"], errors="coerce")
    if {"substrate_class", "dp"}.issubset(df.columns):
        df["substrate_dp"] = df["substrate_class"].astype(str) + " DP" + df["dp"].astype("Int64").astype(str)
    return df


def nice(name: str) -> str:
    name = re.sub(r"^(delta_|mean_|median_)", "", name)
    name = name.replace("occupancy_weighted_", "").replace("_fraction", "").replace("_", " ")
    return name


def finite(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        values = pd.to_numeric(out[col], errors="coerce")
        out = out[values.notna() & np.isfinite(values)]
    return out


def heatmap_size(n_rows: int, n_cols: int, min_w: float = 6.0, min_h: float = 4.0) -> tuple[float, float]:
    return max(min_w, min(18.0, 0.42 * n_cols + 3.2)), max(min_h, min(22.0, 0.22 * n_rows + 2.8))


def annotate_counts(ax: plt.Axes, text: str) -> None:
    ax.text(1.01, 0.98, text, transform=ax.transAxes, ha="left", va="top", fontsize=8)


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


def ec_active_substrates(ec_raw: object) -> set[str]:
    raw = str(ec_raw)
    out: set[str] = set()
    if "1.14.99.53" in raw:
        out.add("chitin")
    if "1.14.99.54" in raw or "1.14.99.56" in raw:
        out.add("cellulose")
    if "1.14.99.55" in raw:
        out.add("starch")
    return out


def ec_tokens(ec_raw: object) -> list[str]:
    return re.findall(r"1\.14\.99\.[0-9-]+", str(ec_raw or ""))


def ec_regio_label(ec_raw: object, family_raw: object) -> str:
    family = str(family_raw or "").upper()
    regio: set[str] = set()
    for token in ec_tokens(ec_raw):
        if token in {"1.14.99.53", "1.14.99.54", "1.14.99.55"}:
            regio.add("C1")
        elif token == "1.14.99.56":
            regio.add("C4")
        elif token == "1.14.99.-" and "AA10" in family:
            regio.add("C4")
        elif token == "1.14.99.-" and "AA14" in family:
            regio.add("C1")
        elif token == "1.14.99.-" and "AA17" in family:
            regio.add("C4")
    if regio == {"C1"}:
        return "C1"
    if regio == {"C4"}:
        return "C4"
    if regio == {"C1", "C4"}:
        return "C1+C4"
    return ""


def activity_table(renderer: Renderer) -> pd.DataFrame:
    rows = []
    meta = renderer.read_tsv("input_data/metadata_final_ec_fixed.tsv", base="analyse")
    require(meta, ["UniProt_ID", "CAZy_family", "EC_Number"])
    for _, row in meta.iterrows():
        active = ec_active_substrates(row.get("EC_Number"))
        rows.append(
            {
                "protein_id": row.get("UniProt_ID"),
                "active_substrates": ",".join(sorted(active)),
                "regio": ec_regio_label(row.get("EC_Number"), row.get("CAZy_family")),
            }
        )
    return pd.DataFrame(rows).dropna(subset=["protein_id"]).drop_duplicates("protein_id")


def add_right_wrong(renderer: Renderer, df: pd.DataFrame) -> pd.DataFrame:
    out = df.merge(activity_table(renderer), on="protein_id", how="left")
    out["active_substrate_set"] = out["active_substrates"].fillna("").map(lambda x: set(filter(None, str(x).split(","))))
    out["ligand_relation"] = [
        "right" if str(s) in active else ("wrong" if active else "unknown")
        for s, active in zip(out["substrate_class"].astype(str), out["active_substrate_set"], strict=False)
    ]
    return out[out["ligand_relation"].isin(["right", "wrong"])].copy()


def paired_stats(df: pd.DataFrame, metric_cols: list[str], group_cols: list[str] | None = None) -> pd.DataFrame:
    group_cols = group_cols or []
    rows = []
    for keys, group in df.groupby(group_cols, observed=True) if group_cols else [((), df)]:
        key_values = keys if isinstance(keys, tuple) else (keys,)
        for metric in metric_cols:
            if metric not in group.columns:
                continue
            sub = group[["protein_id", "ligand_relation", metric]].dropna()
            sub[metric] = pd.to_numeric(sub[metric], errors="coerce")
            pivot = sub.groupby(["protein_id", "ligand_relation"], observed=True)[metric].mean().unstack()
            if not {"right", "wrong"}.issubset(pivot.columns):
                continue
            paired = pivot.dropna(subset=["right", "wrong"])
            delta = paired["right"] - paired["wrong"]
            p_value = np.nan
            test = "not_run"
            if len(delta) >= 8 and np.count_nonzero(delta) > 0:
                try:
                    p_value = stats.wilcoxon(paired["right"], paired["wrong"], zero_method="wilcox").pvalue
                    test = "paired_wilcoxon"
                except Exception:
                    p_value = np.nan
                    test = "paired_wilcoxon_failed"
            elif len(delta) >= 3:
                try:
                    p_value = stats.binomtest(int((delta > 0).sum()), int((delta != 0).sum()), 0.5).pvalue
                    test = "paired_sign_test"
                except Exception:
                    p_value = np.nan
                    test = "paired_sign_test_failed"
            row = {
                "metric": metric,
                "n_pairs": len(paired),
                "median_right": paired["right"].median(),
                "median_wrong": paired["wrong"].median(),
                "median_delta": delta.median(),
                "n_positive_delta": int((delta > 0).sum()),
                "n_negative_delta": int((delta < 0).sum()),
                "test": test,
                "p_value": p_value,
            }
            row.update(dict(zip(group_cols, key_values, strict=False)))
            rows.append(row)
    return pd.DataFrame(rows)


def plot_cbm_qc_clusterable(renderer: Renderer) -> None:
    fid = "cbm_paired/qc_clusterable_proportions_core_vs_full_length.png"
    title = "CBM paired hard-QC and clusterable proportions"
    sources = [
        "postprocess/15_cbm_paired_analysis/cbm_construct_condition_summary.tsv",
        "postprocess/15_cbm_paired_analysis/cbm_paired_comparison_table.tsv",
    ]
    try:
        summary = add_categories(renderer.read_tsv(sources[0]))
        paired = renderer.read_tsv(sources[1])
        require(summary, ["condition_id", "construct_type", "n_generated", "n_qc_pass", "qc_pass_fraction", "any_valid_cluster"])
        require(paired, ["domain_only_condition_id", "full_length_condition_id"])
        ids = pd.concat(
            [
                paired[["domain_only_condition_id"]].rename(columns={"domain_only_condition_id": "condition_id"}).assign(pair_id=paired.index),
                paired[["full_length_condition_id"]].rename(columns={"full_length_condition_id": "condition_id"}).assign(pair_id=paired.index),
            ],
            ignore_index=True,
        )
        df = ids.merge(summary, on="condition_id", how="left")
        df["any_valid_cluster"] = bool_col(df["any_valid_cluster"])
        df = num(df, ["n_generated", "n_qc_pass", "qc_pass_fraction"])
        df["metric_qc"] = df["qc_pass_fraction"]
        df["metric_clusterable"] = df["any_valid_cluster"].astype(float)
        long = df.melt(
            id_vars=["pair_id", "construct_type", "n_generated", "n_qc_pass"],
            value_vars=["metric_qc", "metric_clusterable"],
            var_name="metric",
            value_name="proportion",
        )
        long["metric"] = long["metric"].map({"metric_qc": "Hard-QC pass", "metric_clusterable": "Valid cluster"})
        fig, axes = plt.subplots(1, 2, figsize=(8.5, 4.8), sharey=True)
        for ax, metric in zip(axes, ["Hard-QC pass", "Valid cluster"], strict=False):
            sub = long[long["metric"] == metric]
            for _, pair in sub.groupby("pair_id", observed=True):
                if pair["construct_type"].nunique() == 2:
                    pair = pair.sort_values("construct_type", key=lambda s: s.map({"domain_only": 0, "full_length": 1}))
                    ax.plot(pair["construct_type"], pair["proportion"], color="#999999", alpha=0.22, linewidth=0.8)
            sns.pointplot(data=sub, x="construct_type", y="proportion", order=CONSTRUCT_ORDER, color="#222222", errorbar=("pi", 50), ax=ax)
            ax.set_title(metric)
            ax.set_xlabel("")
            ax.set_ylabel("Proportion")
            ax.set_ylim(-0.03, 1.03)
            n_pairs = paired.shape[0]
            n_gen = int(df.loc[df["construct_type"].astype(str).eq("domain_only"), "n_generated"].sum())
            n_gen_fl = int(df.loc[df["construct_type"].astype(str).eq("full_length"), "n_generated"].sum())
            annotate_counts(ax, f"n pairs={n_pairs}\nposes domain={n_gen}\nposes full={n_gen_fl}")
        fig.suptitle(title)
        renderer.save(fig, fid, title, sources, "matched CBM condition IDs", len(long))
    except Exception as exc:
        renderer.skip(fid, title, sources, str(exc))


def plot_cbm_clustered_matrix(renderer: Renderer) -> None:
    fid = "cbm_paired/clustered_domain_only_vs_full_length_paired_matrix.png"
    title = "Clustered CBM pairs: domain-only versus full-length"
    sources = [
        "postprocess/15_cbm_paired_analysis/cbm_paired_comparison_table.tsv",
        "postprocess/15_cbm_paired_analysis/cbm_construct_condition_summary.tsv",
    ]
    try:
        paired = renderer.read_tsv(sources[0])
        summary = renderer.read_tsv(sources[1])
        paired = paired[bool_col(paired["any_valid_cluster_domain_only"]) & bool_col(paired["any_valid_cluster_full_length"])].copy()
        metrics = ["top_cluster_occupancy", "cluster_entropy", "occupancy_gini", "n_clusters"]
        require(summary, ["condition_id", "construct_type", *metrics])
        rows = []
        for _, row in paired.iterrows():
            for construct, cid in [("domain_only", row["domain_only_condition_id"]), ("full_length", row["full_length_condition_id"])]:
                hit = summary[summary["condition_id"].eq(cid)]
                if hit.empty:
                    continue
                values = hit.iloc[0][metrics].to_dict()
                values.update({"row_label": f"{row['protein_id']} {row['substrate_class']} DP{row['dp']}", "construct": construct})
                rows.append(values)
        df = pd.DataFrame(rows)
        long = df.melt(id_vars=["row_label", "construct"], value_vars=metrics, var_name="metric", value_name="value")
        long["column"] = long["metric"].map(nice) + " | " + long["construct"].astype(str)
        pivot = long.pivot_table(index="row_label", columns="column", values="value", aggfunc="mean")
        pivot = pivot.dropna(how="all")
        width, height = heatmap_size(len(pivot), len(pivot.columns), 9, 5)
        fig, ax = plt.subplots(figsize=(width, height))
        sns.heatmap(pivot, cmap="viridis", ax=ax, cbar_kws={"label": "Metric value"})
        ax.set_title(title)
        ax.set_xlabel("")
        ax.set_ylabel("Protein/substrate/DP")
        renderer.save(fig, fid, title, sources, "pairs where both constructs have valid clusters", len(pivot))
    except Exception as exc:
        renderer.skip(fid, title, sources, str(exc))


def plot_cbm_delta_bias(renderer: Renderer) -> None:
    fid = "cbm_paired/cbm_delta_c4_c1_geometry_bias_clustered.png"
    title = "CBM delta C4-C1 geometry bias, clustered pairs"
    source = "postprocess/15_cbm_paired_analysis/cbm_paired_comparison_table.tsv"
    try:
        df = add_categories(renderer.read_tsv(source))
        require(df, ["protein_id", "substrate_class", "dp", "delta_C4_minus_C1_geometry_bias"])
        if "any_valid_cluster_domain_only" in df.columns and "any_valid_cluster_full_length" in df.columns:
            df = df[bool_col(df["any_valid_cluster_domain_only"]) & bool_col(df["any_valid_cluster_full_length"])].copy()
        df = num(df, ["delta_C4_minus_C1_geometry_bias"])
        df["column"] = df["substrate_class"].astype(str) + " DP" + df["dp"].astype("Int64").astype(str)
        pivot = df.pivot_table(index="protein_id", columns="column", values="delta_C4_minus_C1_geometry_bias", aggfunc="mean")
        pivot = pivot.loc[pivot.mean(axis=1).sort_values(ascending=False).index]
        max_abs = np.nanmax(np.abs(pivot.to_numpy())) if pivot.size else 1
        max_abs = 1 if not math.isfinite(max_abs) or max_abs == 0 else max_abs
        fig, axes = plt.subplots(1, 2, figsize=(13, max(4.5, min(12, 0.22 * len(pivot) + 2.5))), gridspec_kw={"width_ratios": [2.2, 1]})
        sns.heatmap(pivot, cmap="vlag", center=0, vmin=-max_abs, vmax=max_abs, ax=axes[0], cbar_kws={"label": "Delta C4-C1 bias"})
        axes[0].set_title("Per clustered pair")
        axes[0].set_xlabel("Substrate/DP")
        axes[0].set_ylabel("Protein")
        sns.violinplot(data=df, y="delta_C4_minus_C1_geometry_bias", color="#d9d9d9", cut=0, ax=axes[1])
        sns.stripplot(data=df, y="delta_C4_minus_C1_geometry_bias", color="#222222", alpha=0.45, ax=axes[1])
        axes[1].axhline(0, color="#555555", linestyle="--", linewidth=1)
        axes[1].set_title("Delta distribution")
        axes[1].set_ylabel("Full-length - domain-only")
        annotate_counts(axes[1], f"n pairs={len(df)}\nmedian={df['delta_C4_minus_C1_geometry_bias'].median():.3g}")
        fig.suptitle(title)
        renderer.save(fig, fid, title, [source], "clustered pairs only", len(df))
    except Exception as exc:
        renderer.skip(fid, title, [source], str(exc))


def plot_pocket_rmsd_by_ligand(renderer: Renderer) -> None:
    fid = "crystal/pocket_rmsd_by_ligand_type.png"
    title = "Local pocket RMSD by ligand type"
    sources = ["combined/crystal_anchor_table.tsv", "combined/condition_table.tsv"]
    try:
        anchor = renderer.read_tsv(sources[0])
        cond = add_categories(renderer.read_tsv(sources[1]))
        require(anchor, ["condition_id", "local_pocket_rmsd"])
        df = anchor.merge(cond[["condition_id", "substrate_class", "dp", "construct_type"]].drop_duplicates("condition_id"), on="condition_id", how="left")
        df = add_categories(finite(num(df, ["local_pocket_rmsd"]), ["local_pocket_rmsd"]))
        fig, ax = plt.subplots(figsize=(9.5, 5))
        sns.boxplot(data=df, x="substrate_class", y="local_pocket_rmsd", hue="construct_type", palette=CONSTRUCT_PALETTE, showfliers=False, ax=ax)
        sns.stripplot(data=df, x="substrate_class", y="local_pocket_rmsd", hue="construct_type", dodge=True, color="#222222", alpha=0.35, size=3, ax=ax, legend=False)
        ax.set_title(title)
        ax.set_xlabel("Ligand/substrate type")
        ax.set_ylabel("Local pocket RMSD")
        ax.legend(title="Construct", frameon=False)
        annotate_counts(ax, f"n rows={len(df)}\nn proteins={df['protein_id'].nunique() if 'protein_id' in df.columns else 'NA'}")
        renderer.save(fig, fid, title, sources, "non-missing local_pocket_rmsd", len(df))
    except Exception as exc:
        renderer.skip(fid, title, sources, str(exc))


def plot_holo_ifp_vs_rmsd(renderer: Renderer) -> None:
    fid = "crystal/holo_ifp_similarity_vs_local_pocket_rmsd.png"
    title = "Holo-only IFP similarity versus local pocket RMSD"
    sources = ["combined/crystal_anchor_table.tsv", "combined/condition_table.tsv"]
    try:
        anchor = renderer.read_tsv(sources[0])
        cond = renderer.read_tsv(sources[1])
        require(anchor, ["condition_id", "local_pocket_rmsd", "ifp_tanimoto", "ifp_comparison_eligible"])
        df = anchor.merge(cond[["condition_id", "substrate_class", "dp", "construct_type"]].drop_duplicates("condition_id"), on="condition_id", how="left")
        df = df[bool_col(df["ifp_comparison_eligible"])].copy()
        if "crystal_ifp_contact_eligible" in df.columns:
            df = df[bool_col(df["crystal_ifp_contact_eligible"])].copy()
        df = add_categories(finite(num(df, ["local_pocket_rmsd", "ifp_tanimoto"]), ["local_pocket_rmsd", "ifp_tanimoto"]))
        fig, ax = plt.subplots(figsize=(7, 5.5))
        sns.scatterplot(data=df, x="local_pocket_rmsd", y="ifp_tanimoto", hue="substrate_class", style="construct_type", palette=SUBSTRATE_PALETTE, s=70, ax=ax)
        if len(df) >= 3:
            sns.regplot(data=df, x="local_pocket_rmsd", y="ifp_tanimoto", scatter=False, lowess=True, color="#222222", ax=ax)
        ax.set_title(title)
        ax.set_xlabel("Local pocket RMSD")
        ax.set_ylabel("IFP Tanimoto")
        ax.set_ylim(-0.03, 1.03)
        ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
        renderer.save(fig, fid, title, sources, "ifp_comparison_eligible and crystal_ifp_contact_eligible rows", len(df))
    except Exception as exc:
        renderer.skip(fid, title, sources, str(exc))


def plot_his_brace_quality(renderer: Renderer) -> None:
    fid = "crystal/his_brace_quality_by_ligand_protein.png"
    title = "His-brace quality by protein and ligand type"
    sources = ["combined/pose_geometry.tsv", "combined/condition_table.tsv"]
    try:
        geom = renderer.read_tsv(sources[0])
        cond = renderer.read_tsv(sources[1])
        require(geom, ["protein_id", "ligand_id", "brace_integrity_flag", "his_brace_angle_deg"])
        geom = num(geom, ["his_brace_angle_deg"])
        # pose_geometry has no condition_id; ligand_id is enough to add substrate metadata in this run.
        meta = cond[["ligand_id", "substrate_class"]].drop_duplicates("ligand_id")
        df = geom.merge(meta, on="ligand_id", how="left")
        df["quality"] = np.select(
            [
                df["brace_integrity_flag"].astype(str).str.lower().isin(["true", "ok", "1", "yes"]),
                df["his_brace_angle_deg"].notna() & (df["his_brace_angle_deg"] <= 35),
                df["his_brace_angle_deg"].notna() & (df["his_brace_angle_deg"] <= 55),
            ],
            ["flag_ok", "angle_good", "angle_marginal"],
            default="poor_or_missing",
        )
        counts = df.groupby(["protein_id", "substrate_class", "quality"], observed=True).size().rename("n").reset_index()
        counts["total"] = counts.groupby(["protein_id", "substrate_class"], observed=True)["n"].transform("sum")
        counts["fraction"] = counts["n"] / counts["total"]
        # Show the fraction of non-poor poses; keep detailed counts in table.
        renderer.write_table(counts, "tables/his_brace_quality_counts.tsv")
        good = counts[counts["quality"].isin(["flag_ok", "angle_good"])].groupby(["protein_id", "substrate_class"], observed=True)["fraction"].sum().reset_index()
        pivot = good.pivot_table(index="protein_id", columns="substrate_class", values="fraction", aggfunc="mean")
        width, height = heatmap_size(len(pivot), len(pivot.columns))
        fig, ax = plt.subplots(figsize=(width, height))
        sns.heatmap(pivot, cmap="YlGnBu", vmin=0, vmax=1, ax=ax, cbar_kws={"label": "Fraction flag_ok/angle_good"})
        ax.set_title(title)
        ax.set_xlabel("Ligand/substrate")
        ax.set_ylabel("Protein")
        renderer.save(fig, fid, title, sources, "pose-level brace flag/angle categories", len(df))
    except Exception as exc:
        renderer.skip(fid, title, sources, str(exc))


def plot_iptm_vs_qc(renderer: Renderer) -> None:
    fid = "crystal/iptm_vs_qc_clusterable.png"
    title = "ipTM versus QC and clusterability"
    source = "combined/condition_table.tsv"
    try:
        df = add_categories(renderer.read_tsv(source))
        require(df, ["mean_iptm", "n_stage1_pass", "n_generated", "any_valid_cluster", "top_cluster_occupancy"])
        df = num(df, ["mean_iptm", "n_stage1_pass", "n_generated", "top_cluster_occupancy"])
        df["qc_pass_fraction"] = df["n_stage1_pass"] / df["n_generated"].replace(0, np.nan)
        df["valid_cluster"] = np.where(bool_col(df["any_valid_cluster"]), "valid cluster", "no valid cluster")
        df = finite(df, ["mean_iptm", "qc_pass_fraction"])
        if df.empty:
            raise ValueError("no non-missing mean_iptm rows available")
        fig, axes = plt.subplots(1, 2, figsize=(11.5, 5), sharex=True)
        sns.scatterplot(data=df, x="mean_iptm", y="qc_pass_fraction", hue="valid_cluster", style="construct_type", s=55, ax=axes[0])
        sns.scatterplot(data=df, x="mean_iptm", y="top_cluster_occupancy", hue="valid_cluster", style="construct_type", s=55, ax=axes[1], legend=False)
        for ax, ylabel in zip(axes, ["Hard-QC pass fraction", "Top-cluster occupancy"], strict=False):
            ax.set_xlabel("Mean ipTM")
            ax.set_ylabel(ylabel)
            ax.set_ylim(-0.03, 1.03)
        axes[0].legend(frameon=False)
        fig.suptitle(title)
        renderer.save(fig, fid, title, [source], "non-missing mean_iptm and QC fraction", len(df))
    except Exception as exc:
        renderer.skip(fid, title, [source], str(exc))


def plot_right_wrong_c1c4(renderer: Renderer) -> None:
    fid = "descriptive/right_wrong_ligand_c1_c4_plausible_fraction.png"
    title = "C1/C4 plausible fraction for right versus wrong ligand"
    source = "combined/condition_table.tsv"
    try:
        df = add_categories(add_right_wrong(renderer, renderer.read_tsv(source)))
        metrics = ["occupancy_weighted_c1_plausible_fraction", "occupancy_weighted_c4_plausible_fraction"]
        require(df, ["protein_id", "ligand_relation", *metrics])
        df = num(df, metrics)
        long = df.melt(id_vars=["protein_id", "ligand_relation", "substrate_class", "dp"], value_vars=metrics, var_name="site", value_name="plausible_fraction")
        long["site"] = long["site"].map({metrics[0]: "C1", metrics[1]: "C4"})
        long = finite(long, ["plausible_fraction"])
        fig, axes = plt.subplots(1, 2, figsize=(8.5, 5), sharey=True)
        for ax, site in zip(axes, ["C1", "C4"], strict=False):
            sub = long[long["site"] == site]
            sns.boxplot(data=sub, x="ligand_relation", y="plausible_fraction", order=["right", "wrong"], palette=RELATION_PALETTE, showfliers=False, ax=ax)
            sns.stripplot(data=sub, x="ligand_relation", y="plausible_fraction", order=["right", "wrong"], color="#222222", alpha=0.25, size=2, ax=ax)
            ax.set_title(site)
            ax.set_xlabel("")
            ax.set_ylabel("Occupancy-weighted plausible fraction")
            ax.set_ylim(-0.03, 1.03)
            annotate_counts(ax, f"n proteins={sub['protein_id'].nunique()}\nn rows={len(sub)}")
        stats_df = paired_stats(df, metrics)
        renderer.write_table(stats_df, "tables/right_wrong_c1_c4_nonparametric_stats.tsv")
        fig.suptitle(title)
        renderer.save(fig, fid, title, [source, "input_data/metadata_final_ec_fixed.tsv"], "right/wrong ligand relation from EC activity metadata", len(long))
    except Exception as exc:
        renderer.skip(fid, title, [source], str(exc))


def plot_right_wrong_angles_distances(renderer: Renderer) -> None:
    fid = "descriptive/right_wrong_ligand_attack_angle_distance.png"
    title = "Attack angle and Cu distance for right versus wrong ligand"
    source = "combined/condition_table.tsv"
    try:
        df = add_categories(add_right_wrong(renderer, renderer.read_tsv(source)))
        metrics = [
            "occupancy_weighted_Cu_oxyl_H_C1_angle_median",
            "occupancy_weighted_Cu_oxyl_H_C4_angle_median",
            "occupancy_weighted_Cu_C1_distance_median",
            "occupancy_weighted_Cu_C4_distance_median",
        ]
        require(df, ["ligand_relation", *metrics])
        df = num(df, metrics)
        long = df.melt(id_vars=["protein_id", "ligand_relation"], value_vars=metrics, var_name="metric", value_name="value")
        long["metric_label"] = long["metric"].map(nice)
        long = finite(long, ["value"])
        g = sns.catplot(data=long, x="ligand_relation", y="value", col="metric_label", col_wrap=2, order=["right", "wrong"], kind="box", showfliers=False, height=3.5, aspect=1.25, palette=RELATION_PALETTE)
        for ax in g.axes.flat:
            metric = ax.get_title().replace("metric_label = ", "")
            sub = long[long["metric_label"].eq(metric)]
            sns.stripplot(data=sub, x="ligand_relation", y="value", order=["right", "wrong"], color="#222222", alpha=0.22, size=2, ax=ax)
            ax.set_xlabel("")
        g.set_axis_labels("", "Value")
        g.fig.suptitle(title, y=1.04)
        stats_df = paired_stats(df, metrics)
        existing = renderer.output_dir / "tables/right_wrong_c1_c4_nonparametric_stats.tsv"
        if existing.is_file():
            old = pd.read_csv(existing, sep="\t")
            stats_df = pd.concat([old, stats_df], ignore_index=True)
        renderer.write_table(stats_df, "tables/right_wrong_c1_c4_nonparametric_stats.tsv")
        renderer.save(g.fig, fid, title, [source, "input_data/metadata_final_ec_fixed.tsv"], "right/wrong ligand relation from EC activity metadata", len(long))
    except Exception as exc:
        renderer.skip(fid, title, [source], str(exc))


def plot_contact_chemistry_right_wrong(renderer: Renderer) -> None:
    fid = "descriptive/contact_chemistry_right_wrong_by_substrate_dp.png"
    title = "Contact chemistry for right versus wrong ligand"
    source = "combined/condition_patch_summary.tsv"
    try:
        df = add_categories(add_right_wrong(renderer, renderer.read_tsv(source)))
        metrics = ["aromatic_contact_fraction", "polar_contact_fraction", "charged_contact_fraction", "hbond_contact_fraction", "catalytic_surface_contact_fraction", "non_core_contact_fraction"]
        metrics = [m for m in metrics if m in df.columns]
        require(df, ["protein_id", "substrate_class", "dp", "ligand_relation", *metrics])
        df = num(df, metrics)
        long = df.melt(id_vars=["protein_id", "substrate_class", "dp", "substrate_dp", "ligand_relation"], value_vars=metrics, var_name="metric", value_name="fraction")
        long["metric_label"] = long["metric"].map(nice)
        long = finite(long, ["fraction"])
        g = sns.catplot(data=long, x="substrate_dp", y="fraction", hue="ligand_relation", col="metric_label", col_wrap=3, kind="point", errorbar=("pi", 50), palette=RELATION_PALETTE, height=3.2, aspect=1.2, sharey=True)
        for ax in g.axes.flat:
            ax.tick_params(axis="x", rotation=45)
            ax.set_ylim(-0.03, 1.03)
        g.set_axis_labels("", "Fraction")
        g.fig.suptitle(title, y=1.04)
        stats_df = paired_stats(df, metrics, ["substrate_class"])
        renderer.write_table(stats_df, "tables/substrate_right_wrong_nonparametric_stats.tsv")
        renderer.save(g.fig, fid, title, [source, "input_data/metadata_final_ec_fixed.tsv"], "right/wrong ligand relation from EC activity metadata", len(long))
    except Exception as exc:
        renderer.skip(fid, title, [source], str(exc))


def plot_geometry_right_wrong(renderer: Renderer) -> None:
    fid = "descriptive/geometry_plausibility_right_wrong_by_substrate_dp.png"
    title = "Geometry plausibility for right versus wrong ligand"
    source = "combined/condition_table.tsv"
    try:
        df = add_categories(add_right_wrong(renderer, renderer.read_tsv(source)))
        require(df, ["occupancy_weighted_c1_plausible_fraction", "occupancy_weighted_c4_plausible_fraction"])
        df = num(df, ["occupancy_weighted_c1_plausible_fraction", "occupancy_weighted_c4_plausible_fraction"])
        df["occupancy_weighted_any_plausibility"] = df[["occupancy_weighted_c1_plausible_fraction", "occupancy_weighted_c4_plausible_fraction"]].max(axis=1)
        df["occupancy_weighted_C1_minus_C4_plausibility"] = df["occupancy_weighted_c1_plausible_fraction"] - df["occupancy_weighted_c4_plausible_fraction"]
        metrics = ["occupancy_weighted_any_plausibility", "occupancy_weighted_C1_minus_C4_plausibility"]
        long = df.melt(id_vars=["protein_id", "substrate_dp", "substrate_class", "dp", "ligand_relation"], value_vars=metrics, var_name="metric", value_name="value")
        long["metric_label"] = long["metric"].map(nice)
        long = finite(long, ["value"])
        g = sns.catplot(data=long, x="substrate_dp", y="value", hue="ligand_relation", col="metric_label", kind="point", errorbar=("pi", 50), palette=RELATION_PALETTE, height=4, aspect=1.35)
        for ax in g.axes.flat:
            ax.tick_params(axis="x", rotation=45)
            ax.axhline(0, color="#555555", linestyle="--", linewidth=0.8)
        g.set_axis_labels("", "Value")
        g.fig.suptitle(title, y=1.05)
        stats_df = paired_stats(df, metrics, ["substrate_class"])
        old_path = renderer.output_dir / "tables/substrate_right_wrong_nonparametric_stats.tsv"
        if old_path.is_file():
            stats_df = pd.concat([pd.read_csv(old_path, sep="\t"), stats_df], ignore_index=True)
        renderer.write_table(stats_df, "tables/substrate_right_wrong_nonparametric_stats.tsv")
        renderer.save(g.fig, fid, title, [source, "input_data/metadata_final_ec_fixed.tsv"], "right/wrong ligand relation from EC activity metadata", len(long))
    except Exception as exc:
        renderer.skip(fid, title, [source], str(exc))


def plot_predictive_baseline(renderer: Renderer) -> None:
    fid = "predictive/predictive_performance_with_pr_auc_baseline.png"
    title = "Predictive PR-AUC performance with prevalence baseline"
    source = "postprocess/10_predictive/cv_results/*_metrics.tsv"
    try:
        frames = []
        for path in sorted((renderer.results_dir / "postprocess/10_predictive/cv_results").glob("*_metrics.tsv")):
            if path.name.endswith("_fold_metrics.tsv"):
                continue
            frame = pd.read_csv(path, sep="\t")
            frame["metric_file"] = path.name
            frames.append(frame)
        df = pd.concat(frames, ignore_index=True)
        require(df, ["model_name", "n_positive", "n_negative", "mean_pr_auc", "std_pr_auc"])
        df = num(df, ["n_positive", "n_negative", "mean_pr_auc", "std_pr_auc"])
        df["baseline_pr_auc"] = df["n_positive"] / (df["n_positive"] + df["n_negative"])
        long = df.melt(id_vars=["model_name"], value_vars=["baseline_pr_auc", "mean_pr_auc"], var_name="metric", value_name="score")
        long["metric"] = long["metric"].map({"baseline_pr_auc": "Baseline", "mean_pr_auc": "Model"})
        fig, ax = plt.subplots(figsize=(9.5, 5))
        sns.barplot(data=long, x="model_name", y="score", hue="metric", palette={"Baseline": "#bdbdbd", "Model": "#3b6ea8"}, ax=ax)
        model = df.reset_index()
        for _, row in model.iterrows():
            x = row["index"] + 0.2
            ax.errorbar(x, row["mean_pr_auc"], yerr=row["std_pr_auc"], color="#222222", capsize=3, fmt="none")
            ax.text(row["index"], 1.02, f"n+={int(row['n_positive'])}\nn-={int(row['n_negative'])}", ha="center", va="bottom", fontsize=7, transform=ax.get_xaxis_transform())
        ax.set_ylim(0, 1.18)
        ax.set_xlabel("")
        ax.set_ylabel("PR-AUC")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=35, ha="right")
        ax.set_title(title)
        ax.legend(frameon=False)
        renderer.write_table(df.assign(delta_vs_baseline=df["mean_pr_auc"] - df["baseline_pr_auc"]), "tables/predictive_performance_and_predictors.tsv")
        renderer.save(fig, fid, title, [source], "all *_metrics.tsv files", len(df))
    except Exception as exc:
        renderer.skip(fid, title, [source], str(exc))


def plot_top_predictors(renderer: Renderer) -> None:
    fid = "predictive/top_predictors_summary.png"
    title = "Top predictors across exploratory models"
    source = "postprocess/10_predictive/cv_results/*_metrics.tsv"
    try:
        rows = []
        for path in sorted((renderer.results_dir / "postprocess/10_predictive/cv_results").glob("*_metrics.tsv")):
            if path.name.endswith("_fold_metrics.tsv"):
                continue
            df = pd.read_csv(path, sep="\t")
            for _, row in df.iterrows():
                values = json.loads(str(row.get("feature_importances_json", "{}")))
                for feature, importance in values.items():
                    rows.append({"model_name": row["model_name"], "feature": feature, "importance": importance})
        imp = pd.DataFrame(rows)
        require(imp, ["model_name", "feature", "importance"])
        imp = num(imp, ["importance"])
        order = imp.groupby("feature")["importance"].mean().sort_values(ascending=False).index
        fig, ax = plt.subplots(figsize=(9, 5.5))
        sns.barplot(data=imp, y="feature", x="importance", hue="model_name", order=order, ax=ax)
        ax.set_title(title)
        ax.set_xlabel("Absolute standardized coefficient / importance")
        ax.set_ylabel("")
        ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
        perf_path = renderer.output_dir / "tables/predictive_performance_and_predictors.tsv"
        if perf_path.is_file():
            perf = pd.read_csv(perf_path, sep="\t")
            top = imp.sort_values(["model_name", "importance"], ascending=[True, False]).groupby("model_name").head(3)
            top_txt = top.groupby("model_name").apply(lambda x: "; ".join(x["feature"])).rename("top_predictors").reset_index()
            perf = perf.merge(top_txt, on="model_name", how="left")
            renderer.write_table(perf, "tables/predictive_performance_and_predictors.tsv")
        renderer.save(fig, fid, title, [source], "feature_importances_json parsed from metric files", len(imp))
    except Exception as exc:
        renderer.skip(fid, title, [source], str(exc))


def plot_family_heatmaps(renderer: Renderer) -> None:
    title = "Family alignment residue enrichment"
    sources = [
        "postprocess/08_family_residue_enrichment/family_residue_enrichment.tsv",
        "postprocess/08_family_residue_enrichment/family_substrate_residue_enrichment.tsv",
    ]
    try:
        df = renderer.read_tsv(sources[0])
        sub_df = renderer.read_tsv(sources[1])
        df = num(df, ["alignment_column", "n_proteins_observed", "mean_condition_residue_contact_score", "mean_c1_minus_c4_weighted_delta"])
        sub_df = num(sub_df, ["alignment_column", "target_minus_other_substrate_residue_contact_delta"])
        ranked_rows = []
        for family, fam in df.groupby("family_label", observed=True):
            fam = fam[fam["n_proteins_observed"] >= 2].copy()
            top_cols = fam.sort_values("mean_condition_residue_contact_score", ascending=False).head(50)["alignment_column"]
            fam = fam[fam["alignment_column"].isin(top_cols)]
            pivot = fam.pivot_table(index="alignment_column", values=["mean_condition_residue_contact_score", "mean_c1_minus_c4_weighted_delta"], aggfunc="mean")
            if pivot.empty:
                continue
            fig, ax = plt.subplots(figsize=heatmap_size(len(pivot), len(pivot.columns), 5, 5))
            sns.heatmap(pivot, cmap="mako", ax=ax, cbar_kws={"label": "Score / delta"})
            ax.set_title(f"{title}: {family}")
            ax.set_xlabel("Metric")
            ax.set_ylabel("Alignment column")
            fid = f"family_residue/family_{family}_alignment_enrichment_heatmap.png"
            renderer.save(fig, fid, f"{title}: {family}", [sources[0]], "n_proteins_observed >= 2; top 50 columns", len(fam))
            ranked_rows.append(fam.assign(source="alignment_enrichment"))
            fam_sub = sub_df[sub_df["family_label"].eq(family)].copy()
            if not fam_sub.empty:
                top_sub_cols = fam_sub.reindex(fam_sub["target_minus_other_substrate_residue_contact_delta"].abs().sort_values(ascending=False).index).head(70)["alignment_column"]
                fam_sub = fam_sub[fam_sub["alignment_column"].isin(top_sub_cols)]
                pivot2 = fam_sub.pivot_table(index="alignment_column", columns="target_substrate", values="target_minus_other_substrate_residue_contact_delta", aggfunc="mean")
                fig, ax = plt.subplots(figsize=heatmap_size(len(pivot2), len(pivot2.columns), 5, 5))
                max_abs = np.nanmax(np.abs(pivot2.to_numpy())) if pivot2.size else 1
                max_abs = 1 if not math.isfinite(max_abs) or max_abs == 0 else max_abs
                sns.heatmap(pivot2, cmap="vlag", center=0, vmin=-max_abs, vmax=max_abs, ax=ax, cbar_kws={"label": "Target - other contact delta"})
                ax.set_title(f"Substrate enrichment: {family}")
                ax.set_xlabel("Target substrate")
                ax.set_ylabel("Alignment column")
                fid2 = f"family_residue/family_{family}_substrate_enrichment_heatmap.png"
                renderer.save(fig, fid2, f"Substrate enrichment: {family}", [sources[1]], "top absolute substrate deltas", len(fam_sub))
        if ranked_rows:
            ranked = pd.concat(ranked_rows, ignore_index=True)
            ranked.sort_values(["family_label", "mean_condition_residue_contact_score"], ascending=[True, False]).to_csv(
                renderer.output_dir / "tables/family_enriched_residues_ranked.tsv", sep="\t", index=False
            )
    except Exception as exc:
        renderer.skip("family_residue/family_alignment_enrichment_heatmaps.png", title, sources, str(exc))


def plot_top_occupancy_nonzero(renderer: Renderer) -> None:
    fid = "descriptive/top_cluster_occupancy_distribution_nonzero.png"
    title = "Top-cluster occupancy distribution, nonzero values"
    source = "combined/condition_table.tsv"
    try:
        df = add_categories(renderer.read_tsv(source))
        require(df, ["top_cluster_occupancy", "construct_type", "substrate_class"])
        df = finite(num(df, ["top_cluster_occupancy"]), ["top_cluster_occupancy"])
        excluded = int((df["top_cluster_occupancy"] <= 0).sum())
        df = df[df["top_cluster_occupancy"] > 0].copy()
        fig, ax = plt.subplots(figsize=(9.5, 5))
        sns.boxplot(data=df, x="substrate_class", y="top_cluster_occupancy", hue="construct_type", palette=CONSTRUCT_PALETTE, showfliers=False, ax=ax)
        sns.stripplot(data=df, x="substrate_class", y="top_cluster_occupancy", hue="construct_type", dodge=True, color="#222222", alpha=0.25, size=2, ax=ax, legend=False)
        ax.set_ylim(-0.03, 1.03)
        ax.set_title(title)
        ax.set_xlabel("Substrate")
        ax.set_ylabel("Top-cluster occupancy")
        ax.legend(frameon=False, title="Construct")
        annotate_counts(ax, f"n shown={len(df)}\nzero excluded={excluded}")
        renderer.save(fig, fid, title, [source], "top_cluster_occupancy > 0", len(df))
    except Exception as exc:
        renderer.skip(fid, title, [source], str(exc))


def plot_cluster_distribution(renderer: Renderer) -> None:
    fid = "descriptive/cluster_distribution_summary.png"
    title = "Cluster distribution and stability summaries"
    source = "combined/condition_table.tsv"
    try:
        df = add_categories(renderer.read_tsv(source))
        metrics = ["n_clusters", "noise_fraction", "cluster_entropy", "occupancy_gini"]
        require(df, ["construct_type", "substrate_class", *metrics])
        df = num(df, metrics)
        long = df.melt(id_vars=["construct_type", "substrate_class", "dp"], value_vars=metrics, var_name="metric", value_name="value")
        long["metric_label"] = long["metric"].map(nice)
        long = finite(long, ["value"])
        g = sns.displot(data=long, x="value", hue="construct_type", col="metric_label", col_wrap=2, kind="hist", kde=True, palette=CONSTRUCT_PALETTE, height=3.6, aspect=1.25, facet_kws={"sharex": False, "sharey": False})
        g.set_axis_labels("Value", "Count")
        g.fig.suptitle(title, y=1.04)
        summary = df.groupby(["construct_type", "substrate_class", "dp"], observed=True).agg(
            n_conditions=("condition_id", "nunique"),
            n_clusterable=("any_valid_cluster", lambda x: bool_col(x).sum()),
            median_top_occupancy=("top_cluster_occupancy", "median"),
            median_entropy=("cluster_entropy", "median"),
            median_gini=("occupancy_gini", "median"),
            median_ligand_rmsd=("median_ligand_rmsd", "median"),
            iqr_ligand_rmsd=("iqr_ligand_rmsd", "median"),
        ).reset_index()
        renderer.write_table(summary, "tables/cluster_stability_summary.tsv")
        renderer.save(g.fig, fid, title, [source], "all finite cluster stability metrics", len(long))
    except Exception as exc:
        renderer.skip(fid, title, [source], str(exc))


def plot_within_cluster_rmsd(renderer: Renderer) -> None:
    fid = "descriptive/within_cluster_ligand_rmsd_distribution.png"
    title = "Within-condition and cluster ligand RMSD distributions"
    sources = ["combined/condition_table.tsv", "combined/cluster_table.tsv"]
    try:
        cond = add_categories(renderer.read_tsv(sources[0]))
        cluster = add_categories(renderer.read_tsv(sources[1]))
        frames = []
        if "median_ligand_rmsd" in cond.columns:
            frames.append(cond[["construct_type", "substrate_class", "median_ligand_rmsd"]].rename(columns={"median_ligand_rmsd": "rmsd"}).assign(level="condition median"))
        if "median_ligand_rmsd_to_reference" in cluster.columns:
            frames.append(cluster[["construct_type", "substrate_class", "median_ligand_rmsd_to_reference"]].rename(columns={"median_ligand_rmsd_to_reference": "rmsd"}).assign(level="cluster to reference"))
        df = pd.concat(frames, ignore_index=True)
        df = finite(num(df, ["rmsd"]), ["rmsd"])
        fig, ax = plt.subplots(figsize=(9.5, 5))
        sns.boxplot(data=df, x="substrate_class", y="rmsd", hue="level", showfliers=False, ax=ax)
        sns.stripplot(data=df, x="substrate_class", y="rmsd", hue="level", dodge=True, color="#222222", alpha=0.18, size=2, ax=ax, legend=False)
        ax.set_title(title)
        ax.set_xlabel("Substrate")
        ax.set_ylabel("Ligand RMSD")
        ax.legend(frameon=False)
        renderer.save(fig, fid, title, sources, "finite condition and/or cluster ligand RMSD metrics", len(df))
    except Exception as exc:
        renderer.skip(fid, title, sources, str(exc))


def plot_clusterable_vs_pdb(renderer: Renderer) -> None:
    fid = "crystal/clusterable_rate_vs_pre_cutoff_crystal_count.png"
    title = "Clusterable rate versus pre-AF3-cutoff PDB count"
    sources = ["combined/condition_table.tsv", "input_data/metadata_pdb_count.tsv"]
    try:
        cond = renderer.read_tsv(sources[0])
        pdb = renderer.read_tsv(sources[1], base="analyse")
        require(cond, ["protein_id", "any_valid_cluster"])
        require(pdb, ["UniProt_ID", "PDB_Before_AF3_Cutoff_Count"])
        cond["any_valid_cluster_bool"] = bool_col(cond["any_valid_cluster"])
        rates = cond.groupby("protein_id", observed=True).agg(clusterable_rate=("any_valid_cluster_bool", "mean"), n_conditions=("condition_id", "nunique")).reset_index()
        pdb = num(pdb, ["PDB_Before_AF3_Cutoff_Count"]).rename(columns={"UniProt_ID": "protein_id"})
        df = rates.merge(pdb[["protein_id", "PDB_Before_AF3_Cutoff_Count", "PDB_IDs_Before_AF3_Cutoff"]], on="protein_id", how="left")
        df["PDB_Before_AF3_Cutoff_Count"] = df["PDB_Before_AF3_Cutoff_Count"].fillna(0)
        fig, ax = plt.subplots(figsize=(7, 5.2))
        sns.scatterplot(data=df, x="PDB_Before_AF3_Cutoff_Count", y="clusterable_rate", size="n_conditions", sizes=(35, 160), color="#3b6ea8", ax=ax)
        if len(df) >= 3:
            sns.regplot(data=df, x="PDB_Before_AF3_Cutoff_Count", y="clusterable_rate", scatter=False, lowess=True, color="#222222", ax=ax)
        ax.set_title(title)
        ax.set_xlabel("PDB structures before AF3 cutoff")
        ax.set_ylabel("Protein-level valid-cluster rate")
        ax.set_ylim(-0.03, 1.03)
        annotate_counts(ax, f"n proteins={len(df)}")
        renderer.save(fig, fid, title, sources, "condition_table joined to metadata_pdb_count by UniProt/protein_id", len(df))
    except Exception as exc:
        renderer.skip(fid, title, sources, str(exc))


FIGURES: list[Callable[[Renderer], None]] = [
    plot_cbm_qc_clusterable,
    plot_cbm_clustered_matrix,
    plot_cbm_delta_bias,
    plot_pocket_rmsd_by_ligand,
    plot_holo_ifp_vs_rmsd,
    plot_his_brace_quality,
    plot_iptm_vs_qc,
    plot_right_wrong_c1c4,
    plot_right_wrong_angles_distances,
    plot_contact_chemistry_right_wrong,
    plot_geometry_right_wrong,
    plot_predictive_baseline,
    plot_top_predictors,
    plot_family_heatmaps,
    plot_top_occupancy_nonzero,
    plot_cluster_distribution,
    plot_within_cluster_rmsd,
    plot_clusterable_vs_pdb,
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--analyse-dir", type=Path, default=Path("Masteroppgave_clean/analyse"))
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_style()
    output_dir = args.output_dir or args.results_dir / "postprocess" / "figure2"
    renderer = Renderer(args.results_dir, output_dir, args.analyse_dir)
    for func in FIGURES:
        func(renderer)
    renderer.manifest()
    ok = sum(r.status == "ok" for r in renderer.records)
    skipped = sum(r.status != "ok" for r in renderer.records)
    print(f"Rendered {ok} figure records; skipped {skipped}.")
    print(f"Manifest: {renderer.output_dir / 'figure_manifest.tsv'}")
    return 2 if args.strict and skipped else 0


if __name__ == "__main__":
    raise SystemExit(main())
