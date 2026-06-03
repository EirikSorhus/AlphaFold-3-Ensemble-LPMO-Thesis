#!/usr/bin/env python3
"""Render an extra noise-versus-clustered pose fraction plot."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


SUBSTRATE_ORDER = ["chitin", "cellulose", "amylose", "starch"]
CONSTRUCT_ORDER = ["domain_only", "full_length"]
CONSTRUCT_LABELS = {"domain_only": "Catalytic domain", "full_length": "Full-length"}
LIGHT_PALETTE = {"domain_only": "#b8cce4", "full_length": "#e8bd89"}
DARK_PALETTE = {"domain_only": "#2f5f8f", "full_length": "#a85f17"}


def centered_violin_boxplot(ax: plt.Axes, df: pd.DataFrame, *, value_col: str, y_label: str) -> None:
    values: list[np.ndarray] = []
    positions: list[float] = []
    colors: list[tuple[str, str]] = []
    substrates = [s for s in SUBSTRATE_ORDER if s in set(df["substrate_class"].dropna().astype(str))]
    if not substrates:
        substrates = sorted(df["substrate_class"].dropna().astype(str).unique())
    offsets = {"domain_only": -0.18, "full_length": 0.18}
    for ix, substrate in enumerate(substrates):
        for construct in CONSTRUCT_ORDER:
            sub = df[(df["substrate_class"].astype(str) == substrate) & (df["construct_type"].astype(str) == construct)]
            arr = pd.to_numeric(sub[value_col], errors="coerce").dropna().to_numpy()
            if arr.size == 0:
                continue
            values.append(arr)
            positions.append(ix + offsets[construct])
            colors.append((LIGHT_PALETTE[construct], DARK_PALETTE[construct]))

    violin = ax.violinplot(values, positions=positions, widths=0.32, showmeans=False, showmedians=False, showextrema=False)
    for body, (light, _) in zip(violin["bodies"], colors, strict=False):
        body.set_facecolor(light)
        body.set_edgecolor("#666666")
        body.set_alpha(0.9)
        body.set_linewidth(0.9)

    boxes = ax.boxplot(values, positions=positions, widths=0.075, patch_artist=True, showfliers=False, manage_ticks=False)
    for box, (_, dark) in zip(boxes["boxes"], colors, strict=False):
        box.set_facecolor(dark)
        box.set_edgecolor("#303030")
        box.set_linewidth(0.9)
    for key in ["whiskers", "caps", "medians"]:
        for artist in boxes[key]:
            artist.set_color("#303030")
            artist.set_linewidth(0.9)

    ax.set_xticks(range(len(substrates)))
    ax.set_xticklabels(substrates)
    ax.set_xlabel("Substrate")
    ax.set_ylabel(y_label)
    ax.set_ylim(-0.03, 1.03)
    ax.grid(False)
    ax.tick_params(axis="x", which="major", length=0, width=0.8, direction="out", colors="#000000")
    ax.tick_params(axis="y", which="major", length=3.5, width=0.8, direction="out", colors="#000000")
    for spine in ax.spines.values():
        spine.set_color("#000000")
        spine.set_linewidth(0.8)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.08, 1.06, label, transform=ax.transAxes, ha="left", va="top", fontsize=13, weight="bold")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    table_dir = args.output_dir / "tables"
    table_dir.mkdir(exist_ok=True)

    df = pd.read_csv(args.results_dir / "combined" / "condition_table.tsv", sep="\t", low_memory=False)
    for col in ["noise_fraction", "n_noise", "n_ifp_clustered"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[df["construct_type"].isin(CONSTRUCT_ORDER)].copy()
    df = df[df["noise_fraction"].notna()].copy()
    df["clustered_pose_fraction"] = 1.0 - df["noise_fraction"]

    summary = (
        df.groupby(["construct_type", "substrate_class"], observed=True)
        .agg(
            n_conditions=("condition_id", "nunique"),
            median_noise_fraction=("noise_fraction", "median"),
            median_clustered_pose_fraction=("clustered_pose_fraction", "median"),
            mean_noise_fraction=("noise_fraction", "mean"),
            mean_clustered_pose_fraction=("clustered_pose_fraction", "mean"),
        )
        .reset_index()
    )
    summary["construct"] = summary["construct_type"].map(CONSTRUCT_LABELS)
    summary.to_csv(table_dir / "extra_noise_vs_clustered_pose_fraction.tsv", sep="\t", index=False)

    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#000000",
            "axes.labelcolor": "#000000",
            "xtick.color": "#000000",
            "ytick.color": "#000000",
            "text.color": "#000000",
            "axes.grid": False,
            "legend.frameon": False,
            "font.size": 9,
        }
    )
    fig, ax = plt.subplots(figsize=(7.0, 5.3), constrained_layout=True)
    centered_violin_boxplot(ax, df, value_col="noise_fraction", y_label="Noise fraction")
    panel_label(ax, "A")
    handles = [Patch(facecolor=DARK_PALETTE[c], edgecolor="#303030", label=CONSTRUCT_LABELS[c]) for c in CONSTRUCT_ORDER]
    ax.legend(handles=handles, title="Construct", frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0)
    out = args.output_dir / "extra_noise_vs_clustered_pose_fraction.png"
    fig.savefig(out, dpi=240, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
