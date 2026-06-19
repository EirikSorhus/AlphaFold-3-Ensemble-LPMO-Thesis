#!/usr/bin/env python3
"""Render thesis-oriented main figure candidates into postprocess/figures3."""

from __future__ import annotations

import argparse
import itertools
import json
import math
import re
import shlex
import textwrap
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from sklearn.metrics import precision_recall_curve, roc_auc_score, roc_curve


SUBSTRATE_ORDER = ["chitin", "cellulose", "starch", "amylose"]
DEMONSTRATED_ACTIVITY_GROUP_ORDER = ["Chitin", "Cellulose", "Mixed cellulose/chitin", "Starch", "Other"]
CONSTRUCT_ORDER = ["domain_only", "full_length"]
ACTIVITY_ORDER = ["Demonstrated activity", "No demonstrated activity"]
SUBSTRATE_PALETTE = {"chitin": "#2b8cbe", "cellulose": "#41ab5d", "starch": "#e34a33", "amylose": "#e34a33"}
CONSTRUCT_PALETTE = {"domain_only": "#3b6ea8", "full_length": "#d17a22"}
ACTIVITY_PALETTE = {"Demonstrated activity": "#2f7f4f", "No demonstrated activity": "#9b3a3a"}
FAMILY_PALETTE = {"AA9": "#3b6ea8", "AA10": "#d17a22"}
MARKER_ORDER = [
    "Catalytic domain | Demonstrated activity",
    "Catalytic domain | No demonstrated activity",
    "Full-length | Demonstrated activity",
    "Full-length | No demonstrated activity",
]
MARKERS = {
    "Catalytic domain | Demonstrated activity": "o",
    "Catalytic domain | No demonstrated activity": "X",
    "Full-length | Demonstrated activity": "^",
    "Full-length | No demonstrated activity": "P",
}
BASE_FONTSIZE = 12.0
LEGEND_FONTSIZE = 10.5
PANEL_FONTSIZE = 16
CONSTRUCT_DISPLAY = {
    "domain_only": "Catalytic domain",
    "full_length": "Full-length",
    "Domain only": "Catalytic domain",
    "Full length": "Full-length",
    "full length": "Full-length",
    "Full-length": "Full-length",
    "full-length": "Full-length",
    "catalytic domain": "Catalytic domain",
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


class Renderer:
    def __init__(self, results_dir: Path, analyse_dir: Path, output_dir: Path) -> None:
        self.results_dir = results_dir.resolve()
        self.analyse_dir = analyse_dir.resolve()
        self.output_dir = output_dir.resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / "tables").mkdir(exist_ok=True)
        self.records: list[FigureRecord] = []
        self.cache: dict[Path, pd.DataFrame] = {}

    def rel(self, path: Path) -> str:
        try:
            return str(path.resolve().relative_to(self.results_dir))
        except ValueError:
            return str(path)

    def tsv(self, rel_path: str, *, base: str = "results") -> pd.DataFrame:
        root = self.results_dir if base == "results" else self.analyse_dir
        path = root / rel_path
        if path not in self.cache:
            self.cache[path] = pd.read_csv(path, sep="\t", low_memory=False)
        return self.cache[path].copy()

    def save(self, fig: plt.Figure, figure_id: str, title: str, sources: list[str], filters: str, n: int) -> None:
        out = self.output_dir / figure_id
        fig.savefig(out, dpi=600, bbox_inches="tight")
        fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig)
        self.records.append(FigureRecord(figure_id, title, ";".join(sources), filters, int(n), self.rel(out)))

    def skip(self, figure_id: str, title: str, sources: list[str], message: str) -> None:
        self.records.append(FigureRecord(figure_id, title, ";".join(sources), "not_rendered", 0, self.rel(self.output_dir / figure_id), "skipped", message))

    def table(self, df: pd.DataFrame, name: str) -> None:
        df.to_csv(self.output_dir / "tables" / name, sep="\t", index=False)

    def manifest(self) -> None:
        pd.DataFrame([r.__dict__ for r in self.records]).to_csv(self.output_dir / "figure_manifest.tsv", sep="\t", index=False)


def style() -> None:
    sns.set_theme(
        context="paper",
        style="ticks",
        font_scale=1.05,
        rc={
            "font.size": BASE_FONTSIZE,
            "axes.labelsize": BASE_FONTSIZE + 0.6,
            "xtick.labelsize": BASE_FONTSIZE - 0.5,
            "ytick.labelsize": BASE_FONTSIZE - 0.5,
            "legend.fontsize": LEGEND_FONTSIZE,
            "legend.title_fontsize": LEGEND_FONTSIZE,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#000000",
            "axes.labelcolor": "#000000",
            "xtick.color": "#000000",
            "ytick.color": "#000000",
            "text.color": "#000000",
            "axes.titleweight": "bold",
            "axes.grid": False,
            "legend.frameon": False,
        },
    )


def need(df: pd.DataFrame, cols: list[str]) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(", ".join(missing))


def num(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for c in cols:
        if c in out.columns:
            out[c] = pd.to_numeric(out[c], errors="coerce")
    return out


def bools(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().isin(["true", "1", "yes", "ja"])


def cats(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "construct_type" not in out.columns and "_construct_type" in out.columns:
        out["construct_type"] = out["_construct_type"]
    for col, order in [("construct_type", CONSTRUCT_ORDER), ("substrate_class", SUBSTRATE_ORDER)]:
        if col in out.columns:
            values = set(out[col].dropna().astype(str))
            out[col] = pd.Categorical(out[col].astype(str), [x for x in order if x in values] + sorted(values - set(order)), ordered=True)
    if "dp" in out.columns:
        out["dp"] = pd.to_numeric(out["dp"], errors="coerce")
    if {"substrate_class", "dp"}.issubset(out.columns):
        out["substrate_dp"] = out["substrate_class"].astype(str) + " DP" + out["dp"].astype("Int64").astype(str)
    return out


def finite(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for c in cols:
        values = pd.to_numeric(out[c], errors="coerce")
        out = out[values.notna() & np.isfinite(values)].copy()
    return out


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


def ec_active_substrates_with_other(text: object) -> set[str]:
    out: set[str] = set()
    for token in ec_tokens(text):
        if token == "1.14.99.53":
            out.add("chitin")
        elif token in {"1.14.99.54", "1.14.99.56"}:
            out.add("cellulose")
        elif token == "1.14.99.55":
            out.add("starch")
        elif token.startswith("1.14.99."):
            out.add("other")
    return out


def demonstrated_activity_group(active: set[str]) -> str:
    if {"cellulose", "chitin"}.issubset(active):
        return "Mixed cellulose/chitin"
    if "cellulose" in active:
        return "Cellulose"
    if "chitin" in active:
        return "Chitin"
    if "starch" in active:
        return "Starch"
    if "other" in active:
        return "Other"
    return ""


def ec_tokens(text: object) -> list[str]:
    return re.findall(r"1\.14\.99\.[0-9-]+", str(text or ""))


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


def activity_table(r: Renderer) -> pd.DataFrame:
    rows = []
    meta = r.tsv("input_data/metadata_final_ec_fixed.tsv", base="analyse")
    for _, row in meta.iterrows():
        active = ec_active_substrates(row.get("EC_Number"))
        rows.append(
            {
                "protein_id": row.get("UniProt_ID"),
                "family": row.get("CAZy_family"),
                "active_substrates": ",".join(sorted(active)),
                "regio": ec_regio_label(row.get("EC_Number"), row.get("CAZy_family")),
            }
        )
    out = pd.DataFrame(rows).dropna(subset=["protein_id"])
    return out.drop_duplicates("protein_id")


def add_activity(r: Renderer, df: pd.DataFrame) -> pd.DataFrame:
    out = df.merge(activity_table(r), on="protein_id", how="left")
    active_sets = out["active_substrates"].fillna("").map(lambda x: set(filter(None, str(x).split(","))))
    out["activity_evidence"] = [
        ACTIVITY_ORDER[0] if ("starch" if str(substrate) == "amylose" else str(substrate)) in active else (ACTIVITY_ORDER[1] if active else np.nan)
        for substrate, active in zip(out["substrate_class"].astype(str), active_sets, strict=False)
    ]
    return out[out["activity_evidence"].notna()].copy()


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.11, 1.08, label, transform=ax.transAxes, ha="left", va="top", fontsize=PANEL_FONTSIZE, weight="bold")


def nice(name: str) -> str:
    return re.sub(r"^(occupancy_weighted_|delta_|mean_|median_)", "", name).replace("_", " ")


def sentence_case_label(label: object) -> str:
    text = str(label)
    if not text:
        return text
    mapped = CONSTRUCT_DISPLAY.get(text, text)
    mapped = re.sub(r"\bfull[ -]?length\b", "Full-length", mapped, flags=re.IGNORECASE)
    mapped = re.sub(r"\bproven active\b", "demonstrated activity", mapped, flags=re.IGNORECASE)
    mapped = re.sub(r"\bno proven activity\b", "no demonstrated activity", mapped, flags=re.IGNORECASE)
    mapped = re.sub(r"\bproven activity\b", "demonstrated activity", mapped, flags=re.IGNORECASE)
    mapped = re.sub(r"\bproven\b", "demonstrated", mapped, flags=re.IGNORECASE)
    replacements = {
        "construct_type": "Construct type",
        "substrate_class": "Substrate class",
        "model_name": "Model name",
        "activity evidence": "Activity evidence",
        "Symbol_group": "Symbol group",
        "symbol_group": "Symbol group",
        "site": "Site",
    }
    mapped = replacements.get(mapped, mapped)
    return mapped[:1].upper() + mapped[1:]


def wrap_label(label: object, width: int = 18) -> str:
    text = sentence_case_label(label)
    text = text.replace("Chitin-active", "Chitin\nactive")
    text = text.replace("Cellulose-active", "Cellulose\nactive")
    text = text.replace("Mixed C1/C4", "Mixed\nC1/C4")
    text = text.replace(" | ", "\n")
    text = text.replace(" versus ", "\nversus ")
    text = text.replace(" vs ", "\nvs ")
    if len(text) <= width or "\n" in text:
        return text
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False, break_on_hyphens=False))


def wrap_ticklabels(ax: plt.Axes, *, axis: str = "x", width: int = 14) -> None:
    if axis == "x":
        labels = [wrap_label(label.get_text(), width) for label in ax.get_xticklabels()]
        ax.set_xticks(ax.get_xticks())
        ax.set_xticklabels(labels)
    else:
        labels = [wrap_label(label.get_text(), width) for label in ax.get_yticklabels()]
        ax.set_yticks(ax.get_yticks())
        ax.set_yticklabels(labels)


def compact_margins(ax: plt.Axes, *, x: float = 0.01, y: float = 0.04) -> None:
    ax.margins(x=x, y=y)


def substrate_dp_ticklabels(ax: plt.Axes, *, rotation: float = 40, abbreviate: bool = False) -> None:
    labels = []
    abbr = {"amylose": "Amy.", "cellulose": "Cell.", "chitin": "Chit.", "starch": "Starch"}
    for label in ax.get_xticklabels():
        raw = label.get_text()
        if abbreviate:
            parts = raw.split()
            if len(parts) == 2:
                text = f"{abbr.get(parts[0].lower(), parts[0].title())}\n{parts[1].upper()}"
            else:
                text = raw.title().replace("Dp", "DP")
        else:
            text = raw.replace(" ", "\n").title().replace("Dp", "DP")
        labels.append(text)
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(labels, rotation=rotation, ha="right" if rotation else "center")


def title_case_label(name: str) -> str:
    return nice(name).replace("ligand", "substrate").replace("c1", "C1").replace("c4", "C4").replace("ifp", "IFP").title().replace("C1", "C1").replace("C4", "C4").replace("Ifp", "IFP")


def strip_residue_chain(label: object) -> str:
    return re.sub(r"\.A\b", "", str(label))


def centered_violin_boxplot(
    ax: plt.Axes,
    df: pd.DataFrame,
    *,
    value_col: str,
    y_label: str,
    light_palette: dict[str, str],
    dark_palette: dict[str, str],
) -> None:
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
            positions.append(ix + offsets.get(construct, 0.0))
            colors.append((light_palette.get(construct, "#cccccc"), dark_palette.get(construct, "#555555")))

    if values:
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


def average_ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and sorted_values[end] == sorted_values[start]:
            end += 1
        avg_rank = 0.5 * (start + 1 + end)
        ranks[order[start:end]] = avg_rank
        start = end
    return ranks


def kruskal_wallis_h(values: np.ndarray, labels: np.ndarray) -> float:
    if len(values) == 0:
        return np.nan
    unique_labels = pd.unique(labels)
    if len(unique_labels) < 2:
        return np.nan
    ranks = average_ranks(values.astype(float))
    n_total = float(len(values))
    h_value = 0.0
    for label in unique_labels:
        group_ranks = ranks[labels == label]
        if len(group_ranks) == 0:
            continue
        h_value += (float(group_ranks.sum()) ** 2) / float(len(group_ranks))
    h_value = (12.0 / (n_total * (n_total + 1.0))) * h_value - 3.0 * (n_total + 1.0)
    _, counts = np.unique(values, return_counts=True)
    tie_sum = np.sum(counts**3 - counts)
    tie_correction = 1.0 - tie_sum / (n_total**3 - n_total) if n_total > 1 else 1.0
    if tie_correction > 0:
        h_value /= tie_correction
    return float(max(h_value, 0.0))


def permutation_kruskal(
    values: np.ndarray,
    labels: np.ndarray,
    *,
    n_permutations: int = 10000,
    seed: int = 20260604,
) -> tuple[float, float]:
    observed = kruskal_wallis_h(values, labels)
    if not np.isfinite(observed):
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n_permutations):
        permuted = rng.permutation(labels)
        if kruskal_wallis_h(values, permuted) >= observed - 1e-12:
            hits += 1
    p_value = (hits + 1.0) / (n_permutations + 1.0)
    return observed, float(p_value)


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    q_values = [np.nan] * len(p_values)
    valid = [(idx, p) for idx, p in enumerate(p_values) if np.isfinite(p)]
    if not valid:
        return q_values
    order = sorted(valid, key=lambda item: item[1], reverse=True)
    running = 1.0
    m = len(valid)
    for rank_from_end, (idx, p_value) in enumerate(order, start=1):
        rank = m - rank_from_end + 1
        running = min(running, p_value * m / rank)
        q_values[idx] = float(min(running, 1.0))
    return q_values


def annotate(ax: plt.Axes, text: str) -> None:
    ax.text(0.99, 0.98, text, transform=ax.transAxes, ha="right", va="top", fontsize=9.8)


def metric_n_text(df: pd.DataFrame, group_col: str, value_col: str, *, groups: list[str] | None = None, prefix: str = "n") -> str:
    if groups is None:
        groups = [str(group) for group in df[group_col].dropna().unique()]
    rows = []
    for group in groups:
        n_value = int(df[df[group_col].astype(str).eq(str(group))][value_col].notna().sum())
        rows.append(f"{group}: {prefix}={n_value}")
    return "\n".join(rows)


def axis_break_marks(ax: plt.Axes) -> None:
    kwargs = {"transform": ax.transAxes, "color": "#000000", "clip_on": False, "linewidth": 1.0}
    ax.plot((0.055, 0.075), (-0.014, 0.014), **kwargs)
    ax.plot((0.085, 0.105), (-0.014, 0.014), **kwargs)
    ax.plot((-0.014, 0.014), (0.055, 0.075), **kwargs)
    ax.plot((-0.014, 0.014), (0.085, 0.105), **kwargs)


def ticklabels_are_names(labels: list[plt.Text]) -> bool:
    texts = [label.get_text().strip() for label in labels if label.get_text().strip()]
    if not texts:
        return False
    for text in texts:
        cleaned = text.replace(",", "").replace("−", "-").replace("%", "")
        try:
            float(cleaned)
        except ValueError:
            return True
    return False


def polish(fig: plt.Figure) -> None:
    for ax in fig.axes:
        ax.set_xlabel(sentence_case_label(ax.get_xlabel()))
        ax.set_ylabel(sentence_case_label(ax.get_ylabel()))
        legend = ax.get_legend()
        if legend is not None:
            if legend.get_title() is not None:
                legend.get_title().set_text(wrap_label(legend.get_title().get_text(), 18))
                legend.get_title().set_fontsize(LEGEND_FONTSIZE)
            for text in legend.get_texts():
                text.set_text(wrap_label(text.get_text(), 24))
                text.set_fontsize(LEGEND_FONTSIZE)
        ax.grid(False)
        x_length = 0.0 if ticklabels_are_names(ax.get_xticklabels()) else 3.5
        y_length = 0.0 if ticklabels_are_names(ax.get_yticklabels()) else 3.5
        ax.tick_params(axis="x", which="major", length=x_length, width=0.8, direction="out", colors="#000000")
        ax.tick_params(axis="y", which="major", length=y_length, width=0.8, direction="out", colors="#000000")
        ax.tick_params(axis="both", which="minor", length=2.0, width=0.6, direction="out", colors="#000000")
        for spine in ax.spines.values():
            spine.set_color("#000000")
            spine.set_linewidth(0.8)
    for legend in fig.legends:
        if legend.get_title() is not None:
            legend.get_title().set_text(wrap_label(legend.get_title().get_text(), 18))
            legend.get_title().set_fontsize(LEGEND_FONTSIZE)
        for text in legend.get_texts():
            text.set_text(wrap_label(text.get_text(), 24))
            text.set_fontsize(LEGEND_FONTSIZE)


def save_clean(r: Renderer, fig: plt.Figure, figure_id: str, title: str, sources: list[str], filters: str, n: int) -> None:
    for ax in fig.axes:
        ax.set_title("")
    if getattr(fig, "_suptitle", None) is not None:
        fig._suptitle.set_text("")
    try:
        fig.set_constrained_layout_pads(w_pad=0.015, h_pad=0.015, wspace=0.03, hspace=0.035)
    except Exception:
        pass
    polish(fig)
    r.save(fig, figure_id, title, sources, filters, n)


def _kabsch_align(mobile: np.ndarray, reference: np.ndarray, points: np.ndarray) -> np.ndarray:
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (reference - reference_center)
    u, _s, vt = np.linalg.svd(covariance)
    sign = np.sign(np.linalg.det(u @ vt))
    rotation = u @ np.diag([1.0, 1.0, sign]) @ vt
    return (points - mobile_center) @ rotation + reference_center


def _atom_site_arrays(cif_path: str) -> tuple[np.ndarray, np.ndarray]:
    headers: list[str] = []
    in_atom_site = False
    ca_coords: list[list[float]] = []
    ligand_coords: list[list[float]] = []
    with open(cif_path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line == "#":
                continue
            if line == "loop_":
                headers = []
                in_atom_site = False
                continue
            if line.startswith("_"):
                if line.startswith("_atom_site."):
                    headers.append(line.replace("_atom_site.", ""))
                    in_atom_site = True
                elif in_atom_site and headers:
                    break
                continue
            if not in_atom_site or not headers:
                continue
            parts = shlex.split(line)
            if len(parts) < len(headers):
                continue
            row = dict(zip(headers, parts, strict=False))
            atom_name = row.get("label_atom_id") or row.get("auth_atom_id") or ""
            atom_type = (row.get("type_symbol") or "").upper()
            chain = row.get("auth_asym_id") or row.get("label_asym_id") or ""
            try:
                coord = [float(row["Cartn_x"]), float(row["Cartn_y"]), float(row["Cartn_z"])]
            except (KeyError, ValueError):
                continue
            if chain == "A" and atom_name == "CA" and atom_type == "C":
                ca_coords.append(coord)
            elif chain in {"B", "C", "D"} and atom_type != "H":
                ligand_coords.append(coord)
    return np.asarray(ca_coords, dtype=float), np.asarray(ligand_coords, dtype=float)


@lru_cache(maxsize=1024)
def _cached_atom_site_arrays(cif_path: str) -> tuple[np.ndarray, np.ndarray]:
    return _atom_site_arrays(cif_path)


def ligand_rmsd_to_medoid(mobile_cif: str, medoid_cif: str) -> float:
    if not mobile_cif or not medoid_cif or not Path(mobile_cif).exists() or not Path(medoid_cif).exists():
        return np.nan
    if mobile_cif == medoid_cif:
        return 0.0
    reference_ca, reference_ligand = _cached_atom_site_arrays(medoid_cif)
    mobile_ca, mobile_ligand = _atom_site_arrays(mobile_cif)
    if reference_ca.shape != mobile_ca.shape or reference_ca.size == 0:
        return np.nan
    if reference_ligand.shape != mobile_ligand.shape or reference_ligand.size == 0:
        return np.nan
    aligned_ligand = _kabsch_align(mobile_ca, reference_ca, mobile_ligand)
    squared = np.sum((aligned_ligand - reference_ligand) ** 2, axis=1)
    return float(np.sqrt(np.mean(squared)))


def qc_attrition_tables(r: Renderer) -> tuple[pd.DataFrame, pd.DataFrame]:
    files = sorted(r.results_dir.glob("*/shards/shard_*/qc_attrition_table.tsv"))
    frames = []
    for path in files:
        frame = pd.read_csv(path, sep="\t")
        frame["qc_source_file"] = str(path.relative_to(r.results_dir))
        frames.append(frame)
    if not frames:
        raise FileNotFoundError("no shard qc_attrition_table.tsv files found")
    qc = cats(pd.concat(frames, ignore_index=True))
    qc = num(qc, ["n_generated", "n_prep_error", "n_stage1_hard_fail", "n_stage1_pass"])
    category_rows = []
    for _, row in qc.iterrows():
        counts = {"Ligand too far": 0.0, "Privateer": 0.0, "PoseBusters": 0.0, "Other hard QC": 0.0}
        raw = row.get("hard_fail_reason_counts", "{}")
        try:
            parsed = json.loads(raw) if str(raw).strip() else {}
        except json.JSONDecodeError:
            parsed = {}
        for key, value in parsed.items():
            label = str(key).lower()
            if "ligand too far" in label:
                counts["Ligand too far"] += float(value)
            elif "privateer" in label:
                counts["Privateer"] += float(value)
            elif "posebusters" in label:
                counts["PoseBusters"] += float(value)
            else:
                counts["Other hard QC"] += float(value)
        total_reasons = sum(counts.values())
        hard_fail = float(row.get("n_stage1_hard_fail") or 0.0)
        if hard_fail and total_reasons:
            for key, value in counts.items():
                counts[key] = hard_fail * value / total_reasons
        elif hard_fail:
            counts["Other hard QC"] = hard_fail
        for stage in ["Ligand too far", "Privateer", "PoseBusters", "Other hard QC"]:
            category_rows.append(
                {
                    "condition_id": row.get("condition_id"),
                    "construct_type": row.get("construct_type"),
                    "substrate_class": row.get("substrate_class"),
                    "dp": row.get("dp"),
                    "stage": stage,
                    "n": counts[stage],
                }
            )
        category_rows.append({**{k: row.get(k) for k in ["condition_id", "construct_type", "substrate_class", "dp"]}, "stage": "QC pass", "n": row.get("n_stage1_pass", 0)})
        if float(row.get("n_prep_error") or 0.0) > 0:
            category_rows.append({**{k: row.get(k) for k in ["condition_id", "construct_type", "substrate_class", "dp"]}, "stage": "Preparation error", "n": row.get("n_prep_error", 0)})
    return qc, cats(pd.DataFrame(category_rows))


def figure1(r: Renderer) -> None:
    fid, title = "figure1_dataset_workflow.png", "Dataset"
    sources = ["combined/condition_table.tsv", "*/shards/shard_*/qc_attrition_table.tsv"]
    try:
        df = cats(r.tsv(sources[0]))
        df = num(df, ["n_generated", "n_stage1_pass", "n_ifp_success", "n_ifp_clustered", "n_cluster_rows"])
        qc, qc_long = qc_attrition_tables(r)
        fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.4), constrained_layout=True, gridspec_kw={"width_ratios": [0.95, 1.75]})
        ax1, ax2 = axes
        counts_by_construct = df.groupby("construct_type", observed=True).agg(
            proteins=("protein_id", "nunique"),
            conditions=("condition_id", "nunique"),
            poses=("n_generated", "sum"),
        ).reset_index()
        long = counts_by_construct.melt(id_vars=["construct_type"], value_vars=["proteins", "conditions", "poses"], var_name="unit", value_name="n")
        sns.barplot(data=long, x="unit", y="n", hue="construct_type", palette=CONSTRUCT_PALETTE, errorbar=None, ax=ax1)
        ax1.set_title("Dataset size")
        ax1.set_xlabel("")
        ax1.set_ylabel("Count")
        ax1.set_yscale("log")
        ax1.legend(title="Construct", frameon=False, loc="upper left")
        compact_margins(ax1, x=0.02)
        panel_label(ax1, "A")
        fail_order = ["Ligand too far", "Privateer", "PoseBusters"]
        colors = {
            "Ligand too far": "#8c6d31",
            "Privateer": "#756bb1",
            "PoseBusters": "#c44e52",
            "QC pass": "#3a923a",
            "IFP success": "#4c78a8",
            "Clusterable": "#72b7b2",
            "Geometry scorable": "#f58518",
        }
        construct_labels = {"domain_only": "Catalytic domain", "full_length": "Full-length"}
        failure = qc_long[qc_long["stage"].isin(fail_order)].groupby(["construct_type", "stage"], observed=True)["n"].sum().reset_index()
        construct_counts = df.groupby("construct_type", observed=True).agg(
            Generated=("n_generated", "sum"),
            **{
                "QC pass": ("n_stage1_pass", "sum"),
                "IFP success": ("n_ifp_success", "sum"),
                "Clusterable": ("n_ifp_clustered", "sum"),
                "Geometry scorable": ("n_cluster_rows", "sum"),
            },
        ).reset_index()
        rows = []
        for _, row in construct_counts.iterrows():
            construct = str(row["construct_type"])
            generated = float(row["Generated"])
            ligand_too_far = float(failure[(failure["construct_type"].astype(str) == construct) & (failure["stage"] == "Ligand too far")]["n"].sum())
            privateer = float(failure[(failure["construct_type"].astype(str) == construct) & (failure["stage"] == "Privateer")]["n"].sum())
            posebusters = float(failure[(failure["construct_type"].astype(str) == construct) & (failure["stage"] == "PoseBusters")]["n"].sum())
            cumulative = [
                ("Generated", generated),
                ("Distance gate", generated - ligand_too_far),
                ("Privateer", generated - ligand_too_far - privateer),
                ("PoseBusters", generated - ligand_too_far - privateer - posebusters),
                ("QC pass", float(row["QC pass"])),
                ("IFP success", float(row["IFP success"])),
                ("Clusterable", float(row["Clusterable"])),
                ("Geometry scorable", float(row["Geometry scorable"])),
            ]
            previous = math.inf
            for stage, value in cumulative:
                value = max(0.0, min(float(value), previous))
                rows.append({"construct_type": construct, "stage": stage, "n": value})
                previous = value
        attr = pd.DataFrame(rows)
        attr["n"] = attr["n"].round().astype(int)
        stage_order = ["Generated", "Distance gate", "Privateer", "PoseBusters", "QC pass", "IFP success", "Clusterable", "Geometry scorable"]
        stage_display = {
            "Generated": "Generated",
            "Distance gate": "Distance\ngate",
            "Privateer": "Privateer",
            "PoseBusters": "PoseBusters",
            "QC pass": "QC pass",
            "IFP success": "IFP\nsuccess",
            "Clusterable": "Clusterable",
            "Geometry scorable": "Geometry\nscoreable",
        }
        y = np.arange(len(stage_order))
        offsets = {"domain_only": -0.18, "full_length": 0.18}
        for construct in [c for c in CONSTRUCT_ORDER if c in set(attr["construct_type"])]:
            sub = attr[attr["construct_type"].eq(construct)].set_index("stage").reindex(stage_order)
            ax2.barh(y + offsets[construct], sub["n"], height=0.32, label=construct_labels.get(construct, construct), color=CONSTRUCT_PALETTE[construct])
        hard_start = stage_order.index("Distance gate") - 0.48
        hard_end = stage_order.index("PoseBusters") + 0.48
        ax2.axhspan(hard_start, hard_end, color="#c44e52", alpha=0.11, zorder=0)
        ax2.axhline(hard_start, color="#c44e52", linewidth=1.0, alpha=0.7)
        ax2.axhline(hard_end, color="#c44e52", linewidth=1.0, alpha=0.7)
        ax2.set_yticks(y)
        ax2.set_yticklabels([stage_display[stage] for stage in stage_order])
        ax2.invert_yaxis()
        ax2.set_title("Cumulative pass count through hard QC and downstream steps")
        ax2.set_xlabel("Count passing step")
        ax2.set_ylabel("")
        ax2.legend(title="Construct", frameon=False, loc="lower right")
        compact_margins(ax2, x=0.015)
        panel_label(ax2, "B")
        fig.suptitle(title, fontsize=15, weight="bold")
        r.table(counts_by_construct, "figure1_dataset_summary.tsv")
        r.table(attr, "figure1_qc_failure_categories.tsv")
        save_clean(r, fig, fid, title, sources, "all conditions; QC failure categories from shard QC tables", len(qc))
    except Exception as exc:
        r.skip(fid, title, sources, str(exc))


def figure2(r: Renderer) -> None:
    fid, title = "figure2_prediction_quality_clusterability.png", "Clusterability and binding-mode reproducibility"
    source = "combined/condition_table.tsv"
    try:
        df = cats(r.tsv(source))
        metrics = ["top_cluster_occupancy", "n_clusters", "cluster_entropy", "noise_fraction"]
        df = num(df, metrics)
        df["valid_cluster"] = bools(df["any_valid_cluster"])
        fig, axes = plt.subplots(5, 1, figsize=(7.9, 15.6), constrained_layout=True)
        rate = df.groupby(["construct_type", "substrate_dp"], observed=True).agg(fraction=("valid_cluster", "mean"), n=("condition_id", "nunique")).reset_index()
        ax_c, ax_a, ax_b, ax_d, ax_e = axes
        counts = df["n_clusters"].dropna().astype(int)
        cluster_count = counts.value_counts().sort_index()
        cluster_count = cluster_count.reindex(range(0, int(cluster_count.index.max()) + 1), fill_value=0)
        ax_c.bar(cluster_count.index, cluster_count.values, color="#587a9d")
        ax_c.set_yscale("log", base=2)
        ax_c.set_ylim(0.5, max(cluster_count.max() * 1.9, 2))
        for x_value, y_value in cluster_count.items():
            ax_c.text(x_value, max(y_value, 0.55) * 1.08, str(int(y_value)), ha="center", va="bottom", fontsize=10)
        ax_c.set_xlabel("Number of clusters")
        ax_c.set_ylabel("Conditions")
        ax_c.set_xticks(cluster_count.index)
        panel_label(ax_c, "A")
        sns.pointplot(data=rate, x="substrate_dp", y="fraction", hue="construct_type", palette=CONSTRUCT_PALETTE, ax=ax_a)
        ax_a.set_ylim(-0.03, 1.03)
        substrate_dp_ticklabels(ax_a, rotation=35)
        ax_a.set_xlabel("")
        ax_a.set_ylabel("Fractions of conditions clusterable")
        annotate(ax_a, f"n conditions={df['condition_id'].nunique()}")
        panel_label(ax_a, "B")
        occ = df[df["n_clusters"].fillna(0).gt(1) & df["top_cluster_occupancy"].gt(0)].copy()
        sns.histplot(data=occ, x="top_cluster_occupancy", hue="construct_type", bins=25, stat="count", element="step", palette=CONSTRUCT_PALETTE, ax=ax_b)
        ax_b.set_xlabel("Top-cluster occupancy")
        ax_b.set_ylabel("Count")
        if ax_b.get_legend() is not None:
            ax_b.get_legend().remove()
        ax_b.legend(
            handles=[Patch(facecolor=CONSTRUCT_PALETTE[construct], edgecolor=CONSTRUCT_PALETTE[construct], alpha=0.35, label=CONSTRUCT_DISPLAY[construct]) for construct in CONSTRUCT_ORDER],
            title="Construct",
            loc="upper left",
            frameon=False,
        )
        compact_margins(ax_b, x=0.005)
        annotate(ax_b, "multi-cluster conditions only (n clusters > 1)")
        panel_label(ax_b, "C")
        ent = finite(df[df["valid_cluster"]].copy(), ["cluster_entropy"])
        sns.violinplot(data=ent, x="substrate_class", y="cluster_entropy", hue="construct_type", cut=0, inner="quartile", palette=CONSTRUCT_PALETTE, ax=ax_d)
        ax_d.set_xlabel("Substrate")
        ax_d.set_ylabel("Normalized entropy")
        ax_d.legend(title="Construct", loc="upper right", frameon=False)
        compact_margins(ax_d, x=0.02)
        panel_label(ax_d, "D")
        df = num(df, ["n_ifp_clustered"])
        formal_allowed = bools(df["formal_clustering_allowed"]) if "formal_clustering_allowed" in df.columns else pd.Series(False, index=df.index)
        noise = finite(df[formal_allowed | df["n_ifp_clustered"].fillna(0).gt(0)].copy(), ["noise_fraction"])
        centered_violin_boxplot(
            ax_e,
            noise,
            value_col="noise_fraction",
            y_label="Noise fraction",
            light_palette={"domain_only": "#b8cce4", "full_length": "#e8bd89"},
            dark_palette={"domain_only": "#2f5f8f", "full_length": "#a85f17"},
        )
        ax_e.set_ylim(-0.03, 1.03)
        handles = [Patch(facecolor=color, edgecolor="#303030", label=CONSTRUCT_DISPLAY[construct]) for construct, color in {"domain_only": "#2f5f8f", "full_length": "#a85f17"}.items()]
        ax_e.legend(handles=handles, title="Construct", frameon=False, loc="upper right")
        compact_margins(ax_e, x=0.02)
        panel_label(ax_e, "E")
        fig.suptitle(title, fontsize=15, weight="bold")
        r.table(rate, "figure2_valid_cluster_fraction.tsv")
        r.table(cluster_count.rename("conditions").reset_index().rename(columns={"index": "n_clusters"}), "figure2_cluster_count_distribution.tsv")
        save_clean(r, fig, fid, title, [source], "condition-level summaries", len(df))
    except Exception as exc:
        r.skip(fid, title, [source], str(exc))


def figure3(r: Renderer) -> None:
    fid, title = "figure3_c1_c4_geometry.png", "C1/C4 geometry of predicted binding modes"
    sources = ["combined/cluster_table.tsv", "combined/condition_table.tsv"]
    try:
        clusters = cats(add_activity(r, r.tsv(sources[0])))
        cond = cats(add_activity(r, r.tsv(sources[1])))
        cluster_cols = [
            "Cu_C1_distance_median", "Cu_C4_distance_median", "oxyl_H_C1_distance_median",
            "oxyl_H_C4_distance_median", "occupancy", "Cu_oxyl_H_C1_angle_median", "Cu_oxyl_H_C4_angle_median",
        ]
        clusters = finite(num(clusters, cluster_cols), ["Cu_C1_distance_median", "Cu_C4_distance_median"])
        clusters["symbol_group"] = clusters["construct_type"].astype(str).map({"domain_only": "Catalytic domain", "full_length": "Full-length"}) + " | " + clusters["activity_evidence"].map({ACTIVITY_ORDER[0]: "Demonstrated activity", ACTIVITY_ORDER[1]: "No demonstrated activity"})
        fig, axes = plt.subplot_mosaic([["A", "A"], ["B", "B"], ["D1", "D2"], ["C", "C"]], figsize=(13.2, 15.3), constrained_layout=True)
        sns.scatterplot(data=clusters, x="Cu_C1_distance_median", y="Cu_C4_distance_median", hue="substrate_class", style="symbol_group", style_order=MARKER_ORDER, markers=MARKERS, size="occupancy", sizes=(18, 150), alpha=0.75, palette=SUBSTRATE_PALETTE, ax=axes["A"])
        values = clusters[["Cu_C1_distance_median", "Cu_C4_distance_median"]].to_numpy()
        lo = max(3.0, math.floor(float(np.nanpercentile(values, 1))))
        hi = math.ceil(float(np.nanpercentile(values, 99)))
        axes["A"].plot([lo, hi], [lo, hi], "--", color="#555555", linewidth=1)
        axes["A"].set_xlim(lo, hi)
        axes["A"].set_ylim(lo, hi)
        axes["A"].set_title("Cu to reactive carbon")
        axes["A"].set_xlabel("Cu-C1 distance (Å)")
        axes["A"].set_ylabel("Cu-C4 distance (Å)")
        axis_break_marks(axes["A"])
        axes["A"].legend(bbox_to_anchor=(1.005, 1), loc="upper left", fontsize=LEGEND_FONTSIZE, frameon=False, borderaxespad=0.0)
        panel_label(axes["A"], "A")
        ox = finite(clusters, ["oxyl_H_C1_distance_median", "oxyl_H_C4_distance_median"])
        sns.scatterplot(data=ox, x="oxyl_H_C1_distance_median", y="oxyl_H_C4_distance_median", hue="substrate_class", style="symbol_group", style_order=MARKER_ORDER, markers=MARKERS, size="occupancy", sizes=(18, 145), alpha=0.75, palette=SUBSTRATE_PALETTE, ax=axes["B"])
        axes["B"].axvspan(1.5, 4.0, color="#e6e6e6", alpha=0.5, zorder=0)
        axes["B"].axhspan(1.5, 4.0, color="#e6e6e6", alpha=0.5, zorder=0)
        axes["B"].axvline(2.1, color="#555555", linestyle=":", linewidth=1)
        axes["B"].axhline(2.1, color="#555555", linestyle=":", linewidth=1)
        axes["B"].set_title("Oxyl-H distances")
        axes["B"].set_xlabel("Oxyl-H(C1) distance (Å)")
        axes["B"].set_ylabel("Oxyl-H(C4) distance (Å)")
        axis_break_marks(axes["B"])
        axes["B"].legend(bbox_to_anchor=(1.005, 1), loc="upper left", fontsize=LEGEND_FONTSIZE, frameon=False, borderaxespad=0.0)
        panel_label(axes["B"], "B")
        metrics = ["occupancy_weighted_c1_plausible_fraction", "occupancy_weighted_c4_plausible_fraction"]
        cond = num(cond, metrics + ["occupancy_weighted_Cu_oxyl_H_C1_angle_median", "occupancy_weighted_Cu_oxyl_H_C4_angle_median", "occupancy_weighted_Cu_C1_distance_median", "occupancy_weighted_Cu_C4_distance_median"])
        geom_cond = cond[bools(cond["any_valid_cluster"])].copy()
        long = geom_cond.melt(id_vars=["protein_id", "condition_id", "activity_evidence"], value_vars=metrics, var_name="site", value_name="fraction")
        long["site"] = long["site"].map({metrics[0]: "C1", metrics[1]: "C4"})
        long = finite(long, ["fraction"])
        sns.boxplot(data=long, x="activity_evidence", y="fraction", hue="site", order=ACTIVITY_ORDER, showfliers=False, ax=axes["C"])
        sns.stripplot(data=long, x="activity_evidence", y="fraction", hue="site", order=ACTIVITY_ORDER, dodge=True, color="#222222", alpha=0.2, size=2, ax=axes["C"], legend=False)
        axes["C"].set_title("Plausible geometry fraction")
        axes["C"].set_xlabel("")
        axes["C"].set_ylabel("Occupancy-weighted fraction")
        axes["C"].tick_params(axis="x", rotation=15)
        n_rows = (
            long.groupby(["activity_evidence", "site"], observed=True)["condition_id"]
            .nunique()
            .reindex(pd.MultiIndex.from_product([ACTIVITY_ORDER, ["C1", "C4"]]), fill_value=0)
        )
        axes["C"].set_xticks([0, 1])
        axes["C"].set_xticklabels([
            f"Demonstrated\nactivity\nC1/C4 n={int(n_rows[(ACTIVITY_ORDER[0], 'C1')])}/{int(n_rows[(ACTIVITY_ORDER[0], 'C4')])}",
            f"No demonstrated\nactivity\nC1/C4 n={int(n_rows[(ACTIVITY_ORDER[1], 'C1')])}/{int(n_rows[(ACTIVITY_ORDER[1], 'C4')])}",
        ])
        axes["C"].legend(title="Site", frameon=False, loc="upper right")
        panel_label(axes["C"], "E")
        angle_metrics = ["occupancy_weighted_Cu_oxyl_H_C1_angle_median", "occupancy_weighted_Cu_oxyl_H_C4_angle_median"]
        angle = cond.melt(id_vars=["activity_evidence"], value_vars=angle_metrics, var_name="metric", value_name="value")
        angle["metric"] = angle["metric"].map(nice)
        sns.boxplot(data=finite(angle, ["value"]), y="metric", x="value", hue="activity_evidence", hue_order=ACTIVITY_ORDER, showfliers=False, palette=ACTIVITY_PALETTE, ax=axes["D1"])
        axes["D1"].set_yticks(axes["D1"].get_yticks())
        axes["D1"].set_yticklabels(["Cu oxyl H\nC1 angle\nmedian", "Cu oxyl H\nC4 angle\nmedian"])
        axes["D1"].set_title("Attack angle metrics")
        axes["D1"].set_xlabel("Degrees")
        axes["D1"].set_ylabel("")
        handles, labels = axes["D1"].get_legend_handles_labels()
        if axes["D1"].get_legend() is not None:
            axes["D1"].get_legend().remove()
        panel_label(axes["D1"], "C")
        distance_metrics = ["occupancy_weighted_Cu_C1_distance_median", "occupancy_weighted_Cu_C4_distance_median"]
        dist = cond.melt(id_vars=["activity_evidence"], value_vars=distance_metrics, var_name="metric", value_name="value")
        dist["metric"] = dist["metric"].map(nice)
        sns.boxplot(data=finite(dist, ["value"]), y="metric", x="value", hue="activity_evidence", hue_order=ACTIVITY_ORDER, showfliers=False, palette=ACTIVITY_PALETTE, ax=axes["D2"])
        axes["D2"].set_yticks(axes["D2"].get_yticks())
        axes["D2"].set_yticklabels(["Cu C1\ndistance\nmedian", "Cu C4\ndistance\nmedian"])
        axes["D2"].set_title("Cu-distance metrics")
        axes["D2"].set_xlabel("Distance (Å)")
        axes["D2"].set_ylabel("")
        if axes["D2"].get_legend() is not None:
            axes["D2"].get_legend().remove()
        axes["D2"].legend(handles, labels, title="Activity evidence", frameon=False, loc="upper left", bbox_to_anchor=(1.005, 1), fontsize=LEGEND_FONTSIZE, borderaxespad=0.0)
        panel_label(axes["D2"], "D")
        fig.suptitle(title, fontsize=15, weight="bold")
        save_clean(r, fig, fid, title, sources + ["input_data/metadata_final_ec_fixed.tsv"], "cluster medians and EC-activity-annotated conditions", len(clusters))
    except Exception as exc:
        r.skip(fid, title, sources, str(exc))


def figure3_angles_by_substrate(r: Renderer) -> None:
    fid, title = "figure3_angles_by_substrate.png", "Angles by substrate"
    source = "combined/condition_table.tsv"
    try:
        df = cats(r.tsv(source))
        meta = r.tsv("input_data/metadata_final_ec_fixed.tsv", base="analyse")
        protein_col = "UniProt_ID" if "UniProt_ID" in meta.columns else "Entry"
        ec_col = "EC_Number" if "EC_Number" in meta.columns else "EC number"
        meta["protein_id"] = meta[protein_col].astype(str)
        meta["active_substrate_set"] = meta[ec_col].map(ec_active_substrates_with_other)
        meta["experimental_substrate_group"] = meta["active_substrate_set"].map(demonstrated_activity_group)
        meta["active_substrates_with_other"] = meta["active_substrate_set"].map(lambda active: ",".join(sorted(active)))
        df = df.merge(meta[["protein_id", "experimental_substrate_group", "active_substrates_with_other"]], on="protein_id", how="left")
        df = df[bools(df["any_valid_cluster"])].copy()
        active_sets = df["active_substrates_with_other"].fillna("").map(lambda value: set(filter(None, str(value).split(","))))
        condition_substrates = df["substrate_class"].astype(str).replace({"amylose": "starch"})
        df["activity_evidence"] = [
            ACTIVITY_ORDER[0] if substrate in active else ACTIVITY_ORDER[1]
            for substrate, active in zip(condition_substrates, active_sets, strict=False)
        ]
        angle_metrics = ["occupancy_weighted_Cu_oxyl_H_C1_angle_median", "occupancy_weighted_Cu_oxyl_H_C4_angle_median"]
        df = num(df, angle_metrics)
        long = df.melt(
            id_vars=["protein_id", "condition_id", "experimental_substrate_group", "activity_evidence"],
            value_vars=angle_metrics,
            var_name="site",
            value_name="angle_deg",
        )
        long["site"] = long["site"].map({angle_metrics[0]: "C1", angle_metrics[1]: "C4"})
        long = finite(long, ["angle_deg"])
        substrate_order = [group for group in DEMONSTRATED_ACTIVITY_GROUP_ORDER if group in set(long["experimental_substrate_group"].astype(str))]
        label_map = {
            "Chitin": "Chitin",
            "Cellulose": "Cellulose",
            "Mixed cellulose/chitin": "Mixed\ncellulose/\nchitin",
            "Starch": "Starch",
            "Other": "Other",
        }
        fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8), constrained_layout=True, sharey=True)
        for ax, site in zip(axes, ["C1", "C4"], strict=False):
            sub = long[long["site"].eq(site)]
            sns.boxplot(
                data=sub,
                x="experimental_substrate_group",
                y="angle_deg",
                hue="activity_evidence",
                order=substrate_order,
                hue_order=ACTIVITY_ORDER,
                showfliers=False,
                palette=ACTIVITY_PALETTE,
                width=0.56,
                ax=ax,
            )
            sns.stripplot(
                data=sub,
                x="experimental_substrate_group",
                y="angle_deg",
                hue="activity_evidence",
                order=substrate_order,
                hue_order=ACTIVITY_ORDER,
                dodge=True,
                color="#222222",
                alpha=0.16,
                size=2,
                ax=ax,
                legend=False,
            )
            all_summary = sub.groupby("experimental_substrate_group", observed=True)["angle_deg"].agg(mean="mean", std="std", n="count").reindex(substrate_order)
            x = np.arange(len(substrate_order))
            ax.errorbar(
                x,
                all_summary["mean"].to_numpy(),
                yerr=all_summary["std"].fillna(0).to_numpy(),
                fmt="D",
                color="#000000",
                ecolor="#000000",
                markersize=4.8,
                linewidth=1.0,
                capsize=3,
                label="All mean ± SD",
                zorder=5,
            )
            ax.set_xlabel("Experimental substrate activity")
            ax.set_ylabel("Angle (degrees)" if site == "C1" else "")
            ax.set_xticks(x)
            ax.set_xticklabels([label_map.get(label, label) for label in substrate_order])
            ax.set_title(f"{site} angle")
            handles, labels = ax.get_legend_handles_labels()
            if ax.get_legend() is not None:
                ax.get_legend().remove()
            if site == "C4":
                keep = []
                seen = set()
                for handle, label in zip(handles, labels, strict=False):
                    if label not in seen and label in [*ACTIVITY_ORDER, "All mean ± SD"]:
                        keep.append((handle, label))
                        seen.add(label)
                ax.legend([item[0] for item in keep], [item[1] for item in keep], title="Group", frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5), borderaxespad=0.0)
            compact_margins(ax, x=0.03)
            panel_label(ax, site)
        r.table(long, "figure3_angles_by_substrate_values.tsv")
        save_clean(r, fig, fid, title, [source, "input_data/metadata_final_ec_fixed.tsv"], "valid-cluster conditions with finite occupancy-weighted Cu oxyl H angle medians", len(long))
    except Exception as exc:
        r.skip(fid, title, [source], str(exc))


def figure4(r: Renderer) -> None:
    fid, title = "figure4_activity_evidence_paired.png", "Demonstrated activity versus no demonstrated activity substrates"
    source = "combined/condition_table.tsv"
    try:
        df = cats(add_activity(r, r.tsv(source)))
        metrics = [
            "occupancy_weighted_c1_plausible_fraction",
            "occupancy_weighted_c4_plausible_fraction",
            "top_cluster_occupancy",
            "cluster_entropy",
            "occupancy_weighted_oxyl_H_score_C1_median",
            "occupancy_weighted_oxyl_H_score_C4_median",
        ]
        df = num(df, metrics)
        df = df[bools(df["any_valid_cluster"])].copy()
        long = df.melt(id_vars=["protein_id", "condition_id", "activity_evidence"], value_vars=metrics, var_name="metric", value_name="value")
        long["metric"] = long["metric"].map(nice)
        long = finite(long, ["value"])
        metric_order = list(long["metric"].drop_duplicates())
        fig, axes = plt.subplots(2, 3, figsize=(11.4, 7.2), constrained_layout=True)
        for idx, (ax, metric) in enumerate(zip(axes.flat, metric_order, strict=False)):
            sub = long[long["metric"].eq(metric)]
            for _, protein in sub.groupby("protein_id", observed=True):
                if protein["activity_evidence"].nunique() == 2:
                    protein = protein.sort_values("activity_evidence", key=lambda s: s.map({ACTIVITY_ORDER[0]: 0, ACTIVITY_ORDER[1]: 1}))
                    ax.plot(protein["activity_evidence"], protein["value"], color="#777777", alpha=0.12, linewidth=0.7)
            sns.boxplot(data=sub, x="activity_evidence", y="value", hue="activity_evidence", order=ACTIVITY_ORDER, hue_order=ACTIVITY_ORDER, showfliers=False, width=0.42, palette=ACTIVITY_PALETTE, ax=ax, legend=False)
            sns.stripplot(data=sub, x="activity_evidence", y="value", order=ACTIVITY_ORDER, color="#222222", alpha=0.2, size=2, ax=ax)
            ax.tick_params(axis="x", rotation=15)
            ax.set_xlabel("")
            ax.set_ylabel(metric)
            ax.set_xticks([0, 1])
            counts = sub.groupby("activity_evidence", observed=True)["condition_id"].nunique()
            ax.set_xticklabels([
                f"Demonstrated activity\nn={int(counts.get(ACTIVITY_ORDER[0], 0))}",
                f"No demonstrated activity\nn={int(counts.get(ACTIVITY_ORDER[1], 0))}",
            ])
            wrap_ticklabels(ax, axis="x", width=14)
            compact_margins(ax, x=0.02)
            panel_label(ax, chr(ord("A") + idx))
        fig.suptitle(title, y=1.03, fontsize=15, weight="bold")
        summary = long.groupby(["metric", "activity_evidence"], observed=True).agg(n=("value", "size"), median=("value", "median"), proteins=("protein_id", "nunique")).reset_index()
        r.table(summary, "figure4_activity_metric_summary.tsv")
        save_clean(r, fig, fid, title, [source, "input_data/metadata_final_ec_fixed.tsv"], "EC-activity-annotated conditions", len(long))
    except Exception as exc:
        r.skip(fid, title, [source], str(exc))


def figure4_demonstrated_activity_substrate_boxplots(r: Renderer) -> None:
    fid = "figure4_demonstrated_activity_substrate_boxplots.png"
    title = "Demonstrated activity substrate distributions"
    source = "combined/condition_table.tsv"
    meta_source = "input_data/metadata_final_ec_fixed.tsv"
    try:
        df = cats(r.tsv(source))
        meta = r.tsv(meta_source, base="analyse")
        protein_col = "UniProt_ID" if "UniProt_ID" in meta.columns else "Entry"
        ec_col = "EC_Number" if "EC_Number" in meta.columns else "EC number"
        meta["protein_id"] = meta[protein_col].astype(str)
        meta["active_substrate_set"] = meta[ec_col].map(ec_active_substrates_with_other)
        meta["demonstrated_substrate_group"] = meta["active_substrate_set"].map(demonstrated_activity_group)
        meta["active_substrates_with_other"] = meta["active_substrate_set"].map(lambda active: ",".join(sorted(active)))
        df = df.merge(
            meta[["protein_id", "demonstrated_substrate_group", "active_substrates_with_other"]],
            on="protein_id",
            how="left",
        )
        metrics = [
            "occupancy_weighted_c1_plausible_fraction",
            "occupancy_weighted_c4_plausible_fraction",
            "top_cluster_occupancy",
            "cluster_entropy",
            "occupancy_weighted_oxyl_H_score_C1_median",
            "occupancy_weighted_oxyl_H_score_C4_median",
        ]
        df = num(df, metrics)
        df = df[bools(df["any_valid_cluster"])].copy()
        df["condition_substrate"] = df["substrate_class"].astype(str).replace({"amylose": "starch"})
        group_match = {
            "Chitin": df["condition_substrate"].eq("chitin"),
            "Cellulose": df["condition_substrate"].eq("cellulose"),
            "Mixed cellulose/chitin": df["condition_substrate"].isin(["cellulose", "chitin"]),
            "Starch": df["condition_substrate"].eq("starch"),
            "Other": df["demonstrated_substrate_group"].eq("Other"),
        }
        demonstrated_mask = pd.Series(False, index=df.index)
        for group, mask in group_match.items():
            demonstrated_mask |= df["demonstrated_substrate_group"].eq(group) & mask
        df = df[demonstrated_mask].copy()
        group_counts = df.groupby("demonstrated_substrate_group", observed=True)["condition_id"].nunique()
        group_order = [group for group in DEMONSTRATED_ACTIVITY_GROUP_ORDER if int(group_counts.get(group, 0)) >= 3]
        df = df[df["demonstrated_substrate_group"].isin(group_order)].copy()
        group_counts = df.groupby("demonstrated_substrate_group", observed=True)["condition_id"].nunique()
        group_display = {
            "Chitin": "Chitin",
            "Cellulose": "Cellulose",
            "Mixed cellulose/chitin": "Mixed\ncellulose/\nchitin",
            "Starch": "Starch",
            "Other": "Other",
        }
        long = df.melt(
            id_vars=["protein_id", "condition_id", "demonstrated_substrate_group", "active_substrates_with_other"],
            value_vars=metrics,
            var_name="metric",
            value_name="value",
        )
        long["metric"] = long["metric"].map(nice)
        long = finite(long, ["value"])
        metric_order = list(long["metric"].drop_duplicates())
        test_rows: list[dict[str, object]] = []
        auroc_rows: list[dict[str, object]] = []
        for metric_idx, metric in enumerate(metric_order):
            sub = long[long["metric"].eq(metric)].copy()
            finite_counts = sub.groupby("demonstrated_substrate_group", observed=True)["value"].size()
            tested_groups = [group for group in group_order if int(finite_counts.get(group, 0)) >= 3]
            tested = sub[sub["demonstrated_substrate_group"].isin(tested_groups)].copy()
            if len(tested_groups) >= 2:
                values = tested["value"].astype(float).to_numpy()
                labels = tested["demonstrated_substrate_group"].astype(str).to_numpy()
                h_value, p_value = permutation_kruskal(values, labels, n_permutations=10000, seed=20260604 + metric_idx)
            else:
                h_value, p_value = np.nan, np.nan
            test_rows.append(
                {
                    "metric": metric,
                    "n_groups_tested": len(tested_groups),
                    "groups_tested": ",".join(tested_groups),
                    "kruskal_wallis_h": h_value,
                    "permutation_p": p_value,
                    "n_permutations": 10000 if len(tested_groups) >= 2 else 0,
                }
            )
            for group_a, group_b in itertools.combinations(tested_groups, 2):
                values_a = sub.loc[sub["demonstrated_substrate_group"].eq(group_a), "value"].astype(float).to_numpy()
                values_b = sub.loc[sub["demonstrated_substrate_group"].eq(group_b), "value"].astype(float).to_numpy()
                combined = np.concatenate([values_a, values_b])
                labels = np.concatenate([np.ones(len(values_a), dtype=int), np.zeros(len(values_b), dtype=int)])
                auroc = roc_auc_score(labels, combined) if len(np.unique(labels)) == 2 else np.nan
                auroc_rows.append(
                    {
                        "metric": metric,
                        "group_a_positive": group_a,
                        "group_b_negative": group_b,
                        "n_group_a": len(values_a),
                        "n_group_b": len(values_b),
                        "auroc_group_a_higher": auroc,
                        "auc_centered": auroc - 0.5 if np.isfinite(auroc) else np.nan,
                    }
                )
        q_values = benjamini_hochberg([float(row["permutation_p"]) for row in test_rows])
        for row, q_value in zip(test_rows, q_values, strict=False):
            row["fdr_q"] = q_value
        test_df = pd.DataFrame(test_rows)
        q_lookup = dict(zip(test_df["metric"], test_df["fdr_q"], strict=False)) if len(test_df) else {}
        fig, axes = plt.subplots(2, 3, figsize=(11.4, 7.2), constrained_layout=True)
        for idx, (ax, metric) in enumerate(zip(axes.flat, metric_order, strict=False)):
            sub = long[long["metric"].eq(metric)]
            sns.boxplot(data=sub, x="demonstrated_substrate_group", y="value", order=group_order, showfliers=False, width=0.46, color="#6f91b6", ax=ax)
            sns.stripplot(data=sub, x="demonstrated_substrate_group", y="value", order=group_order, color="#222222", alpha=0.18, size=2, ax=ax)
            ax.set_xlabel("")
            ax.set_ylabel(metric)
            ax.set_xticks(range(len(group_order)))
            counts = sub.groupby("demonstrated_substrate_group", observed=True)["condition_id"].nunique()
            ax.set_xticklabels([f"{group_display.get(group, group)}\nn={int(counts.get(group, 0))}" for group in group_order])
            if idx in {0, 1, 2, 4, 5}:
                y_min, y_max = ax.get_ylim()
                ax.set_ylim(y_min, y_max + 0.12 * (y_max - y_min))
            q_value = q_lookup.get(metric, np.nan)
            annotate(ax, f"perm. q={q_value:.3g}" if np.isfinite(q_value) else "perm. q=NA")
            wrap_ticklabels(ax, axis="x", width=13)
            compact_margins(ax, x=0.02)
            panel_label(ax, chr(ord("A") + idx))
        fig.suptitle(title, y=1.03, fontsize=15, weight="bold")
        summary = (
            long.groupby(["metric", "demonstrated_substrate_group"], observed=True)
            .agg(
                n_values=("value", "size"),
                median=("value", "median"),
                q1=("value", lambda values: values.quantile(0.25)),
                q3=("value", lambda values: values.quantile(0.75)),
                proteins=("protein_id", "nunique"),
                conditions=("condition_id", "nunique"),
            )
            .reset_index()
        )
        summary["x_axis_label"] = summary.apply(
            lambda row: f"{group_display.get(row['demonstrated_substrate_group'], row['demonstrated_substrate_group'])}\nn={int(row['conditions'])}",
            axis=1,
        )
        r.table(summary, "figure4_demonstrated_activity_substrate_boxplot_summary.tsv")
        r.table(test_df, "figure4_demonstrated_activity_substrate_distribution_tests.tsv")
        r.table(pd.DataFrame(auroc_rows), "figure4_demonstrated_activity_substrate_pairwise_auroc.tsv")
        save_clean(r, fig, fid, title, [source, meta_source], "EC-demonstrated substrate groups; groups with fewer than three unique condition_id rows removed", len(long))
    except Exception as exc:
        r.skip(fid, title, [source], str(exc))


def contact_fraction_difference_by_substrate(df: pd.DataFrame, contact_metrics: list[str]) -> pd.DataFrame:
    rows = []
    group_specs = [(s, df[df["substrate_class"].astype(str).eq(s)]) for s in SUBSTRATE_ORDER if s in set(df["substrate_class"].astype(str))]
    group_specs.append(("All", df))
    for substrate, sub in group_specs:
        grouped = sub.groupby("activity_evidence", observed=True)[contact_metrics].agg(["mean", "std", "count"])
        for metric in contact_metrics:
            active_mean = grouped.loc[ACTIVITY_ORDER[0], (metric, "mean")] if ACTIVITY_ORDER[0] in grouped.index else np.nan
            inactive_mean = grouped.loc[ACTIVITY_ORDER[1], (metric, "mean")] if ACTIVITY_ORDER[1] in grouped.index else np.nan
            active_sd = grouped.loc[ACTIVITY_ORDER[0], (metric, "std")] if ACTIVITY_ORDER[0] in grouped.index else np.nan
            inactive_sd = grouped.loc[ACTIVITY_ORDER[1], (metric, "std")] if ACTIVITY_ORDER[1] in grouped.index else np.nan
            active_n = grouped.loc[ACTIVITY_ORDER[0], (metric, "count")] if ACTIVITY_ORDER[0] in grouped.index else 0
            inactive_n = grouped.loc[ACTIVITY_ORDER[1], (metric, "count")] if ACTIVITY_ORDER[1] in grouped.index else 0
            if np.isfinite(active_mean) and np.isfinite(inactive_mean):
                rows.append(
                    {
                        "substrate_class": substrate,
                        "contact_type": metric,
                        "delta": active_mean - inactive_mean,
                        "sd_demonstrated": active_sd,
                        "sd_no_demonstrated": inactive_sd,
                        "n_demonstrated": int(active_n),
                        "n_no_demonstrated": int(inactive_n),
                    }
                )
    diff = pd.DataFrame(rows)
    if diff.empty:
        return diff
    label_map = {
        "aromatic_contact_fraction": "Aromatic",
        "polar_contact_fraction": "Polar",
        "charged_contact_fraction": "Charged",
        "hbond_contact_fraction": "Hydrogen bond",
    }
    diff["contact_type_label"] = diff["contact_type"].map(label_map)
    return diff


def plot_contact_fraction_difference_by_substrate(ax: plt.Axes, diff: pd.DataFrame, *, legend: bool) -> None:
    contact_order = ["Hydrogen bond", "Charged", "Polar", "Aromatic"]
    substrate_order = [s for s in SUBSTRATE_ORDER if s in set(diff["substrate_class"].astype(str))]
    if "All" in set(diff["substrate_class"].astype(str)):
        substrate_order.append("All")
    y_base = {label: idx for idx, label in enumerate(contact_order)}
    offsets = np.linspace(-0.3, 0.3, max(len(substrate_order), 1))
    for offset, substrate in zip(offsets, substrate_order, strict=False):
        sub = diff[diff["substrate_class"].astype(str).eq(substrate)]
        color = "#000000" if substrate == "All" else SUBSTRATE_PALETTE.get(substrate, "#555555")
        y = sub["contact_type_label"].map(y_base).astype(float).to_numpy() + offset
        xerr = np.vstack([
            sub["sd_no_demonstrated"].fillna(0).to_numpy(),
            sub["sd_demonstrated"].fillna(0).to_numpy(),
        ])
        ax.errorbar(
            sub["delta"],
            y,
            xerr=xerr,
            fmt="o",
            label=substrate,
            color=color,
            ecolor=color,
            alpha=0.95 if substrate == "All" else 0.78,
            elinewidth=1.2 if substrate == "All" else 1.0,
            capsize=3,
            markersize=5.5 if substrate == "All" else 4.8,
        )
    ax.axvline(0, color="#777777", linestyle="--", linewidth=1.0)
    ax.set_yticks([y_base[label] for label in contact_order])
    ax.set_yticklabels(contact_order)
    ax.set_xlabel("Mean fraction difference with group SD")
    ax.set_ylabel("Contact type")
    ax.set_title("")
    if legend:
        ax.legend(title="Substrate", frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.015), ncol=len(substrate_order), borderaxespad=0.0, columnspacing=1.0, handletextpad=0.4)


def figure5(r: Renderer) -> None:
    fid, title = "figure5_contact_chemistry.png", "Contact chemistry by substrate class, DP, and activity evidence"
    sources = ["combined/condition_table.tsv", "combined/protein_condition_residue_scores.tsv"]
    try:
        df = cats(add_activity(r, r.tsv(sources[0])))
        df = df[bools(df["any_valid_cluster"])].copy()
        contact_metrics = ["aromatic_contact_fraction", "polar_contact_fraction", "charged_contact_fraction", "hbond_contact_fraction"]
        df = num(df, contact_metrics)
        long = df.melt(id_vars=["protein_id", "substrate_dp", "substrate_class", "activity_evidence"], value_vars=contact_metrics, var_name="contact_type", value_name="fraction")
        long["contact_type"] = long["contact_type"].str.replace("_contact_fraction", "", regex=False).str.replace("hbond", "hydrogen bond")
        substrate_dp_order = [
            substrate + " DP" + str(dp)
            for substrate in ["chitin", "amylose", "cellulose", "starch"]
            for dp in [4, 6, 8]
            if substrate + " DP" + str(dp) in set(long["substrate_dp"].astype(str))
        ]
        mean = long.groupby(["activity_evidence", "substrate_dp", "contact_type"], observed=True)["fraction"].mean().reset_index()
        fig = plt.figure(figsize=(10.2, 7.25), constrained_layout=True)
        gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.0, 0.78], height_ratios=[1.0, 0.44])
        ax_a = fig.add_subplot(gs[0, 0:2])
        ax_b = fig.add_subplot(gs[:, 2])
        ax_c = fig.add_subplot(gs[1, 0:2])
        sns.barplot(data=mean, x="substrate_dp", y="fraction", hue="contact_type", order=substrate_dp_order, ax=ax_a)
        ax_a.set_title("Mean contact fractions")
        ax_a.set_xlabel("")
        ax_a.set_ylabel("Occupancy-weighted fraction")
        ax_a.set_ylim(0, 0.55)
        substrate_dp_ticklabels(ax_a, rotation=0, abbreviate=True)
        ax_a.legend(title="Contact type", loc="lower left", bbox_to_anchor=(0.0, 1.015), ncol=4, frameon=False, borderaxespad=0.0)
        compact_margins(ax_a, x=0.02)
        panel_label(ax_a, "A")
        diff = contact_fraction_difference_by_substrate(df, contact_metrics)
        plot_contact_fraction_difference_by_substrate(ax_c, diff, legend=True)
        compact_margins(ax_c, x=0.02, y=0.01)
        panel_label(ax_c, "C")
        valid_condition_ids = set(df["condition_id"].astype(str))
        res = cats(r.tsv(sources[1]))
        res = res[res["condition_id"].astype(str).isin(valid_condition_ids)].copy()
        res = num(res, ["residue_contact_score"])
        top = res.groupby("residue_label", observed=True)["residue_contact_score"].mean().sort_values(ascending=False).head(35).index
        res = res[res["residue_label"].isin(top)].copy()
        pivot = res.pivot_table(index="residue_label", columns="substrate_class", values="residue_contact_score", aggfunc="mean")
        pivot.index = [strip_residue_chain(label) for label in pivot.index]
        sns.heatmap(pivot, cmap="mako", ax=ax_b, cbar_kws={"label": "Contact score", "shrink": 0.78, "pad": 0.01})
        ax_b.set_title("Top 35 residue contact scores")
        ax_b.set_xlabel("Substrate")
        ax_b.set_ylabel("Residue")
        ax_b.tick_params(axis="y", labelsize=8.2)
        wrap_ticklabels(ax_b, axis="y", width=10)
        panel_label(ax_b, "B")
        fig.suptitle(title, fontsize=15, weight="bold")
        r.table(diff, "figure5_contact_fraction_difference.tsv")
        save_clean(r, fig, fid, title, sources, "non-VdW contact fractions; top residue scores", len(long))
    except Exception as exc:
        r.skip(fid, title, sources, str(exc))


def figure5_full_residue_heatmap(r: Renderer) -> None:
    fid, title = "figure5c_full_residue_contact_heatmap.png", "Full residue contact heatmap"
    source = "combined/protein_condition_residue_scores.tsv"
    try:
        cond = r.tsv("combined/condition_table.tsv")
        valid_condition_ids = set(cond.loc[bools(cond["any_valid_cluster"]), "condition_id"].astype(str))
        res = cats(r.tsv(source))
        res = res[res["condition_id"].astype(str).isin(valid_condition_ids)].copy()
        res = finite(num(res, ["residue_contact_score"]), ["residue_contact_score"])
        pivot_all = res.pivot_table(index="residue_label", columns="substrate_class", values="residue_contact_score", aggfunc="mean")
        selected: set[str] = set()
        for substrate in [s for s in SUBSTRATE_ORDER if s in pivot_all.columns]:
            selected.update(pivot_all[substrate].dropna().sort_values(ascending=False).head(20).index.astype(str))
        pivot = pivot_all.loc[pivot_all.index.astype(str).isin(selected)].copy()
        pivot = pivot.loc[pivot.max(axis=1).sort_values(ascending=False).index]
        pivot.index = [strip_residue_chain(label) for label in pivot.index]
        height = max(8.0, min(20.0, 0.16 * len(pivot) + 3.0))
        fig, ax = plt.subplots(figsize=(5.6, height), constrained_layout=True)
        sns.heatmap(pivot, cmap="mako", ax=ax, yticklabels=True, cbar_kws={"label": "Mean contact score", "shrink": 0.9, "pad": 0.012})
        ax.set_title("Top 20 residue contacts per substrate")
        ax.set_xlabel("Substrate")
        ax.set_ylabel("Residue")
        ax.tick_params(axis="y", labelsize=7.4 if len(pivot) > 45 else 8.6)
        wrap_ticklabels(ax, axis="y", width=11)
        r.table(pivot.reset_index(), "figure5c_full_residue_contact_heatmap.tsv")
        save_clean(r, fig, fid, title, [source, "combined/condition_table.tsv"], "valid-cluster conditions; union of top 20 residue labels by mean contact score within each substrate", len(pivot))
    except Exception as exc:
        r.skip(fid, title, [source], str(exc))


def figure5_all_residue_heatmap(r: Renderer) -> None:
    fid, title = "figure5c_all_residue_contact_heatmap.png", "All residue contact heatmap"
    source = "combined/protein_condition_residue_scores.tsv"
    try:
        cond = r.tsv("combined/condition_table.tsv")
        valid_condition_ids = set(cond.loc[bools(cond["any_valid_cluster"]), "condition_id"].astype(str))
        res = cats(r.tsv(source))
        res = res[res["condition_id"].astype(str).isin(valid_condition_ids)].copy()
        res = finite(num(res, ["residue_contact_score"]), ["residue_contact_score"])
        pivot = res.pivot_table(index="residue_label", columns="substrate_class", values="residue_contact_score", aggfunc="mean")
        pivot = pivot.loc[pivot.max(axis=1).sort_values(ascending=False).index]
        pivot.index = [strip_residue_chain(label) for label in pivot.index]
        height = max(18.0, min(96.0, 0.055 * len(pivot) + 4.0))
        fig, ax = plt.subplots(figsize=(5.4, height), constrained_layout=True)
        sns.heatmap(pivot, cmap="mako", ax=ax, yticklabels=True, cbar_kws={"label": "Mean contact score", "shrink": 0.95, "pad": 0.01})
        ax.set_title(title)
        ax.set_xlabel("Substrate")
        ax.set_ylabel("Residue")
        ax.tick_params(axis="y", labelsize=2.8)
        r.table(pivot.reset_index(), "figure5c_all_residue_contact_heatmap.tsv")
        save_clean(r, fig, fid, title, [source, "combined/condition_table.tsv"], "valid-cluster conditions; all residue labels by mean contact score within each substrate", len(pivot))
    except Exception as exc:
        r.skip(fid, title, [source], str(exc))


def figure5_residue_contact_activity_difference_heatmap(r: Renderer) -> None:
    fid, title = "figure5_residue_contact_activity_difference_heatmap.png", "Residue contact difference for demonstrated versus non-demonstrated substrate conditions"
    sources = ["combined/protein_condition_residue_scores.tsv", "combined/condition_table.tsv"]
    try:
        cond = cats(add_activity(r, r.tsv(sources[1])))
        cond = cond[bools(cond["any_valid_cluster"])].copy()
        cond = cond[["condition_id", "activity_evidence"]].drop_duplicates("condition_id")
        res = cats(r.tsv(sources[0]))
        res = res.merge(cond, on="condition_id", how="inner")
        res = finite(num(res, ["residue_contact_score"]), ["residue_contact_score"])
        means = (
            res.groupby(["residue_label", "substrate_class", "activity_evidence"], observed=True)["residue_contact_score"]
            .mean()
            .reset_index()
        )
        pivot = means.pivot_table(
            index="residue_label",
            columns=["substrate_class", "activity_evidence"],
            values="residue_contact_score",
            aggfunc="mean",
        )
        substrate_cols = [substrate for substrate in SUBSTRATE_ORDER if substrate in pivot.columns.get_level_values(0)]
        diff = pd.DataFrame(index=pivot.index)
        for substrate in substrate_cols:
            active = pivot[(substrate, ACTIVITY_ORDER[0])] if (substrate, ACTIVITY_ORDER[0]) in pivot.columns else pd.Series(np.nan, index=pivot.index)
            inactive = pivot[(substrate, ACTIVITY_ORDER[1])] if (substrate, ACTIVITY_ORDER[1]) in pivot.columns else pd.Series(np.nan, index=pivot.index)
            diff[substrate] = active - inactive
        diff = diff.dropna(how="all")
        selected = diff.abs().max(axis=1).sort_values(ascending=False).head(35).index
        plot = diff.loc[selected, substrate_cols].copy()
        plot.index = [strip_residue_chain(label) for label in plot.index]
        max_abs = float(np.nanmax(np.abs(plot.to_numpy()))) if plot.size else 1.0
        if not math.isfinite(max_abs) or max_abs == 0:
            max_abs = 1.0
        fig, ax = plt.subplots(figsize=(5.8, max(8.0, 0.22 * len(plot) + 2.8)), constrained_layout=True)
        sns.heatmap(
            plot,
            cmap="vlag",
            center=0,
            vmin=-max_abs,
            vmax=max_abs,
            ax=ax,
            cbar_kws={"label": "Demonstrated - no demonstrated contact score", "shrink": 0.88, "pad": 0.012},
        )
        ax.set_title(title)
        ax.set_xlabel("Substrate")
        ax.set_ylabel("Residue")
        ax.tick_params(axis="y", labelsize=8.0)
        wrap_ticklabels(ax, axis="y", width=10)
        r.table(plot.reset_index().rename(columns={"index": "residue"}), "figure5_residue_contact_activity_difference_heatmap.tsv")
        save_clean(r, fig, fid, title, sources, "valid-cluster conditions; mean residue contact score difference between demonstrated and no demonstrated activity within each substrate", len(plot))
    except Exception as exc:
        r.skip(fid, title, sources, str(exc))


def figure5_contact_fraction_difference_by_substrate(r: Renderer) -> None:
    fid, title = "figure5_contact_fraction_difference_by_substrate.png", "Contact fraction difference by demonstrated substrate"
    source = "combined/condition_table.tsv"
    try:
        df = cats(add_activity(r, r.tsv(source)))
        df = df[bools(df["any_valid_cluster"])].copy()
        contact_metrics = ["aromatic_contact_fraction", "polar_contact_fraction", "charged_contact_fraction", "hbond_contact_fraction"]
        df = num(df, contact_metrics)
        diff = contact_fraction_difference_by_substrate(df, contact_metrics)
        if diff.empty:
            raise ValueError("no substrate classes had both demonstrated and no demonstrated activity rows")
        fig, ax = plt.subplots(figsize=(8.8, 5.0), constrained_layout=True)
        plot_contact_fraction_difference_by_substrate(ax, diff, legend=True)
        compact_margins(ax, x=0.03, y=0.08)
        r.table(diff, "figure5_contact_fraction_difference_by_substrate.tsv")
        save_clean(r, fig, fid, title, [source], "valid-cluster conditions; demonstrated minus no demonstrated contact fraction difference per substrate class", len(diff))
    except Exception as exc:
        r.skip(fid, title, [source], str(exc))


def _render_cbm_paired_analysis(
    r: Renderer,
    *,
    figure_id: str,
    title: str,
    contact_summary_name: str,
    bridge_fraction_threshold: float | None = None,
) -> None:
    sources = ["postprocess/15_cbm_paired_analysis/cbm_construct_condition_summary.tsv", "postprocess/15_cbm_paired_analysis/cbm_paired_comparison_table.tsv"]
    try:
        summary = cats(r.tsv(sources[0]))
        paired = cats(r.tsv(sources[1]))
        need(paired, ["protein_id", "domain_only_condition_id", "full_length_condition_id"])
        filters = "paired CBM comparison tables"
        paired_for_panel_a = paired.copy()
        if bridge_fraction_threshold is not None:
            need(paired, ["bridge_fraction_full_length"])
            paired = finite(num(paired, ["bridge_fraction_full_length"]), ["bridge_fraction_full_length"])
            paired = paired[paired["bridge_fraction_full_length"] > bridge_fraction_threshold].copy()
            filters = f"panel A: all paired CBM comparison rows; panels B-E: bridge_fraction_full_length > {bridge_fraction_threshold:g}"
            if paired.empty:
                raise ValueError(filters)
        metrics = [
            "top_cluster_occupancy",
            "cluster_entropy",
            "geometry_plausible_fraction",
            "C4_minus_C1_geometry_bias",
            "catalytic_surface_contact_fraction",
            "aromatic_contact_fraction",
            "polar_contact_fraction",
        ]
        metrics = [metric for metric in metrics if metric in summary.columns][:4]
        panel_a_pair_map = paired_for_panel_a[["protein_id", "domain_only_condition_id", "full_length_condition_id"]].copy()
        panel_a_pair_map["pair_id"] = np.arange(len(panel_a_pair_map))
        panel_a_long_ids = pd.concat(
            [
                panel_a_pair_map[["pair_id", "protein_id", "domain_only_condition_id"]].rename(columns={"domain_only_condition_id": "condition_id"}).assign(expected_construct="domain_only"),
                panel_a_pair_map[["pair_id", "protein_id", "full_length_condition_id"]].rename(columns={"full_length_condition_id": "condition_id"}).assign(expected_construct="full_length"),
            ],
            ignore_index=True,
        )
        panel_a_matched = panel_a_long_ids.merge(summary[["condition_id", "construct_type", "any_valid_cluster"]], on="condition_id", how="left")
        panel_a_matched["valid_cluster"] = bools(panel_a_matched["any_valid_cluster"])
        pair_map = paired[["protein_id", "domain_only_condition_id", "full_length_condition_id"]].copy()
        pair_map["pair_id"] = np.arange(len(pair_map))
        long_ids = pd.concat(
            [
                pair_map[["pair_id", "protein_id", "domain_only_condition_id"]].rename(columns={"domain_only_condition_id": "condition_id"}).assign(expected_construct="domain_only"),
                pair_map[["pair_id", "protein_id", "full_length_condition_id"]].rename(columns={"full_length_condition_id": "condition_id"}).assign(expected_construct="full_length"),
            ],
            ignore_index=True,
        )
        matched = long_ids.merge(summary[["condition_id", "construct_type", "any_valid_cluster", *metrics]], on="condition_id", how="left")
        matched["valid_cluster"] = bools(matched["any_valid_cluster"])
        fig = plt.figure(figsize=(11.6, 7.1), constrained_layout=True)
        gs = fig.add_gridspec(2, 6, height_ratios=[1.0, 1.0])
        axes = [
            fig.add_subplot(gs[0, 0:2]),
            fig.add_subplot(gs[0, 2:4]),
            fig.add_subplot(gs[0, 4:6]),
            fig.add_subplot(gs[1, 0:3]),
            fig.add_subplot(gs[1, 3:6]),
        ]
        ax0 = axes[0]
        for _, pair in panel_a_matched.groupby("pair_id", observed=True):
            if pair["expected_construct"].nunique() == 2:
                pair = pair.sort_values("expected_construct", key=lambda s: s.map({"domain_only": 0, "full_length": 1}))
                ax0.plot(pair["expected_construct"], pair["valid_cluster"].astype(float), color="#999999", alpha=0.16, linewidth=0.8)
        sns.pointplot(data=panel_a_matched, x="expected_construct", y="valid_cluster", order=CONSTRUCT_ORDER, color="#000000", errorbar=("pi", 50), ax=ax0)
        ax0.set_ylim(-0.03, 1.18)
        ax0.set_xlabel("")
        ax0.set_ylabel("Valid-Cluster Fraction")
        ax0.set_xticks([0, 1])
        ax0.set_xticklabels(["Catalytic domain", "Full-length"])
        ax0.set_xlim(-0.08, 1.08)
        annotate(ax0, f"n paired conditions={len(paired_for_panel_a)}\nn proteins={paired_for_panel_a['protein_id'].nunique()}")
        panel_label(ax0, "A")
        value_long = matched.melt(
            id_vars=["pair_id", "protein_id", "condition_id", "expected_construct"],
            value_vars=metrics,
            var_name="metric",
            value_name="value",
        )
        value_long = finite(num(value_long, ["value"]), ["value"])
        value_long["metric_label"] = value_long["metric"].map(nice)
        xmap = {"domain_only": 0, "full_length": 1}
        for idx, (ax, metric) in enumerate(zip(axes[1:5], value_long["metric_label"].drop_duplicates(), strict=False)):
            metric_df = value_long[value_long["metric_label"] == metric]
            for _, pair in metric_df.groupby("pair_id", observed=True):
                if pair["expected_construct"].nunique() < 2:
                    continue
                pair = pair.sort_values("expected_construct", key=lambda s: s.map(xmap))
                ax.plot(pair["expected_construct"].map(xmap), pair["value"], color="#999999", alpha=0.18, linewidth=0.8)
            sns.pointplot(data=metric_df, x="expected_construct", y="value", order=CONSTRUCT_ORDER, errorbar=("pi", 50), color="#000000", ax=ax)
            ax.set_xlabel("")
            ax.set_ylabel(title_case_label(metric))
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["Catalytic domain", "Full-length"])
            ax.set_xlim(-0.08, 1.08)
            if idx in {0, 2, 3}:
                y_min, y_max = ax.get_ylim()
                ax.set_ylim(y_min, y_max + 0.16 * (y_max - y_min))
            wrap_ticklabels(ax, axis="x", width=12)
            annotate(ax, f"n conditions={metric_df['condition_id'].nunique()}\nn proteins={metric_df['protein_id'].nunique()}")
            panel_label(ax, chr(ord("B") + idx))
        if "non_core_ligand_contact_fraction_full_length" in paired.columns:
            contact = finite(num(paired, ["non_core_ligand_contact_fraction_full_length", "bridge_fraction_full_length"]), ["non_core_ligand_contact_fraction_full_length"])
            summary_rows = contact.groupby("substrate_class", observed=True).agg(
                n_pairs=("protein_id", "size"),
                n_nonzero_non_core=("non_core_ligand_contact_fraction_full_length", lambda s: int((s > 0).sum())),
                median_non_core=("non_core_ligand_contact_fraction_full_length", "median"),
                max_non_core=("non_core_ligand_contact_fraction_full_length", "max"),
                n_nonzero_bridge=("bridge_fraction_full_length", lambda s: int((pd.to_numeric(s, errors="coerce") > 0).sum()) if "bridge_fraction_full_length" in contact.columns else 0),
            ).reset_index()
            r.table(summary_rows, contact_summary_name)
        fig.suptitle(title, fontsize=15, weight="bold")
        save_clean(r, fig, figure_id, title, sources, filters, len(paired))
    except Exception as exc:
        r.skip(figure_id, title, sources, str(exc))


def figure6(r: Renderer) -> None:
    _render_cbm_paired_analysis(
        r,
        figure_id="figure6_cbm_paired_analysis.png",
        title="CBM paired analysis: catalytic domain versus full-length",
        contact_summary_name="figure6_non_core_contact_summary.tsv",
    )


def figure6b_cbm_bridge_high_fraction(r: Renderer) -> None:
    _render_cbm_paired_analysis(
        r,
        figure_id="figure6b_cbm_bridge_high_fraction.png",
        title="CBM paired analysis for high bridge fraction",
        contact_summary_name="figure6b_bridge_high_fraction_contact_summary.tsv",
        bridge_fraction_threshold=0.5,
    )


def figure6c_cbm_delta_heatmap(r: Renderer) -> None:
    fid, title = "figure6c_cbm_delta_heatmap.png", "CBM delta geometry/support heatmap"
    source = "postprocess/15_cbm_paired_analysis/cbm_paired_comparison_table.tsv"
    try:
        paired = cats(r.tsv(source))
        delta_cols = [c for c in ["delta_C1_compatible_fraction", "delta_C4_compatible_fraction", "delta_geometry_plausible_fraction", "delta_top_cluster_occupancy", "delta_cluster_entropy", "delta_C4_minus_C1_geometry_bias"] if c in paired.columns]
        heat = num(paired, delta_cols)
        heat["row"] = heat["protein_id"].astype(str) + " " + heat["substrate_class"].astype(str) + " DP" + heat["dp"].astype("Int64").astype(str)
        selected_rows: set[str] = set()
        for col in delta_cols:
            ranked = heat[["row", col]].dropna().groupby("row", observed=True)[col].mean()
            selected_rows.update(ranked[ranked > 0].sort_values(ascending=False).head(10).index.astype(str))
            selected_rows.update(ranked[ranked < 0].sort_values(ascending=True).head(10).index.astype(str))
        pivot = heat.pivot_table(index="row", values=delta_cols, aggfunc="mean")
        pivot = pivot.loc[pivot.index.astype(str).isin(selected_rows)].copy()
        if len(pivot):
            order = pivot.abs().max(axis=1).sort_values(ascending=False).index
            pivot = pivot.loc[order]
        max_abs = np.nanmax(np.abs(pivot.to_numpy())) if pivot.size else 1
        if not math.isfinite(max_abs) or max_abs == 0:
            max_abs = 1
        fig, ax = plt.subplots(figsize=(7.4, max(7.2, min(50.0, 0.12 * len(pivot) + 2.8))), constrained_layout=True)
        heat_pivot = pivot.rename(columns=nice)
        sns.heatmap(heat_pivot, cmap="vlag", center=0, vmin=-max_abs, vmax=max_abs, yticklabels=True, ax=ax, cbar_kws={"label": "Full-length - catalytic domain", "shrink": 0.9, "pad": 0.012})
        ax.set_title(title)
        ax.set_xlabel("")
        ax.set_ylabel("Paired condition")
        ax.tick_params(axis="y", labelsize=5.0 if len(pivot) > 200 else 7.0)
        wrap_ticklabels(ax, axis="x", width=13)
        r.table(pivot.reset_index(), "figure6c_cbm_delta_heatmap.tsv")
        save_clean(r, fig, fid, title, [source], "union of top 10 positive and top 10 negative full-length minus domain-only changes per metric", len(pivot))
    except Exception as exc:
        r.skip(fid, title, [source], str(exc))


def figure6d_cbm_delta_heatmap_high_bridge(r: Renderer) -> None:
    fid, title = "figure6d_cbm_delta_heatmap_high_bridge.png", "CBM delta geometry/support heatmap for high bridge fraction"
    source = "postprocess/15_cbm_paired_analysis/cbm_paired_comparison_table.tsv"
    try:
        paired = cats(r.tsv(source))
        need(paired, ["bridge_fraction_full_length"])
        paired = finite(num(paired, ["bridge_fraction_full_length"]), ["bridge_fraction_full_length"])
        paired = paired[paired["bridge_fraction_full_length"] > 0.5].copy()
        if paired.empty:
            raise ValueError("no paired rows with bridge_fraction_full_length > 0.5")
        delta_cols = [c for c in ["delta_C1_compatible_fraction", "delta_C4_compatible_fraction", "delta_geometry_plausible_fraction", "delta_top_cluster_occupancy", "delta_cluster_entropy", "delta_C4_minus_C1_geometry_bias"] if c in paired.columns]
        heat = num(paired, delta_cols)
        heat["row"] = heat["protein_id"].astype(str) + " " + heat["substrate_class"].astype(str) + " DP" + heat["dp"].astype("Int64").astype(str)
        selected_rows: set[str] = set()
        for col in delta_cols:
            ranked = heat[["row", col]].dropna().groupby("row", observed=True)[col].mean()
            selected_rows.update(ranked[ranked > 0].sort_values(ascending=False).head(10).index.astype(str))
            selected_rows.update(ranked[ranked < 0].sort_values(ascending=True).head(10).index.astype(str))
        pivot = heat.pivot_table(index="row", values=delta_cols, aggfunc="mean")
        pivot = pivot.loc[pivot.index.astype(str).isin(selected_rows)].copy()
        if len(pivot):
            order = pivot.abs().max(axis=1).sort_values(ascending=False).index
            pivot = pivot.loc[order]
        max_abs = np.nanmax(np.abs(pivot.to_numpy())) if pivot.size else 1
        if not math.isfinite(max_abs) or max_abs == 0:
            max_abs = 1
        fig, ax = plt.subplots(figsize=(7.4, max(7.2, min(50.0, 0.12 * len(pivot) + 2.8))), constrained_layout=True)
        sns.heatmap(pivot.rename(columns=nice), cmap="vlag", center=0, vmin=-max_abs, vmax=max_abs, yticklabels=True, ax=ax, cbar_kws={"label": "Full-length - catalytic domain", "shrink": 0.9, "pad": 0.012})
        ax.set_title(title)
        ax.set_xlabel("")
        ax.set_ylabel("Paired condition")
        ax.tick_params(axis="y", labelsize=5.0 if len(pivot) > 200 else 7.0)
        wrap_ticklabels(ax, axis="x", width=13)
        r.table(pivot.reset_index(), "figure6d_cbm_delta_heatmap_high_bridge.tsv")
        save_clean(r, fig, fid, title, [source], "bridge_fraction_full_length > 0.5; union of top 10 positive and top 10 negative full-length minus domain-only changes per metric", len(pivot))
    except Exception as exc:
        r.skip(fid, title, [source], str(exc))


def figure6c_all_cbm_delta_heatmap(r: Renderer) -> None:
    fid, title = "figure6c_all_cbm_delta_heatmap.png", "All CBM delta geometry/support heatmap"
    source = "postprocess/15_cbm_paired_analysis/cbm_paired_comparison_table.tsv"
    try:
        paired = cats(r.tsv(source))
        delta_cols = [c for c in ["delta_C1_compatible_fraction", "delta_C4_compatible_fraction", "delta_geometry_plausible_fraction", "delta_top_cluster_occupancy", "delta_cluster_entropy", "delta_C4_minus_C1_geometry_bias"] if c in paired.columns]
        heat = num(paired, delta_cols)
        heat["row"] = heat["protein_id"].astype(str) + " " + heat["substrate_class"].astype(str) + " DP" + heat["dp"].astype("Int64").astype(str)
        pivot = heat.pivot_table(index="row", values=delta_cols, aggfunc="mean")
        if len(pivot):
            pivot = pivot.loc[pivot.abs().max(axis=1).sort_values(ascending=False).index]
        max_abs = np.nanmax(np.abs(pivot.to_numpy())) if pivot.size else 1
        if not math.isfinite(max_abs) or max_abs == 0:
            max_abs = 1
        fig, ax = plt.subplots(figsize=(7.0, max(12.0, min(46.0, 0.085 * len(pivot) + 3.5))), constrained_layout=True)
        sns.heatmap(pivot.rename(columns=nice), cmap="vlag", center=0, vmin=-max_abs, vmax=max_abs, yticklabels=True, ax=ax, cbar_kws={"label": "Full-length - catalytic domain", "shrink": 0.95, "pad": 0.01})
        ax.set_title(title)
        ax.set_xlabel("")
        ax.set_ylabel("Paired condition")
        ax.tick_params(axis="y", labelsize=4.6)
        wrap_ticklabels(ax, axis="x", width=13)
        r.table(pivot.reset_index(), "figure6c_all_cbm_delta_heatmap.tsv")
        save_clean(r, fig, fid, title, [source], "all full-length minus catalytic-domain paired-condition changes per metric", len(pivot))
    except Exception as exc:
        r.skip(fid, title, [source], str(exc))


def figure7(r: Renderer) -> None:
    fid, title = "figure7_crystal_anchoring.png", "Reference structure comparison"
    sources = ["combined/crystal_anchor_table.tsv", "combined/condition_table.tsv", "input_data/metadata_pdb_count.tsv"]
    try:
        anchor = r.tsv(sources[0])
        cond = cats(r.tsv(sources[1]))
        df = cats(anchor.merge(cond[["condition_id", "substrate_class", "dp", "construct_type"]].drop_duplicates("condition_id"), on="condition_id", how="left"))
        fig = plt.figure(figsize=(10.2, 7.4), constrained_layout=True)
        gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 0.68], width_ratios=[0.35, 1.0, 1.0, 0.35])
        axes = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[1, 1:3])]
        rmsd = finite(num(df, ["local_pocket_rmsd"]), ["local_pocket_rmsd"])
        sns.boxplot(data=rmsd, x="substrate_class", y="local_pocket_rmsd", hue="construct_type", showfliers=False, palette=CONSTRUCT_PALETTE, ax=axes[0])
        sns.stripplot(data=rmsd, x="substrate_class", y="local_pocket_rmsd", hue="construct_type", dodge=True, color="#111111", alpha=0.45, size=2.2, ax=axes[0], legend=False)
        axes[0].set_title("Local pocket RMSD")
        axes[0].set_xlabel("Substrate")
        axes[0].set_ylabel("RMSD (Å)")
        axes[0].legend(title="Construct", loc="upper right", frameon=False)
        compact_margins(axes[0], x=0.02)
        panel_label(axes[0], "A")
        sim = finite(num(df, ["local_pocket_rmsd", "ifp_tanimoto"]), ["local_pocket_rmsd", "ifp_tanimoto"])
        sns.scatterplot(data=sim, x="local_pocket_rmsd", y="ifp_tanimoto", hue="substrate_class", palette=SUBSTRATE_PALETTE, ax=axes[1])
        axes[1].set_ylim(-0.03, 1.03)
        axes[1].set_title("IFP similarity versus RMSD")
        axes[1].set_xlabel("Local pocket RMSD (Å)")
        axes[1].set_ylabel("IFP Tanimoto")
        axes[1].legend(fontsize=LEGEND_FONTSIZE, loc="upper right", frameon=False)
        compact_margins(axes[1], x=0.02)
        panel_label(axes[1], "B")
        pdb = num(r.tsv("input_data/metadata_pdb_count.tsv", base="analyse"), ["PDB_Before_AF3_Cutoff_Count"]).rename(columns={"UniProt_ID": "protein_id"})
        rates = cond.assign(valid=bools(cond["any_valid_cluster"])).groupby("protein_id", observed=True)["valid"].mean().reset_index()
        rates = rates.merge(pdb[["protein_id", "PDB_Before_AF3_Cutoff_Count"]], on="protein_id", how="left").fillna({"PDB_Before_AF3_Cutoff_Count": 0})
        sns.scatterplot(data=rates, x="PDB_Before_AF3_Cutoff_Count", y="valid", color="#3b6ea8", ax=axes[2])
        x = rates["PDB_Before_AF3_Cutoff_Count"].astype(float).to_numpy()
        y = rates["valid"].astype(float).to_numpy()
        r2 = np.nan
        if len(rates) >= 2 and np.nanstd(x) > 0:
            slope, intercept = np.polyfit(x, y, 1)
            xs = np.linspace(np.nanmin(x), np.nanmax(x), 100)
            axes[2].plot(xs, slope * xs + intercept, color="#555555", linestyle="--", linewidth=1.2)
            pred = slope * x + intercept
            ss_res = np.nansum((y - pred) ** 2)
            ss_tot = np.nansum((y - np.nanmean(y)) ** 2)
            r2 = 1.0 - ss_res / ss_tot if ss_tot else np.nan
        axes[2].set_ylim(-0.03, 1.03)
        axes[2].set_title("Clusterability versus crystal reference count")
        axes[2].set_xlabel("X-ray crystal references before AF3 cutoff")
        axes[2].set_ylabel("Valid-cluster fraction")
        axes[2].xaxis.set_major_locator(MaxNLocator(integer=True))
        if len(rates):
            xmin = float(np.nanmin(x))
            xmax = float(np.nanmax(x))
            axes[2].set_xlim(xmin - 0.35, xmax + 0.35)
        annotate(axes[2], f"linear R²={r2:.2f}" if np.isfinite(r2) else "linear R²=NA")
        panel_label(axes[2], "C")
        fig.suptitle(title, fontsize=15, weight="bold")
        save_clean(r, fig, fid, title, sources, "crystal anchor rows and protein-level clusterability", len(df))
    except Exception as exc:
        r.skip(fid, title, sources, str(exc))


def figure7c_clusterability_bubble(r: Renderer) -> None:
    fid, title = "figure7c_clusterability_bubble.png", "Clusterability versus crystal reference count with duplicate-value counts"
    sources = ["combined/condition_table.tsv", "input_data/metadata_pdb_count.tsv"]
    try:
        cond = cats(r.tsv(sources[0]))
        pdb = num(r.tsv(sources[1], base="analyse"), ["PDB_Before_AF3_Cutoff_Count"]).rename(columns={"UniProt_ID": "protein_id"})
        rates = cond.assign(valid=bools(cond["any_valid_cluster"])).groupby("protein_id", observed=True)["valid"].mean().reset_index()
        rates = rates.merge(pdb[["protein_id", "PDB_Before_AF3_Cutoff_Count"]], on="protein_id", how="left").fillna({"PDB_Before_AF3_Cutoff_Count": 0})
        rates["PDB_Before_AF3_Cutoff_Count"] = pd.to_numeric(rates["PDB_Before_AF3_Cutoff_Count"], errors="coerce").fillna(0).astype(int)
        grouped = rates.groupby(["PDB_Before_AF3_Cutoff_Count", "valid"], observed=True).agg(n_proteins=("protein_id", "nunique")).reset_index()
        fig, ax = plt.subplots(figsize=(6.6, 4.6), constrained_layout=True)
        sns.scatterplot(
            data=grouped,
            x="PDB_Before_AF3_Cutoff_Count",
            y="valid",
            size="n_proteins",
            sizes=(35, 360),
            color="#3b6ea8",
            alpha=0.78,
            edgecolor="#1f3552",
            linewidth=0.7,
            legend="brief",
            ax=ax,
        )
        x = rates["PDB_Before_AF3_Cutoff_Count"].astype(float).to_numpy()
        y = rates["valid"].astype(float).to_numpy()
        r2 = np.nan
        if len(rates) >= 2 and np.nanstd(x) > 0:
            slope, intercept = np.polyfit(x, y, 1)
            xs = np.linspace(np.nanmin(x), np.nanmax(x), 100)
            ax.plot(xs, slope * xs + intercept, color="#555555", linestyle="--", linewidth=1.2)
            pred = slope * x + intercept
            ss_res = np.nansum((y - pred) ** 2)
            ss_tot = np.nansum((y - np.nanmean(y)) ** 2)
            r2 = 1.0 - ss_res / ss_tot if ss_tot else np.nan
        ax.set_ylim(-0.03, 1.03)
        if len(rates):
            xmin = float(np.nanmin(x))
            xmax = float(np.nanmax(x))
            ax.set_xlim(xmin - 0.35, xmax + 0.35)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_xlabel("X-ray crystal references before AF3 cutoff")
        ax.set_ylabel("Valid-cluster fraction")
        annotate(ax, f"linear R²={r2:.2f}" if np.isfinite(r2) else "linear R²=NA")
        legend = ax.get_legend()
        if legend is not None:
            legend.set_title("Proteins\nper dot")
            legend.set_frame_on(False)
            legend.set_bbox_to_anchor((1.01, 0.5))
            legend._loc = 6
            legend._legend_box.sep = 10
            for handle in legend.legend_handles:
                if hasattr(handle, "set_alpha"):
                    handle.set_alpha(0.78)
        r.table(grouped, "figure7c_clusterability_bubble_counts.tsv")
        save_clean(r, fig, fid, title, sources, "protein-level valid-cluster fraction aggregated by identical PDB-count and clusterability values", len(grouped))
    except Exception as exc:
        r.skip(fid, title, sources, str(exc))


def figure8(r: Renderer) -> None:
    fid, title = "figure8_family_aligned_residue_enrichment.png", "Residues associated with experimental regioselectivity and substrate activity"
    sources = ["postprocess/08_family_residue_enrichment/family_aligned_residue_table.tsv", "input_data/metadata_final_ec_fixed.tsv"]
    try:
        aligned = num(
            r.tsv(sources[0]),
            ["alignment_column", "mean_condition_residue_contact_score", "max_condition_residue_contact_score"],
        )
        meta = activity_table(r)[["protein_id", "active_substrates", "regio"]].drop_duplicates("protein_id")
        df = aligned.merge(meta, on="protein_id", how="left")
        df["score"] = pd.to_numeric(df["mean_condition_residue_contact_score"], errors="coerce").fillna(0.0)
        df = df[df["score"] > 0].copy()
        df["active_set"] = df["active_substrates"].fillna("").map(lambda value: set(filter(None, str(value).split(","))))
        df["regio_norm"] = df["regio"].fillna("").astype(str).str.upper().str.replace(" ", "", regex=False)
        group_defs = {
            "C1": df["regio_norm"].eq("C1"),
            "C4": df["regio_norm"].eq("C4"),
            "Mixed C1/C4": df["regio_norm"].isin(["C1+C4", "C4+C1", "C1/C4", "MIXED"]),
            "Chitin-active": df["active_set"].map(lambda active: "chitin" in active),
            "Cellulose-active": df["active_set"].map(lambda active: "cellulose" in active),
        }
        identity_groups = ["C1", "C4", "Mixed C1/C4", "Chitin-active", "Cellulose-active"]
        focused_groups = ["C1", "C4", "Chitin-active", "Cellulose-active"]
        identity_rows = []
        for label, mask in group_defs.items():
            sub = df[mask].copy()
            total = sub["score"].sum()
            for residue, score in sub.groupby("alignment_residue", observed=True)["score"].sum().items():
                identity_rows.append({"residue": residue, "group": label, "share_percent": 100.0 * score / total if total else 0.0})
        identity = pd.DataFrame(identity_rows)
        identity_pivot = identity.pivot_table(index="residue", columns="group", values="share_percent", aggfunc="sum").fillna(0.0)
        identity_pivot = identity_pivot.reindex(columns=identity_groups)
        identity_pivot = identity_pivot.loc[identity_pivot.max(axis=1).sort_values(ascending=False).index]
        consensus = (
            df.groupby(["family_label", "alignment_column", "alignment_residue"], observed=True)["score"]
            .sum()
            .reset_index()
            .sort_values("score", ascending=False)
            .drop_duplicates(["family_label", "alignment_column"])
        )
        consensus["hotspot_label"] = consensus["family_label"].astype(str) + " col" + consensus["alignment_column"].astype("Int64").astype(str) + " " + consensus["alignment_residue"].astype(str)
        hotspot_rows = []
        for label in focused_groups:
            sub = df[group_defs[label]].copy()
            grouped = sub.groupby(["family_label", "alignment_column"], observed=True)["score"].mean().reset_index()
            grouped["group"] = label
            hotspot_rows.append(grouped)
        hotspot_long = pd.concat(hotspot_rows, ignore_index=True).merge(
            consensus[["family_label", "alignment_column", "alignment_residue", "hotspot_label"]],
            on=["family_label", "alignment_column"],
            how="left",
        )
        hotspot_pivot_all = hotspot_long.pivot_table(index="hotspot_label", columns="group", values="score", aggfunc="mean").fillna(0.0)
        top_labels = hotspot_pivot_all.max(axis=1).sort_values(ascending=False).head(34).index
        hotspot_pivot = hotspot_pivot_all.loc[top_labels, focused_groups]
        enrich_rows = []
        eps = 1e-9
        comparison_specs = [
            ("C1", "C1-only vs C4-only", group_defs["C4"]),
            ("Mixed C1/C4", "Mixed C1/C4 vs rest", ~group_defs["Mixed C1/C4"]),
            ("Chitin-active", "Chitin-active vs non-chitin", ~group_defs["Chitin-active"]),
            ("Cellulose-active", "Cellulose-active vs non-cellulose", ~group_defs["Cellulose-active"]),
        ]
        for label, comparison, out_mask in comparison_specs:
            in_sub = df[group_defs[label]]
            out_sub = df[out_mask]
            in_total = in_sub["score"].sum()
            out_total = out_sub["score"].sum()
            residues = sorted(set(in_sub["alignment_residue"].dropna().astype(str)) | set(out_sub["alignment_residue"].dropna().astype(str)))
            for residue in residues:
                in_share = in_sub.loc[in_sub["alignment_residue"].astype(str).eq(residue), "score"].sum() / in_total if in_total else 0.0
                out_share = out_sub.loc[out_sub["alignment_residue"].astype(str).eq(residue), "score"].sum() / out_total if out_total else 0.0
                enrich_rows.append({"residue": residue, "comparison": comparison, "log2_enrichment": math.log2((in_share + eps) / (out_share + eps))})
        enrich = pd.DataFrame(enrich_rows)
        enrich_pivot = enrich.pivot_table(index="residue", columns="comparison", values="log2_enrichment", aggfunc="mean").fillna(0.0)
        enrich_pivot = enrich_pivot.reindex(columns=[spec[1] for spec in comparison_specs])
        enrich_pivot = enrich_pivot.loc[identity_pivot.index.intersection(enrich_pivot.index)]
        class_map = {
            **{residue: "Aromatic" for residue in ["F", "W", "Y"]},
            **{residue: "Polar uncharged" for residue in ["S", "T", "N", "Q", "C"]},
            **{residue: "Acidic" for residue in ["D", "E"]},
            **{residue: "Basic" for residue in ["K", "R", "H"]},
            **{residue: "Hydrophobic" for residue in ["A", "V", "I", "L", "M"]},
            **{residue: "Gly/Pro/special" for residue in ["G", "P"]},
        }
        class_order = ["Aromatic", "Polar uncharged", "Acidic", "Basic", "Hydrophobic", "Gly/Pro/special"]
        class_rows = []
        df["residue_class"] = df["alignment_residue"].astype(str).map(class_map).fillna("Gly/Pro/special")
        for label in focused_groups:
            sub = df[group_defs[label]].copy()
            total = sub["score"].sum()
            for residue_class, score in sub.groupby("residue_class", observed=True)["score"].sum().items():
                class_rows.append({"residue_class": residue_class, "group": label, "share_percent": 100.0 * score / total if total else 0.0})
        class_pivot = pd.DataFrame(class_rows).pivot_table(index="residue_class", columns="group", values="share_percent", aggfunc="sum").fillna(0.0)
        class_pivot = class_pivot.reindex(index=class_order, columns=focused_groups).fillna(0.0)
        fig, axes = plt.subplots(2, 2, figsize=(9.2, 12.6), constrained_layout=True)
        ax_a, ax_b, ax_c, ax_d = axes.flat
        heat_kws = {"shrink": 0.82, "pad": 0.012}
        sns.heatmap(identity_pivot, cmap="YlGnBu", ax=ax_a, cbar_kws={**heat_kws, "label": "Contact-weighted residue share (%)"})
        ax_a.set_title("Experimental residue identity by label")
        ax_a.set_xlabel("")
        ax_a.set_ylabel("Residue type")
        wrap_ticklabels(ax_a, axis="x", width=11)
        panel_label(ax_a, "A")
        max_abs = max(1.0, float(np.nanmax(np.abs(enrich_pivot.to_numpy()))) if enrich_pivot.size else 1.0)
        sns.heatmap(enrich_pivot, cmap="vlag", center=0, vmin=-max_abs, vmax=max_abs, ax=ax_b, cbar_kws={**heat_kws, "label": "log2 enrichment"})
        ax_b.set_title("Experimental residue identity log2 enrichment")
        ax_b.set_xlabel("")
        ax_b.set_ylabel("Residue type")
        wrap_ticklabels(ax_b, axis="x", width=13)
        panel_label(ax_b, "B")
        sns.heatmap(hotspot_pivot, cmap="mako", ax=ax_c, yticklabels=True, cbar_kws={**heat_kws, "label": "Contact-weighted score"})
        ax_c.set_title("Experimental aligned hotspot signal")
        ax_c.set_xlabel("")
        ax_c.set_ylabel("Family + consensus residue + alignment column")
        ax_c.set_yticks(np.arange(len(hotspot_pivot.index)) + 0.5)
        ax_c.set_yticklabels([wrap_label(label, 18) for label in hotspot_pivot.index], fontsize=6.9)
        ax_c.tick_params(axis="y", labelsize=5.9)
        wrap_ticklabels(ax_c, axis="x", width=11)
        panel_label(ax_c, "C")
        sns.heatmap(class_pivot, cmap="YlOrBr", ax=ax_d, cbar_kws={**heat_kws, "label": "Contact-weighted residue class share (%)"})
        ax_d.set_title("Residue class by experimental label")
        ax_d.set_xlabel("")
        ax_d.set_ylabel("Residue class")
        wrap_ticklabels(ax_d, axis="x", width=11)
        wrap_ticklabels(ax_d, axis="y", width=14)
        panel_label(ax_d, "D")
        fig.suptitle(title, fontsize=15, weight="bold")
        r.table(identity_pivot.reset_index(), "figure8_residue_identity_share.tsv")
        r.table(hotspot_pivot.reset_index(), "figure8_top_family_residue_enrichment.tsv")
        r.table(enrich_pivot.reset_index(), "figure8_residue_log2_enrichment.tsv")
        r.table(class_pivot.reset_index(), "figure8_residue_class_share.tsv")
        save_clean(r, fig, fid, title, sources, "activity/regio-labeled residue and hotspot summaries", len(df))
    except Exception as exc:
        r.skip(fid, title, sources, str(exc))


def figure9(r: Renderer) -> None:
    fid, title = "figure9_predictive_performance.png", "Predictive performance and most informative features"
    source = "postprocess/10_predictive/cv_results"
    try:
        cv_dir = r.results_dir / source
        metrics = []
        predictions = []
        for path in sorted(cv_dir.glob("*_metrics.tsv")):
            if path.name.endswith("_fold_metrics.tsv"):
                continue
            m = pd.read_csv(path, sep="\t")
            m["model_name"] = path.name.replace("_metrics.tsv", "")
            metrics.append(m)
        for path in sorted(cv_dir.glob("*_predictions.tsv")):
            p = pd.read_csv(path, sep="\t")
            p["model_name"] = path.name.replace("_predictions.tsv", "")
            predictions.append(p)
        met = num(pd.concat(metrics, ignore_index=True), ["n_positive", "n_negative", "mean_pr_auc", "std_pr_auc", "mean_balanced_accuracy", "std_balanced_accuracy"])
        pred = num(pd.concat(predictions, ignore_index=True), ["y_true", "y_score"])
        roc_auc_rows = []
        for name, sub in pred.groupby("model_name", observed=True):
            y_true = pd.to_numeric(sub["y_true"], errors="coerce")
            y_score = pd.to_numeric(sub["y_score"], errors="coerce")
            ok = y_true.notna() & y_score.notna()
            if y_true[ok].nunique() < 2:
                roc_auc = np.nan
            else:
                roc_auc = roc_auc_score(y_true[ok], y_score[ok])
            roc_auc_rows.append({"model_name": name, "roc_auc": roc_auc})
        met = met.merge(pd.DataFrame(roc_auc_rows), on="model_name", how="left")

        fig = plt.figure(figsize=(12.6, 13.2), constrained_layout=True)
        gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.12, 1.25])
        ax_perf = fig.add_subplot(gs[0, :])
        ax_pr = fig.add_subplot(gs[1, 0])
        ax_roc = fig.add_subplot(gs[1, 1])
        ax_imp = fig.add_subplot(gs[2, :])
        for name, sub in pred.groupby("model_name", observed=True):
            precision, recall, _ = precision_recall_curve(sub["y_true"], sub["y_score"])
            line, = ax_pr.plot(recall, precision, label=name)
            baseline = sub["y_true"].mean()
            ax_pr.axhline(baseline, color=line.get_color(), linestyle="--", linewidth=1.0, alpha=0.65, label=f"{name} baseline ({baseline:.2f})")
            if sub["y_true"].nunique() > 1:
                fpr, tpr, _ = roc_curve(sub["y_true"], sub["y_score"])
                auc = float(met.loc[met["model_name"] == name, "roc_auc"].iloc[0])
                model_label = wrap_label(name.replace("_", " "), width=17)
                ax_roc.plot(fpr, tpr, color=line.get_color(), label=f"{model_label}\nROC-AUC {auc:.2f}")
        ax_pr.set_title("Precision-recall curves")
        ax_pr.set_xlabel("Recall")
        ax_pr.set_ylabel("Precision")
        ax_pr.legend(fontsize=LEGEND_FONTSIZE, loc="upper right", bbox_to_anchor=(-0.24, 1.0), frameon=False, borderaxespad=0.0, labelspacing=0.25)
        panel_label(ax_pr, "B")
        ax_roc.plot([0, 1], [0, 1], color="#777777", linestyle="--", linewidth=1.0, label="Random baseline")
        ax_roc.set_xlim(0, 1)
        ax_roc.set_ylim(0, 1)
        ax_roc.set_title("ROC curves")
        ax_roc.set_xlabel("False positive rate")
        ax_roc.set_ylabel("True positive rate")
        ax_roc.legend(fontsize=LEGEND_FONTSIZE, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, borderaxespad=0.0, labelspacing=0.35)
        panel_label(ax_roc, "C")
        met["baseline_pr_auc"] = met["n_positive"] / (met["n_positive"] + met["n_negative"])
        perf = met.melt(id_vars=["model_name"], value_vars=["mean_pr_auc", "baseline_pr_auc", "roc_auc", "mean_balanced_accuracy"], var_name="metric", value_name="value")
        perf["metric"] = perf["metric"].map({"mean_pr_auc": "AUPRC", "baseline_pr_auc": "Baseline AUPRC", "roc_auc": "ROC-AUC", "mean_balanced_accuracy": "Balanced accuracy"})
        sns.barplot(data=perf, y="model_name", x="value", hue="metric", ax=ax_perf)
        ax_perf.set_xlim(0, 1)
        ax_perf.set_title("Cross-validated performance")
        ax_perf.set_xlabel("Score")
        ax_perf.set_ylabel("")
        ax_perf.legend(title="Metric", frameon=False, loc="upper right", bbox_to_anchor=(-0.16, 1.0), borderaxespad=0.0)
        panel_label(ax_perf, "A")
        rows = []
        for _, row in met.iterrows():
            for feature, value in json.loads(str(row.get("feature_importances_json", "{}"))).items():
                rows.append({"model_name": row["model_name"], "feature": feature, "importance": value})
        imp = num(pd.DataFrame(rows), ["importance"])
        top = imp.groupby("feature", observed=True)["importance"].mean().sort_values(ascending=False).head(12).index
        sns.barplot(data=imp[imp["feature"].isin(top)], y="feature", x="importance", hue="model_name", ax=ax_imp)
        ax_imp.set_title("Top feature importances")
        ax_imp.set_xlabel("Importance")
        ax_imp.set_ylabel("")
        ax_imp.legend(title="Model name", frameon=False, loc="lower right")
        panel_label(ax_imp, "D")
        fig.suptitle(title, fontsize=15, weight="bold")
        r.table(met.assign(delta_vs_baseline=met["mean_pr_auc"] - met["baseline_pr_auc"], roc_auc_minus_random=met["roc_auc"] - 0.5), "figure9_predictive_metrics.tsv")
        save_clean(r, fig, fid, title, [source + "/*_metrics.tsv", source + "/*_predictions.tsv"], "grouped CV result files", len(met))
    except Exception as exc:
        r.skip(fid, title, [source], str(exc))


def figure10_within_cluster_rmsd(r: Renderer) -> None:
    fid, title = "figure10_within_cluster_substrate_rmsd.png", "Within-cluster substrate and IFP medoid distances"
    sources = ["combined/cluster_assignments.tsv", "combined/cluster_table.tsv", "combined/medoid_manifest.tsv", "*/shards/shard_*/pose_manifest.tsv"]
    try:
        cache_path = r.output_dir / "tables" / "figure10_within_cluster_medoid_pose_metrics.tsv"
        summary_cache_path = r.output_dir / "tables" / "figure10_within_cluster_medoid_cluster_metrics.tsv"
        cluster = r.tsv("combined/cluster_table.tsv")
        expected_constructs = set(cluster["construct_type"].dropna().astype(str).unique())
        use_cache = False
        if cache_path.exists() and summary_cache_path.exists():
            pose_metrics = pd.read_csv(cache_path, sep="\t")
            cluster_metrics = pd.read_csv(summary_cache_path, sep="\t")
            use_cache = expected_constructs.issubset(set(pose_metrics.get("construct_type", pd.Series(dtype=str)).dropna().astype(str).unique()))
        if not use_cache:
            assignments = num(r.tsv("combined/cluster_assignments.tsv"), ["cluster_id", "distance_to_cluster_representative"])
            assignments = assignments[bools(assignments["cluster_member_flag"]) & assignments["cluster_id"].ge(0)].copy()
            medoids = r.tsv("combined/medoid_manifest.tsv")
            pose_files = sorted(r.results_dir.glob("*/shards/shard_*/pose_manifest.tsv"))
            pose_manifest = pd.concat(
                [pd.read_csv(path, sep="\t", usecols=["condition_id", "pose_id", "normalized_cif"]) for path in pose_files],
                ignore_index=True,
            ).drop_duplicates(["condition_id", "pose_id"])
            pose_metrics = (
                assignments.merge(
                    cluster[["condition_id", "cluster_id", "construct_type", "substrate_class", "dp", "medoid_pose_id"]],
                    on=["condition_id", "cluster_id"],
                    how="left",
                )
                .merge(medoids[["condition_id", "cluster_id", "medoid_structure_path"]], on=["condition_id", "cluster_id"], how="left")
                .merge(pose_manifest.rename(columns={"normalized_cif": "pose_structure_path"}), on=["condition_id", "pose_id"], how="left")
            )
            rmsd_values = []
            for row in pose_metrics.itertuples(index=False):
                rmsd_values.append(ligand_rmsd_to_medoid(str(row.pose_structure_path), str(row.medoid_structure_path)))
            pose_metrics["ligand_rmsd_to_cluster_medoid"] = rmsd_values
            pose_metrics = finite(pose_metrics, ["ligand_rmsd_to_cluster_medoid", "distance_to_cluster_representative"])
            cluster_metrics = (
                pose_metrics.groupby(["condition_id", "cluster_id", "construct_type", "substrate_class", "dp"], observed=True)
                .agg(
                    n_member_poses=("pose_id", "size"),
                    ligand_rmsd_to_medoid_median=("ligand_rmsd_to_cluster_medoid", "median"),
                    medoid_ifp_distance_mae=("distance_to_cluster_representative", lambda s: float(np.mean(np.abs(s)))),
                    medoid_ifp_distance_rmsd=("distance_to_cluster_representative", lambda s: float(np.sqrt(np.mean(np.square(s))))),
                )
                .reset_index()
            )
            r.table(pose_metrics, "figure10_within_cluster_medoid_pose_metrics.tsv")
            r.table(cluster_metrics, "figure10_within_cluster_medoid_cluster_metrics.tsv")
        fig, axes = plt.subplots(1, 2, figsize=(10.8, 5.2), constrained_layout=True)
        light_palette = {"domain_only": "#b8cce4", "full_length": "#e8bd89"}
        dark_palette = {"domain_only": "#2f5f8f", "full_length": "#a85f17"}
        centered_violin_boxplot(
            axes[0],
            pose_metrics,
            value_col="ligand_rmsd_to_cluster_medoid",
            y_label="Substrate RMSD (Å)",
            light_palette=light_palette,
            dark_palette=dark_palette,
        )
        axes[0].set_title("Substrate RMSD to cluster medoid")
        panel_label(axes[0], "A")
        centered_violin_boxplot(
            axes[1],
            cluster_metrics,
            value_col="medoid_ifp_distance_mae",
            y_label="Mean absolute IFP distance to medoid (Jaccard)",
            light_palette=light_palette,
            dark_palette=dark_palette,
        )
        axes[1].set_title("Medoid IFP distance MAE (Jaccard)")
        axes[1].set_ylim(0, 1.0)
        panel_label(axes[1], "B")
        legend_handles = [Patch(facecolor=dark_palette[c], edgecolor="#303030", label=CONSTRUCT_DISPLAY[c]) for c in CONSTRUCT_ORDER]
        axes[0].legend(handles=legend_handles, title="Construct", frameon=False, fontsize=LEGEND_FONTSIZE, loc="upper right")
        axes[1].legend(handles=legend_handles, title="Construct", frameon=False, fontsize=LEGEND_FONTSIZE, loc="upper center")
        for ax in axes:
            compact_margins(ax, x=0.02)
        fig.suptitle(title, fontsize=15, weight="bold")
        save_clean(r, fig, fid, title, sources, "cluster-member poses with finite substrate RMSD to cluster medoid and Jaccard IFP distance to medoid", len(pose_metrics))
    except Exception as exc:
        r.skip(fid, title, sources, str(exc))


FIGURES = [
    figure1,
    figure2,
    figure3,
    figure3_angles_by_substrate,
    figure4,
    figure4_demonstrated_activity_substrate_boxplots,
    figure5,
    figure5_full_residue_heatmap,
    figure5_all_residue_heatmap,
    figure5_residue_contact_activity_difference_heatmap,
    figure5_contact_fraction_difference_by_substrate,
    figure6,
    figure6b_cbm_bridge_high_fraction,
    figure6c_cbm_delta_heatmap,
    figure6d_cbm_delta_heatmap_high_bridge,
    figure6c_all_cbm_delta_heatmap,
    figure7,
    figure7c_clusterability_bubble,
    figure8,
    figure9,
    figure10_within_cluster_rmsd,
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
    style()
    out = args.output_dir or args.results_dir / "postprocess" / "figures3"
    renderer = Renderer(args.results_dir, args.analyse_dir, out)
    for func in FIGURES:
        func(renderer)
    renderer.manifest()
    ok = sum(r.status == "ok" for r in renderer.records)
    skipped = sum(r.status != "ok" for r in renderer.records)
    print(f"Rendered {ok} figures; skipped {skipped}.")
    print(f"Manifest: {renderer.output_dir / 'figure_manifest.tsv'}")
    return 2 if args.strict and skipped else 0


if __name__ == "__main__":
    raise SystemExit(main())
