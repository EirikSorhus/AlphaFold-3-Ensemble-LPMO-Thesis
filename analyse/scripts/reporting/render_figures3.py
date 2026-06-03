#!/usr/bin/env python3
"""Render thesis-oriented main figure candidates into postprocess/figures3."""

from __future__ import annotations

import argparse
import json
import math
import re
import shlex
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
from sklearn.metrics import precision_recall_curve


SUBSTRATE_ORDER = ["chitin", "cellulose", "starch", "amylose"]
CONSTRUCT_ORDER = ["domain_only", "full_length"]
ACTIVITY_ORDER = ["Demonstrated activity substrate", "No demonstrated activity substrate"]
SUBSTRATE_PALETTE = {"chitin": "#2b8cbe", "cellulose": "#41ab5d", "starch": "#e34a33", "amylose": "#e34a33"}
CONSTRUCT_PALETTE = {"domain_only": "#3b6ea8", "full_length": "#d17a22"}
ACTIVITY_PALETTE = {"Demonstrated activity substrate": "#2f7f4f", "No demonstrated activity substrate": "#9b3a3a"}
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
        fig.savefig(out, dpi=240, bbox_inches="tight")
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
        font_scale=0.88,
        rc={
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
    ax.text(-0.11, 1.08, label, transform=ax.transAxes, ha="left", va="top", fontsize=13, weight="bold")


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


def annotate(ax: plt.Axes, text: str) -> None:
    ax.text(0.99, 0.98, text, transform=ax.transAxes, ha="right", va="top", fontsize=7.5)


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
                legend.get_title().set_text(sentence_case_label(legend.get_title().get_text()))
            for text in legend.get_texts():
                text.set_text(sentence_case_label(text.get_text()))
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
            legend.get_title().set_text(sentence_case_label(legend.get_title().get_text()))
        for text in legend.get_texts():
            text.set_text(sentence_case_label(text.get_text()))


def save_clean(r: Renderer, fig: plt.Figure, figure_id: str, title: str, sources: list[str], filters: str, n: int) -> None:
    for ax in fig.axes:
        ax.set_title("")
    if getattr(fig, "_suptitle", None) is not None:
        fig._suptitle.set_text("")
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
        fig, axes = plt.subplots(1, 2, figsize=(13.8, 5.7), constrained_layout=True, gridspec_kw={"width_ratios": [1.0, 1.85]})
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
        ax2.set_yticklabels(stage_order)
        ax2.invert_yaxis()
        ax2.set_title("Cumulative pass count through hard QC and downstream steps")
        ax2.set_xlabel("Count passing step")
        ax2.set_ylabel("")
        ax2.legend(title="Construct", frameon=False, loc="lower right")
        panel_label(ax2, "B")
        fig.suptitle(title, fontsize=14, weight="bold")
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
        fig, axes = plt.subplots(5, 1, figsize=(8.2, 17.8), constrained_layout=True)
        rate = df.groupby(["construct_type", "substrate_dp"], observed=True).agg(fraction=("valid_cluster", "mean"), n=("condition_id", "nunique")).reset_index()
        ax_c, ax_a, ax_b, ax_d, ax_e = axes
        counts = df["n_clusters"].dropna().astype(int)
        cluster_count = counts.value_counts().sort_index()
        cluster_count = cluster_count.reindex(range(0, int(cluster_count.index.max()) + 1), fill_value=0)
        ax_c.bar(cluster_count.index, cluster_count.values, color="#587a9d")
        ax_c.set_yscale("log", base=2)
        ax_c.set_ylim(0.5, max(cluster_count.max() * 1.9, 2))
        for x_value, y_value in cluster_count.items():
            ax_c.text(x_value, max(y_value, 0.55) * 1.08, str(int(y_value)), ha="center", va="bottom", fontsize=8)
        ax_c.set_xlabel("Number of clusters")
        ax_c.set_ylabel("Conditions")
        ax_c.set_xticks(cluster_count.index)
        panel_label(ax_c, "A")
        sns.pointplot(data=rate, x="substrate_dp", y="fraction", hue="construct_type", palette=CONSTRUCT_PALETTE, ax=ax_a)
        ax_a.set_ylim(-0.03, 1.03)
        ax_a.tick_params(axis="x", rotation=35)
        ax_a.set_xlabel("")
        ax_a.set_ylabel("Fractions of conditions clusterable")
        annotate(ax_a, f"n conditions={df['condition_id'].nunique()}")
        panel_label(ax_a, "B")
        occ = df[df["n_clusters"].fillna(0).gt(1) & df["top_cluster_occupancy"].gt(0)].copy()
        sns.histplot(data=occ, x="top_cluster_occupancy", hue="construct_type", bins=25, stat="count", element="step", palette=CONSTRUCT_PALETTE, ax=ax_b)
        ax_b.set_xlabel("Top-cluster occupancy")
        ax_b.set_ylabel("Count")
        annotate(ax_b, "multi-cluster conditions only (n clusters > 1)")
        panel_label(ax_b, "C")
        ent = finite(df, ["cluster_entropy"])
        sns.violinplot(data=ent, x="substrate_class", y="cluster_entropy", hue="construct_type", cut=0, inner="quartile", palette=CONSTRUCT_PALETTE, ax=ax_d)
        ax_d.set_xlabel("Substrate")
        ax_d.set_ylabel("Normalized entropy")
        panel_label(ax_d, "D")
        noise = finite(df, ["noise_fraction"])
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
        panel_label(ax_e, "E")
        fig.suptitle(title, fontsize=14, weight="bold")
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
        fig, axes = plt.subplot_mosaic([["A", "A"], ["B", "B"], ["D1", "D2"], ["C", "C"]], figsize=(13.8, 16.5), constrained_layout=True)
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
        axes["A"].legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
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
        axes["B"].legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
        panel_label(axes["B"], "B")
        metrics = ["occupancy_weighted_c1_plausible_fraction", "occupancy_weighted_c4_plausible_fraction"]
        cond = num(cond, metrics + ["occupancy_weighted_Cu_oxyl_H_C1_angle_median", "occupancy_weighted_Cu_oxyl_H_C4_angle_median", "occupancy_weighted_Cu_C1_distance_median", "occupancy_weighted_Cu_C4_distance_median"])
        long = cond.melt(id_vars=["protein_id", "activity_evidence"], value_vars=metrics, var_name="site", value_name="fraction")
        long["site"] = long["site"].map({metrics[0]: "C1", metrics[1]: "C4"})
        sns.boxplot(data=long, x="activity_evidence", y="fraction", hue="site", order=ACTIVITY_ORDER, showfliers=False, ax=axes["C"])
        sns.stripplot(data=long, x="activity_evidence", y="fraction", hue="site", order=ACTIVITY_ORDER, dodge=True, color="#222222", alpha=0.2, size=2, ax=axes["C"], legend=False)
        axes["C"].set_title("Plausible geometry fraction")
        axes["C"].set_xlabel("")
        axes["C"].set_ylabel("Occupancy-weighted fraction")
        axes["C"].tick_params(axis="x", rotation=15)
        axes["C"].legend(title="Site", frameon=False, loc="upper right")
        panel_label(axes["C"], "E")
        angle_metrics = ["occupancy_weighted_Cu_oxyl_H_C1_angle_median", "occupancy_weighted_Cu_oxyl_H_C4_angle_median"]
        angle = cond.melt(id_vars=["activity_evidence"], value_vars=angle_metrics, var_name="metric", value_name="value")
        angle["metric"] = angle["metric"].map(nice)
        sns.boxplot(data=finite(angle, ["value"]), y="metric", x="value", hue="activity_evidence", hue_order=ACTIVITY_ORDER, showfliers=False, palette=ACTIVITY_PALETTE, ax=axes["D1"])
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
        axes["D2"].set_title("Cu-distance metrics")
        axes["D2"].set_xlabel("Distance (Å)")
        axes["D2"].set_ylabel("")
        if axes["D2"].get_legend() is not None:
            axes["D2"].get_legend().remove()
        axes["D2"].legend(handles, labels, title="Activity evidence", frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=7)
        panel_label(axes["D2"], "D")
        fig.suptitle(title, fontsize=14, weight="bold")
        save_clean(r, fig, fid, title, sources + ["input_data/metadata_final_ec_fixed.tsv"], "cluster medians and EC-activity-annotated conditions", len(clusters))
    except Exception as exc:
        r.skip(fid, title, sources, str(exc))


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
        long = df.melt(id_vars=["protein_id", "activity_evidence"], value_vars=metrics, var_name="metric", value_name="value")
        long["metric"] = long["metric"].map(nice)
        long = finite(long, ["value"])
        metric_order = list(long["metric"].drop_duplicates())
        fig, axes = plt.subplots(2, 3, figsize=(12.0, 7.8), constrained_layout=True)
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
            ax.set_xticklabels(["Demonstrated activity", "No demonstrated activity"])
            panel_label(ax, chr(ord("A") + idx))
        fig.suptitle(title, y=1.03, fontsize=14, weight="bold")
        summary = long.groupby(["metric", "activity_evidence"], observed=True).agg(n=("value", "size"), median=("value", "median"), proteins=("protein_id", "nunique")).reset_index()
        r.table(summary, "figure4_activity_metric_summary.tsv")
        save_clean(r, fig, fid, title, [source, "input_data/metadata_final_ec_fixed.tsv"], "EC-activity-annotated conditions", len(long))
    except Exception as exc:
        r.skip(fid, title, [source], str(exc))


def figure5(r: Renderer) -> None:
    fid, title = "figure5_contact_chemistry.png", "Contact chemistry by substrate class, DP, and activity evidence"
    sources = ["combined/condition_table.tsv", "combined/protein_condition_residue_scores.tsv"]
    try:
        df = cats(add_activity(r, r.tsv(sources[0])))
        contact_metrics = ["aromatic_contact_fraction", "polar_contact_fraction", "charged_contact_fraction", "hbond_contact_fraction"]
        df = num(df, contact_metrics)
        long = df.melt(id_vars=["protein_id", "substrate_dp", "substrate_class", "activity_evidence"], value_vars=contact_metrics, var_name="contact_type", value_name="fraction")
        long["contact_type"] = long["contact_type"].str.replace("_contact_fraction", "", regex=False).str.replace("hbond", "hydrogen bond")
        mean = long.groupby(["activity_evidence", "substrate_dp", "contact_type"], observed=True)["fraction"].mean().reset_index()
        fig = plt.figure(figsize=(12.8, 10.4), constrained_layout=True)
        gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.25])
        ax_a = fig.add_subplot(gs[0, 0])
        ax_b = fig.add_subplot(gs[0, 1])
        ax_c = fig.add_subplot(gs[1, :])
        sns.barplot(data=mean, x="substrate_dp", y="fraction", hue="contact_type", ax=ax_a)
        ax_a.set_title("Mean contact fractions")
        ax_a.set_xlabel("")
        ax_a.set_ylabel("Occupancy-weighted fraction")
        ax_a.tick_params(axis="x", rotation=35)
        panel_label(ax_a, "A")
        grouped = df.groupby("activity_evidence", observed=True)[contact_metrics].agg(["mean", "std"])
        diff_rows = []
        for metric in contact_metrics:
            active_mean = grouped.loc[ACTIVITY_ORDER[0], (metric, "mean")] if ACTIVITY_ORDER[0] in grouped.index else np.nan
            inactive_mean = grouped.loc[ACTIVITY_ORDER[1], (metric, "mean")] if ACTIVITY_ORDER[1] in grouped.index else np.nan
            active_sd = grouped.loc[ACTIVITY_ORDER[0], (metric, "std")] if ACTIVITY_ORDER[0] in grouped.index else np.nan
            inactive_sd = grouped.loc[ACTIVITY_ORDER[1], (metric, "std")] if ACTIVITY_ORDER[1] in grouped.index else np.nan
            diff_rows.append({"contact_type": metric, "delta": active_mean - inactive_mean, "sd_active": active_sd, "sd_no_activity": inactive_sd})
        diff = pd.DataFrame(diff_rows)
        diff["contact_type"] = diff["contact_type"].str.replace("_contact_fraction", "", regex=False).str.replace("hbond", "hydrogen bond")
        y_pos = np.arange(len(diff))
        ax_b.errorbar(
            diff["delta"],
            y_pos,
            xerr=np.vstack([diff["sd_no_activity"].fillna(0).to_numpy(), diff["sd_active"].fillna(0).to_numpy()]),
            fmt="o",
            color="#222222",
            ecolor="#555555",
            elinewidth=1.0,
            capsize=3,
        )
        ax_b.set_yticks(y_pos)
        ax_b.set_yticklabels(diff["contact_type"])
        ax_b.axvline(0, color="#777777", linestyle="--")
        ax_b.set_title("Demonstrated activity - no demonstrated activity")
        ax_b.set_xlabel("Mean fraction difference with group SD")
        ax_b.set_ylabel("")
        panel_label(ax_b, "B")
        res = cats(r.tsv(sources[1]))
        res = num(res, ["residue_contact_score"])
        top = res.groupby("residue_label", observed=True)["residue_contact_score"].mean().sort_values(ascending=False).head(35).index
        res = res[res["residue_label"].isin(top)].copy()
        pivot = res.pivot_table(index="residue_label", columns="substrate_class", values="residue_contact_score", aggfunc="mean")
        pivot.index = [strip_residue_chain(label) for label in pivot.index]
        sns.heatmap(pivot, cmap="mako", ax=ax_c, cbar_kws={"label": "Contact score"})
        ax_c.set_title("Top 35 residue contact scores")
        ax_c.set_xlabel("Substrate")
        ax_c.set_ylabel("Residue")
        panel_label(ax_c, "C")
        fig.suptitle(title, fontsize=14, weight="bold")
        r.table(diff, "figure5_contact_fraction_difference.tsv")
        save_clean(r, fig, fid, title, sources, "non-VdW contact fractions; top residue scores", len(long))
    except Exception as exc:
        r.skip(fid, title, sources, str(exc))


def figure5_full_residue_heatmap(r: Renderer) -> None:
    fid, title = "figure5c_full_residue_contact_heatmap.png", "Full residue contact heatmap"
    source = "combined/protein_condition_residue_scores.tsv"
    try:
        res = cats(r.tsv(source))
        res = finite(num(res, ["residue_contact_score"]), ["residue_contact_score"])
        pivot_all = res.pivot_table(index="residue_label", columns="substrate_class", values="residue_contact_score", aggfunc="mean")
        selected: set[str] = set()
        for substrate in [s for s in SUBSTRATE_ORDER if s in pivot_all.columns]:
            selected.update(pivot_all[substrate].dropna().sort_values(ascending=False).head(20).index.astype(str))
        pivot = pivot_all.loc[pivot_all.index.astype(str).isin(selected)].copy()
        pivot = pivot.loc[pivot.max(axis=1).sort_values(ascending=False).index]
        pivot.index = [strip_residue_chain(label) for label in pivot.index]
        height = max(9.0, min(24.0, 0.18 * len(pivot) + 3.5))
        fig, ax = plt.subplots(figsize=(7.5, height), constrained_layout=True)
        sns.heatmap(pivot, cmap="mako", ax=ax, yticklabels=True, cbar_kws={"label": "Mean contact score"})
        ax.set_title("Top 20 residue contacts per substrate")
        ax.set_xlabel("Substrate")
        ax.set_ylabel("Residue")
        ax.tick_params(axis="y", labelsize=6.5 if len(pivot) > 45 else 7.5)
        r.table(pivot.reset_index(), "figure5c_full_residue_contact_heatmap.tsv")
        save_clean(r, fig, fid, title, [source], "union of top 20 residue labels by mean contact score within each substrate", len(pivot))
    except Exception as exc:
        r.skip(fid, title, [source], str(exc))


def figure6(r: Renderer) -> None:
    fid, title = "figure6_cbm_paired_analysis.png", "CBM paired analysis: catalytic domain versus full-length"
    sources = ["postprocess/15_cbm_paired_analysis/cbm_construct_condition_summary.tsv", "postprocess/15_cbm_paired_analysis/cbm_paired_comparison_table.tsv"]
    try:
        summary = cats(r.tsv(sources[0]))
        paired = cats(r.tsv(sources[1]))
        need(paired, ["protein_id", "domain_only_condition_id", "full_length_condition_id"])
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
        fig, axes_grid = plt.subplots(2, 3, figsize=(13.2, 7.8), constrained_layout=True)
        axes = axes_grid.flat
        ax0 = axes[0]
        for _, pair in matched.groupby("pair_id", observed=True):
            if pair["expected_construct"].nunique() == 2:
                pair = pair.sort_values("expected_construct", key=lambda s: s.map({"domain_only": 0, "full_length": 1}))
                ax0.plot(pair["expected_construct"], pair["valid_cluster"].astype(float), color="#999999", alpha=0.16, linewidth=0.8)
        sns.pointplot(data=matched, x="expected_construct", y="valid_cluster", order=CONSTRUCT_ORDER, color="#000000", errorbar=("pi", 50), ax=ax0)
        ax0.set_ylim(-0.03, 1.03)
        ax0.set_xlabel("")
        ax0.set_ylabel("Valid-Cluster Fraction")
        ax0.set_xticks([0, 1])
        ax0.set_xticklabels(["Catalytic domain", "Full-length"])
        annotate(ax0, f"n proteins={paired['protein_id'].nunique()}\nn paired conditions={len(paired)}")
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
            panel_label(ax, chr(ord("B") + idx))
        axes_grid[1, 2].set_axis_off()
        if "non_core_ligand_contact_fraction_full_length" in paired.columns:
            contact = finite(num(paired, ["non_core_ligand_contact_fraction_full_length", "bridge_fraction_full_length"]), ["non_core_ligand_contact_fraction_full_length"])
            summary_rows = contact.groupby("substrate_class", observed=True).agg(
                n_pairs=("protein_id", "size"),
                n_nonzero_non_core=("non_core_ligand_contact_fraction_full_length", lambda s: int((s > 0).sum())),
                median_non_core=("non_core_ligand_contact_fraction_full_length", "median"),
                max_non_core=("non_core_ligand_contact_fraction_full_length", "max"),
                n_nonzero_bridge=("bridge_fraction_full_length", lambda s: int((pd.to_numeric(s, errors="coerce") > 0).sum()) if "bridge_fraction_full_length" in contact.columns else 0),
            ).reset_index()
            r.table(summary_rows, "figure6_non_core_contact_summary.tsv")
        fig.suptitle(title, fontsize=14, weight="bold")
        save_clean(r, fig, fid, title, sources, "paired CBM comparison tables", len(paired))
    except Exception as exc:
        r.skip(fid, title, sources, str(exc))


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
        fig, ax = plt.subplots(figsize=(9.5, max(8.0, min(62.0, 0.14 * len(pivot) + 3.0))), constrained_layout=True)
        sns.heatmap(pivot.rename(columns=nice), cmap="vlag", center=0, vmin=-max_abs, vmax=max_abs, yticklabels=True, ax=ax, cbar_kws={"label": "Full-length - catalytic domain"})
        ax.set_title(title)
        ax.set_xlabel("")
        ax.set_ylabel("Paired condition")
        ax.tick_params(axis="y", labelsize=3.5 if len(pivot) > 200 else 5.5)
        r.table(pivot.reset_index(), "figure6c_cbm_delta_heatmap.tsv")
        save_clean(r, fig, fid, title, [source], "union of top 10 positive and top 10 negative full-length minus domain-only changes per metric", len(pivot))
    except Exception as exc:
        r.skip(fid, title, [source], str(exc))


def figure7(r: Renderer) -> None:
    fid, title = "figure7_crystal_anchoring.png", "Reference structure comparison"
    sources = ["combined/crystal_anchor_table.tsv", "combined/condition_table.tsv", "input_data/metadata_pdb_count.tsv"]
    try:
        anchor = r.tsv(sources[0])
        cond = cats(r.tsv(sources[1]))
        df = cats(anchor.merge(cond[["condition_id", "substrate_class", "dp", "construct_type"]].drop_duplicates("condition_id"), on="condition_id", how="left"))
        fig = plt.figure(figsize=(12.0, 9.2), constrained_layout=True)
        gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05])
        axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])]
        rmsd = finite(num(df, ["local_pocket_rmsd"]), ["local_pocket_rmsd"])
        sns.boxplot(data=rmsd, x="substrate_class", y="local_pocket_rmsd", hue="construct_type", showfliers=False, palette=CONSTRUCT_PALETTE, ax=axes[0])
        sns.stripplot(data=rmsd, x="substrate_class", y="local_pocket_rmsd", hue="construct_type", dodge=True, color="#111111", alpha=0.45, size=2.2, ax=axes[0], legend=False)
        axes[0].set_title("Local pocket RMSD")
        axes[0].set_xlabel("Substrate")
        axes[0].set_ylabel("RMSD (Å)")
        panel_label(axes[0], "A")
        sim = finite(num(df, ["local_pocket_rmsd", "ifp_tanimoto"]), ["local_pocket_rmsd", "ifp_tanimoto"])
        sns.scatterplot(data=sim, x="local_pocket_rmsd", y="ifp_tanimoto", hue="substrate_class", palette=SUBSTRATE_PALETTE, ax=axes[1])
        axes[1].set_ylim(-0.03, 1.03)
        axes[1].set_title("IFP similarity versus RMSD")
        axes[1].set_xlabel("Local pocket RMSD (Å)")
        axes[1].set_ylabel("IFP Tanimoto")
        axes[1].legend(fontsize=7, loc="upper right", frameon=False)
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
        annotate(axes[2], f"linear R²={r2:.2f}" if np.isfinite(r2) else "linear R²=NA")
        panel_label(axes[2], "C")
        fig.suptitle(title, fontsize=14, weight="bold")
        save_clean(r, fig, fid, title, sources, "crystal anchor rows and protein-level clusterability", len(df))
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
        fig, axes = plt.subplots(4, 1, figsize=(9.8, 18.5), constrained_layout=True)
        ax_a, ax_b, ax_c, ax_d = axes
        sns.heatmap(identity_pivot, cmap="YlGnBu", ax=ax_a, cbar_kws={"label": "Contact-weighted residue share (%)"})
        ax_a.set_title("Experimental residue identity by label")
        ax_a.set_xlabel("")
        ax_a.set_ylabel("Residue type")
        panel_label(ax_a, "A")
        max_abs = max(1.0, float(np.nanmax(np.abs(enrich_pivot.to_numpy()))) if enrich_pivot.size else 1.0)
        sns.heatmap(enrich_pivot, cmap="vlag", center=0, vmin=-max_abs, vmax=max_abs, ax=ax_b, cbar_kws={"label": "log2 enrichment"})
        ax_b.set_title("Experimental residue identity log2 enrichment")
        ax_b.set_xlabel("")
        ax_b.set_ylabel("Residue type")
        panel_label(ax_b, "B")
        sns.heatmap(hotspot_pivot, cmap="mako", ax=ax_c, cbar_kws={"label": "Contact-weighted score"})
        ax_c.set_title("Experimental aligned hotspot signal")
        ax_c.set_xlabel("")
        ax_c.set_ylabel("Family + consensus residue + alignment column")
        ax_c.tick_params(axis="y", labelsize=6.5)
        panel_label(ax_c, "C")
        sns.heatmap(class_pivot, cmap="YlOrBr", ax=ax_d, cbar_kws={"label": "Contact-weighted residue class share (%)"})
        ax_d.set_title("Residue class by experimental label")
        ax_d.set_xlabel("")
        ax_d.set_ylabel("Residue class")
        panel_label(ax_d, "D")
        fig.suptitle(title, fontsize=14, weight="bold")
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
        fig, axes = plt.subplots(3, 1, figsize=(9.4, 14.2), constrained_layout=True)
        ax_perf, ax_pr, ax_imp = axes
        for name, sub in pred.groupby("model_name", observed=True):
            precision, recall, _ = precision_recall_curve(sub["y_true"], sub["y_score"])
            line, = ax_pr.plot(recall, precision, label=name)
            baseline = sub["y_true"].mean()
            ax_pr.axhline(baseline, color=line.get_color(), linestyle="--", linewidth=1.0, alpha=0.65, label=f"{name} baseline ({baseline:.2f})")
        ax_pr.set_title("Precision-recall curves")
        ax_pr.set_xlabel("Recall")
        ax_pr.set_ylabel("Precision")
        ax_pr.legend(fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
        panel_label(ax_pr, "B")
        met["baseline_pr_auc"] = met["n_positive"] / (met["n_positive"] + met["n_negative"])
        perf = met.melt(id_vars=["model_name"], value_vars=["mean_pr_auc", "baseline_pr_auc", "mean_balanced_accuracy"], var_name="metric", value_name="value")
        perf["metric"] = perf["metric"].map({"mean_pr_auc": "AUPRC", "baseline_pr_auc": "Baseline AUPRC", "mean_balanced_accuracy": "Balanced accuracy"})
        sns.barplot(data=perf, y="model_name", x="value", hue="metric", ax=ax_perf)
        ax_perf.set_xlim(0, 1)
        ax_perf.set_title("Cross-validated performance")
        ax_perf.set_xlabel("Score")
        ax_perf.set_ylabel("")
        ax_perf.legend(title="Metric", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
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
        panel_label(ax_imp, "C")
        fig.suptitle(title, fontsize=14, weight="bold")
        r.table(met.assign(delta_vs_baseline=met["mean_pr_auc"] - met["baseline_pr_auc"]), "figure9_predictive_metrics.tsv")
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
        fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.5), constrained_layout=True)
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
        panel_label(axes[1], "B")
        legend_handles = [Patch(facecolor=dark_palette[c], edgecolor="#303030", label=c) for c in CONSTRUCT_ORDER]
        for ax in axes:
            ax.legend(handles=legend_handles, title="Construct", frameon=False, fontsize=7)
        fig.suptitle(title, fontsize=14, weight="bold")
        save_clean(r, fig, fid, title, sources, "cluster-member poses with finite substrate RMSD to cluster medoid and Jaccard IFP distance to medoid", len(pose_metrics))
    except Exception as exc:
        r.skip(fid, title, sources, str(exc))


FIGURES = [
    figure1,
    figure2,
    figure3,
    figure4,
    figure5,
    figure5_full_residue_heatmap,
    figure6,
    figure6c_cbm_delta_heatmap,
    figure7,
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
