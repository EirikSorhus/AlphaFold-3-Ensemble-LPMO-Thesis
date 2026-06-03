#!/usr/bin/env python3
"""Render selected ligand-focused PyMOL panels for full-length medoid examples."""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from PIL import Image, ImageDraw, ImageFont


RESULTS_DEFAULT = Path("results/2full_pipeline_array_rerun_downstream_20260528_164450")


@dataclass(frozen=True)
class RenderCase:
    protein_id: str
    substrate_class: str
    dp: int
    cluster_id: int
    pose_id: str
    structure_path: Path
    cbm_start: int
    cbm_end: int

    @property
    def condition_label(self) -> str:
        return f"{self.protein_id} {self.substrate_class} DP{self.dp}"

    @property
    def condition_id(self) -> str:
        return f"{self.protein_id}__full_length__{self.substrate_class}_DP{self.dp}"

    @property
    def short_label(self) -> str:
        return f"{self.condition_label} cluster {self.cluster_id}"


CASES = [
    RenderCase(
        protein_id="A0A6C0PI79",
        substrate_class="cellulose",
        dp=6,
        cluster_id=0,
        pose_id="A0A6C0PI79_CEL6_seed-3_sample-4_model",
        structure_path=Path("/cluster/work/projects/nn1003k/eirik/Masteroppgave_clean/analyse/results/2full_pipeline_array/full_length/shards/shard_002/cases/0465_A0A6C0PI79_CEL6_seed-3_sample-4_model/normalize/normalized.cif"),
        cbm_start=288,
        cbm_end=316,
    ),
    RenderCase(
        protein_id="A0A6C0PI79",
        substrate_class="cellulose",
        dp=6,
        cluster_id=1,
        pose_id="A0A6C0PI79_CEL6_seed-3_sample-2_model",
        structure_path=Path("/cluster/work/projects/nn1003k/eirik/Masteroppgave_clean/analyse/results/2full_pipeline_array/full_length/shards/shard_002/cases/0463_A0A6C0PI79_CEL6_seed-3_sample-2_model/normalize/normalized.cif"),
        cbm_start=288,
        cbm_end=316,
    ),
    RenderCase(
        protein_id="M2QEW2",
        substrate_class="cellulose",
        dp=8,
        cluster_id=0,
        pose_id="M2QEW2_CEL8_seed-15_sample-2_model",
        structure_path=Path("/cluster/work/projects/nn1003k/eirik/Masteroppgave_clean/analyse/results/2full_pipeline_array/full_length/shards/shard_006/cases/1123_M2QEW2_CEL8_seed-15_sample-2_model/normalize/normalized.cif"),
        cbm_start=311,
        cbm_end=339,
    ),
    RenderCase(
        protein_id="M2QEW2",
        substrate_class="cellulose",
        dp=8,
        cluster_id=1,
        pose_id="M2QEW2_CEL8_seed-14_sample-2_model",
        structure_path=Path("/cluster/work/projects/nn1003k/eirik/Masteroppgave_clean/analyse/results/2full_pipeline_array/full_length/shards/shard_006/cases/1118_M2QEW2_CEL8_seed-14_sample-2_model/normalize/normalized.cif"),
        cbm_start=311,
        cbm_end=339,
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DEFAULT)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Defaults to <results-dir>/postprocess/figures/structures/pymol_selected_examples.",
    )
    parser.add_argument("--width", type=int, default=1800)
    parser.add_argument("--height", type=int, default=1400)
    parser.add_argument("--ray", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument(
        "--whole-molecule",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Use a very wide whole-molecule camera with a large effective slab.",
    )
    return parser.parse_args()


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def read_ifp_contacts(results_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    contacts = pd.read_csv(results_dir / "combined/pose_residue_contact_table.tsv", sep="\t", low_memory=False)
    ifp = pd.read_csv(results_dir / "combined/pose_ifp_table.tsv", sep="\t", low_memory=False)
    contacts = contacts[contacts["pose_id"].isin([case.pose_id for case in CASES])].copy()
    contacts = contacts[contacts["_construct_type"].astype(str).eq("full_length")].copy()
    if "contact_present" in contacts.columns:
        contacts = contacts[contacts["contact_present"].astype(str).str.lower().isin(["true", "1", "yes"])]
    ifp = ifp[ifp["pose_id"].isin([case.pose_id for case in CASES])].copy()
    ifp = ifp[ifp["_construct_type"].astype(str).eq("full_length")].copy()
    return contacts, ifp


def ifp_residue_set(contacts: pd.DataFrame, pose_id: str) -> set[tuple[str, int]]:
    sub = contacts[contacts["pose_id"] == pose_id].copy()
    return {
        (str(row.residue_chain), int(row.residue_number))
        for row in sub.itertuples()
        if pd.notna(row.residue_chain) and pd.notna(row.residue_number)
    }


def ifp_interaction_summary(ifp: pd.DataFrame, pose_id: str) -> dict[str, int]:
    row = ifp[ifp["pose_id"] == pose_id]
    if row.empty:
        return {}
    raw = row.iloc[0].get("ifp_interaction_occurrence_counts", "")
    if not isinstance(raw, str) or not raw.strip():
        return {}
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError:
        return {}
    return {str(key): int(value) for key, value in payload.items()}


def launch_pymol():
    import pymol

    pymol.finish_launching(["pymol", "-cq"])
    from pymol import cmd

    return cmd


def get_pymol_contact_residues(cmd, obj_name: str) -> set[tuple[str, int]]:
    model = cmd.get_model(f"{obj_name} and polymer.protein and byres (polymer.protein within 4.0 of ligand)")
    residues: set[tuple[str, int]] = set()
    for atom in model.atom:
        try:
            residues.add((atom.chain, int(atom.resi)))
        except ValueError:
            continue
    return residues


def selection_from_residues(obj_name: str, residues: set[tuple[str, int]]) -> str:
    if not residues:
        return "none"
    parts = [f"(chain {chain} and resi {resi})" for chain, resi in sorted(residues)]
    return f"{obj_name} and (" + " or ".join(parts) + ")"


def draw_case(
    cmd,
    case: RenderCase,
    contacts: pd.DataFrame,
    image_path: Path,
    *,
    surface: bool,
    width: int,
    height: int,
    ray: bool,
    whole_molecule: bool,
) -> dict:
    obj = "pose_obj"
    cmd.reinitialize()
    cmd.load(str(case.structure_path), obj)
    cmd.remove("solvent")
    cmd.remove("hydrogens")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("antialias", 2)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("stick_radius", 0.18)
    cmd.set("sphere_scale", 0.45)
    cmd.set("dash_gap", 0.25)
    cmd.set("dash_radius", 0.045)
    cmd.set("label_size", 18)
    cmd.set("label_color", "black")
    cmd.hide("everything", obj)

    cmd.select("protein_sel", f"{obj} and polymer.protein")
    cmd.select("cbm_domain", f"{obj} and polymer.protein and chain A and resi {case.cbm_start}-{case.cbm_end}")
    cmd.select("core_domain", f"{obj} and polymer.protein and chain A and resi 1-{case.cbm_start - 1}")
    cmd.select("non_cbm_other", f"{obj} and polymer.protein and not cbm_domain and not core_domain")
    cmd.select("ligand", f"{obj} and organic and not polymer.protein")
    cmd.select("cu", f"{obj} and resn CU")
    cmd.select("ligand_heavy", "ligand and not hydro")

    cmd.show("cartoon", "protein_sel")
    cmd.color("slate", "core_domain")
    cmd.color("lightblue", "non_cbm_other")
    cmd.color("teal", "cbm_domain")
    cmd.show("surface", "protein_sel")
    cmd.set("transparency", 0.68, "protein_sel")
    cmd.color("gray70", "core_domain")
    cmd.color("palecyan", "non_cbm_other")
    cmd.color("cyan", "cbm_domain")

    cmd.show("sticks", "ligand")
    cmd.color("orange", "ligand")
    cmd.show("spheres", "cu")
    cmd.color("copper", "cu")

    ifp_residues = ifp_residue_set(contacts, case.pose_id)
    contact_sel = selection_from_residues(obj, ifp_residues)
    cmd.select("ifp_contact_residues", contact_sel)
    cmd.show("sticks", "ifp_contact_residues")
    cmd.color("hotpink", "ifp_contact_residues")

    # Highlight PyMOL-detected nearby residues that were not in IFP.
    pymol_contacts = get_pymol_contact_residues(cmd, obj)
    extra_pymol = pymol_contacts - ifp_residues
    if extra_pymol:
        cmd.select("pymol_only_contacts", selection_from_residues(obj, extra_pymol))
        cmd.show("sticks", "pymol_only_contacts")
        cmd.color("yelloworange", "pymol_only_contacts")

    cmd.select("polar_contacts", "(ligand_heavy and elem O+N) within 3.5 of (protein and elem O+N)")
    cmd.distance("polar_distances", "ligand_heavy and elem O+N", "protein and elem O+N within 3.5 of ligand_heavy", cutoff=3.5, mode=2)
    cmd.color("gold", "polar_distances")
    cmd.hide("labels", "polar_distances")

    for atom_name, color in [("C1", "magenta"), ("C4", "purple")]:
        if cmd.count_atoms(f"ligand and name {atom_name}") and cmd.count_atoms("cu"):
            dist_name = f"Cu_{atom_name}"
            cmd.distance(dist_name, "cu", f"ligand and name {atom_name}", cutoff=12.0)
            cmd.color(color, dist_name)

    if whole_molecule:
        view_selection = "all"
        cmd.orient(view_selection)
        cmd.zoom(view_selection, 18, complete=1)
        cmd.clip("slab", 5000)
    else:
        view_selection = "protein_sel or ligand or cu"
        cmd.orient(view_selection)
        cmd.zoom(view_selection, 6, complete=1)
        cmd.clip("slab", 1000)
    cmd.png(str(image_path), width=width, height=height, dpi=220, ray=1 if ray else 0)
    return {
        "pose_id": case.pose_id,
        "condition_id": case.condition_id,
        "cluster_id": case.cluster_id,
        "ifp_contact_residue_count": len(ifp_residues),
        "pymol_4A_contact_residue_count": len(pymol_contacts),
        "intersection_count": len(ifp_residues & pymol_contacts),
        "ifp_only": residue_list(ifp_residues - pymol_contacts),
        "pymol_only": residue_list(extra_pymol),
        "shared": residue_list(ifp_residues & pymol_contacts),
    }


def residue_list(residues: set[tuple[str, int]]) -> str:
    return ",".join(f"{chain}:{resi}" for chain, resi in sorted(residues))


def add_label(image: Image.Image, label: str) -> Image.Image:
    image = image.convert("RGB")
    draw = ImageDraw.Draw(image)
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", 34)
    except OSError:
        font = ImageFont.load_default()
    box_h = 64
    draw.rectangle([(0, 0), (image.width, box_h)], fill=(255, 255, 255))
    draw.text((22, 16), label, fill=(20, 20, 20), font=font)
    return image


def make_panel(image_paths: list[Path], labels: list[str], output_path: Path, columns: int) -> None:
    labeled = [add_label(Image.open(path), label) for path, label in zip(image_paths, labels, strict=True)]
    cell_w = max(image.width for image in labeled)
    cell_h = max(image.height for image in labeled)
    rows = math.ceil(len(labeled) / columns)
    panel = Image.new("RGB", (columns * cell_w, rows * cell_h), "white")
    for idx, image in enumerate(labeled):
        row, col = divmod(idx, columns)
        panel.paste(image, (col * cell_w, row * cell_h))
    panel.save(output_path)


def main() -> int:
    args = parse_args()
    results_dir = args.results_dir.resolve()
    output_dir = (
        args.output_dir
        or results_dir / "postprocess/figures/structures/pymol_selected_examples"
    ).resolve()
    singles_dir = output_dir / "single_views"
    output_dir.mkdir(parents=True, exist_ok=True)
    singles_dir.mkdir(parents=True, exist_ok=True)

    contacts, ifp = read_ifp_contacts(results_dir)
    cmd = launch_pymol()
    report_rows = []
    single_images: dict[str, Path] = {}
    for case in CASES:
        path = singles_dir / f"{safe_name(case.condition_label)}_cluster{case.cluster_id}.png"
        report_rows.append(
            draw_case(
                cmd,
                case,
                contacts,
                path,
                surface=True,
                width=args.width,
                height=args.height,
                ray=args.ray,
                whole_molecule=args.whole_molecule,
            )
        )
        single_images[f"{case.condition_id}:{case.cluster_id}"] = path

    for condition_id in sorted({case.condition_id for case in CASES}):
        cases = [case for case in CASES if case.condition_id == condition_id]
        images = [single_images[f"{case.condition_id}:{case.cluster_id}"] for case in cases]
        labels = [f"{case.condition_label} cluster {case.cluster_id}" for case in cases]
        make_panel(images, labels, output_dir / f"{safe_name(condition_id)}_clusters_panel.png", columns=2)

    for protein_id in sorted({case.protein_id for case in CASES}):
        cases = [case for case in CASES if case.protein_id == protein_id]
        images = [single_images[f"{case.condition_id}:{case.cluster_id}"] for case in cases]
        labels = [f"{case.substrate_class} DP{case.dp} cluster {case.cluster_id}" for case in cases]
        make_panel(images, labels, output_dir / f"{protein_id}_all_selected_clusters_panel.png", columns=2)

    report = pd.DataFrame(report_rows)
    report["ifp_interaction_counts"] = [
        json.dumps(ifp_interaction_summary(ifp, pose_id), sort_keys=True)
        for pose_id in report["pose_id"]
    ]
    report.to_csv(output_dir / "pymol_vs_ifp_contact_check.tsv", sep="\t", index=False)
    for image_path in single_images.values():
        image_path.unlink(missing_ok=True)
    singles_dir.rmdir()
    print(f"Wrote PyMOL images and contact report to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
