#!/usr/bin/env python3
"""Split a clustering-pilot selection manifest into independent Slurm shards."""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from dataclasses import replace
from pathlib import Path
from typing import Any

from lpmo_pipeline.analysis.clustering_pilot_real_case import (
    PilotSelection,
    PilotSelectionManifest,
    discover_selected_pose_inputs,
    load_selection_manifest,
    serialize_selection_manifest,
    write_selection_manifest,
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Split one real-case clustering-pilot selection manifest into smaller "
            "selection manifests. Shards preserve the original construct/config and "
            "only partition the protein-target selections."
        ),
    )
    parser.add_argument("--input-manifest", required=True, help="Source selection manifest YAML.")
    parser.add_argument("--output-dir", required=True, help="Directory for shard manifests.")
    parser.add_argument(
        "--shard-size",
        type=int,
        default=1,
        help="Number of selections per shard. Kept as a compatibility alias for --proteins-per-shard.",
    )
    parser.add_argument(
        "--proteins-per-shard",
        type=int,
        default=0,
        help="Maximum number of protein selections per shard. Defaults to --shard-size when omitted.",
    )
    parser.add_argument(
        "--balance-by",
        choices=("protein_count", "estimated_pose_count"),
        default="protein_count",
        help="Shard packing strategy. estimated_pose_count performs one manifest discovery pass and balances shards by expected pose load.",
    )
    parser.add_argument(
        "--output-root-template",
        default="",
        help=(
            "Optional output-root template recorded in the shard index, for example "
            "'/path/to/run/{construct_type}_shards/{shard_id}/production_output'."
        ),
    )
    parser.add_argument("--summary-path", default="", help="Optional JSON summary path.")
    parser.add_argument("--index-path", default="", help="Optional TSV index path.")
    return parser.parse_args()


def _selection_sort_key(selection: PilotSelection) -> tuple[str, str]:
    return (selection.protein_id, selection.target or "")


def _selection_weight_map(
    manifest: PilotSelectionManifest,
    *,
    balance_by: str,
) -> dict[tuple[str, str | None], int]:
    if balance_by == "protein_count":
        return {
            (selection.protein_id, selection.target): 1
            for selection in manifest.selections
        }

    discovery_payload = discover_selected_pose_inputs(manifest)
    pose_inputs = discovery_payload.get("pose_inputs") or []
    protein_counts: Counter[str] = Counter()
    pair_counts: Counter[tuple[str, str]] = Counter()
    for pose in pose_inputs:
        protein_id = str(pose.get("protein_id") or "")
        ligand_id = str(pose.get("ligand_id") or "")
        if not protein_id:
            continue
        protein_counts[protein_id] += 1
        if ligand_id:
            pair_counts[(protein_id, ligand_id)] += 1

    weights: dict[tuple[str, str | None], int] = {}
    for selection in manifest.selections:
        if selection.target:
            weight = pair_counts.get((selection.protein_id, selection.target), 0)
        else:
            weight = protein_counts.get(selection.protein_id, 0)
        weights[(selection.protein_id, selection.target)] = max(1, int(weight))
    return weights


def _build_shards(
    manifest: PilotSelectionManifest,
    *,
    proteins_per_shard: int,
    balance_by: str,
) -> list[dict[str, Any]]:
    proteins_per_shard = max(1, int(proteins_per_shard))
    selection_weights = _selection_weight_map(manifest, balance_by=balance_by)
    indexed_selections = [
        {
            "original_index": index,
            "selection": selection,
            "estimated_pose_count": int(selection_weights[(selection.protein_id, selection.target)]),
        }
        for index, selection in enumerate(manifest.selections)
    ]
    indexed_selections.sort(
        key=lambda entry: (
            -int(entry["estimated_pose_count"]),
            *_selection_sort_key(entry["selection"]),
        )
    )

    if not indexed_selections:
        return []

    n_shards = max(1, int(math.ceil(len(indexed_selections) / float(proteins_per_shard))))
    shards: list[dict[str, Any]] = [
        {"entries": [], "estimated_pose_count": 0}
        for _ in range(n_shards)
    ]
    for entry in indexed_selections:
        candidate = None
        for shard in shards:
            if len(shard["entries"]) >= proteins_per_shard:
                continue
            if candidate is None or int(shard["estimated_pose_count"]) < int(candidate["estimated_pose_count"]):
                candidate = shard
        if candidate is None:
            raise RuntimeError("No shard with free capacity while assigning selections")
        candidate["entries"].append(entry)
        candidate["estimated_pose_count"] = int(candidate["estimated_pose_count"]) + int(entry["estimated_pose_count"])

    normalized_shards: list[dict[str, Any]] = []
    for shard in shards:
        entries = sorted(shard["entries"], key=lambda entry: int(entry["original_index"]))
        normalized_shards.append(
            {
                "selections": tuple(entry["selection"] for entry in entries),
                "estimated_pose_count": int(shard["estimated_pose_count"]),
            }
        )
    return normalized_shards


def main() -> int:
    args = _parse_args()
    input_manifest = Path(args.input_manifest).resolve()
    output_dir = Path(args.output_dir).resolve()
    proteins_per_shard = max(1, int(args.proteins_per_shard or args.shard_size))

    manifest = load_selection_manifest(input_manifest)
    output_dir.mkdir(parents=True, exist_ok=True)

    index_path = Path(args.index_path).resolve() if args.index_path else output_dir / "shard_index.tsv"
    summary_path = Path(args.summary_path).resolve() if args.summary_path else output_dir / "shard_summary.json"

    rows: list[dict[str, Any]] = []
    shard_specs = _build_shards(
        manifest,
        proteins_per_shard=proteins_per_shard,
        balance_by=args.balance_by,
    )
    for array_task_id, shard_spec in enumerate(shard_specs, start=1):
        shard_id = f"shard_{array_task_id:04d}"
        selections = shard_spec["selections"]
        shard_manifest = replace(manifest, selections=selections)
        shard_manifest_path = output_dir / f"{shard_id}_selection.yaml"
        write_selection_manifest(shard_manifest, shard_manifest_path)
        protein_ids = sorted({selection.protein_id for selection in selections})
        targets = sorted({selection.target or "" for selection in selections})
        output_root = ""
        if args.output_root_template:
            output_root = args.output_root_template.format(
                construct_type=manifest.construct_type,
                shard_id=shard_id,
            )
        rows.append(
            {
                "array_task_id": array_task_id,
                "shard_id": shard_id,
                "construct_type": manifest.construct_type,
                "n_selections": len(selections),
                "estimated_pose_count": int(shard_spec["estimated_pose_count"]),
                "manifest_path": str(shard_manifest_path),
                "protein_ids": ",".join(protein_ids),
                "targets": ",".join(targets),
                "output_root": output_root,
            }
        )

    index_path.parent.mkdir(parents=True, exist_ok=True)
    with index_path.open("w", newline="") as handle:
        fieldnames = [
            "array_task_id",
            "shard_id",
            "construct_type",
            "n_selections",
            "estimated_pose_count",
            "manifest_path",
            "protein_ids",
            "targets",
            "output_root",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "input_manifest": str(input_manifest),
        "output_dir": str(output_dir),
        "index_path": str(index_path),
        "construct_type": manifest.construct_type,
        "proteins_per_shard": proteins_per_shard,
        "balance_by": args.balance_by,
        "n_input_selections": len(manifest.selections),
        "n_shards": len(rows),
        "total_estimated_pose_count": sum(int(row["estimated_pose_count"]) for row in rows),
        "source_manifest": serialize_selection_manifest(manifest),
        "shards": rows,
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
