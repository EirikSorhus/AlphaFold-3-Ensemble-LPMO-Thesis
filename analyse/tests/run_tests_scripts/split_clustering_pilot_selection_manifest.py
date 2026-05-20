#!/usr/bin/env python3
"""Split a clustering-pilot selection manifest into independent Slurm shards."""
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import replace
from pathlib import Path
from typing import Any

from lpmo_pipeline.analysis.clustering_pilot_real_case import (
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
        help="Number of selections per shard. Use 1 for maximum Slurm parallelism.",
    )
    parser.add_argument("--summary-path", default="", help="Optional JSON summary path.")
    parser.add_argument("--index-path", default="", help="Optional TSV index path.")
    return parser.parse_args()


def _chunks(values: tuple[Any, ...], size: int) -> list[tuple[Any, ...]]:
    return [values[index : index + size] for index in range(0, len(values), size)]


def main() -> int:
    args = _parse_args()
    input_manifest = Path(args.input_manifest).resolve()
    output_dir = Path(args.output_dir).resolve()
    shard_size = max(1, int(args.shard_size))

    manifest = load_selection_manifest(input_manifest)
    output_dir.mkdir(parents=True, exist_ok=True)

    index_path = Path(args.index_path).resolve() if args.index_path else output_dir / "shard_index.tsv"
    summary_path = Path(args.summary_path).resolve() if args.summary_path else output_dir / "shard_summary.json"

    rows: list[dict[str, Any]] = []
    for array_task_id, selections in enumerate(_chunks(manifest.selections, shard_size), start=1):
        shard_id = f"shard_{array_task_id:04d}"
        shard_manifest = replace(manifest, selections=selections)
        shard_manifest_path = output_dir / f"{shard_id}_selection.yaml"
        write_selection_manifest(shard_manifest, shard_manifest_path)
        protein_ids = sorted({selection.protein_id for selection in selections})
        targets = sorted({selection.target or "" for selection in selections})
        rows.append(
            {
                "array_task_id": array_task_id,
                "shard_id": shard_id,
                "construct_type": manifest.construct_type,
                "n_selections": len(selections),
                "manifest_path": str(shard_manifest_path),
                "protein_ids": ",".join(protein_ids),
                "targets": ",".join(targets),
            }
        )

    index_path.parent.mkdir(parents=True, exist_ok=True)
    with index_path.open("w", newline="") as handle:
        fieldnames = [
            "array_task_id",
            "shard_id",
            "construct_type",
            "n_selections",
            "manifest_path",
            "protein_ids",
            "targets",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "input_manifest": str(input_manifest),
        "output_dir": str(output_dir),
        "index_path": str(index_path),
        "construct_type": manifest.construct_type,
        "shard_size": shard_size,
        "n_input_selections": len(manifest.selections),
        "n_shards": len(rows),
        "source_manifest": serialize_selection_manifest(manifest),
        "shards": rows,
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
