#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import traceback
from pathlib import Path
from typing import Any

from lpmo_pipeline.analysis.family_enrichment_postprocess import (
    read_tsv,
    run_family_enrichment_postprocess,
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge Stage 16b outputs from staged production shards and run the optional "
            "AA9/AA10 family enrichment postprocess."
        )
    )
    parser.add_argument(
        "--production-output-root",
        required=True,
        help="Production output directory or staged shard root to validate.",
    )
    parser.add_argument(
        "--protein-metadata",
        default="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/input_data/metadata_final_ec_fixed.tsv",
        help="Protein metadata TSV used to resolve family labels and catalytic sequence groups.",
    )
    parser.add_argument(
        "--core-fasta",
        default=(
            "/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/input_data/"
            "lpmo_core_domain_2026-03-14_06-52-15_deduplicated.fasta"
        ),
        help="Deduplicated catalytic-core FASTA used as the canonical sequence source.",
    )
    parser.add_argument(
        "--output-dir",
        default="",
        help=(
            "Optional output directory for merged Stage 16b inputs and the family enrichment "
            "artifacts. Defaults to <production-output-root>/family_enrichment_validation."
        ),
    )
    parser.add_argument(
        "--alignment-dir",
        default="",
        help="Optional directory with precomputed AA9.aligned.fasta / AA10.aligned.fasta files.",
    )
    parser.add_argument(
        "--mafft-executable",
        default="mafft",
        help="MAFFT executable used when precomputed family alignments are absent.",
    )
    return parser.parse_args()


def _resolve_production_outputs(path: Path) -> list[Path]:
    if (path / "protein_condition_residue_scores.tsv").exists():
        return [path]

    shard_outputs = sorted(
        child / "production_output"
        for child in path.iterdir()
        if child.is_dir() and child.name.startswith("shard_") and (child / "production_output").exists()
    )
    if shard_outputs:
        return shard_outputs

    raise FileNotFoundError(
        "Expected either a production output directory or a staged shard root with shard_*/production_output children: "
        f"{path}"
    )


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key in seen:
                continue
            seen.add(key)
            fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def main() -> int:
    args = _parse_args()
    production_output_root = Path(args.production_output_root).resolve()
    protein_metadata_path = Path(args.protein_metadata).resolve()
    core_fasta_path = Path(args.core_fasta).resolve()
    alignment_dir = Path(args.alignment_dir).resolve() if args.alignment_dir else None
    output_dir = (
        Path(args.output_dir).resolve()
        if args.output_dir
        else (production_output_root / "family_enrichment_validation").resolve()
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    summary: dict[str, Any] = {
        "production_output_root": str(production_output_root),
        "protein_metadata": str(protein_metadata_path),
        "core_fasta": str(core_fasta_path),
        "alignment_dir": str(alignment_dir) if alignment_dir is not None else None,
        "output_dir": str(output_dir),
    }

    try:
        production_outputs = _resolve_production_outputs(production_output_root)
        merged_scores_rows: list[dict[str, Any]] = []
        merged_delta_rows: list[dict[str, Any]] = []
        loaded_outputs: list[str] = []
        skipped_outputs: list[dict[str, str]] = []

        for production_output in production_outputs:
            scores_path = production_output / "protein_condition_residue_scores.tsv"
            delta_path = production_output / "protein_residue_regio_delta.tsv"
            if not scores_path.exists() or not delta_path.exists():
                skipped_outputs.append(
                    {
                        "path": str(production_output),
                        "reason": "missing_stage16b_tables",
                    }
                )
                continue
            merged_scores_rows.extend(read_tsv(scores_path))
            merged_delta_rows.extend(read_tsv(delta_path))
            loaded_outputs.append(str(production_output))

        if not loaded_outputs:
            raise FileNotFoundError(
                "No production outputs with both protein_condition_residue_scores.tsv and "
                "protein_residue_regio_delta.tsv were available under "
                f"{production_output_root}"
            )

        merged_input_dir = output_dir / "merged_stage16b"
        merged_scores_path = merged_input_dir / "protein_condition_residue_scores.tsv"
        merged_delta_path = merged_input_dir / "protein_residue_regio_delta.tsv"
        _write_tsv(merged_scores_path, merged_scores_rows)
        _write_tsv(merged_delta_path, merged_delta_rows)

        result = run_family_enrichment_postprocess(
            protein_condition_residue_scores_path=merged_scores_path,
            protein_residue_regio_delta_path=merged_delta_path,
            protein_metadata_path=protein_metadata_path,
            core_fasta_path=core_fasta_path,
            output_dir=output_dir,
            alignment_dir=alignment_dir,
            mafft_executable=args.mafft_executable,
        )

        validation_status = "ok" if result.processed_families else "skipped"
        summary.update(
            {
                "merged_scores_path": str(merged_scores_path),
                "merged_delta_path": str(merged_delta_path),
                "n_production_outputs": len(production_outputs),
                "n_loaded_production_outputs": len(loaded_outputs),
                "n_skipped_production_outputs": len(skipped_outputs),
                "loaded_production_outputs": loaded_outputs,
                "skipped_production_outputs": skipped_outputs,
                "n_merged_condition_rows": len(merged_scores_rows),
                "n_merged_delta_rows": len(merged_delta_rows),
                "family_enrichment_summary_path": str(result.summary_path),
                "family_substrate_residue_enrichment_path": str(
                    result.family_substrate_residue_enrichment_path
                ),
                "family_wrong_ligand_residue_enrichment_path": str(
                    result.family_wrong_ligand_residue_enrichment_path
                ),
                "processed_families": result.processed_families,
                "skipped_families": result.skipped_families,
                "validation_status": validation_status,
            }
        )
    except Exception as exc:
        summary.update(
            {
                "validation_status": "error",
                "error": str(exc),
                "traceback": traceback.format_exc(),
            }
        )
        summary_path = output_dir / "family_enrichment_validation_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return 1

    summary_path = output_dir / "family_enrichment_validation_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
