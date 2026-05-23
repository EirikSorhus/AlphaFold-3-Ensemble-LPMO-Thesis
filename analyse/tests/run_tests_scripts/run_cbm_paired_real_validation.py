#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import traceback
from pathlib import Path
from typing import Any

from lpmo_pipeline.analysis.cbm_comparison import run_cbm_paired_analysis
from lpmo_pipeline.analysis.cluster_signatures import write_cluster_table_tsv
from lpmo_pipeline.analysis.condition_summary import (
    build_condition_table_rows,
    read_tsv,
    write_condition_table,
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate domain-only/full-length summary tables from staged real outputs, "
            "merge them, and validate the condition-level CBM paired analysis on real rows."
        )
    )
    parser.add_argument(
        "--domain-only-root",
        default=(
            "/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/tests_results/"
            "clustering_pilot_staged/domain_only_shards"
        ),
        help="Staged domain-only shard root or one production output directory.",
    )
    parser.add_argument(
        "--full-length-root",
        default=(
            "/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/tests_results/"
            "clustering_pilot_staged/full_length_shards"
        ),
        help="Staged full-length shard root or one production output directory.",
    )
    parser.add_argument(
        "--protein-metadata",
        default=(
            "/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/input_data/"
            "metadata_final_ec_fixed.tsv"
        ),
        help="Protein metadata TSV used for family and CBM labels.",
    )
    parser.add_argument(
        "--output-dir",
        default="",
        help=(
            "Optional output directory. Defaults to "
            "<analyse>/tests/tests_results/cbm_paired_real_validation_<timestamp-less manual run>."
        ),
    )
    return parser.parse_args()


def _resolve_production_outputs(path: Path) -> list[Path]:
    if (path / "qc_attrition_table.tsv").exists() or (path / "cluster_signatures.json").exists():
        return [path]

    shard_outputs = sorted(
        child / "production_output"
        for child in path.iterdir()
        if child.is_dir() and child.name.startswith("shard_") and (child / "production_output").exists()
    )
    if shard_outputs:
        return shard_outputs

    raise FileNotFoundError(
        "Expected a production output directory or a staged shard root with shard_*/production_output children: "
        f"{path}"
    )


def _load_cluster_rows(production_output: Path) -> list[dict[str, Any]]:
    cluster_table_path = production_output / "cluster_table.tsv"
    if cluster_table_path.exists():
        return read_tsv(cluster_table_path)

    cluster_signatures_path = production_output / "cluster_signatures.json"
    if not cluster_signatures_path.exists():
        raise FileNotFoundError(
            f"Neither cluster_table.tsv nor cluster_signatures.json exists under {production_output}"
        )

    payload = json.loads(cluster_signatures_path.read_text())
    if isinstance(payload, dict):
        clusters = payload.get("clusters") or []
    elif isinstance(payload, list):
        clusters = payload
    else:
        raise TypeError("cluster_signatures.json must contain a list or {'clusters': [...]} payload")
    return list(clusters)


def _write_merged_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
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


def _build_summary_tables(production_root: Path, output_dir: Path) -> dict[str, Any]:
    qc_attrition_rows: list[dict[str, Any]] = []
    condition_cluster_summary_rows: list[dict[str, Any]] = []
    condition_convergence_summary_rows: list[dict[str, Any]] = []
    condition_patch_summary_rows: list[dict[str, Any]] = []
    pose_confidence_rows: list[dict[str, Any]] = []
    cluster_rows: list[dict[str, Any]] = []
    loaded_outputs: list[str] = []
    skipped_outputs: list[dict[str, str]] = []

    for production_output in _resolve_production_outputs(production_root):
        try:
            cluster_rows_for_output = _load_cluster_rows(production_output)
        except FileNotFoundError as exc:
            skipped_outputs.append({"path": str(production_output), "reason": str(exc)})
            continue

        qc_attrition_rows.extend(read_tsv(production_output / "qc_attrition_table.tsv"))
        condition_cluster_summary_rows.extend(read_tsv(production_output / "condition_cluster_summary.tsv"))
        condition_convergence_summary_rows.extend(read_tsv(production_output / "condition_convergence_summary.tsv"))
        condition_patch_summary_rows.extend(read_tsv(production_output / "condition_patch_summary.tsv"))
        pose_confidence_rows.extend(read_tsv(production_output / "pose_confidence.tsv"))
        cluster_rows.extend(cluster_rows_for_output)
        loaded_outputs.append(str(production_output))

    if not loaded_outputs:
        raise FileNotFoundError(f"No complete production outputs were available under {production_root}")

    cluster_table_path = output_dir / "generated_cluster_table.tsv"
    condition_table_path = output_dir / "generated_condition_table.tsv"
    write_cluster_table_tsv(cluster_rows, cluster_table_path)
    condition_rows = build_condition_table_rows(
        qc_attrition_rows=qc_attrition_rows,
        condition_cluster_summary_rows=condition_cluster_summary_rows,
        condition_convergence_summary_rows=condition_convergence_summary_rows,
        condition_patch_summary_rows=condition_patch_summary_rows,
        pose_confidence_rows=pose_confidence_rows,
        cluster_table_rows=cluster_rows,
    )
    write_condition_table(condition_rows, condition_table_path)

    return {
        "production_root": str(production_root),
        "output_dir": str(output_dir),
        "generated_condition_table": str(condition_table_path),
        "generated_cluster_table": str(cluster_table_path),
        "n_condition_rows": len(condition_rows),
        "n_cluster_rows": len(cluster_rows),
        "loaded_production_outputs": loaded_outputs,
        "skipped_production_outputs": skipped_outputs,
    }


def main() -> int:
    args = _parse_args()
    domain_root = Path(args.domain_only_root).resolve()
    full_root = Path(args.full_length_root).resolve()
    protein_metadata_path = Path(args.protein_metadata).resolve()
    project_root = Path(__file__).resolve().parents[2]
    output_dir = (
        Path(args.output_dir).resolve()
        if args.output_dir
        else (project_root / "tests/tests_results/cbm_paired_real_validation").resolve()
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    summary: dict[str, Any] = {
        "domain_only_root": str(domain_root),
        "full_length_root": str(full_root),
        "protein_metadata": str(protein_metadata_path),
        "output_dir": str(output_dir),
    }

    try:
        domain_summary = _build_summary_tables(domain_root, output_dir / "domain_only_summary_validation")
        full_summary = _build_summary_tables(full_root, output_dir / "full_length_summary_validation")

        domain_condition_rows = read_tsv(Path(domain_summary["generated_condition_table"]))
        full_condition_rows = read_tsv(Path(full_summary["generated_condition_table"]))
        domain_cluster_rows = read_tsv(Path(domain_summary["generated_cluster_table"]))
        full_cluster_rows = read_tsv(Path(full_summary["generated_cluster_table"]))

        overlap_keys = {
            (str(row.get("protein_id", "")), str(row.get("substrate_class", "")), str(row.get("dp", "")))
            for row in domain_condition_rows
            if row.get("protein_id")
        } & {
            (str(row.get("protein_id", "")), str(row.get("substrate_class", "")), str(row.get("dp", "")))
            for row in full_condition_rows
            if row.get("protein_id")
        }

        merged_condition_path = output_dir / "merged_inputs" / "merged_condition_table.tsv"
        merged_cluster_path = output_dir / "merged_inputs" / "merged_cluster_table.tsv"
        _write_merged_tsv(merged_condition_path, domain_condition_rows + full_condition_rows)
        _write_merged_tsv(merged_cluster_path, domain_cluster_rows + full_cluster_rows)

        cbm_result = run_cbm_paired_analysis(
            condition_table_path=merged_condition_path,
            cluster_table_path=merged_cluster_path,
            protein_metadata_path=protein_metadata_path,
            output_dir=output_dir,
        )

        paired_rows = read_tsv(cbm_result.table_paths["paired_comparison"])
        primary_rows = read_tsv(cbm_result.table_paths["primary_metric_summary"])
        errors: list[str] = []
        if not overlap_keys:
            errors.append("No overlapping protein_id x substrate_class x dp keys between domain-only and full-length tables")
        if len(paired_rows) != len(overlap_keys):
            errors.append(
                f"Paired-row count mismatch: cbm_paired_comparison_table.tsv has {len(paired_rows)} rows, expected {len(overlap_keys)} overlap keys"
            )
        if any(not row.get("domain_only_condition_id") for row in paired_rows):
            errors.append("Encountered paired CBM rows with empty domain_only_condition_id")
        if any(not row.get("full_length_condition_id") for row in paired_rows):
            errors.append("Encountered paired CBM rows with empty full_length_condition_id")
        if any(not row.get("family_label") for row in paired_rows):
            errors.append("Encountered paired CBM rows with empty family_label after metadata join")
        if any(not row.get("cbm_type") for row in paired_rows):
            errors.append("Encountered paired CBM rows with empty cbm_type after metadata join")
        if len({(row.get("protein_id"), row.get("substrate_class"), row.get("dp")) for row in paired_rows}) != len(paired_rows):
            errors.append("Duplicate protein_id x substrate_class x dp keys detected in paired CBM rows")

        primary_by_metric = {row.get("metric_name", ""): row for row in primary_rows}
        delta_qc_summary = primary_by_metric.get("delta_qc_pass_fraction")
        bridge_summary = primary_by_metric.get("bridge_fraction")
        if delta_qc_summary is None:
            errors.append("delta_qc_pass_fraction primary summary row is missing")
        elif int(float(delta_qc_summary.get("n_pairs") or 0)) != len(paired_rows):
            errors.append("delta_qc_pass_fraction primary summary row has unexpected n_pairs")
        elif delta_qc_summary.get("statistical_test") not in {"sign_test", "paired_wilcoxon_signed_rank"}:
            errors.append(
                "delta_qc_pass_fraction primary summary row did not report sign_test or paired_wilcoxon_signed_rank"
            )
        if bridge_summary is None:
            errors.append("bridge_fraction primary summary row is missing")
        elif bridge_summary.get("statistical_test") != "descriptive_full_length_only":
            errors.append("bridge_fraction primary summary row did not remain descriptive_full_length_only")

        summary.update(
            {
                "domain_only_summary": domain_summary,
                "full_length_summary": full_summary,
                "merged_condition_table": str(merged_condition_path),
                "merged_cluster_table": str(merged_cluster_path),
                "overlap_pair_key_count": len(overlap_keys),
                "overlap_pair_key_sample": [list(key) for key in sorted(overlap_keys)[:10]],
                "cbm_paired_analysis_summary": str(cbm_result.summary_path),
                "cbm_tables": {name: str(path) for name, path in sorted(cbm_result.table_paths.items())},
                "n_paired_rows": len(paired_rows),
                "n_nonempty_family_labels": sum(bool(row.get("family_label")) for row in paired_rows),
                "n_nonempty_cbm_types": sum(bool(row.get("cbm_type")) for row in paired_rows),
                "validation_passed": not errors,
                "errors": errors,
            }
        )
    except Exception as exc:
        summary.update(
            {
                "validation_passed": False,
                "error": f"{exc.__class__.__name__}: {exc}",
                "traceback": traceback.format_exc(),
            }
        )
        summary_path = output_dir / "cbm_paired_real_validation_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True))
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 1

    summary_path = output_dir / "cbm_paired_real_validation_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True))
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary.get("validation_passed") else 1


if __name__ == "__main__":
    raise SystemExit(main())