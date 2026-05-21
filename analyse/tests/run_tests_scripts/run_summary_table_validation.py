#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import traceback
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from lpmo_pipeline.analysis.cluster_signatures import write_cluster_table_tsv
from lpmo_pipeline.analysis.condition_summary import (
    build_condition_table_rows,
    build_protein_summary_rows,
    read_tsv,
    write_condition_table,
    write_protein_summary_table,
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate cluster/condition/protein summary-table consistency from an "
            "analysis-core production output directory."
        )
    )
    parser.add_argument(
        "--production-output",
        required=True,
        help="Production output directory to validate.",
    )
    parser.add_argument(
        "--output-dir",
        default="",
        help=(
            "Optional directory for generated cluster/condition/protein summary tables and "
            "the validation summary JSON. Defaults to <production-output>/summary_validation."
        ),
    )
    return parser.parse_args()


def _as_float(value: Any) -> float | None:
    if value in {None, "", "None"}:
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(numeric):
        return None
    return numeric


def _as_int(value: Any) -> int:
    numeric = _as_float(value)
    return int(numeric) if numeric is not None else 0


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


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
        "Expected either a production output directory or a staged shard root with shard_*/production_output children: "
        f"{path}"
    )


def _load_cluster_table_rows(production_output: Path) -> tuple[list[dict[str, Any]], str]:
    cluster_table_path = production_output / "cluster_table.tsv"
    if cluster_table_path.exists():
        return read_tsv(cluster_table_path), "cluster_table.tsv"

    cluster_signatures_path = production_output / "cluster_signatures.json"
    if not cluster_signatures_path.exists():
        raise FileNotFoundError(
            "Neither cluster_table.tsv nor cluster_signatures.json exists under "
            f"{production_output}"
        )

    payload = json.loads(cluster_signatures_path.read_text())
    if isinstance(payload, dict):
        clusters = payload.get("clusters") or []
    elif isinstance(payload, list):
        clusters = payload
    else:
        raise TypeError("cluster_signatures.json must contain a list or {'clusters': [...]} payload")
    return list(clusters), "cluster_signatures.json"


def _validate_condition_rows(
    condition_rows: list[dict[str, Any]],
    cluster_rows: list[dict[str, Any]],
) -> list[str]:
    errors: list[str] = []
    clusters_by_condition: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in cluster_rows:
        condition_id = str(row.get("condition_id", ""))
        if condition_id:
            clusters_by_condition[condition_id].append(row)

    seen_condition_ids: set[str] = set()
    for row in condition_rows:
        condition_id = str(row.get("condition_id", ""))
        if not condition_id:
            errors.append("Encountered condition row with empty condition_id")
            continue
        if condition_id in seen_condition_ids:
            errors.append(f"Duplicate condition_id in condition table: {condition_id}")
            continue
        seen_condition_ids.add(condition_id)

        expected_cluster_rows = clusters_by_condition.get(condition_id, [])
        expected_n_cluster_rows = len(expected_cluster_rows)
        actual_n_cluster_rows = _as_int(row.get("n_cluster_rows"))
        if actual_n_cluster_rows != expected_n_cluster_rows:
            errors.append(
                f"Condition {condition_id} has n_cluster_rows={actual_n_cluster_rows}, "
                f"expected {expected_n_cluster_rows}"
            )

        expected_total_occupancy = sum(
            _as_float(cluster_row.get("occupancy")) or 0.0 for cluster_row in expected_cluster_rows
        )
        actual_total_occupancy = _as_float(row.get("cluster_total_occupancy"))
        if expected_cluster_rows:
            if actual_total_occupancy is None or not math.isclose(
                actual_total_occupancy,
                expected_total_occupancy,
                rel_tol=1e-9,
                abs_tol=1e-9,
            ):
                errors.append(
                    f"Condition {condition_id} has cluster_total_occupancy={actual_total_occupancy}, "
                    f"expected {expected_total_occupancy}"
                )
        elif actual_total_occupancy is not None:
            errors.append(
                f"Condition {condition_id} has cluster_total_occupancy={actual_total_occupancy} despite no cluster rows"
            )

    return errors


def _validate_protein_rows(
    protein_rows: list[dict[str, Any]],
    condition_rows: list[dict[str, Any]],
) -> list[str]:
    errors: list[str] = []
    n_conditions_by_protein = Counter(str(row.get("protein_id", "")) for row in condition_rows)
    n_cluster_conditions_by_protein = Counter(
        str(row.get("protein_id", "")) for row in condition_rows if _as_int(row.get("n_clusters")) > 0
    )

    seen_protein_ids: set[str] = set()
    for row in protein_rows:
        protein_id = str(row.get("protein_id", ""))
        if not protein_id:
            errors.append("Encountered protein summary row with empty protein_id")
            continue
        if protein_id in seen_protein_ids:
            errors.append(f"Duplicate protein_id in protein summary table: {protein_id}")
            continue
        seen_protein_ids.add(protein_id)

        expected_n_conditions = n_conditions_by_protein.get(protein_id, 0)
        actual_n_conditions = _as_int(row.get("n_conditions"))
        if actual_n_conditions != expected_n_conditions:
            errors.append(
                f"Protein {protein_id} has n_conditions={actual_n_conditions}, expected {expected_n_conditions}"
            )

        expected_n_cluster_conditions = n_cluster_conditions_by_protein.get(protein_id, 0)
        actual_n_cluster_conditions = _as_int(row.get("n_conditions_with_clusters"))
        if actual_n_cluster_conditions != expected_n_cluster_conditions:
            errors.append(
                f"Protein {protein_id} has n_conditions_with_clusters={actual_n_cluster_conditions}, "
                f"expected {expected_n_cluster_conditions}"
            )

    missing_proteins = sorted(set(n_conditions_by_protein) - seen_protein_ids)
    for protein_id in missing_proteins:
        if protein_id:
            errors.append(f"Protein {protein_id} is present in condition rows but missing from protein summary")
    return errors


def main() -> int:
    args = _parse_args()
    production_output = Path(args.production_output).resolve()
    output_dir = (
        Path(args.output_dir).resolve()
        if args.output_dir
        else (production_output / "summary_validation").resolve()
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    summary: dict[str, Any] = {
        "production_output": str(production_output),
        "output_dir": str(output_dir),
    }

    try:
        production_outputs = _resolve_production_outputs(production_output)
        qc_attrition_rows: list[dict[str, Any]] = []
        condition_cluster_summary_rows: list[dict[str, Any]] = []
        condition_convergence_summary_rows: list[dict[str, Any]] = []
        condition_patch_summary_rows: list[dict[str, Any]] = []
        pose_confidence_rows: list[dict[str, Any]] = []
        cluster_rows: list[dict[str, Any]] = []
        cluster_sources: Counter[str] = Counter()
        skipped_production_outputs: list[dict[str, str]] = []
        loaded_production_outputs: list[Path] = []
        for output_root in production_outputs:
            try:
                cluster_rows_for_output, cluster_source = _load_cluster_table_rows(output_root)
            except FileNotFoundError as exc:
                skipped_production_outputs.append({"path": str(output_root), "reason": str(exc)})
                continue

            qc_attrition_rows.extend(read_tsv(output_root / "qc_attrition_table.tsv"))
            condition_cluster_summary_rows.extend(read_tsv(output_root / "condition_cluster_summary.tsv"))
            condition_convergence_summary_rows.extend(read_tsv(output_root / "condition_convergence_summary.tsv"))
            condition_patch_summary_rows.extend(read_tsv(output_root / "condition_patch_summary.tsv"))
            pose_confidence_rows.extend(read_tsv(output_root / "pose_confidence.tsv"))
            cluster_rows.extend(cluster_rows_for_output)
            cluster_sources[cluster_source] += 1
            loaded_production_outputs.append(output_root)

        if not loaded_production_outputs:
            raise FileNotFoundError(
                "No complete production outputs with cluster summary surfaces were available under "
                f"{production_output}"
            )

        generated_cluster_table_path = output_dir / "generated_cluster_table.tsv"
        generated_condition_table_path = output_dir / "generated_condition_table.tsv"
        generated_protein_summary_path = output_dir / "generated_protein_summary_table.tsv"

        write_cluster_table_tsv(cluster_rows, generated_cluster_table_path)
        condition_rows = build_condition_table_rows(
            qc_attrition_rows=qc_attrition_rows,
            condition_cluster_summary_rows=condition_cluster_summary_rows,
            condition_convergence_summary_rows=condition_convergence_summary_rows,
            condition_patch_summary_rows=condition_patch_summary_rows,
            pose_confidence_rows=pose_confidence_rows,
            cluster_table_rows=cluster_rows,
        )
        write_condition_table(condition_rows, generated_condition_table_path)
        protein_rows = build_protein_summary_rows(condition_rows)
        write_protein_summary_table(protein_rows, generated_protein_summary_path)

        qc_condition_ids = {str(row.get("condition_id", "")) for row in qc_attrition_rows if row.get("condition_id")}
        condition_condition_ids = {str(row.get("condition_id", "")) for row in condition_rows if row.get("condition_id")}
        cluster_condition_ids = {str(row.get("condition_id", "")) for row in cluster_rows if row.get("condition_id")}
        cluster_only_condition_ids = sorted(cluster_condition_ids - qc_condition_ids)

        errors: list[str] = []
        if qc_condition_ids != condition_condition_ids:
            errors.append(
                "Condition-table condition_id set differs from qc_attrition_table.tsv master set"
            )
        if cluster_only_condition_ids:
            errors.append(
                f"Cluster rows reference condition_ids not present in qc_attrition_table.tsv: {cluster_only_condition_ids[:10]}"
            )
        errors.extend(_validate_condition_rows(condition_rows, cluster_rows))
        errors.extend(_validate_protein_rows(protein_rows, condition_rows))

        n_conditions_with_cluster_rows = sum(1 for row in condition_rows if _as_int(row.get("n_cluster_rows")) > 0)
        n_conditions_with_valid_clusters = sum(1 for row in condition_rows if _as_bool(row.get("any_valid_cluster")))
        n_conditions_without_valid_clusters = len(condition_rows) - n_conditions_with_valid_clusters
        n_conditions_zero_contact_eligible = sum(
            1 for row in condition_rows if _as_int(row.get("n_contact_eligible")) == 0
        )
        n_conditions_with_null_ifp_signal = sum(
            1 for row in condition_rows if (_as_float(row.get("null_ifp_fraction")) or 0.0) > 0.0
        )

        if n_conditions_with_cluster_rows == 0:
            errors.append("Validation target does not contain any retained-cluster conditions")
        if n_conditions_without_valid_clusters == 0:
            errors.append("Validation target does not contain any no-valid-cluster conditions")
        if n_conditions_zero_contact_eligible == 0 and n_conditions_with_null_ifp_signal == 0:
            errors.append(
                "Validation target does not contain any zero-contact-eligible or null-IFP-signal conditions"
            )

        summary.update(
            {
                "n_production_outputs": len(production_outputs),
                "production_outputs": [str(path) for path in loaded_production_outputs],
                "n_skipped_production_outputs": len(skipped_production_outputs),
                "skipped_production_outputs": skipped_production_outputs,
                "cluster_table_sources": dict(cluster_sources),
                "generated_cluster_table": str(generated_cluster_table_path),
                "generated_condition_table": str(generated_condition_table_path),
                "generated_protein_summary_table": str(generated_protein_summary_path),
                "n_qc_attrition_conditions": len(qc_attrition_rows),
                "n_cluster_rows": len(cluster_rows),
                "n_condition_rows": len(condition_rows),
                "n_protein_rows": len(protein_rows),
                "n_conditions_with_cluster_rows": n_conditions_with_cluster_rows,
                "n_conditions_with_valid_clusters": n_conditions_with_valid_clusters,
                "n_conditions_without_valid_clusters": n_conditions_without_valid_clusters,
                "n_conditions_zero_contact_eligible": n_conditions_zero_contact_eligible,
                "n_conditions_with_null_ifp_signal": n_conditions_with_null_ifp_signal,
                "validation_passed": not errors,
                "errors": errors,
            }
        )

        summary_path = output_dir / "summary_table_validation_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True))
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0 if not errors else 1
    except Exception as exc:
        summary.update(
            {
                "validation_passed": False,
                "error": f"{exc.__class__.__name__}: {exc}",
                "traceback": traceback.format_exc(),
            }
        )
        summary_path = output_dir / "summary_table_validation_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True))
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())