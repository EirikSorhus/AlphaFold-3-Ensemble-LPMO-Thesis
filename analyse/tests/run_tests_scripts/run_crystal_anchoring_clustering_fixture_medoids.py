#!/usr/bin/env python3
"""Run crystal anchoring for medoids in the clustering-stage fixture."""

from __future__ import annotations

import argparse
import csv
import json
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from lpmo_pipeline.analysis.crystal_anchoring import (
    DEFAULT_CRYSTAL_ROOT,
    DEFAULT_REFERENCE_INDEX_CSV,
    CRYSTAL_GEOMETRY_COLUMNS,
    load_crystal_reference_records,
    run_crystal_reference_screen,
    write_crystal_reference_screen_report,
)


DEFAULT_FIXTURE_ROOT = Path("tests/fixtures/clustering_stage_outputs")
DEFAULT_RESULTS_ROOT = Path("tests/tests_results")


@dataclass(frozen=True)
class FixtureMedoid:
    condition_id: str
    cluster_id: str
    medoid_pose_id: str
    medoid_structure_path: Path
    medoid_ifp_distance_sum: str
    protein_id: str
    construct_type: str
    target: str
    ligand_id: str


def _parse_condition_id(condition_id: str) -> tuple[str, str, str]:
    parts = condition_id.split("__")
    if len(parts) < 3:
        return condition_id, "", ""
    return parts[0], parts[1], parts[2]


def _parse_ligand_id(medoid_pose_id: str, target: str) -> str:
    prefix = medoid_pose_id.split("_seed-", 1)[0]
    tokens = prefix.split("_")
    if len(tokens) >= 2:
        return tokens[-1]
    target_prefix = target.split("_", 1)[0].lower()
    dp = "".join(ch for ch in target if ch.isdigit())
    if target_prefix.startswith("cellulose"):
        return f"CEL{dp}"
    if target_prefix.startswith("chitin"):
        return f"NAG{dp}"
    if target_prefix.startswith("amylose"):
        return f"STA{dp}"
    return target


def _split_filter(value: str) -> set[str]:
    return {item.strip() for item in value.split(",") if item.strip()}


def _load_fixture_medoids(fixture_root: Path) -> list[FixtureMedoid]:
    medoid_manifest = fixture_root / "medoid_manifest.tsv"
    if not medoid_manifest.exists():
        raise FileNotFoundError(f"Missing fixture medoid manifest: {medoid_manifest}")

    medoids: list[FixtureMedoid] = []
    with medoid_manifest.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            condition_id = row.get("condition_id", "")
            protein_id, construct_type, target = _parse_condition_id(condition_id)
            medoid_pose_id = row.get("medoid_pose_id", "")
            structure_path = row.get("medoid_structure_path", "")
            medoids.append(
                FixtureMedoid(
                    condition_id=condition_id,
                    cluster_id=row.get("cluster_id", ""),
                    medoid_pose_id=medoid_pose_id,
                    medoid_structure_path=Path(structure_path) if structure_path else Path(),
                    medoid_ifp_distance_sum=row.get("medoid_ifp_distance_sum", ""),
                    protein_id=protein_id,
                    construct_type=construct_type,
                    target=target,
                    ligand_id=_parse_ligand_id(medoid_pose_id, target),
                )
            )
    return medoids


def _has_crystal_references(
    protein_id: str,
    *,
    reference_index_csv: Path,
    crystal_root: Path,
    cache: dict[str, bool],
) -> bool:
    if protein_id in cache:
        return cache[protein_id]
    try:
        cache[protein_id] = bool(
            load_crystal_reference_records(
                protein_id,
                reference_index_csv=reference_index_csv,
                crystal_root=crystal_root,
            )
        )
    except Exception:
        cache[protein_id] = False
    return cache[protein_id]


def _comparison_rows(
    *,
    medoid: FixtureMedoid,
    report_path: Path,
    report: Any,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for comparison in report.comparisons:
        crystal_geometry = comparison.crystal_geometry or {}
        rows.append(
            {
                "condition_id": medoid.condition_id,
                "protein_id": medoid.protein_id,
                "construct_type": medoid.construct_type,
                "target": medoid.target,
                "ligand_id": medoid.ligand_id,
                "cluster_id": medoid.cluster_id,
                "representative_pose_id": medoid.medoid_pose_id,
                "representative_role": "cluster_medoid",
                "medoid_pose_id": medoid.medoid_pose_id,
                "medoid_structure_path": str(medoid.medoid_structure_path),
                "medoid_ifp_distance_sum": medoid.medoid_ifp_distance_sum,
                "report_path": str(report_path),
                "pdb_code": comparison.pdb_code,
                "source_cif": comparison.source_cif,
                "prepared_subset_cif": comparison.prepared_subset_cif,
                "selected_protein_chain": comparison.selected_protein_chain,
                "ligand_chain_ids": ",".join(comparison.ligand_chain_ids),
                "copper_chain_ids": ",".join(comparison.copper_chain_ids),
                "preferred_protein_chain": comparison.preferred_protein_chain,
                "preferred_chain_has_ligand": comparison.preferred_chain_has_ligand,
                "used_fallback_protein_chain": comparison.used_fallback_protein_chain,
                "representative_pose_ifp_status": comparison.representative_pose_ifp_status,
                "crystal_ifp_status": comparison.crystal_ifp_status,
                "crystal_ifp_contact_eligible": comparison.crystal_ifp_contact_eligible,
                "crystal_ifp_exclusion_class": comparison.crystal_ifp_exclusion_class,
                "crystal_n_vdw_interactions": comparison.crystal_n_vdw_interactions,
                "crystal_n_non_vdw_interactions": comparison.crystal_n_non_vdw_interactions,
                "crystal_n_non_vdw_contact_residues": comparison.crystal_n_non_vdw_contact_residues,
                "ifp_comparison_eligible": comparison.ifp_comparison_eligible,
                "ifp_tanimoto": comparison.ifp_tanimoto,
                "pocket_rmsd": comparison.pocket_rmsd,
                "pocket_rmsd_below_threshold": comparison.pocket_rmsd_below_threshold,
                "status": comparison.status,
                "error": comparison.error or "",
                "crystal_Cu_C1_distance": crystal_geometry.get("Cu_C1_distance", ""),
                "crystal_Cu_C4_distance": crystal_geometry.get("Cu_C4_distance", ""),
                "crystal_geometry_status_C1": crystal_geometry.get("geometry_status_C1", ""),
                "crystal_geometry_status_C4": crystal_geometry.get("geometry_status_C4", ""),
            }
        )
    return rows


def _geometry_rows(*, medoid: FixtureMedoid, report: Any) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for comparison in report.comparisons:
        if not comparison.crystal_geometry:
            continue
        row = {column: "" for column in CRYSTAL_GEOMETRY_COLUMNS}
        row.update(comparison.crystal_geometry)
        row.update(
            {
                "condition_id": medoid.condition_id,
                "cluster_id": medoid.cluster_id,
                "representative_pose_id": medoid.medoid_pose_id,
                "representative_role": "cluster_medoid",
                "medoid_pose_id": medoid.medoid_pose_id,
                "protein_id": medoid.protein_id,
                "pdb_code": comparison.pdb_code,
                "source_cif": comparison.source_cif,
                "prepared_subset_cif": comparison.prepared_subset_cif,
                "selected_protein_chain": comparison.selected_protein_chain,
                "ligand_chain_ids": ",".join(comparison.ligand_chain_ids),
                "copper_chain_ids": ",".join(comparison.copper_chain_ids),
            }
        )
        rows.append(row)
    return rows


def _write_tsv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def _diagnostic_summary_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    def summarize(scope: str, selected_rows: list[dict[str, Any]]) -> dict[str, Any]:
        total = len(selected_rows)
        eligible = sum(1 for row in selected_rows if str(row.get("crystal_ifp_contact_eligible")) == "True")
        vdw_only = sum(1 for row in selected_rows if row.get("crystal_ifp_exclusion_class") == "vdw_only")
        zero_or_null = sum(
            1
            for row in selected_rows
            if row.get("crystal_ifp_exclusion_class") in {"zero_contacts", "null_ifp"}
        )
        low_specific = sum(
            1 for row in selected_rows if row.get("crystal_ifp_exclusion_class") == "low_specific_contact"
        )
        return {
            "scope": scope,
            "n_total": total,
            "n_contact_eligible": eligible,
            "n_not_contact_eligible": total - eligible,
            "n_vdw_only": vdw_only,
            "n_zero_or_null_contacts": zero_or_null,
            "n_low_specific_contact": low_specific,
            "contact_eligible_fraction": eligible / total if total else 0.0,
            "not_contact_eligible_fraction": (total - eligible) / total if total else 0.0,
            "vdw_only_fraction": vdw_only / total if total else 0.0,
            "low_specific_contact_fraction": low_specific / total if total else 0.0,
        }

    unique_by_reference: dict[tuple[str, str, str], dict[str, Any]] = {}
    for row in rows:
        key = (
            str(row.get("protein_id", "")),
            str(row.get("pdb_code", "")),
            str(row.get("selected_protein_chain", "")),
        )
        unique_by_reference.setdefault(key, row)

    return [
        summarize("unique_crystal_reference", list(unique_by_reference.values())),
        summarize("medoid_vs_crystal_comparison", rows),
    ]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run crystal anchoring for medoids in tests/fixtures/clustering_stage_outputs. "
            "By default, all fixture medoids whose protein has a crystal reference are screened."
        ),
    )
    parser.add_argument("--fixture-root", default=str(DEFAULT_FIXTURE_ROOT))
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--run-id", default="")
    parser.add_argument("--reference-index-csv", default=str(DEFAULT_REFERENCE_INDEX_CSV))
    parser.add_argument("--crystal-root", default=str(DEFAULT_CRYSTAL_ROOT))
    parser.add_argument("--include-proteins", default="", help="Comma-separated protein IDs to include.")
    parser.add_argument("--include-conditions", default="", help="Comma-separated condition IDs to include.")
    parser.add_argument("--max-medoids", type=int, default=0, help="Optional cap after filtering; 0 means all.")
    parser.add_argument(
        "--allow-empty",
        action="store_true",
        help="Return success even if no medoid-vs-crystal comparisons were produced.",
    )
    parser.add_argument(
        "--fail-on-screen-error",
        action="store_true",
        help="Return failure if any medoid screen raises an exception.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    fixture_root = Path(args.fixture_root).resolve()
    run_dir = Path(args.run_dir).resolve()
    run_id = args.run_id or run_dir.name
    reference_index_csv = Path(args.reference_index_csv).resolve()
    crystal_root = Path(args.crystal_root).resolve()
    include_proteins = _split_filter(args.include_proteins)
    include_conditions = _split_filter(args.include_conditions)

    run_dir.mkdir(parents=True, exist_ok=True)
    output_root = run_dir / "crystal_anchoring"

    summary: dict[str, Any] = {
        "run_id": run_id,
        "run_dir": str(run_dir),
        "fixture_root": str(fixture_root),
        "reference_index_csv": str(reference_index_csv),
        "crystal_root": str(crystal_root),
        "include_proteins": sorted(include_proteins),
        "include_conditions": sorted(include_conditions),
        "max_medoids": args.max_medoids,
    }

    medoids = _load_fixture_medoids(fixture_root)
    reference_cache: dict[str, bool] = {}
    filtered: list[FixtureMedoid] = []
    skipped_missing_structure: list[dict[str, str]] = []
    skipped_no_reference: dict[str, int] = {}

    for medoid in medoids:
        if include_proteins and medoid.protein_id not in include_proteins:
            continue
        if include_conditions and medoid.condition_id not in include_conditions:
            continue
        if not medoid.medoid_structure_path.exists():
            skipped_missing_structure.append(
                {
                    "condition_id": medoid.condition_id,
                    "cluster_id": medoid.cluster_id,
                    "medoid_pose_id": medoid.medoid_pose_id,
                    "medoid_structure_path": str(medoid.medoid_structure_path),
                }
            )
            continue
        if not _has_crystal_references(
            medoid.protein_id,
            reference_index_csv=reference_index_csv,
            crystal_root=crystal_root,
            cache=reference_cache,
        ):
            skipped_no_reference[medoid.protein_id] = skipped_no_reference.get(medoid.protein_id, 0) + 1
            continue
        filtered.append(medoid)

    if args.max_medoids > 0:
        filtered = filtered[: args.max_medoids]

    summary.update(
        {
            "n_fixture_medoids": len(medoids),
            "n_selected_medoids": len(filtered),
            "skipped_missing_structure_count": len(skipped_missing_structure),
            "skipped_no_reference_by_protein": skipped_no_reference,
            "skipped_missing_structure_preview": skipped_missing_structure[:20],
        }
    )

    anchor_rows: list[dict[str, Any]] = []
    geometry_rows: list[dict[str, Any]] = []
    medoid_runs: list[dict[str, Any]] = []
    screen_errors: list[dict[str, str]] = []

    for index, medoid in enumerate(filtered, start=1):
        medoid_dir = output_root / medoid.condition_id / medoid.medoid_pose_id
        report_path = medoid_dir / "crystal_reference_screen.json"
        medoid_summary: dict[str, Any] = {
            "index": index,
            "condition_id": medoid.condition_id,
            "cluster_id": medoid.cluster_id,
            "protein_id": medoid.protein_id,
            "ligand_id": medoid.ligand_id,
            "medoid_pose_id": medoid.medoid_pose_id,
            "medoid_structure_path": str(medoid.medoid_structure_path),
            "report_path": str(report_path),
        }
        try:
            report = run_crystal_reference_screen(
                medoid.medoid_structure_path,
                protein_id=medoid.protein_id,
                ligand_id=medoid.ligand_id,
                representative_pose_id=medoid.medoid_pose_id,
                output_dir=medoid_dir,
                reference_index_csv=reference_index_csv,
                crystal_root=crystal_root,
            )
            write_crystal_reference_screen_report(report, report_path)
            rows = _comparison_rows(medoid=medoid, report_path=report_path, report=report)
            anchor_rows.extend(rows)
            geometry_rows.extend(_geometry_rows(medoid=medoid, report=report))
            status_counts: dict[str, int] = {}
            for row in rows:
                status = str(row.get("status", ""))
                status_counts[status] = status_counts.get(status, 0) + 1
            medoid_summary.update(
                {
                    "status": "ok",
                    "comparison_count": len(rows),
                    "status_counts": status_counts,
                    "pocket_rmsd_count": sum(1 for row in rows if row.get("pocket_rmsd") not in {"", None}),
                    "ifp_comparison_eligible_count": sum(
                        1 for row in rows if str(row.get("ifp_comparison_eligible")) == "True"
                    ),
                }
            )
        except Exception as exc:
            error = {
                "condition_id": medoid.condition_id,
                "cluster_id": medoid.cluster_id,
                "medoid_pose_id": medoid.medoid_pose_id,
                "error": f"{exc.__class__.__name__}: {exc}",
                "traceback": traceback.format_exc(),
            }
            screen_errors.append(error)
            medoid_summary.update({"status": "error", "error": error["error"]})
        medoid_runs.append(medoid_summary)

    anchor_columns = [
        "condition_id",
        "protein_id",
        "construct_type",
        "target",
        "ligand_id",
        "cluster_id",
        "representative_pose_id",
        "representative_role",
        "medoid_pose_id",
        "medoid_structure_path",
        "medoid_ifp_distance_sum",
        "report_path",
        "pdb_code",
        "source_cif",
        "prepared_subset_cif",
        "selected_protein_chain",
        "ligand_chain_ids",
        "copper_chain_ids",
        "preferred_protein_chain",
        "preferred_chain_has_ligand",
        "used_fallback_protein_chain",
        "representative_pose_ifp_status",
        "crystal_ifp_status",
        "crystal_ifp_contact_eligible",
        "crystal_ifp_exclusion_class",
        "crystal_n_vdw_interactions",
        "crystal_n_non_vdw_interactions",
        "crystal_n_non_vdw_contact_residues",
        "ifp_comparison_eligible",
        "ifp_tanimoto",
        "pocket_rmsd",
        "pocket_rmsd_below_threshold",
        "crystal_Cu_C1_distance",
        "crystal_Cu_C4_distance",
        "crystal_geometry_status_C1",
        "crystal_geometry_status_C4",
        "status",
        "error",
    ]
    geometry_columns = [
        "condition_id",
        "cluster_id",
        "representative_pose_id",
        "representative_role",
        "medoid_pose_id",
        *CRYSTAL_GEOMETRY_COLUMNS,
    ]
    diagnostic_columns = [
        "scope",
        "n_total",
        "n_contact_eligible",
        "n_not_contact_eligible",
        "n_vdw_only",
        "n_zero_or_null_contacts",
        "n_low_specific_contact",
        "contact_eligible_fraction",
        "not_contact_eligible_fraction",
        "vdw_only_fraction",
        "low_specific_contact_fraction",
    ]

    anchor_table = run_dir / "fixture_crystal_anchor_table.tsv"
    geometry_table = run_dir / "fixture_crystal_geometry_table.tsv"
    diagnostic_table = run_dir / "fixture_crystal_ifp_diagnostic_summary.tsv"
    medoid_run_table = run_dir / "fixture_crystal_medoid_run_table.tsv"
    summary_path = run_dir / "fixture_crystal_anchoring_summary.json"

    _write_tsv(anchor_table, anchor_rows, anchor_columns)
    _write_tsv(geometry_table, geometry_rows, geometry_columns)
    _write_tsv(diagnostic_table, _diagnostic_summary_rows(anchor_rows), diagnostic_columns)
    _write_tsv(
        medoid_run_table,
        medoid_runs,
        [
            "index",
            "condition_id",
            "cluster_id",
            "protein_id",
            "ligand_id",
            "medoid_pose_id",
            "medoid_structure_path",
            "report_path",
            "status",
            "comparison_count",
            "pocket_rmsd_count",
            "ifp_comparison_eligible_count",
            "status_counts",
            "error",
        ],
    )

    status_counts: dict[str, int] = {}
    for row in anchor_rows:
        status = str(row.get("status", ""))
        status_counts[status] = status_counts.get(status, 0) + 1

    summary.update(
        {
            "comparison_count": len(anchor_rows),
            "geometry_row_count": len(geometry_rows),
            "status_counts": status_counts,
            "screen_error_count": len(screen_errors),
            "screen_errors": screen_errors,
            "anchor_table": str(anchor_table),
            "geometry_table": str(geometry_table),
            "diagnostic_table": str(diagnostic_table),
            "medoid_run_table": str(medoid_run_table),
        }
    )
    summary_path.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    if not args.allow_empty and not anchor_rows:
        return 1
    if args.fail_on_screen_error and screen_errors:
        return 1
    if filtered and not any(row.get("pocket_rmsd") not in {"", None} for row in anchor_rows):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
