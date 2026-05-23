#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import traceback
from pathlib import Path
from typing import Any

from lpmo_pipeline.analysis.crystal_anchoring import (
    CrystalReferencePoseComparison,
    CrystalReferenceRecord,
    filter_crystal_reference_records,
    load_crystal_reference_records,
    run_crystal_reference_screen,
    write_crystal_reference_screen_report,
)


DEFAULT_ANALYSE_ROOT = Path("/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse")
DEFAULT_REFERENCE_INDEX_CSV = DEFAULT_ANALYSE_ROOT / "input_data" / "pdb_structure_data.csv"
DEFAULT_CRYSTAL_ROOT = DEFAULT_ANALYSE_ROOT / "crystal_structures"
DEFAULT_WORK_CORE_ROOT = DEFAULT_ANALYSE_ROOT.parent / "structure_pipeline" / "work_core"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Find holo crystal references with matching AF3 ligand conditions, prepare both sides "
            "for ProLIF, and run crystal-vs-AF3 IFP comparisons."
        ),
    )
    parser.add_argument("--run-dir", required=True, help="Output directory for artifacts and summaries.")
    parser.add_argument("--run-id", default="", help="Optional run identifier written into the summary.")
    parser.add_argument(
        "--reference-index-csv",
        default=str(DEFAULT_REFERENCE_INDEX_CSV),
        help="Crystal reference CSV used to resolve holo structures.",
    )
    parser.add_argument(
        "--crystal-root",
        default=str(DEFAULT_CRYSTAL_ROOT),
        help="Root directory containing crystal mmCIF subsets.",
    )
    parser.add_argument(
        "--work-core-root",
        default=str(DEFAULT_WORK_CORE_ROOT),
        help="AF3 work_core root used to select the best-ranked pose for each ligand condition.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Optional limit on the number of matched protein-ligand conditions to run.",
    )
    return parser.parse_args()


def _normalize_reference_protein_id(value: str) -> str:
    return "__".join(token.strip() for token in str(value).split() if token.strip())


def _iter_candidate_proteins(reference_index_csv: Path) -> list[str]:
    proteins: set[str] = set()
    with reference_index_csv.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=";")
        for row in reader:
            if not str(row.get("Carbohydrate_Ligands", "") or "").strip():
                continue
            if not str(row.get("DP", "") or "").strip():
                continue
            proteins.add(_normalize_reference_protein_id(row.get("Uniprot_ID", "")))
    return sorted(protein for protein in proteins if protein)


def _select_best_ranked_af3_pose(
    protein_id: str,
    ligand_id: str,
    *,
    work_core_root: Path,
) -> tuple[Path, str, float]:
    runs_root = work_core_root / ligand_id / "af3" / "runs"
    if not runs_root.is_dir():
        raise FileNotFoundError(f"Missing AF3 runs root for {protein_id} {ligand_id}: {runs_root}")

    run_dirs = sorted(
        (path for path in runs_root.iterdir() if path.is_dir() and path.name.isdigit()),
        key=lambda path: int(path.name),
        reverse=True,
    )
    for run_dir in run_dirs:
        case_dir = run_dir / f"{protein_id}_{ligand_id}"
        if not case_dir.is_dir():
            continue

        best_score: float | None = None
        best_cif: Path | None = None
        for confidence_path in sorted(case_dir.glob("seed-*_sample-*/*_summary_confidences.json")):
            confidence_data = json.loads(confidence_path.read_text())
            ranking_score = float(confidence_data.get("ranking_score", 0.0))
            cif_path = confidence_path.with_name(
                f"{confidence_path.name.removesuffix('_summary_confidences.json')}_model.cif"
            )
            if not cif_path.exists():
                continue
            if best_score is None or ranking_score > best_score:
                best_score = ranking_score
                best_cif = cif_path.resolve()

        if best_cif is not None and best_score is not None:
            return best_cif, run_dir.name, best_score

    raise FileNotFoundError(
        f"No ranked AF3 pose found for {protein_id} {ligand_id} under {runs_root}"
    )


def _artifact_path_if_exists(path: Path) -> str:
    return str(path) if path.exists() else ""


def _representative_artifacts(output_dir: Path) -> dict[str, str]:
    representative_dir = output_dir / "representative_pose"
    return {
        "representative_normalized_cif": _artifact_path_if_exists(
            representative_dir / "normalize" / "normalized.cif"
        ),
        "representative_complex_for_prolif_pdb": _artifact_path_if_exists(
            representative_dir / "analysis_export" / "complex_for_prolif.pdb"
        ),
        "representative_ligand_pdb": _artifact_path_if_exists(
            representative_dir / "analysis_export" / "ligand_only_for_prolif.pdb"
        ),
    }


def _comparison_artifacts(comparison: CrystalReferencePoseComparison) -> dict[str, str]:
    normalized_cif = Path(comparison.prepared_normalized_cif) if comparison.prepared_normalized_cif else None
    reference_dir = normalized_cif.parent.parent if normalized_cif is not None else None
    return {
        "crystal_complex_for_prolif_pdb": (
            _artifact_path_if_exists(reference_dir / "analysis_export" / "complex_for_prolif.pdb")
            if reference_dir is not None
            else ""
        ),
        "crystal_ligand_pdb": (
            _artifact_path_if_exists(reference_dir / "analysis_export" / "ligand_only_for_prolif.pdb")
            if reference_dir is not None
            else ""
        ),
    }


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _discover_condition_records(
    *,
    reference_index_csv: Path,
    crystal_root: Path,
) -> list[tuple[str, str, list[CrystalReferenceRecord]]]:
    discovered: list[tuple[str, str, list[CrystalReferenceRecord]]] = []
    for protein_id in _iter_candidate_proteins(reference_index_csv):
        records = load_crystal_reference_records(
            protein_id,
            reference_index_csv=reference_index_csv,
            crystal_root=crystal_root,
        )
        ligand_conditions = sorted(
            {
                f"{str(record.carbohydrate_ligands).strip().upper()}{record.dp}"
                for record in records
                if record.expected_ligand
                and str(record.carbohydrate_ligands).strip()
                and record.dp is not None
            }
        )
        for ligand_id in ligand_conditions:
            matched_records = filter_crystal_reference_records(records, ligand_id=ligand_id)
            if matched_records:
                discovered.append((protein_id, ligand_id, matched_records))
    discovered.sort(key=lambda item: (item[1], item[0]))
    return discovered


def main() -> int:
    args = _parse_args()
    run_dir = Path(args.run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    reference_index_csv = Path(args.reference_index_csv).resolve()
    crystal_root = Path(args.crystal_root).resolve()
    work_core_root = Path(args.work_core_root).resolve()
    run_id = args.run_id or run_dir.name

    summary: dict[str, Any] = {
        "run_dir": str(run_dir),
        "run_id": run_id,
        "reference_index_csv": str(reference_index_csv),
        "crystal_root": str(crystal_root),
        "work_core_root": str(work_core_root),
        "cases": [],
        "unmatched_conditions": [],
    }

    try:
        discovered = _discover_condition_records(
            reference_index_csv=reference_index_csv,
            crystal_root=crystal_root,
        )
        if args.limit > 0:
            discovered = discovered[: args.limit]
        summary["n_discovered_conditions"] = len(discovered)

        matched_rows: list[dict[str, Any]] = []
        for index, (protein_id, ligand_id, records) in enumerate(discovered, start=1):
            pdb_codes = [record.pdb_code for record in records]
            case_dir = run_dir / "cases" / f"{index:04d}_{protein_id}_{ligand_id}"
            case_dir.mkdir(parents=True, exist_ok=True)

            case_summary: dict[str, Any] = {
                "index": index,
                "protein_id": protein_id,
                "ligand_id": ligand_id,
                "pdb_codes": pdb_codes,
                "n_crystal_references": len(records),
                "case_dir": str(case_dir),
            }

            try:
                representative_cif, representative_run_id, representative_ranking_score = _select_best_ranked_af3_pose(
                    protein_id,
                    ligand_id,
                    work_core_root=work_core_root,
                )
            except FileNotFoundError as exc:
                case_summary.update(
                    {
                        "status": "af3_missing",
                        "error": str(exc),
                    }
                )
                summary["unmatched_conditions"].append(case_summary)
                continue

            output_dir = case_dir / "crystal_anchoring"
            report = run_crystal_reference_screen(
                representative_cif,
                protein_id=protein_id,
                ligand_id=ligand_id,
                representative_pose_id=representative_cif.stem,
                output_dir=output_dir,
                reference_index_csv=reference_index_csv,
                crystal_root=crystal_root,
                records=records,
            )

            report_path = output_dir / "crystal_reference_screen.json"
            write_crystal_reference_screen_report(report, report_path)

            representative_artifacts = _representative_artifacts(output_dir)
            comparison_status_counts: dict[str, int] = {}
            comparisons_payload: list[dict[str, Any]] = []

            for comparison in report.comparisons:
                comparison_status_counts[comparison.status] = comparison_status_counts.get(comparison.status, 0) + 1
                comparison_artifacts = _comparison_artifacts(comparison)
                comparison_payload = {
                    "pdb_code": comparison.pdb_code,
                    "source_cif": comparison.source_cif,
                    "prepared_subset_cif": comparison.prepared_subset_cif,
                    "prepared_normalized_cif": comparison.prepared_normalized_cif,
                    "comparison_status": comparison.status,
                    "comparison_error": comparison.error,
                    "ifp_comparison_eligible": comparison.ifp_comparison_eligible,
                    "ifp_tanimoto": comparison.ifp_tanimoto,
                    "pocket_rmsd": comparison.pocket_rmsd,
                    "crystal_ifp_contact_eligible": comparison.crystal_ifp_contact_eligible,
                    "crystal_ifp_exclusion_class": comparison.crystal_ifp_exclusion_class,
                    **comparison_artifacts,
                }
                comparisons_payload.append(comparison_payload)
                matched_rows.append(
                    {
                        "protein_id": protein_id,
                        "ligand_id": ligand_id,
                        "pdb_code": comparison.pdb_code,
                        "representative_cif": str(representative_cif),
                        "representative_run_id": representative_run_id,
                        "representative_ranking_score": representative_ranking_score,
                        "source_cif": comparison.source_cif,
                        "prepared_subset_cif": comparison.prepared_subset_cif,
                        "prepared_normalized_cif": comparison.prepared_normalized_cif,
                        "representative_normalized_cif": representative_artifacts["representative_normalized_cif"],
                        "representative_complex_for_prolif_pdb": representative_artifacts["representative_complex_for_prolif_pdb"],
                        "representative_ligand_pdb": representative_artifacts["representative_ligand_pdb"],
                        "crystal_complex_for_prolif_pdb": comparison_artifacts["crystal_complex_for_prolif_pdb"],
                        "crystal_ligand_pdb": comparison_artifacts["crystal_ligand_pdb"],
                        "comparison_status": comparison.status,
                        "ifp_comparison_eligible": comparison.ifp_comparison_eligible,
                        "ifp_tanimoto": comparison.ifp_tanimoto,
                        "pocket_rmsd": comparison.pocket_rmsd,
                        "crystal_ifp_contact_eligible": comparison.crystal_ifp_contact_eligible,
                        "crystal_ifp_exclusion_class": comparison.crystal_ifp_exclusion_class,
                        "report_path": str(report_path),
                    }
                )

            case_summary.update(
                {
                    "status": "completed",
                    "representative_cif": str(representative_cif),
                    "representative_run_id": representative_run_id,
                    "representative_ranking_score": representative_ranking_score,
                    "report_path": str(report_path),
                    **representative_artifacts,
                    "comparison_status_counts": comparison_status_counts,
                    "comparisons": comparisons_payload,
                }
            )
            summary["cases"].append(case_summary)

        matched_tsv = run_dir / "matched_pairs.tsv"
        _write_tsv(matched_tsv, matched_rows)

        summary.update(
            {
                "matched_pairs_tsv": str(matched_tsv),
                "n_completed_conditions": len(summary["cases"]),
                "n_unmatched_conditions": len(summary["unmatched_conditions"]),
                "n_total_crystal_references_compared": sum(
                    len(case.get("comparisons", [])) for case in summary["cases"]
                ),
                "n_ifp_comparison_eligible": sum(
                    1
                    for case in summary["cases"]
                    for comparison in case.get("comparisons", [])
                    if comparison.get("ifp_comparison_eligible")
                ),
                "n_non_ok_comparisons": sum(
                    1
                    for case in summary["cases"]
                    for comparison in case.get("comparisons", [])
                    if comparison.get("comparison_status") != "ok"
                ),
            }
        )

        summary_path = run_dir / "crystal_fil_test_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))

        if summary["n_completed_conditions"] == 0:
            return 1
        return 0
    except Exception as exc:
        summary["error"] = f"{exc.__class__.__name__}: {exc}"
        summary["traceback"] = traceback.format_exc()
        summary_path = run_dir / "crystal_fil_test_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
