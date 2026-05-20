#!/usr/bin/env python3

import argparse
import json
import traceback
from pathlib import Path
from typing import Any

from lpmo_pipeline.analysis.crystal_anchoring import (
    run_crystal_reference_screen,
    write_crystal_reference_screen_report,
)


DEFAULT_PROTEIN_ID = "A0A0S2GKZ1"
DEFAULT_LIGAND_ID = "CEL4"
DEFAULT_WORK_CORE_ROOT = Path(
    "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core"
)


def _select_best_ranked_af3_pose(
    protein_id: str,
    ligand_id: str,
    *,
    work_core_root: Path = DEFAULT_WORK_CORE_ROOT,
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


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run crystal reference preparation and comparison on one real AF3 representative pose.",
    )
    parser.add_argument("--run-dir", required=True, help="Output directory for crystal anchoring artifacts.")
    parser.add_argument("--run-id", default="", help="Optional run identifier written to the summary.")
    parser.add_argument(
        "--representative-cif",
        default="",
        help="Optional explicit representative AF3 CIF. If omitted, the best-ranked AF3 pose is selected automatically.",
    )
    parser.add_argument(
        "--protein-id",
        default=DEFAULT_PROTEIN_ID,
        help="Protein identifier to resolve under crystal_structures/ and the reference CSV.",
    )
    parser.add_argument(
        "--ligand-id",
        default=DEFAULT_LIGAND_ID,
        help="Ligand identifier written into the report metadata.",
    )
    parser.add_argument(
        "--work-core-root",
        default=str(DEFAULT_WORK_CORE_ROOT),
        help="AF3 work_core root used when auto-selecting the best-ranked representative pose.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    run_dir = Path(args.run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    run_id = args.run_id or run_dir.name
    work_core_root = Path(args.work_core_root).resolve()

    representative_selection_mode = "explicit" if args.representative_cif else "auto_best_ranked"
    representative_run_id = ""
    representative_ranking_score: float | None = None
    if args.representative_cif:
        representative_cif = Path(args.representative_cif).resolve()
    else:
        representative_cif, representative_run_id, representative_ranking_score = _select_best_ranked_af3_pose(
            args.protein_id,
            args.ligand_id,
            work_core_root=work_core_root,
        )

    summary: dict[str, Any] = {
        "run_dir": str(run_dir),
        "run_id": run_id,
        "protein_id": args.protein_id,
        "ligand_id": args.ligand_id,
        "representative_cif": str(representative_cif),
        "representative_selection_mode": representative_selection_mode,
        "representative_selection_run_id": representative_run_id,
        "representative_selection_ranking_score": representative_ranking_score,
    }

    if not representative_cif.exists():
        summary["error"] = f"Missing representative CIF: {representative_cif}"
        print(json.dumps(summary, indent=2))
        return 1

    try:
        output_dir = run_dir / "crystal_anchoring_output"
        report = run_crystal_reference_screen(
            representative_cif,
            protein_id=args.protein_id,
            ligand_id=args.ligand_id,
            output_dir=output_dir,
        )

        report_path = output_dir / "crystal_reference_screen.json"
        write_crystal_reference_screen_report(report, report_path)

        comparison_status_counts: dict[str, int] = {}
        fallback_count = 0
        ok_count = 0
        ifp_scored_count = 0
        crystal_ifp_contact_eligible_count = 0
        pocket_rmsd_count = 0
        prepared_subset_paths: list[str] = []
        pdb_codes: list[str] = []
        comparison_ifp_artifacts: list[dict[str, str]] = []
        for comparison in report.comparisons:
            comparison_status_counts[comparison.status] = comparison_status_counts.get(comparison.status, 0) + 1
            if comparison.used_fallback_protein_chain:
                fallback_count += 1
            if comparison.status == "ok":
                ok_count += 1
            if comparison.ifp_tanimoto is not None:
                ifp_scored_count += 1
            if comparison.crystal_ifp_contact_eligible:
                crystal_ifp_contact_eligible_count += 1
            if comparison.pocket_rmsd is not None:
                pocket_rmsd_count += 1
            prepared_subset_paths.append(comparison.prepared_subset_cif)
            pdb_codes.append(comparison.pdb_code)
            comparison_ifp_artifacts.append(
                {
                    "pdb_code": comparison.pdb_code,
                    "crystal_ifp_result_json": comparison.crystal_ifp_result_json,
                    "crystal_pose_ifp_table_tsv": comparison.crystal_pose_ifp_table_tsv,
                    "crystal_ifp_matrix_csv": comparison.crystal_ifp_matrix_csv,
                }
            )

        summary.update(
            {
                "report_path": str(report_path),
                "comparison_count": len(report.comparisons),
                "comparison_status_counts": comparison_status_counts,
                "ok_count": ok_count,
                "ifp_scored_count": ifp_scored_count,
                "crystal_ifp_contact_eligible_count": crystal_ifp_contact_eligible_count,
                "pocket_rmsd_count": pocket_rmsd_count,
                "fallback_count": fallback_count,
                "pdb_codes": pdb_codes,
                "prepared_subset_paths": prepared_subset_paths,
                "representative_pose_ifp_result_json": report.representative_pose_ifp_result_json,
                "representative_pose_pose_ifp_table_tsv": report.representative_pose_pose_ifp_table_tsv,
                "representative_pose_ifp_matrix_csv": report.representative_pose_ifp_matrix_csv,
                "comparison_ifp_artifacts": comparison_ifp_artifacts,
            }
        )

        summary_path = run_dir / "crystal_anchoring_real_cifs_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))

        expected_codes = {"5ACI", "7PXW"} if args.protein_id == "A0A0S2GKZ1" else set(pdb_codes)
        if len(report.comparisons) == 0:
            return 1
        if expected_codes and set(pdb_codes) != expected_codes:
            return 1
        if any(comparison.status in {"pose_ifp_failed", "protonation_failed", "ifp_failed"} for comparison in report.comparisons):
            return 1
        if pocket_rmsd_count == 0:
            return 1
        if any(not Path(path).exists() for path in prepared_subset_paths):
            return 1
        return 0
    except Exception as exc:
        summary["error"] = f"{exc.__class__.__name__}: {exc}"
        summary["traceback"] = traceback.format_exc()
        summary_path = run_dir / "crystal_anchoring_real_cifs_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
