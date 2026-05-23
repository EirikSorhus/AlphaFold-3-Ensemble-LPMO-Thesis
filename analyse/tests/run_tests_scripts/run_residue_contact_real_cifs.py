#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import traceback
from pathlib import Path
from typing import Any

from lpmo_pipeline.analysis.prolif_ifp import IFPResult, compute_ifp_single, write_pose_ifp_table
from lpmo_pipeline.analysis.residue_contact_extraction import (
    build_pose_residue_contact_rows,
    write_pose_residue_contact_table,
)
from lpmo_pipeline.io.analysis_export import export_analysis_artifacts
from lpmo_pipeline.io.normalize_mmcif import NormalizeMMCIFRunner


DEFAULT_CIF_PATHS: tuple[Path, ...] = (
    Path(
        "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/"
        "NAG4/af3/latest/Q7SCE9_NAG4/seed-4_sample-1/Q7SCE9_NAG4_seed-4_sample-1_model.cif"
    ),
    Path(
        "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/"
        "STA6/af3/latest/Q59930_STA6/seed-2_sample-0/Q59930_STA6_seed-2_sample-0_model.cif"
    ),
    Path(
        "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/"
        "STA4/af3/latest/A0A0S2GKZ1_STA4/seed-2_sample-2/A0A0S2GKZ1_STA4_seed-2_sample-2_model.cif"
    ),
)

_FEATURE_SEPARATOR = "|"
_TARGET_PREFIX_TO_SUBSTRATE = {
    "NAG": "chitin",
    "CEL": "cellulose",
    "BGC": "cellulose",
    "STA": "amylose",
    "GLC": "amylose",
}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run normalization, protonation, ProLIF, and residue-contact extraction on real AF3 CIF cases.",
    )
    parser.add_argument(
        "--run-dir",
        required=True,
        help="Output directory for per-case artifacts and residue-contact summaries.",
    )
    parser.add_argument(
        "--run-id",
        default="",
        help="Optional run identifier written into the summary output.",
    )
    parser.add_argument(
        "cifs",
        nargs="*",
        help="Optional explicit CIF paths. Defaults to the three standard real-CIF cases.",
    )
    return parser.parse_args()


def _parse_pose_id(pose_id: str) -> dict[str, Any]:
    match = re.match(
        r"^(?P<protein>[^_]+)_(?P<ligand>[^_]+)_seed-(?P<seed>\d+)_sample-(?P<sample>\d+)_model$",
        pose_id,
    )
    if not match:
        return {"protein_id": "", "ligand_id": "", "seed": None, "sample": None}
    return {
        "protein_id": match.group("protein"),
        "ligand_id": match.group("ligand"),
        "seed": int(match.group("seed")),
        "sample": int(match.group("sample")),
    }


def _condition_id(protein_id: str, ligand_id: str, construct_type: str = "domain_only") -> str:
    match = re.match(r"^(?P<prefix>[A-Za-z]+)(?P<dp>\d+)$", ligand_id)
    if match is None:
        return f"{protein_id}__{construct_type}__unknown_DP0"
    substrate = _TARGET_PREFIX_TO_SUBSTRATE.get(match.group("prefix").upper(), "unknown")
    dp = int(match.group("dp"))
    return f"{protein_id}__{construct_type}__{substrate}_DP{dp}"


def _count_table_rows(path: Path) -> int:
    if not path.exists():
        return 0
    lines = path.read_text().splitlines()
    return max(len(lines) - 1, 0)


def _summarize_feature_contract(result: IFPResult) -> dict[str, Any]:
    ligand_residue_labels: set[str] = set()
    protein_residue_labels: set[str] = set()
    observed_interaction_types: set[str] = set()
    invalid_feature_names: list[str] = []

    for feature_name in result.feature_names:
        parts = feature_name.split(_FEATURE_SEPARATOR)
        if len(parts) != 3:
            invalid_feature_names.append(feature_name)
            continue
        ligand_label, protein_label, interaction_name = parts
        ligand_residue_labels.add(ligand_label)
        protein_residue_labels.add(protein_label)
        observed_interaction_types.add(interaction_name)

    return {
        "feature_name_format": "ligand_residue|protein_residue|interaction",
        "invalid_feature_names": invalid_feature_names,
        "ligand_residue_labels": sorted(ligand_residue_labels),
        "protein_residue_labels": sorted(protein_residue_labels),
        "observed_interaction_types": sorted(observed_interaction_types),
        "has_non_unl_ligand_residues": bool(ligand_residue_labels)
        and all(not label.startswith("UNL") for label in ligand_residue_labels),
    }


def _rebuild_feature_name(row: dict[str, Any]) -> str:
    protein_residue = f"{row['residue_name']}{row['residue_number']}.{row['residue_chain']}"
    return (
        f"{row['ligand_residue_label']}"
        f"{_FEATURE_SEPARATOR}{protein_residue}"
        f"{_FEATURE_SEPARATOR}{row['interaction_type']}"
    )


def main() -> int:
    args = _parse_args()
    run_dir = Path(args.run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    cif_paths = [Path(path).resolve() for path in args.cifs] if args.cifs else list(DEFAULT_CIF_PATHS)
    run_id = args.run_id or run_dir.name

    summary: dict[str, Any] = {
        "run_dir": str(run_dir),
        "run_id": run_id,
        "requested_cifs": [str(path) for path in cif_paths],
        "n_requested": len(cif_paths),
        "cases": [],
    }

    for cif_path in cif_paths:
        if not cif_path.exists():
            summary["error"] = f"Missing CIF: {cif_path}"
            summary_path = run_dir / "residue_contact_real_cifs_summary.json"
            summary_path.write_text(json.dumps(summary, indent=2))
            print(json.dumps(summary, indent=2))
            return 1

    try:
        pose_metadata_by_id: dict[str, dict[str, str]] = {}
        results: list[IFPResult] = []
        all_rows: list[dict[str, Any]] = []

        for index, cif_path in enumerate(cif_paths, start=1):
            pose_id = cif_path.stem
            metadata = _parse_pose_id(pose_id)
            protein_id = str(metadata.get("protein_id") or "")
            ligand_id = str(metadata.get("ligand_id") or "")
            condition_id = _condition_id(protein_id, ligand_id)
            case_dir = run_dir / "cases" / f"{index:04d}_{pose_id}"
            case_dir.mkdir(parents=True, exist_ok=True)

            case_summary: dict[str, Any] = {
                "index": index,
                "pose_id": pose_id,
                "source_cif": str(cif_path),
                "case_dir": str(case_dir),
                "condition_id": condition_id,
                **metadata,
                "status": "preparing",
            }

            normalize_dir = case_dir / "normalize"
            ok_norm, normalized_path = NormalizeMMCIFRunner(cif_path, normalize_dir).run()
            case_summary["normalize_ok"] = bool(ok_norm)
            case_summary["normalize_report_path"] = str(normalize_dir / "normalize_report.json")
            if not ok_norm or normalized_path is None:
                case_summary["status"] = "normalize_failed"
                summary["cases"].append(case_summary)
                continue

            normalized_path = Path(normalized_path).resolve()
            case_summary["normalized_cif"] = str(normalized_path)

            analysis_export_dir = case_dir / "analysis_export"
            ok_export, report = export_analysis_artifacts(normalized_path, analysis_export_dir)
            case_summary["analysis_export_ok"] = bool(ok_export)
            case_summary["analysis_export_report_path"] = str(analysis_export_dir / "analysis_export_report.json")
            case_summary["analysis_export_blockers"] = (report or {}).get("blockers", []) if report else []
            case_summary["analysis_export_warnings"] = (report or {}).get("warnings", []) if report else []

            if not ok_export:
                case_summary["status"] = "analysis_export_failed"
                summary["cases"].append(case_summary)
                continue

            complex_pdb = Path(str((report or {}).get("complex_for_prolif_pdb") or analysis_export_dir / "complex_for_prolif.pdb")).resolve()
            ligand_pdb = Path(str((report or {}).get("ligand_pdb") or analysis_export_dir / "ligand_only_for_prolif.pdb")).resolve()
            case_summary["complex_for_prolif_pdb"] = str(complex_pdb)
            case_summary["ligand_pdb"] = str(ligand_pdb)

            result = compute_ifp_single(
                complex_pdb=complex_pdb,
                ligand_pdb=ligand_pdb,
                pose_id=pose_id,
            )
            results.append(result)
            pose_metadata_by_id[pose_id] = {
                "protein_id": protein_id,
                "condition_id": condition_id,
            }
            case_summary["ifp_status"] = result.status
            case_summary["ifp_error"] = result.error
            case_summary["n_ifp_features"] = len(result.feature_names)
            case_summary["n_total_contacts"] = result.n_total_contacts
            case_summary["feature_contract"] = _summarize_feature_contract(result)

            if result.status != "ok":
                case_summary["status"] = "ifp_failed"
                summary["cases"].append(case_summary)
                continue

            rows = build_pose_residue_contact_rows([result], pose_metadata_by_id)
            all_rows.extend(rows)

            rebuilt_feature_names = [_rebuild_feature_name(row) for row in rows]
            case_summary["residue_contact_row_count"] = len(rows)
            case_summary["residue_contact_present_sum"] = sum(int(row["contact_present"]) for row in rows)
            case_summary["residue_contact_feature_names_match"] = rebuilt_feature_names == result.feature_names
            case_summary["residue_contact_row_count_matches_features"] = len(rows) == len(result.feature_names)
            case_summary["residue_contact_present_sum_matches_ifp"] = (
                case_summary["residue_contact_present_sum"] == result.n_total_contacts
            )
            case_summary["status"] = "ok"
            summary["cases"].append(case_summary)

        pose_ifp_table_path = run_dir / "pose_ifp_table.tsv"
        write_pose_ifp_table(results, pose_ifp_table_path)

        residue_contact_table_path = run_dir / "pose_residue_contact_table.tsv"
        write_pose_residue_contact_table(all_rows, residue_contact_table_path)

        summary.update(
            {
                "pose_ifp_table_tsv": str(pose_ifp_table_path),
                "pose_residue_contact_table_tsv": str(residue_contact_table_path),
                "pose_ifp_rows": _count_table_rows(pose_ifp_table_path),
                "pose_residue_contact_rows": _count_table_rows(residue_contact_table_path),
                "n_ok": sum(1 for case in summary["cases"] if case.get("status") == "ok"),
                "n_errors": sum(1 for case in summary["cases"] if case.get("status") not in {"ok"}),
            }
        )

        summary_path = run_dir / "residue_contact_real_cifs_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))

        if summary["n_ok"] == 0:
            return 1
        if summary["n_ok"] != summary["n_requested"]:
            return 1
        if not residue_contact_table_path.exists() or summary["pose_residue_contact_rows"] == 0:
            return 1
        if any(
            not case.get("feature_contract", {}).get("has_non_unl_ligand_residues", False)
            or case.get("feature_contract", {}).get("invalid_feature_names")
            or not case.get("residue_contact_feature_names_match", False)
            or not case.get("residue_contact_row_count_matches_features", False)
            or not case.get("residue_contact_present_sum_matches_ifp", False)
            for case in summary["cases"]
            if case.get("status") == "ok"
        ):
            return 1
        return 0
    except Exception as exc:
        summary["error"] = f"{exc.__class__.__name__}: {exc}"
        summary["traceback"] = traceback.format_exc()
        summary_path = run_dir / "residue_contact_real_cifs_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
