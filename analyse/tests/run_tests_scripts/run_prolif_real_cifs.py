#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import traceback
from pathlib import Path
from typing import Any

from lpmo_pipeline.analysis.prolif_ifp import (
    IFPBatch,
    IFPResult,
    compute_ifp_single,
    load_prolif_features_config,
    write_ifp_matrix,
    write_pose_ifp_table,
)
from lpmo_pipeline.io.normalize_mmcif import NormalizeMMCIFRunner
from lpmo_pipeline.io.protonate_export import protonate_and_export


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


def _expected_interaction_types() -> list[str]:
    config = load_prolif_features_config()
    configured = [
        str(name)
        for name in config.get("active_interaction_types", ())
        if str(name).strip()
    ]
    return configured


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run normalization, protonation, and ProLIF on three real AF3 CIF cases.",
    )
    parser.add_argument(
        "--run-dir",
        required=True,
        help="Output directory for per-case artifacts and ProLIF summaries.",
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


def _align_results(results: list[IFPResult]) -> IFPBatch:
    feature_names = sorted(
        {
            feature_name
            for result in results
            for feature_name in result.feature_names
        }
    )
    matrix: list[list[int]] = []
    for result in results:
        feature_map = {
            feature_name: int(value)
            for feature_name, value in zip(result.feature_names, result.flat_bitvector)
        }
        matrix.append([int(feature_map.get(feature_name, 0)) for feature_name in feature_names])
    return IFPBatch(results=results, matrix=matrix, feature_names=feature_names)


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

    has_non_unl_ligand_residues = bool(ligand_residue_labels) and all(
        not label.startswith("UNL") for label in ligand_residue_labels
    )

    return {
        "feature_name_format": "ligand_residue|protein_residue|interaction",
        "invalid_feature_names": invalid_feature_names,
        "ligand_residue_labels": sorted(ligand_residue_labels),
        "protein_residue_labels": sorted(protein_residue_labels),
        "observed_interaction_types": sorted(observed_interaction_types),
        "has_non_unl_ligand_residues": has_non_unl_ligand_residues,
    }


def main() -> int:
    args = _parse_args()
    run_dir = Path(args.run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    cif_paths = [Path(path).resolve() for path in args.cifs] if args.cifs else list(DEFAULT_CIF_PATHS)
    run_id = args.run_id or run_dir.name
    expected_interaction_types = _expected_interaction_types()

    summary: dict[str, Any] = {
        "run_dir": str(run_dir),
        "run_id": run_id,
        "expected_interaction_types": expected_interaction_types,
        "requested_cifs": [str(path) for path in cif_paths],
        "n_requested": len(cif_paths),
        "cases": [],
    }

    for cif_path in cif_paths:
        if not cif_path.exists():
            summary["error"] = f"Missing CIF: {cif_path}"
            print(json.dumps(summary, indent=2))
            return 1

    try:
        results: list[IFPResult] = []
        for index, cif_path in enumerate(cif_paths, start=1):
            pose_id = cif_path.stem
            metadata = _parse_pose_id(pose_id)
            case_dir = run_dir / "cases" / f"{index:04d}_{pose_id}"
            case_dir.mkdir(parents=True, exist_ok=True)

            case_summary: dict[str, Any] = {
                "index": index,
                "pose_id": pose_id,
                "source_cif": str(cif_path),
                "case_dir": str(case_dir),
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

            protonation_dir = case_dir / "protonated"
            ok_prot, report = protonate_and_export(normalized_path, protonation_dir)
            case_summary["protonate_ok"] = bool(ok_prot)
            case_summary["protonation_report_path"] = str(protonation_dir / "protonation_report.json")
            case_summary["protonation_blockers"] = (report or {}).get("blockers", []) if report else []
            case_summary["protonation_warnings"] = (report or {}).get("warnings", []) if report else []
            case_summary["complex_h_pdb"] = str(protonation_dir / "complex_H.pdb")
            case_summary["ligand_mol2"] = str(protonation_dir / "ligand_for_prolif.mol2")
            if not ok_prot:
                case_summary["status"] = "protonation_failed"
                summary["cases"].append(case_summary)
                continue

            result = compute_ifp_single(
                complex_pdb=protonation_dir / "complex_H.pdb",
                ligand_mol2=protonation_dir / "ligand_for_prolif.mol2",
                pose_id=pose_id,
            )
            results.append(result)
            feature_contract = _summarize_feature_contract(result)
            missing_interaction_types = [
                interaction_name
                for interaction_name in expected_interaction_types
                if interaction_name not in result.interaction_counts
            ]
            case_summary.update(
                {
                    "status": result.status,
                    "n_total_contacts": result.n_total_contacts,
                    "n_features": len(result.feature_names),
                    "interaction_counts": result.interaction_counts,
                    "missing_interaction_types": missing_interaction_types,
                    **feature_contract,
                    "ifp_error": result.error,
                }
            )
            summary["cases"].append(case_summary)

        prolif_output_dir = run_dir / "prolif_output"
        prolif_output_dir.mkdir(parents=True, exist_ok=True)
        pose_ifp_table_path = prolif_output_dir / "pose_ifp_table.tsv"
        ifp_matrix_path = prolif_output_dir / "ifp_matrix.csv"
        write_pose_ifp_table(results, pose_ifp_table_path)
        write_ifp_matrix(_align_results(results), ifp_matrix_path)

        status_counts: dict[str, int] = {}
        for result in results:
            status_counts[result.status] = status_counts.get(result.status, 0) + 1

        summary.update(
            {
                "pose_ifp_table": str(pose_ifp_table_path),
                "ifp_matrix": str(ifp_matrix_path),
                "pose_ifp_rows": _count_table_rows(pose_ifp_table_path),
                "ifp_matrix_rows": _count_table_rows(ifp_matrix_path),
                "n_normalize_ok": sum(1 for case in summary["cases"] if case.get("normalize_ok")),
                "n_protonate_ok": sum(1 for case in summary["cases"] if case.get("protonate_ok")),
                "n_ifp_results": len(results),
                "n_with_contacts": sum(1 for result in results if result.n_total_contacts > 0),
                "n_with_ligand_resolved_features": sum(
                    1
                    for case in summary["cases"]
                    if case.get("status") == "ok" and case.get("has_non_unl_ligand_residues")
                ),
                "status_counts": status_counts,
            }
        )

        summary_path = run_dir / "prolif_real_cifs_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))

        if summary["n_normalize_ok"] != len(cif_paths):
            return 1
        if summary["n_protonate_ok"] != len(cif_paths):
            return 1
        if len(results) != len(cif_paths):
            return 1
        if any(result.status not in {"ok", "zero_contacts"} for result in results):
            return 1
        if any(case.get("missing_interaction_types") for case in summary["cases"] if case.get("status") in {"ok", "zero_contacts"}):
            return 1
        if any(case.get("invalid_feature_names") for case in summary["cases"] if case.get("status") == "ok"):
            return 1
        if any(
            case.get("n_total_contacts", 0) > 0 and not case.get("has_non_unl_ligand_residues")
            for case in summary["cases"]
        ):
            return 1
        if summary["pose_ifp_rows"] != len(cif_paths):
            return 1
        if summary["ifp_matrix_rows"] != len(cif_paths):
            return 1
        if summary["n_with_contacts"] < 1:
            return 1
        return 0
    except Exception as exc:
        summary["error"] = f"{exc.__class__.__name__}: {exc}"
        summary["traceback"] = traceback.format_exc()
        summary_path = run_dir / "prolif_real_cifs_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())