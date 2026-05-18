#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import traceback
from collections import defaultdict
from pathlib import Path
from typing import Any

from lpmo_pipeline.analysis.convergence_metrics import (
    ConvergencePoseInput,
    compute_condition_convergence,
    write_condition_convergence_summary_tsv,
    write_pose_convergence_tsv,
)
from lpmo_pipeline.io.normalize_mmcif import NormalizeMMCIFRunner
from lpmo_pipeline.io.protonate_export import protonate_and_export


DEFAULT_CIF_PATHS: tuple[Path, ...] = (
    Path(
        "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/"
        "NAG4/af3/latest/Q7SCE9_NAG4/seed-1_sample-0/Q7SCE9_NAG4_seed-1_sample-0_model.cif"
    ),
    Path(
        "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/"
        "NAG4/af3/latest/Q7SCE9_NAG4/seed-1_sample-1/Q7SCE9_NAG4_seed-1_sample-1_model.cif"
    ),
    Path(
        "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/"
        "NAG4/af3/latest/Q7SCE9_NAG4/seed-4_sample-1/Q7SCE9_NAG4_seed-4_sample-1_model.cif"
    ),
    Path(
        "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/"
        "STA6/af3/latest/Q59930_STA6/seed-1_sample-0/Q59930_STA6_seed-1_sample-0_model.cif"
    ),
    Path(
        "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/"
        "STA6/af3/latest/Q59930_STA6/seed-1_sample-1/Q59930_STA6_seed-1_sample-1_model.cif"
    ),
    Path(
        "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/"
        "STA6/af3/latest/Q59930_STA6/seed-2_sample-0/Q59930_STA6_seed-2_sample-0_model.cif"
    ),
)

_TARGET_PREFIX_TO_SUBSTRATE = {
    "NAG": "chitin",
    "CEL": "cellulose",
    "BGC": "cellulose",
    "STA": "amylose",
    "GLC": "amylose",
}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run normalization, protonation, and convergence metrics on a small set of real AF3 CIFs.",
    )
    parser.add_argument(
        "--run-dir",
        required=True,
        help="Output directory for per-case artifacts and convergence summaries.",
    )
    parser.add_argument(
        "--run-id",
        default="",
        help="Optional run identifier written into the summary output.",
    )
    parser.add_argument(
        "cifs",
        nargs="*",
        help="Optional explicit CIF paths. Defaults to a small multi-pose real-CIF set.",
    )
    return parser.parse_args()


def _parse_pose_id(pose_id: str) -> dict[str, Any]:
    match = re.match(
        r"^(?P<protein>[^_]+)_(?P<ligand>[^_]+)_seed-(?P<seed>\d+)_sample-(?P<sample>\d+)_model$",
        pose_id,
    )
    if not match:
        raise ValueError(f"Unexpected pose_id format: {pose_id}")
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
            summary_path = run_dir / "convergence_real_cifs_summary.json"
            summary_path.write_text(json.dumps(summary, indent=2))
            print(json.dumps(summary, indent=2))
            return 1

    try:
        condition_inputs: dict[str, list[ConvergencePoseInput]] = defaultdict(list)
        case_by_pose_id: dict[str, dict[str, Any]] = {}

        for index, cif_path in enumerate(cif_paths, start=1):
            pose_id = cif_path.stem
            metadata = _parse_pose_id(pose_id)
            protein_id = str(metadata["protein_id"])
            ligand_id = str(metadata["ligand_id"])
            seed = int(metadata["seed"])
            sample = int(metadata["sample"])
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
                case_by_pose_id[pose_id] = case_summary
                continue

            normalized_path = Path(normalized_path).resolve()
            case_summary["normalized_cif"] = str(normalized_path)

            protonation_dir = case_dir / "protonated"
            ok_prot, report = protonate_and_export(normalized_path, protonation_dir)
            case_summary["protonate_ok"] = bool(ok_prot)
            case_summary["protonation_report_path"] = str(protonation_dir / "protonation_report.json")
            case_summary["protonation_blockers"] = (report or {}).get("blockers", []) if report else []
            case_summary["protonation_warnings"] = (report or {}).get("warnings", []) if report else []
            if not ok_prot:
                case_summary["status"] = "protonation_failed"
                summary["cases"].append(case_summary)
                case_by_pose_id[pose_id] = case_summary
                continue

            complex_pdb = Path(str((report or {}).get("complex_h_pdb") or protonation_dir / "complex_H.pdb")).resolve()
            case_summary["complex_h_pdb"] = str(complex_pdb)
            case_summary["status"] = "prepared"

            condition_inputs[condition_id].append(
                ConvergencePoseInput(
                    pose_id=pose_id,
                    condition_id=condition_id,
                    complex_pdb=complex_pdb,
                    seed=seed,
                    sample=sample,
                )
            )
            summary["cases"].append(case_summary)
            case_by_pose_id[pose_id] = case_summary

        pose_metrics = []
        condition_summaries = []
        for condition_id in sorted(condition_inputs):
            if len(condition_inputs[condition_id]) < 2:
                raise ValueError(
                    f"Convergence real-data harness requires at least two poses per condition: {condition_id}"
                )

            metrics, condition_summary = compute_condition_convergence(condition_inputs[condition_id])
            pose_metrics.extend(metrics)
            condition_summaries.append(condition_summary)

            for metric in metrics:
                case_summary = case_by_pose_id[metric.pose_id]
                case_summary["reference_pose_id"] = metric.reference_pose_id
                case_summary["ligand_rmsd_to_reference"] = metric.ligand_rmsd_to_reference
                case_summary["convergent_flag"] = metric.convergent_flag
                case_summary["status"] = "ok"

        for condition_summary in condition_summaries:
            summary.setdefault("condition_summaries", []).append(condition_summary.to_row())

        pose_convergence_path = run_dir / "pose_convergence.tsv"
        condition_convergence_path = run_dir / "condition_convergence_summary.tsv"
        write_pose_convergence_tsv(pose_metrics, pose_convergence_path)
        write_condition_convergence_summary_tsv(condition_summaries, condition_convergence_path)

        summary.update(
            {
                "pose_convergence_tsv": str(pose_convergence_path),
                "condition_convergence_summary_tsv": str(condition_convergence_path),
                "pose_convergence_rows": _count_table_rows(pose_convergence_path),
                "condition_convergence_rows": _count_table_rows(condition_convergence_path),
                "n_ok": sum(1 for case in summary["cases"] if case.get("status") == "ok"),
                "n_errors": sum(1 for case in summary["cases"] if case.get("status") not in {"ok"}),
            }
        )

        summary_path = run_dir / "convergence_real_cifs_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))

        if summary["n_ok"] == 0:
            return 1
        if summary["n_ok"] != summary["n_requested"]:
            return 1
        if not pose_convergence_path.exists() or not condition_convergence_path.exists():
            return 1
        if any(case.get("reference_pose_id") != f"{case['protein_id']}_{case['ligand_id']}_seed-1_sample-0_model" for case in summary["cases"] if case.get("status") == "ok"):
            return 1
        if any(case.get("ligand_rmsd_to_reference") is None for case in summary["cases"] if case.get("status") == "ok"):
            return 1
        if any(not (0.0 <= float(condition_summary["convergence_fraction"]) <= 1.0) for condition_summary in summary.get("condition_summaries", [])):
            return 1
        return 0
    except Exception as exc:
        summary["error"] = f"{exc.__class__.__name__}: {exc}"
        summary["traceback"] = traceback.format_exc()
        summary_path = run_dir / "convergence_real_cifs_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())