#!/usr/bin/env python3
"""Small real-CIF timing/profiling harness for the clustering-pilot path.

This runner intentionally does not use the full pilot manifest. It discovers a
small, explicit set of real AF3 CIF poses and runs controlled timing modes so
worker scaling can be checked before the full pilot is submitted again.
"""
from __future__ import annotations

import argparse
import cProfile
import csv
import io
import json
import pstats
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import Any

import yaml

from lpmo_pipeline.analysis.analysis_orchestrator import (
    ProductionRunOptions,
    _discover_pose_inputs,
)
from lpmo_pipeline.analysis.clustering_pilot_real_case import (
    DEFAULT_DOMAIN_ONLY_WORK_ROOT,
    DEFAULT_FULL_LENGTH_WORK_ROOT,
    load_pilot_protein_overview_tsv,
    stage_selected_pose_inputs,
)
from lpmo_pipeline.cli import cmd_run
from lpmo_pipeline.io.cif_to_pdb import convert_cif_to_pdb
from lpmo_pipeline.io.normalize_mmcif import NormalizeMMCIFRunner
from lpmo_pipeline.qc.privateer_runner import prepare_privateer_input


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OVERVIEW_TSV = REPO_ROOT / "input_data" / "clustering_pilot_protein_overview.tsv"
DEFAULT_RUN_DIR = REPO_ROOT / "tests" / "tests_results" / "clustering_pilot_timing_probe"
TIMING_EVENT_FIELDS = [
    "mode",
    "construct_type",
    "n_jobs",
    "step",
    "pose_id",
    "condition_id",
    "protein_id",
    "ligand_id",
    "worker_count",
    "wall_time_s",
    "status",
    "output_dir",
    "detail",
]


@dataclass(frozen=True)
class TimingPose:
    construct_type: str
    pose_id: str
    protein_id: str
    ligand_id: str
    model: str
    source_run_id: str
    discovered_run_id: str | None
    cif_path: Path
    confidence_json_path: Path | None
    seed: int | None
    sample: int | None

    def to_stage_payload(self) -> dict[str, Any]:
        return {
            "pose_id": self.pose_id,
            "protein_id": self.protein_id,
            "ligand_id": self.ligand_id,
            "model": self.model,
            "source_run_id": self.source_run_id,
            "discovered_run_id": self.discovered_run_id,
            "confidence_json_path": str(self.confidence_json_path) if self.confidence_json_path else None,
            "run_status": "succeeded",
            "seed": self.seed,
            "sample": self.sample,
            "cif_path": str(self.cif_path),
        }

    def to_summary_row(self) -> dict[str, Any]:
        return {
            "construct_type": self.construct_type,
            "pose_id": self.pose_id,
            "protein_id": self.protein_id,
            "ligand_id": self.ligand_id,
            "model": self.model,
            "source_run_id": self.source_run_id,
            "discovered_run_id": self.discovered_run_id,
            "seed": self.seed,
            "sample": self.sample,
            "cif_path": str(self.cif_path),
            "confidence_json_path": str(self.confidence_json_path or ""),
        }


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a 12-pose timing/profiling probe for clustering-pilot runtime work.",
    )
    parser.add_argument("--run-dir", default=str(DEFAULT_RUN_DIR), help="Output directory for timing artifacts.")
    parser.add_argument("--overview-tsv", default=str(DEFAULT_OVERVIEW_TSV), help="Pilot overview TSV used only to prioritize proteins.")
    parser.add_argument("--domain-work-root", default=str(DEFAULT_DOMAIN_ONLY_WORK_ROOT), help="Domain-only AF3 work root.")
    parser.add_argument("--full-work-root", default=str(DEFAULT_FULL_LENGTH_WORK_ROOT), help="Full-length AF3 work root.")
    parser.add_argument("--n-domain-poses", type=int, default=8, help="Preferred number of domain-only poses.")
    parser.add_argument("--n-full-poses", type=int, default=4, help="Preferred number of full-length poses.")
    parser.add_argument("--total-poses", type=int, default=12, help="Total number of poses required for the probe.")
    parser.add_argument("--n-jobs", type=int, nargs="+", default=[1, 2, 4], help="Worker counts to test.")
    parser.add_argument(
        "--modes",
        nargs="+",
        choices=["prepare_only", "pilot_like", "qc_full"],
        default=["prepare_only", "pilot_like", "qc_full"],
        help="Timing modes to run.",
    )
    parser.add_argument("--del-branch", default="del_a", choices=["del_a", "del_b"], help="DEL branch for production runs.")
    parser.add_argument("--run-id", default="", help="Optional run id prefix.")
    parser.add_argument(
        "--skip-normalization-profile",
        action="store_true",
        help="Skip cProfile normalization profiling after timing modes.",
    )
    return parser.parse_args()


def _event(
    *,
    mode: str,
    construct_type: str,
    n_jobs: int,
    step: str,
    wall_time_s: float,
    status: str,
    pose: TimingPose | None = None,
    condition_id: str = "",
    output_dir: Path | None = None,
    detail: str = "",
    worker_count: int | None = None,
) -> dict[str, Any]:
    return {
        "mode": mode,
        "construct_type": construct_type,
        "n_jobs": int(n_jobs),
        "step": step,
        "pose_id": pose.pose_id if pose else "",
        "condition_id": condition_id,
        "protein_id": pose.protein_id if pose else "",
        "ligand_id": pose.ligand_id if pose else "",
        "worker_count": int(worker_count if worker_count is not None else n_jobs),
        "wall_time_s": round(float(wall_time_s), 6),
        "status": status,
        "output_dir": str(output_dir or ""),
        "detail": detail,
    }


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TIMING_EVENT_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in TIMING_EVENT_FIELDS})


def _load_priority_proteins(overview_tsv: Path, construct_type: str) -> tuple[str, ...]:
    rows = load_pilot_protein_overview_tsv(overview_tsv)
    if construct_type == "domain_only":
        return tuple(row.uniprot_id for row in rows if row.run_domain_only)
    if construct_type == "full_length":
        return tuple(row.uniprot_id for row in rows if row.run_full_length)
    raise ValueError(f"Unsupported construct_type: {construct_type}")


def _discover_candidates(
    *,
    work_root: Path,
    construct_type: str,
    protein_ids: tuple[str, ...],
) -> tuple[list[TimingPose], list[str], dict[str, Any]]:
    options = ProductionRunOptions(
        run_id=f"timing_probe_{construct_type}",
        output_dir=Path("/tmp") / f"timing_probe_{construct_type}",
        del_variant="a",
        work_root=work_root,
        af3_only=True,
        latest_only=True,
        include_proteins=protein_ids,
        construct_type=construct_type,
        run_posebusters=False,
        run_privateer=False,
        n_jobs=1,
    )
    errors, pose_inputs, summary = _discover_pose_inputs(options)
    candidates = [
        TimingPose(
            construct_type=construct_type,
            pose_id=pose.pose_id,
            protein_id=pose.protein_id,
            ligand_id=pose.ligand_id,
            model=pose.model,
            source_run_id=pose.source_run_id,
            discovered_run_id=pose.discovered_run_id,
            cif_path=pose.cif_path,
            confidence_json_path=pose.confidence_json_path,
            seed=pose.seed,
            sample=pose.sample,
        )
        for pose in pose_inputs
    ]
    candidates.sort(
        key=lambda pose: (
            pose.protein_id,
            pose.ligand_id,
            -1 if pose.seed is None else pose.seed,
            -1 if pose.sample is None else pose.sample,
            pose.pose_id,
        )
    )
    return candidates, errors, summary


def _balanced_select(candidates: list[TimingPose], count: int, *, exclude: set[tuple[str, str]]) -> list[TimingPose]:
    buckets: dict[tuple[str, str], list[TimingPose]] = {}
    for pose in candidates:
        key = (pose.protein_id, pose.ligand_id)
        if (pose.construct_type, pose.pose_id) in exclude:
            continue
        buckets.setdefault(key, []).append(pose)

    selected: list[TimingPose] = []
    bucket_keys = sorted(buckets)
    while len(selected) < count and bucket_keys:
        next_keys: list[tuple[str, str]] = []
        for key in bucket_keys:
            bucket = buckets[key]
            if bucket and len(selected) < count:
                pose = bucket.pop(0)
                selected.append(pose)
                exclude.add((pose.construct_type, pose.pose_id))
            if bucket:
                next_keys.append(key)
        bucket_keys = next_keys
    return selected


def select_timing_poses(
    *,
    overview_tsv: Path,
    domain_work_root: Path,
    full_work_root: Path,
    n_domain_poses: int,
    n_full_poses: int,
    total_poses: int,
) -> tuple[list[TimingPose], dict[str, Any]]:
    domain_candidates, domain_errors, domain_summary = _discover_candidates(
        work_root=domain_work_root,
        construct_type="domain_only",
        protein_ids=_load_priority_proteins(overview_tsv, "domain_only"),
    )
    full_candidates, full_errors, full_summary = _discover_candidates(
        work_root=full_work_root,
        construct_type="full_length",
        protein_ids=_load_priority_proteins(overview_tsv, "full_length"),
    )

    seen: set[tuple[str, str]] = set()
    selected_domain = _balanced_select(domain_candidates, n_domain_poses, exclude=seen)
    selected_full = _balanced_select(full_candidates, n_full_poses, exclude=seen)
    selected = selected_domain + selected_full

    if len(selected) < total_poses:
        selected.extend(
            _balanced_select(
                domain_candidates + full_candidates,
                total_poses - len(selected),
                exclude=seen,
            )
        )

    if len(selected) != total_poses:
        raise RuntimeError(
            f"Timing probe requires {total_poses} poses, but only selected {len(selected)}"
        )

    summary = {
        "domain_work_root": str(domain_work_root),
        "full_work_root": str(full_work_root),
        "n_domain_candidates": len(domain_candidates),
        "n_full_candidates": len(full_candidates),
        "n_selected_domain": sum(1 for pose in selected if pose.construct_type == "domain_only"),
        "n_selected_full_length": sum(1 for pose in selected if pose.construct_type == "full_length"),
        "domain_discovery_errors": domain_errors,
        "full_length_discovery_errors": full_errors,
        "domain_discovery_summary": domain_summary,
        "full_length_discovery_summary": full_summary,
    }
    return selected, summary


def _prepare_one_pose(args: tuple[TimingPose, Path, int]) -> tuple[TimingPose, list[dict[str, Any]], dict[str, Any]]:
    pose, output_root, n_jobs = args
    case_dir = output_root / pose.construct_type / "cases" / pose.pose_id
    case_dir.mkdir(parents=True, exist_ok=True)
    events: list[dict[str, Any]] = []
    case: dict[str, Any] = {
        "construct_type": pose.construct_type,
        "pose_id": pose.pose_id,
        "protein_id": pose.protein_id,
        "ligand_id": pose.ligand_id,
        "case_dir": str(case_dir),
        "status": "preparing",
    }
    total_start = perf_counter()
    try:
        normalize_dir = case_dir / "normalize"
        start = perf_counter()
        normalize_ok, normalized_path = NormalizeMMCIFRunner(pose.cif_path, normalize_dir).run()
        events.append(
            _event(
                mode="prepare_only",
                construct_type=pose.construct_type,
                n_jobs=n_jobs,
                step="prepare.normalize",
                pose=pose,
                wall_time_s=perf_counter() - start,
                status="ok" if normalize_ok else "failed",
                output_dir=normalize_dir,
            )
        )
        if not normalize_ok or normalized_path is None:
            raise RuntimeError("Normalization failed")

        normalized_path = Path(normalized_path).resolve()
        pdb_dir = case_dir / "posebusters_input"
        start = perf_counter()
        pdb_ok, pdb_path = convert_cif_to_pdb(normalized_path, pdb_dir)
        events.append(
            _event(
                mode="prepare_only",
                construct_type=pose.construct_type,
                n_jobs=n_jobs,
                step="prepare.cif_to_pdb",
                pose=pose,
                wall_time_s=perf_counter() - start,
                status="ok" if pdb_ok else "failed",
                output_dir=pdb_dir,
            )
        )
        if not pdb_ok or pdb_path is None:
            raise RuntimeError("CIF to PDB conversion failed")

        privateer_path = case_dir / "privateer_input.cif"
        start = perf_counter()
        prepared_privateer = prepare_privateer_input(normalized_path, privateer_path)
        events.append(
            _event(
                mode="prepare_only",
                construct_type=pose.construct_type,
                n_jobs=n_jobs,
                step="prepare.privateer_input",
                pose=pose,
                wall_time_s=perf_counter() - start,
                status="ok",
                output_dir=prepared_privateer.parent,
            )
        )
        case.update(
            {
                "status": "prepared",
                "normalized_cif": str(normalized_path),
                "posebusters_pdb": str(Path(pdb_path).resolve()),
                "privateer_input_cif": str(prepared_privateer),
            }
        )
    except Exception as exc:
        case.update(
            {
                "status": "prep_error",
                "error": f"{exc.__class__.__name__}: {exc}",
                "traceback": traceback.format_exc(),
            }
        )
    events.append(
        _event(
            mode="prepare_only",
            construct_type=pose.construct_type,
            n_jobs=n_jobs,
            step="prepare.pose_total",
            pose=pose,
            wall_time_s=perf_counter() - total_start,
            status=case["status"],
            output_dir=case_dir,
        )
    )
    return pose, events, case


def run_prepare_only(selected_poses: list[TimingPose], *, run_dir: Path, n_jobs: int) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    output_root = run_dir / f"prepare_only_njobs_{n_jobs}"
    output_root.mkdir(parents=True, exist_ok=True)
    events: list[dict[str, Any]] = []
    cases: list[dict[str, Any]] = []
    start = perf_counter()
    if n_jobs <= 1:
        results = [_prepare_one_pose((pose, output_root, n_jobs)) for pose in selected_poses]
    else:
        with ProcessPoolExecutor(max_workers=min(n_jobs, len(selected_poses))) as executor:
            future_to_pose = {
                executor.submit(_prepare_one_pose, (pose, output_root, n_jobs)): pose
                for pose in selected_poses
            }
            results = [future.result() for future in as_completed(future_to_pose)]
    for _pose, pose_events, case in sorted(results, key=lambda item: item[0].pose_id):
        events.extend(pose_events)
        cases.append(case)
    events.append(
        _event(
            mode="prepare_only",
            construct_type="all",
            n_jobs=n_jobs,
            step="prepare.batch_total",
            wall_time_s=perf_counter() - start,
            status="ok" if all(case["status"] == "prepared" for case in cases) else "error",
            output_dir=output_root,
            detail=f"n_poses={len(selected_poses)}",
        )
    )
    return events, cases


def _write_production_config(
    *,
    config_path: Path,
    run_id: str,
    staged_work_root: Path,
    construct_type: str,
    selected_poses: list[TimingPose],
    run_posebusters: bool,
    run_privateer: bool,
    n_jobs: int,
) -> None:
    payload = {
        "pipeline_version": "2.1",
        "production": {
            "run_id": run_id,
            "work_root": str(staged_work_root),
            "construct_type": construct_type,
            "af3_only": True,
            "latest_only": True,
            "include_targets": sorted({pose.ligand_id for pose in selected_poses}),
            "include_proteins": sorted({pose.protein_id for pose in selected_poses}),
            "run_posebusters": run_posebusters,
            "run_privateer": run_privateer,
            "n_jobs": n_jobs,
            "collect_timing_events": True,
            "clustering_pilot": {
                "enabled": True,
                "label": run_id,
                "minimum_clusterable_n": 10,
                "agglomerative_linkage": "average",
                "agglomerative_distance_threshold": 0.5,
                "agglomerative_min_cluster_size": 10,
            },
        },
    }
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(yaml.safe_dump(payload, sort_keys=False))


def _normalize_analysis_timing_event(
    raw: dict[str, Any],
    *,
    mode: str,
    construct_type: str,
    n_jobs: int,
    output_dir: Path,
) -> dict[str, Any]:
    return {
        "mode": mode,
        "construct_type": construct_type,
        "n_jobs": n_jobs,
        "step": raw.get("step", ""),
        "pose_id": raw.get("pose_id", ""),
        "condition_id": raw.get("condition_id", ""),
        "protein_id": "",
        "ligand_id": "",
        "worker_count": raw.get("worker_count", n_jobs),
        "wall_time_s": raw.get("wall_time_s", 0.0),
        "status": raw.get("status", ""),
        "output_dir": str(output_dir),
        "detail": raw.get("detail", ""),
    }


def run_production_mode(
    selected_poses: list[TimingPose],
    *,
    run_dir: Path,
    mode: str,
    n_jobs: int,
    del_branch: str,
    run_id_prefix: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if mode not in {"pilot_like", "qc_full"}:
        raise ValueError(f"Unsupported production timing mode: {mode}")
    run_posebusters = mode == "qc_full"
    run_privateer = mode == "qc_full"
    events: list[dict[str, Any]] = []
    construct_summaries: list[dict[str, Any]] = []
    for construct_type in ("domain_only", "full_length"):
        construct_poses = [pose for pose in selected_poses if pose.construct_type == construct_type]
        if not construct_poses:
            continue
        construct_dir = run_dir / f"{mode}_n{n_jobs}" / construct_type
        staged_work_root = construct_dir / "staged_work_root"
        staging_start = perf_counter()
        staging_payload = stage_selected_pose_inputs(
            [pose.to_stage_payload() for pose in construct_poses],
            staged_work_root,
        )
        events.append(
            _event(
                mode=mode,
                construct_type=construct_type,
                n_jobs=n_jobs,
                step="production.stage_selected_poses",
                wall_time_s=perf_counter() - staging_start,
                status="ok",
                output_dir=staged_work_root,
                detail=f"n_staged={staging_payload.get('n_staged_poses', 0)}",
            )
        )
        config_path = construct_dir / "production.timing_probe.yaml"
        run_id = f"{run_id_prefix}_{mode}_n{n_jobs}_{construct_type}"
        _write_production_config(
            config_path=config_path,
            run_id=run_id,
            staged_work_root=Path(str(staging_payload["staged_work_root"])),
            construct_type=construct_type,
            selected_poses=construct_poses,
            run_posebusters=run_posebusters,
            run_privateer=run_privateer,
            n_jobs=n_jobs,
        )
        production_output = construct_dir / "production_output"
        production_start = perf_counter()
        exit_code = cmd_run(
            argparse.Namespace(
                mode="production",
                config=config_path,
                output=production_output,
                del_branch=del_branch,
                n_jobs=n_jobs,
            )
        )
        production_wall_time = perf_counter() - production_start
        analysis_summary_path = production_output / "analysis_core_summary.json"
        analysis_summary = json.loads(analysis_summary_path.read_text()) if analysis_summary_path.exists() else {}
        for raw_event in analysis_summary.get("timing_events", []):
            events.append(
                _normalize_analysis_timing_event(
                    raw_event,
                    mode=mode,
                    construct_type=construct_type,
                    n_jobs=n_jobs,
                    output_dir=production_output,
                )
            )
        events.append(
            _event(
                mode=mode,
                construct_type=construct_type,
                n_jobs=n_jobs,
                step="production.run_total",
                wall_time_s=production_wall_time,
                status="ok" if exit_code == 0 else "error",
                output_dir=production_output,
                detail=f"exit_code={exit_code}",
            )
        )
        construct_summaries.append(
            {
                "mode": mode,
                "construct_type": construct_type,
                "n_jobs": n_jobs,
                "run_id": run_id,
                "config_path": str(config_path),
                "production_output": str(production_output),
                "analysis_core_summary_path": str(analysis_summary_path),
                "exit_code": exit_code,
                "n_selected_poses": len(construct_poses),
                "n_prepared": analysis_summary.get("n_prepared", 0),
                "n_analyzed": analysis_summary.get("n_analyzed", 0),
                "pilot_outputs_ready": bool(
                    (analysis_summary.get("clustering_pilot") or {}).get(
                        "pilot_clustering_method_summary_tsv"
                    )
                ),
            }
        )
    return events, construct_summaries


def profile_normalization(selected_poses: list[TimingPose], *, run_dir: Path) -> dict[str, Any]:
    profile_dir = run_dir / "normalization_profile"
    profile_dir.mkdir(parents=True, exist_ok=True)
    profiler = cProfile.Profile()
    pose_summaries: list[dict[str, Any]] = []
    start = perf_counter()
    profiler.enable()
    for pose in selected_poses:
        pose_dir = profile_dir / "cases" / pose.construct_type / pose.pose_id
        pose_start = perf_counter()
        try:
            ok, normalized_path = NormalizeMMCIFRunner(pose.cif_path, pose_dir).run()
            pose_summaries.append(
                {
                    "construct_type": pose.construct_type,
                    "pose_id": pose.pose_id,
                    "status": "ok" if ok else "failed",
                    "wall_time_s": round(perf_counter() - pose_start, 6),
                    "normalized_cif": str(normalized_path or ""),
                }
            )
        except Exception as exc:
            pose_summaries.append(
                {
                    "construct_type": pose.construct_type,
                    "pose_id": pose.pose_id,
                    "status": "error",
                    "wall_time_s": round(perf_counter() - pose_start, 6),
                    "error": f"{exc.__class__.__name__}: {exc}",
                }
            )
    profiler.disable()
    total_wall_time = perf_counter() - start

    stats_path = profile_dir / "normalization_profile.pstats"
    profiler.dump_stats(str(stats_path))
    stats_stream = io.StringIO()
    stats = pstats.Stats(profiler, stream=stats_stream).strip_dirs().sort_stats("cumulative")
    stats.print_stats(40)
    profile_text_path = profile_dir / "normalization_profile_top40.txt"
    profile_text_path.write_text(stats_stream.getvalue())

    stat_rows: list[dict[str, Any]] = []
    for (filename, line_no, function_name), values in stats.stats.items():
        primitive_calls, total_calls, total_time, cumulative_time, _callers = values
        stat_rows.append(
            {
                "filename": filename,
                "line_no": line_no,
                "function_name": function_name,
                "primitive_calls": primitive_calls,
                "total_calls": total_calls,
                "total_time_s": total_time,
                "cumulative_time_s": cumulative_time,
            }
        )
    stat_rows.sort(key=lambda row: float(row["cumulative_time_s"]), reverse=True)

    profile_tsv_path = profile_dir / "normalization_profile_top_cumulative.tsv"
    with profile_tsv_path.open("w", newline="") as handle:
        fieldnames = [
            "filename",
            "line_no",
            "function_name",
            "primitive_calls",
            "total_calls",
            "total_time_s",
            "cumulative_time_s",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in stat_rows[:40]:
            writer.writerow(row)

    remap_rows = [
        row
        for row in stat_rows
        if "remap" in str(row["function_name"]).lower()
        or "chain" in str(row["function_name"]).lower()
        or "remap" in str(row["filename"]).lower()
    ]
    remap_cumulative = sum(float(row["cumulative_time_s"]) for row in remap_rows)
    remap_fraction = remap_cumulative / total_wall_time if total_wall_time > 0 else 0.0
    summary = {
        "profile_dir": str(profile_dir),
        "stats_path": str(stats_path),
        "profile_text_path": str(profile_text_path),
        "profile_tsv_path": str(profile_tsv_path),
        "n_profiled_poses": len(selected_poses),
        "total_wall_time_s": round(total_wall_time, 6),
        "pose_summaries": pose_summaries,
        "top_cumulative_functions": stat_rows[:20],
        "chain_remap_candidate_cumulative_s": round(remap_cumulative, 6),
        "chain_remap_candidate_fraction_of_wall": round(remap_fraction, 6),
        "recommend_chain_schema_change_for_speed": remap_fraction >= 0.10,
        "decision_rule": "Only change AF3-native chain handling for speed if chain/remap cumulative time is >= 10% of normalization wall time.",
    }
    summary_path = profile_dir / "normalization_profile_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    return summary


def _summarize_events(events: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, int, str], dict[str, Any]] = {}
    for event in events:
        key = (
            str(event.get("mode", "")),
            str(event.get("construct_type", "")),
            int(event.get("n_jobs", 0) or 0),
            str(event.get("step", "")),
        )
        row = grouped.setdefault(
            key,
            {
                "mode": key[0],
                "construct_type": key[1],
                "n_jobs": key[2],
                "step": key[3],
                "event_count": 0,
                "wall_time_s_sum": 0.0,
                "error_count": 0,
            },
        )
        row["event_count"] += 1
        row["wall_time_s_sum"] += float(event.get("wall_time_s", 0.0) or 0.0)
        if str(event.get("status", "")) in {"error", "failed", "prep_error"}:
            row["error_count"] += 1
    return [
        {**row, "wall_time_s_sum": round(float(row["wall_time_s_sum"]), 6)}
        for row in sorted(grouped.values(), key=lambda row: (row["mode"], row["n_jobs"], row["construct_type"], row["step"]))
    ]


def main() -> int:
    args = _parse_args()
    run_dir = Path(args.run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    run_id_prefix = args.run_id or run_dir.name
    n_jobs_values = sorted({max(1, int(value)) for value in args.n_jobs})
    events: list[dict[str, Any]] = []
    prepare_cases: list[dict[str, Any]] = []
    production_summaries: list[dict[str, Any]] = []
    profile_summary: dict[str, Any] | None = None

    summary: dict[str, Any] = {
        "run_dir": str(run_dir),
        "run_id_prefix": run_id_prefix,
        "overview_tsv": str(Path(args.overview_tsv).resolve()),
        "n_jobs_values": n_jobs_values,
        "modes": list(args.modes),
        "requested_domain_poses": args.n_domain_poses,
        "requested_full_length_poses": args.n_full_poses,
        "requested_total_poses": args.total_poses,
    }

    try:
        selected_poses, selection_summary = select_timing_poses(
            overview_tsv=Path(args.overview_tsv).resolve(),
            domain_work_root=Path(args.domain_work_root).resolve(),
            full_work_root=Path(args.full_work_root).resolve(),
            n_domain_poses=args.n_domain_poses,
            n_full_poses=args.n_full_poses,
            total_poses=args.total_poses,
        )
        summary["selection"] = selection_summary
        summary["selected_poses"] = [pose.to_summary_row() for pose in selected_poses]

        for n_jobs in n_jobs_values:
            if "prepare_only" in args.modes:
                mode_events, mode_cases = run_prepare_only(selected_poses, run_dir=run_dir, n_jobs=n_jobs)
                events.extend(mode_events)
                prepare_cases.extend(mode_cases)
            for mode in ("pilot_like", "qc_full"):
                if mode not in args.modes:
                    continue
                mode_events, mode_summaries = run_production_mode(
                    selected_poses,
                    run_dir=run_dir,
                    mode=mode,
                    n_jobs=n_jobs,
                    del_branch=args.del_branch,
                    run_id_prefix=run_id_prefix,
                )
                events.extend(mode_events)
                production_summaries.extend(mode_summaries)

        if not args.skip_normalization_profile:
            profile_summary = profile_normalization(selected_poses, run_dir=run_dir)

        timing_events_path = run_dir / "timing_events.tsv"
        _write_tsv(timing_events_path, events)
        summary.update(
            {
                "timing_events_tsv": str(timing_events_path),
                "event_summary": _summarize_events(events),
                "prepare_cases": prepare_cases,
                "production_summaries": production_summaries,
                "normalization_profile": profile_summary,
            }
        )
        summary_path = run_dir / "timing_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        if any(row.get("error_count", 0) for row in summary["event_summary"]):
            return 1
        if any(item.get("exit_code", 0) != 0 for item in production_summaries):
            return 1
        return 0
    except Exception as exc:
        summary["error"] = f"{exc.__class__.__name__}: {exc}"
        summary["traceback"] = traceback.format_exc()
        summary_path = run_dir / "timing_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
