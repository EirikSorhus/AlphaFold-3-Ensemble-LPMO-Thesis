"""Helpers for preparing and resuming real-case clustering-pilot runs."""
from __future__ import annotations

import csv
import json
import shlex
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Collection

import yaml

from lpmo_pipeline.analysis.analysis_orchestrator import ProductionRunOptions, _discover_pose_inputs


DEFAULT_DOMAIN_ONLY_WORK_ROOT = Path(
    "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core"
)
DEFAULT_FULL_LENGTH_WORK_ROOT = Path(
    "/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_full_length"
)
PILOT_KJORING_DOMAIN_ONLY = "Kun katalytisk domene"
PILOT_KJORING_DOMAIN_AND_FULL_LENGTH = "Katalytisk + full-lengde"


@dataclass(frozen=True)
class PilotProteinOverviewRow:
    """One row from the curated pilot protein overview table."""

    pilot_rank: int
    uniprot_id: str
    family: str
    protein_name: str
    organism: str
    ec_numbers: str
    regio: str
    ligand_category: str
    oligo_active_classes: str
    oligo_active_substrates: str
    has_cbm: str
    binding_modules: str
    pilot_kjoring: str
    run_domain_only: bool
    run_full_length: bool
    rationale: str

    def to_row(self) -> dict[str, str | int | bool]:
        return {
            "pilot_rank": self.pilot_rank,
            "uniprot_id": self.uniprot_id,
            "family": self.family,
            "protein_name": self.protein_name,
            "organism": self.organism,
            "ec_numbers": self.ec_numbers,
            "regio": self.regio,
            "ligand_category": self.ligand_category,
            "oligo_active_classes": self.oligo_active_classes,
            "oligo_active_substrates": self.oligo_active_substrates,
            "has_cbm": self.has_cbm,
            "binding_modules": self.binding_modules,
            "pilot_kjoring": self.pilot_kjoring,
            "run_domain_only": self.run_domain_only,
            "run_full_length": self.run_full_length,
            "rationale": self.rationale,
        }


@dataclass(frozen=True)
class PilotSelection:
    """One protein-target selection for a real-case pilot run."""

    protein_id: str
    target: str | None = None


@dataclass(frozen=True)
class PilotSelectionManifest:
    """Manifest describing a resumable real-case clustering pilot run."""

    work_root: Path
    construct_type: str = "domain_only"
    af3_only: bool = True
    latest_only: bool = True
    run_posebusters: bool = False
    run_privateer: bool = False
    clustering_pilot: dict[str, Any] = field(default_factory=dict)
    selections: tuple[PilotSelection, ...] = ()


def _parse_bool(value: str | None, *, default: bool = False) -> bool:
    if value is None:
        return default
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "ja", "y"}:
        return True
    if normalized in {"0", "false", "no", "nei", "n"}:
        return False
    return default


def load_pilot_protein_overview_tsv(overview_path: Path) -> list[PilotProteinOverviewRow]:
    """Load the curated pilot protein overview table."""
    required_columns = {
        "pilot_rank",
        "uniprot_id",
        "family",
        "protein_name",
        "organism",
        "ec_numbers",
        "regio",
        "ligand_category",
        "oligo_active_classes",
        "oligo_active_substrates",
        "has_cbm",
        "binding_modules",
        "pilot_kjoring",
        "rationale",
    }

    with overview_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Pilot overview TSV must contain a header row")
        missing_columns = sorted(required_columns - set(reader.fieldnames))
        if missing_columns:
            raise ValueError(f"Pilot overview TSV missing required columns: {missing_columns}")

        rows: list[PilotProteinOverviewRow] = []
        seen_uniprot_ids: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            uniprot_id = str(row.get("uniprot_id") or "").strip()
            if not uniprot_id:
                raise ValueError(f"Pilot overview TSV line {line_number} is missing uniprot_id")
            if uniprot_id in seen_uniprot_ids:
                raise ValueError(f"Duplicate uniprot_id in pilot overview TSV: {uniprot_id}")
            seen_uniprot_ids.add(uniprot_id)

            pilot_kjoring = str(row.get("pilot_kjoring") or "").strip()
            run_domain_only = _parse_bool(
                row.get("run_domain_only"),
                default=pilot_kjoring in {PILOT_KJORING_DOMAIN_ONLY, PILOT_KJORING_DOMAIN_AND_FULL_LENGTH},
            )
            run_full_length = _parse_bool(
                row.get("run_full_length"),
                default=pilot_kjoring == PILOT_KJORING_DOMAIN_AND_FULL_LENGTH,
            )

            rows.append(
                PilotProteinOverviewRow(
                    pilot_rank=int(str(row.get("pilot_rank") or "0")),
                    uniprot_id=uniprot_id,
                    family=str(row.get("family") or "").strip(),
                    protein_name=str(row.get("protein_name") or "").strip(),
                    organism=str(row.get("organism") or "").strip(),
                    ec_numbers=str(row.get("ec_numbers") or "").strip(),
                    regio=str(row.get("regio") or "").strip(),
                    ligand_category=str(row.get("ligand_category") or "").strip(),
                    oligo_active_classes=str(row.get("oligo_active_classes") or "").strip(),
                    oligo_active_substrates=str(row.get("oligo_active_substrates") or "").strip(),
                    has_cbm=str(row.get("has_cbm") or "").strip(),
                    binding_modules=str(row.get("binding_modules") or "").strip(),
                    pilot_kjoring=pilot_kjoring,
                    run_domain_only=run_domain_only,
                    run_full_length=run_full_length,
                    rationale=str(row.get("rationale") or "").strip(),
                )
            )

    rows.sort(key=lambda row: (row.pilot_rank, row.uniprot_id))
    return rows


def build_selection_manifest_from_overview(
    overview_rows: Collection[PilotProteinOverviewRow],
    *,
    construct_type: str,
    work_root: Path | None = None,
    clustering_pilot: dict[str, Any] | None = None,
    latest_only: bool = True,
    run_posebusters: bool = False,
    run_privateer: bool = False,
) -> PilotSelectionManifest:
    """Build a selection manifest from the curated protein overview."""
    if construct_type not in {"domain_only", "full_length"}:
        raise ValueError(f"Unsupported construct_type: {construct_type}")

    if construct_type == "domain_only":
        selected_rows = [row for row in overview_rows if row.run_domain_only]
        resolved_work_root = (work_root or DEFAULT_DOMAIN_ONLY_WORK_ROOT).resolve()
    else:
        selected_rows = [row for row in overview_rows if row.run_full_length]
        resolved_work_root = (work_root or DEFAULT_FULL_LENGTH_WORK_ROOT).resolve()

    selections = tuple(PilotSelection(protein_id=row.uniprot_id) for row in selected_rows)
    return PilotSelectionManifest(
        work_root=resolved_work_root,
        construct_type=construct_type,
        latest_only=latest_only,
        run_posebusters=run_posebusters,
        run_privateer=run_privateer,
        clustering_pilot=dict(clustering_pilot or {}),
        selections=selections,
    )


def write_selection_manifest(manifest: PilotSelectionManifest, output_path: Path) -> None:
    """Write a pilot selection manifest to YAML."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(yaml.safe_dump(serialize_selection_manifest(manifest), sort_keys=False))


def load_selection_manifest(manifest_path: Path) -> PilotSelectionManifest:
    """Load and validate the real-case pilot selection manifest."""
    payload = yaml.safe_load(manifest_path.read_text()) or {}
    if not isinstance(payload, dict):
        raise TypeError("Pilot selection manifest must be a mapping")

    raw_selections = payload.get("selections") or []
    if not isinstance(raw_selections, list):
        raise TypeError("Pilot selection manifest field 'selections' must be a list")

    selections: list[PilotSelection] = []
    seen_pairs: set[tuple[str, str]] = set()
    for index, selection in enumerate(raw_selections, start=1):
        if not isinstance(selection, dict):
            raise TypeError(f"Selection #{index} must be a mapping")
        protein_id = str(selection.get("protein_id") or "").strip()
        target_raw = selection.get("target")
        target = None if target_raw in {None, ""} else str(target_raw).strip()
        if not protein_id:
            raise ValueError(f"Selection #{index} must define protein_id")
        pair = (target, protein_id)
        if pair in seen_pairs:
            raise ValueError(
                f"Duplicate pilot selection: {protein_id}/{target or 'all_targets'}"
            )
        seen_pairs.add(pair)
        selections.append(PilotSelection(protein_id=protein_id, target=target))

    work_root = Path(str(payload.get("work_root") or "")).resolve()
    if not str(work_root):
        raise ValueError("Pilot selection manifest must define work_root")

    clustering_pilot = payload.get("clustering_pilot") or {}
    if not isinstance(clustering_pilot, dict):
        raise TypeError("Pilot selection manifest field 'clustering_pilot' must be a mapping")

    return PilotSelectionManifest(
        work_root=work_root,
        construct_type=str(payload.get("construct_type") or "domain_only"),
        af3_only=bool(payload.get("af3_only", True)),
        latest_only=bool(payload.get("latest_only", True)),
        run_posebusters=bool(payload.get("run_posebusters", False)),
        run_privateer=bool(payload.get("run_privateer", False)),
        clustering_pilot={str(key): value for key, value in clustering_pilot.items()},
        selections=tuple(selections),
    )


def serialize_selection_manifest(manifest: PilotSelectionManifest) -> dict[str, Any]:
    """Convert the manifest to a JSON-friendly snapshot."""
    payload = {
        "work_root": str(manifest.work_root),
        "construct_type": manifest.construct_type,
        "af3_only": manifest.af3_only,
        "latest_only": manifest.latest_only,
        "run_posebusters": manifest.run_posebusters,
        "run_privateer": manifest.run_privateer,
        "clustering_pilot": manifest.clustering_pilot,
        "selections": [],
    }

    for selection in manifest.selections:
        row = {"protein_id": selection.protein_id}
        if selection.target:
            row["target"] = selection.target
        payload["selections"].append(row)

    return payload


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text())


def _paths_exist(paths: Collection[Path]) -> bool:
    return all(path.exists() or path.is_symlink() for path in paths)


def _load_or_build_json_checkpoint(
    checkpoint_path: Path,
    *,
    force: bool,
    required_paths: Collection[Path],
    build: callable,
) -> tuple[Any, bool]:
    if checkpoint_path.exists() and not force and _paths_exist(required_paths):
        return _read_json(checkpoint_path), True

    payload = build()
    _write_json(checkpoint_path, payload)
    return payload, False


def _serialize_pose_input(pose: Any) -> dict[str, Any]:
    return {
        "pose_id": pose.pose_id,
        "protein_id": pose.protein_id,
        "ligand_id": pose.ligand_id,
        "model": pose.model,
        "source_run_id": pose.source_run_id,
        "discovered_run_id": pose.discovered_run_id,
        "confidence_json_path": str(pose.confidence_json_path) if pose.confidence_json_path else None,
        "run_status": pose.run_status,
        "seed": pose.seed,
        "sample": pose.sample,
        "cif_path": str(pose.cif_path),
    }


def discover_selected_pose_inputs(manifest: PilotSelectionManifest) -> dict[str, Any]:
    """Resolve live work_root selections into concrete pose inputs."""
    selections_by_target: dict[str, set[str]] = {}
    protein_only_selections: set[str] = set()
    for selection in manifest.selections:
        if selection.target:
            selections_by_target.setdefault(selection.target, set()).add(selection.protein_id)
        else:
            protein_only_selections.add(selection.protein_id)

    selection_entries: list[dict[str, Any]] = []
    target_discovery_summaries: list[dict[str, Any]] = []
    all_pose_inputs_by_key: dict[tuple[str, str], dict[str, Any]] = {}
    discovery_errors: list[str] = []
    missing_selections: list[dict[str, str]] = []

    if protein_only_selections:
        requested_proteins = tuple(sorted(protein_only_selections))
        options = ProductionRunOptions(
            run_id="clustering-pilot-real-case-discovery",
            output_dir=manifest.work_root,
            del_variant="a",
            work_root=manifest.work_root,
            af3_only=manifest.af3_only,
            latest_only=manifest.latest_only,
            include_targets=(),
            include_proteins=requested_proteins,
            construct_type=manifest.construct_type,
            run_posebusters=manifest.run_posebusters,
            run_privateer=manifest.run_privateer,
        )
        target_errors, pose_inputs, discovery_summary = _discover_pose_inputs(options)
        discovery_errors.extend(target_errors)
        target_discovery_summaries.append(
            {
                "selection_scope": "protein_only",
                "requested_proteins": list(requested_proteins),
                "n_discovered": len(pose_inputs),
                "discovery_summary": discovery_summary,
                "discovery_errors": list(target_errors),
            }
        )

        for protein_id in requested_proteins:
            selected_poses = [pose for pose in pose_inputs if pose.protein_id == protein_id]
            if not selected_poses:
                missing_selections.append({"protein_id": protein_id})
            serialized_poses = [_serialize_pose_input(pose) for pose in selected_poses]
            selection_entries.append(
                {
                    "protein_id": protein_id,
                    "target": None,
                    "n_discovered": len(serialized_poses),
                    "discovered_targets": sorted({pose["ligand_id"] for pose in serialized_poses}),
                    "run_ids_seen": sorted(
                        {
                            str(pose.get("discovered_run_id") or pose.get("source_run_id") or "")
                            for pose in serialized_poses
                        }
                    ),
                    "first_pose_ids": [pose["pose_id"] for pose in serialized_poses[:10]],
                    "first_source_paths": [pose["cif_path"] for pose in serialized_poses[:10]],
                }
            )
            for pose in serialized_poses:
                all_pose_inputs_by_key[(str(pose["pose_id"]), str(pose["cif_path"]))] = pose

    for target in sorted(selections_by_target):
        requested_proteins = tuple(sorted(selections_by_target[target]))
        options = ProductionRunOptions(
            run_id="clustering-pilot-real-case-discovery",
            output_dir=manifest.work_root,
            del_variant="a",
            work_root=manifest.work_root,
            af3_only=manifest.af3_only,
            latest_only=manifest.latest_only,
            include_targets=(target,),
            include_proteins=requested_proteins,
            construct_type=manifest.construct_type,
            run_posebusters=manifest.run_posebusters,
            run_privateer=manifest.run_privateer,
        )
        target_errors, pose_inputs, discovery_summary = _discover_pose_inputs(options)
        discovery_errors.extend(target_errors)
        target_discovery_summaries.append(
            {
                "target": target,
                "requested_proteins": list(requested_proteins),
                "n_discovered": len(pose_inputs),
                "discovery_summary": discovery_summary,
                "discovery_errors": list(target_errors),
            }
        )

        for protein_id in requested_proteins:
            selected_poses = [
                pose
                for pose in pose_inputs
                if pose.protein_id == protein_id and pose.ligand_id == target
            ]
            if not selected_poses:
                missing_selections.append({"target": target, "protein_id": protein_id})
            serialized_poses = [_serialize_pose_input(pose) for pose in selected_poses]
            selection_entries.append(
                {
                    "target": target,
                    "protein_id": protein_id,
                    "n_discovered": len(serialized_poses),
                    "discovered_targets": [target] if serialized_poses else [],
                    "run_ids_seen": sorted(
                        {
                            str(pose.get("discovered_run_id") or pose.get("source_run_id") or "")
                            for pose in serialized_poses
                        }
                    ),
                    "first_pose_ids": [pose["pose_id"] for pose in serialized_poses[:10]],
                    "first_source_paths": [pose["cif_path"] for pose in serialized_poses[:10]],
                }
            )
            for pose in serialized_poses:
                all_pose_inputs_by_key[(str(pose["pose_id"]), str(pose["cif_path"]))] = pose

    all_pose_inputs = list(all_pose_inputs_by_key.values())

    all_pose_inputs.sort(
        key=lambda pose: (
            str(pose["ligand_id"]),
            str(pose["protein_id"]),
            str(pose["model"]),
            str(pose.get("source_run_id") or ""),
            -1 if pose.get("seed") is None else int(pose["seed"]),
            -1 if pose.get("sample") is None else int(pose["sample"]),
            str(pose["cif_path"]),
        )
    )

    return {
        "work_root": str(manifest.work_root),
        "construct_type": manifest.construct_type,
        "latest_only": manifest.latest_only,
        "n_selected_pairs": len(manifest.selections),
        "n_pose_inputs": len(all_pose_inputs),
        "selection_entries": selection_entries,
        "target_discovery_summaries": target_discovery_summaries,
        "missing_selections": missing_selections,
        "discovery_errors": discovery_errors,
        "pose_inputs": all_pose_inputs,
    }


def _replace_symlink(destination: Path, source: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        destination.unlink()
    destination.symlink_to(source)


def _select_latest_run_id(run_ids: Collection[str]) -> str:
    numeric_run_ids = [int(run_id) for run_id in run_ids if str(run_id).isdigit()]
    if numeric_run_ids:
        return str(max(numeric_run_ids))
    return max(str(run_id) for run_id in run_ids)


def stage_selected_pose_inputs(
    pose_inputs: list[dict[str, Any]],
    staged_work_root: Path,
) -> dict[str, Any]:
    """Stage only the selected poses into a resumable fake work_root."""
    staged_cases: list[dict[str, Any]] = []
    latest_runs_by_target_model: dict[tuple[str, str], set[str]] = {}
    missing_confidence_json: list[str] = []

    for index, pose in enumerate(pose_inputs, start=1):
        ligand_id = str(pose["ligand_id"])
        model = str(pose["model"])
        protein_id = str(pose["protein_id"])
        run_id = str(pose.get("discovered_run_id") or pose.get("source_run_id") or "000001")
        model_dir = staged_work_root / ligand_id / model
        run_dir = model_dir / "runs" / run_id
        ut_dir = run_dir / f"{protein_id}_{ligand_id}"

        seed = pose.get("seed")
        sample = pose.get("sample")
        if seed is not None and sample is not None:
            sample_dir = ut_dir / f"seed-{int(seed)}_sample-{int(sample)}"
        else:
            sample_dir = ut_dir

        source_cif = Path(str(pose["cif_path"]))
        staged_cif = sample_dir / source_cif.name
        _replace_symlink(staged_cif, source_cif)

        confidence_json_path = pose.get("confidence_json_path")
        staged_confidence_json: str | None = None
        if confidence_json_path:
            source_confidence_json = Path(str(confidence_json_path))
            if source_confidence_json.exists():
                staged_confidence_path = sample_dir / source_confidence_json.name
                _replace_symlink(staged_confidence_path, source_confidence_json)
                staged_confidence_json = str(staged_confidence_path)
            else:
                missing_confidence_json.append(str(source_confidence_json))

        latest_runs_by_target_model.setdefault((ligand_id, model), set()).add(run_id)
        staged_cases.append(
            {
                "index": index,
                "pose_id": pose["pose_id"],
                "protein_id": protein_id,
                "ligand_id": ligand_id,
                "model": model,
                "run_id": run_id,
                "seed": seed,
                "sample": sample,
                "source_cif": str(source_cif),
                "staged_cif": str(staged_cif),
                "source_confidence_json": confidence_json_path,
                "staged_confidence_json": staged_confidence_json,
            }
        )

    latest_links: list[dict[str, str]] = []
    for (ligand_id, model), run_ids in sorted(latest_runs_by_target_model.items()):
        latest_run_id = _select_latest_run_id(run_ids)
        latest_link = staged_work_root / ligand_id / model / "latest"
        if latest_link.exists() or latest_link.is_symlink():
            latest_link.unlink()
        latest_link.parent.mkdir(parents=True, exist_ok=True)
        latest_link.symlink_to(Path("runs") / latest_run_id)
        latest_links.append(
            {
                "target": ligand_id,
                "model": model,
                "latest_run_id": latest_run_id,
                "latest_link": str(latest_link),
            }
        )

    return {
        "staged_work_root": str(staged_work_root),
        "n_staged_poses": len(staged_cases),
        "staged_cases": staged_cases,
        "latest_links": latest_links,
        "missing_confidence_json": sorted(set(missing_confidence_json)),
    }


def _build_production_config(
    manifest: PilotSelectionManifest,
    *,
    run_id: str,
    staged_work_root: Path,
) -> dict[str, Any]:
    clustering_pilot = {"enabled": True, **manifest.clustering_pilot}
    clustering_pilot.setdefault("label", run_id)
    selected_targets = sorted({selection.target for selection in manifest.selections if selection.target})

    return {
        "pipeline_version": "2.1",
        "production": {
            "run_id": run_id,
            "work_root": str(staged_work_root),
            "af3_only": manifest.af3_only,
            "latest_only": manifest.latest_only,
            "construct_type": manifest.construct_type,
            "include_targets": selected_targets,
            "include_proteins": sorted({selection.protein_id for selection in manifest.selections}),
            "run_posebusters": manifest.run_posebusters,
            "run_privateer": manifest.run_privateer,
            "clustering_pilot": clustering_pilot,
        },
    }


def write_production_config(
    manifest: PilotSelectionManifest,
    *,
    run_dir: Path,
    run_id: str,
    staged_work_root: Path,
    force: bool,
) -> tuple[Path, dict[str, Any], bool, str]:
    """Write the generated production config for the staged pilot run."""
    config_path = run_dir / "production.clustering_pilot.real_case.yaml"
    config_payload = _build_production_config(manifest, run_id=run_id, staged_work_root=staged_work_root)

    previous_payload = None
    if config_path.exists():
        previous_payload = yaml.safe_load(config_path.read_text()) or {}

    changed = force or previous_payload != config_payload
    if changed:
        config_path.write_text(yaml.safe_dump(config_payload, sort_keys=False))

    if not config_path.exists():
        config_path.write_text(yaml.safe_dump(config_payload, sort_keys=False))
        return config_path, config_payload, True, "written"
    if previous_payload is None:
        return config_path, config_payload, True, "written"
    if changed and force and previous_payload == config_payload:
        return config_path, config_payload, True, "rewritten"
    if changed:
        return config_path, config_payload, True, "updated"
    return config_path, config_payload, False, "verified"


def build_execute_command(
    *,
    python_executable: Path,
    script_path: Path,
    run_dir: Path,
    manifest_path: Path,
    del_branch: str,
    n_jobs: int = 1,
) -> str:
    """Build the exact command needed to execute the prepared pilot run."""
    parts = [
        shlex.quote(str(python_executable)),
        shlex.quote(str(script_path)),
        "--run-dir",
        shlex.quote(str(run_dir)),
        "--selection-manifest",
        shlex.quote(str(manifest_path)),
        "--del-branch",
        shlex.quote(del_branch),
        "--execute",
        "--n-jobs",
        shlex.quote(str(max(1, int(n_jobs)))),
    ]
    return " ".join(parts)


def prepare_real_case_pilot_run(
    *,
    run_dir: Path,
    manifest_path: Path,
    run_id: str,
    del_branch: str,
    python_executable: Path,
    script_path: Path,
    force_steps: Collection[str] = (),
    n_jobs: int = 1,
) -> dict[str, Any]:
    """Prepare a real-case clustering pilot run with resumable checkpoints."""
    checkpoints_dir = run_dir / "tmp" / "clustering_pilot_real_case"
    checkpoints_dir.mkdir(parents=True, exist_ok=True)

    manifest = load_selection_manifest(manifest_path)
    manifest_snapshot = serialize_selection_manifest(manifest)
    manifest_checkpoint_path = checkpoints_dir / "01_selection_manifest.snapshot.json"
    previous_manifest_snapshot = _read_json(manifest_checkpoint_path) if manifest_checkpoint_path.exists() else None
    manifest_changed = (
        previous_manifest_snapshot is not None
        and previous_manifest_snapshot != manifest_snapshot
    )
    _write_json(manifest_checkpoint_path, manifest_snapshot)

    effective_force_steps = set(force_steps)
    if manifest_changed:
        effective_force_steps.update({"discovery", "staging", "config"})

    step_statuses: dict[str, str] = {
        "manifest": "written"
        if previous_manifest_snapshot is None
        else ("updated" if manifest_changed else "verified")
    }

    discovery_checkpoint_path = checkpoints_dir / "02_discovery.json"
    discovery_payload, discovery_reused = _load_or_build_json_checkpoint(
        discovery_checkpoint_path,
        force="discovery" in effective_force_steps,
        required_paths=(),
        build=lambda: discover_selected_pose_inputs(manifest),
    )
    step_statuses["discovery"] = "reused" if discovery_reused else "written"

    staged_work_root = run_dir / "staged_work_root"
    staging_checkpoint_path = checkpoints_dir / "03_staging.json"
    staging_payload, staging_reused = _load_or_build_json_checkpoint(
        staging_checkpoint_path,
        force="staging" in effective_force_steps,
        required_paths=(staged_work_root,),
        build=lambda: stage_selected_pose_inputs(discovery_payload.get("pose_inputs", []), staged_work_root),
    )
    step_statuses["staging"] = "reused" if staging_reused else "written"

    config_path, config_payload, config_changed, config_status = write_production_config(
        manifest,
        run_dir=run_dir,
        run_id=run_id,
        staged_work_root=Path(str(staging_payload["staged_work_root"])),
        force="config" in effective_force_steps,
    )
    step_statuses["config"] = config_status

    config_checkpoint_path = checkpoints_dir / "04_config.json"
    _write_json(
        config_checkpoint_path,
        {
            "config_path": str(config_path),
            "config_payload": config_payload,
        },
    )

    ready_to_run = bool(manifest.selections) and not discovery_payload.get("missing_selections") and bool(
        discovery_payload.get("n_pose_inputs", 0)
    )
    discovery_summary = {
        "work_root": discovery_payload.get("work_root"),
        "construct_type": discovery_payload.get("construct_type"),
        "latest_only": discovery_payload.get("latest_only"),
        "n_selected_pairs": discovery_payload.get("n_selected_pairs", 0),
        "n_pose_inputs": discovery_payload.get("n_pose_inputs", 0),
        "n_selection_entries": len(discovery_payload.get("selection_entries", [])),
        "n_target_discovery_summaries": len(discovery_payload.get("target_discovery_summaries", [])),
        "n_missing_selections": len(discovery_payload.get("missing_selections", [])),
        "n_discovery_errors": len(discovery_payload.get("discovery_errors", [])),
        "selection_entries": discovery_payload.get("selection_entries", []),
        "missing_selections": discovery_payload.get("missing_selections", []),
        "discovery_errors": discovery_payload.get("discovery_errors", []),
    }
    staging_summary = {
        "staged_work_root": staging_payload.get("staged_work_root"),
        "n_staged_poses": staging_payload.get("n_staged_poses", 0),
        "n_latest_links": len(staging_payload.get("latest_links", [])),
        "n_missing_confidence_json": len(staging_payload.get("missing_confidence_json", [])),
        "missing_confidence_json": staging_payload.get("missing_confidence_json", []),
    }

    return {
        "run_id": run_id,
        "del_branch": del_branch,
        "selection_manifest_path": str(manifest_path),
        "checkpoints_dir": str(checkpoints_dir),
        "checkpoints": {
            "manifest": str(manifest_checkpoint_path),
            "discovery": str(discovery_checkpoint_path),
            "staging": str(staging_checkpoint_path),
            "config": str(config_checkpoint_path),
        },
        "manifest": manifest_snapshot,
        "manifest_changed": manifest_changed,
        "config_changed": config_changed,
        "step_statuses": step_statuses,
        "n_selected_pairs": len(manifest.selections),
        "selected_pairs": manifest_snapshot["selections"],
        "n_discovered_poses": discovery_payload.get("n_pose_inputs", 0),
        "discovery": discovery_summary,
        "discovery_summary": discovery_summary,
        "staging": staging_summary,
        "staging_summary": staging_summary,
        "config_path": str(config_path),
        "production_output": str(run_dir / "production_output"),
        "ready_to_run": ready_to_run,
        "next_command": build_execute_command(
            python_executable=python_executable,
            script_path=script_path,
            run_dir=run_dir,
            manifest_path=manifest_path,
            del_branch=del_branch,
            n_jobs=n_jobs,
        ),
    }


def summarize_pilot_execution(
    *,
    run_dir: Path,
    production_output: Path,
    force: bool = False,
) -> tuple[dict[str, Any], bool]:
    """Summarize an executed pilot run from production output files."""
    checkpoints_dir = run_dir / "tmp" / "clustering_pilot_real_case"
    summary_checkpoint_path = checkpoints_dir / "05_execution_summary.json"
    analysis_summary_path = production_output / "analysis_core_summary.json"

    should_rebuild = force or not summary_checkpoint_path.exists()
    if analysis_summary_path.exists() and summary_checkpoint_path.exists():
        should_rebuild = should_rebuild or (
            analysis_summary_path.stat().st_mtime > summary_checkpoint_path.stat().st_mtime
        )
    if summary_checkpoint_path.exists() and not should_rebuild:
        return _read_json(summary_checkpoint_path), True

    payload: dict[str, Any] = {
        "analysis_core_summary_path": str(analysis_summary_path),
        "pilot_outputs_ready": False,
    }

    if analysis_summary_path.exists():
        analysis_summary = json.loads(analysis_summary_path.read_text())
        pilot_summary = analysis_summary.get("clustering_pilot") or {}
        expected_paths = [
            Path(str(pilot_summary[path_key]))
            for path_key in (
                "pilot_ifp_interaction_type_prevalence_tsv",
                "pilot_ifp_feature_prevalence_tsv",
                "pilot_feature_selection_json",
                "pilot_clustering_method_summary_tsv",
            )
            if pilot_summary.get(path_key)
        ]
        payload.update(
            {
                "analysis_core_n_discovered": analysis_summary.get("n_discovered", 0),
                "analysis_core_n_prepared": analysis_summary.get("n_prepared", 0),
                "analysis_core_n_analyzed": analysis_summary.get("n_analyzed", 0),
                "pilot_output_dir": pilot_summary.get("output_dir", ""),
                "pilot_method_summary_tsv": pilot_summary.get(
                    "pilot_clustering_method_summary_tsv", ""
                ),
                "pilot_n_conditions": pilot_summary.get("n_conditions", 0),
                "pilot_n_conditions_formal_clustering_allowed": pilot_summary.get(
                    "n_conditions_formal_clustering_allowed", 0
                ),
                "pilot_n_conditions_insufficient_clusterable_signal": pilot_summary.get(
                    "n_conditions_insufficient_clusterable_signal", 0
                ),
                "pilot_condition_matrix_entries": pilot_summary.get("condition_matrix_entries", []),
                "pilot_outputs_ready": bool(pilot_summary) and _paths_exist(expected_paths),
            }
        )

    _write_json(summary_checkpoint_path, payload)
    return payload, False
