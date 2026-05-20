from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import yaml


SCRIPT_PATH = Path(__file__).parent / "run_tests_scripts" / "run_clustering_pilot_timing_probe.py"


def _load_probe_module():
    spec = importlib.util.spec_from_file_location("run_clustering_pilot_timing_probe", SCRIPT_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _pose(module, construct_type: str, protein_id: str, ligand_id: str, index: int):
    return module.TimingPose(
        construct_type=construct_type,
        pose_id=f"{protein_id}_{ligand_id}_seed-1_sample-{index}_model",
        protein_id=protein_id,
        ligand_id=ligand_id,
        model="af3",
        source_run_id="100001",
        discovered_run_id="100001",
        cif_path=Path(f"/tmp/{protein_id}_{ligand_id}_{index}.cif"),
        confidence_json_path=None,
        seed=1,
        sample=index,
    )


def test_balanced_select_round_robins_conditions() -> None:
    module = _load_probe_module()
    candidates = [
        _pose(module, "domain_only", "P1", "NAG4", 0),
        _pose(module, "domain_only", "P1", "NAG4", 1),
        _pose(module, "domain_only", "P2", "STA6", 0),
        _pose(module, "domain_only", "P2", "STA6", 1),
    ]

    selected = module._balanced_select(candidates, 3, exclude=set())

    assert [pose.pose_id for pose in selected] == [
        "P1_NAG4_seed-1_sample-0_model",
        "P2_STA6_seed-1_sample-0_model",
        "P1_NAG4_seed-1_sample-1_model",
    ]


def test_write_production_config_records_probe_runtime_options(tmp_path: Path) -> None:
    module = _load_probe_module()
    poses = [
        _pose(module, "full_length", "P1", "NAG4", 0),
        _pose(module, "full_length", "P2", "STA6", 0),
    ]
    config_path = tmp_path / "production.timing_probe.yaml"

    module._write_production_config(
        config_path=config_path,
        run_id="timing-test",
        staged_work_root=tmp_path / "staged_work_root",
        construct_type="full_length",
        selected_poses=poses,
        run_posebusters=True,
        run_privateer=True,
        n_jobs=4,
    )

    payload = yaml.safe_load(config_path.read_text())
    production = payload["production"]
    assert production["run_id"] == "timing-test"
    assert production["construct_type"] == "full_length"
    assert production["run_posebusters"] is True
    assert production["run_privateer"] is True
    assert production["n_jobs"] == 4
    assert production["include_targets"] == ["NAG4", "STA6"]
    assert production["include_proteins"] == ["P1", "P2"]
    assert production["clustering_pilot"]["enabled"] is True


def test_summarize_events_groups_by_mode_construct_workers_and_step() -> None:
    module = _load_probe_module()
    events = [
        {
            "mode": "prepare_only",
            "construct_type": "domain_only",
            "n_jobs": 2,
            "step": "prepare.normalize",
            "wall_time_s": 1.5,
            "status": "ok",
        },
        {
            "mode": "prepare_only",
            "construct_type": "domain_only",
            "n_jobs": 2,
            "step": "prepare.normalize",
            "wall_time_s": 2.0,
            "status": "failed",
        },
    ]

    summary = module._summarize_events(events)

    assert summary == [
        {
            "mode": "prepare_only",
            "construct_type": "domain_only",
            "n_jobs": 2,
            "step": "prepare.normalize",
            "event_count": 2,
            "wall_time_s_sum": 3.5,
            "error_count": 1,
        }
    ]
