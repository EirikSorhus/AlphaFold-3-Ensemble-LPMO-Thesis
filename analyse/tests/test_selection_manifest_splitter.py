from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import sys

import yaml


SCRIPT_PATH = (
    Path(__file__).parent / "run_tests_scripts" / "split_clustering_pilot_selection_manifest.py"
)
SRC_PATH = Path(__file__).resolve().parents[1] / "src"


def _load_module():
    if str(SRC_PATH) not in sys.path:
        sys.path.insert(0, str(SRC_PATH))
    spec = importlib.util.spec_from_file_location("split_selection_manifest", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _touch(path: Path, content: str = "test\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return path


def _write_manifest(manifest_path: Path, work_root: Path) -> None:
    manifest_path.write_text(
        yaml.safe_dump(
            {
                "work_root": str(work_root),
                "construct_type": "domain_only",
                "latest_only": False,
                "selections": [
                    {"protein_id": "Q7SCE9"},
                    {"protein_id": "B2ADG1"},
                    {"protein_id": "F4H6A3"},
                ],
            },
            sort_keys=False,
        )
    )


def _make_work_root(tmp_path: Path) -> Path:
    work_root = tmp_path / "work_core"
    _touch(
        work_root / "NAG4" / "af3" / "runs" / "100001" / "Q7SCE9_NAG4" / "seed-1_sample-0" / "Q7SCE9_NAG4_seed-1_sample-0_model.cif"
    )
    _touch(
        work_root / "NAG4" / "af3" / "runs" / "100001" / "Q7SCE9_NAG4" / "seed-1_sample-1" / "Q7SCE9_NAG4_seed-1_sample-1_model.cif"
    )
    _touch(
        work_root / "STA6" / "af3" / "runs" / "100002" / "Q7SCE9_STA6" / "seed-2_sample-0" / "Q7SCE9_STA6_seed-2_sample-0_model.cif"
    )
    _touch(
        work_root / "NAG4" / "af3" / "runs" / "100003" / "B2ADG1_NAG4" / "seed-1_sample-0" / "B2ADG1_NAG4_seed-1_sample-0_model.cif"
    )
    _touch(
        work_root / "NAG4" / "af3" / "runs" / "100004" / "F4H6A3_NAG4" / "seed-1_sample-0" / "F4H6A3_NAG4_seed-1_sample-0_model.cif"
    )
    return work_root


def test_build_shards_balances_by_estimated_pose_count(tmp_path: Path) -> None:
    module = _load_module()
    work_root = _make_work_root(tmp_path)
    manifest_path = tmp_path / "selection.yaml"
    _write_manifest(manifest_path, work_root)
    manifest = module.load_selection_manifest(manifest_path)

    shards = module._build_shards(
        manifest,
        proteins_per_shard=2,
        balance_by="estimated_pose_count",
    )

    assert len(shards) == 2
    shard_weights = sorted(int(shard["estimated_pose_count"]) for shard in shards)
    assert shard_weights == [2, 3]
    shard_proteins = [
        [selection.protein_id for selection in shard["selections"]]
        for shard in shards
    ]
    assert sorted(shard_proteins) == [["B2ADG1", "F4H6A3"], ["Q7SCE9"]]


def test_main_writes_shard_index_with_output_roots(tmp_path: Path) -> None:
    module = _load_module()
    work_root = _make_work_root(tmp_path)
    manifest_path = tmp_path / "selection.yaml"
    output_dir = tmp_path / "shards"
    _write_manifest(manifest_path, work_root)

    argv = sys.argv[:]
    sys.argv = [
        str(SCRIPT_PATH),
        "--input-manifest",
        str(manifest_path),
        "--output-dir",
        str(output_dir),
        "--proteins-per-shard",
        "2",
        "--balance-by",
        "estimated_pose_count",
        "--output-root-template",
        str(tmp_path / "run_root" / "{construct_type}_shards" / "{shard_id}" / "production_output"),
    ]
    try:
        assert module.main() == 0
    finally:
        sys.argv = argv

    summary = json.loads((output_dir / "shard_summary.json").read_text())
    assert summary["balance_by"] == "estimated_pose_count"
    assert summary["proteins_per_shard"] == 2
    assert summary["total_estimated_pose_count"] == 5

    with (output_dir / "shard_index.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert [row["estimated_pose_count"] for row in rows] == ["3", "2"]
    assert rows[0]["output_root"].endswith("domain_only_shards/shard_0001/production_output")
    assert rows[1]["output_root"].endswith("domain_only_shards/shard_0002/production_output")