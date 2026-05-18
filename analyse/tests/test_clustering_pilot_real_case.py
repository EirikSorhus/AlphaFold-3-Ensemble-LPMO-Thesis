from __future__ import annotations

import json
from pathlib import Path

import yaml

from lpmo_pipeline.analysis.clustering_pilot_real_case import (
    build_selection_manifest_from_overview,
    load_pilot_protein_overview_tsv,
    load_selection_manifest,
    prepare_real_case_pilot_run,
    serialize_selection_manifest,
    stage_selected_pose_inputs,
)


def _touch(path: Path, content: str = "test\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return path


def _make_work_root(tmp_path: Path) -> Path:
    work_root = tmp_path / "work_core"
    ut_dir = work_root / "NAG4" / "af3" / "runs" / "100001" / "Q7SCE9_NAG4"
    _touch(ut_dir / "seed-1_sample-0" / "Q7SCE9_NAG4_seed-1_sample-0_model.cif")
    _touch(ut_dir / "seed-1_sample-1" / "Q7SCE9_NAG4_seed-1_sample-1_model.cif")
    sta_dir = work_root / "STA6" / "af3" / "runs" / "100002" / "Q7SCE9_STA6"
    _touch(sta_dir / "seed-2_sample-0" / "Q7SCE9_STA6_seed-2_sample-0_model.cif")
    return work_root


def _write_manifest(manifest_path: Path, work_root: Path) -> None:
    manifest_path.write_text(
        yaml.safe_dump(
            {
                "work_root": str(work_root),
                "construct_type": "domain_only",
                "latest_only": False,
                "clustering_pilot": {
                    "label": "unit-test-pilot",
                    "minimum_clusterable_n": 10,
                },
                "selections": [
                    {"target": "NAG4", "protein_id": "Q7SCE9"},
                ],
            },
            sort_keys=False,
        )
    )


def test_load_selection_manifest_serializes_expected_fields(tmp_path: Path) -> None:
    work_root = _make_work_root(tmp_path)
    manifest_path = tmp_path / "selection.yaml"
    _write_manifest(manifest_path, work_root)

    manifest = load_selection_manifest(manifest_path)
    snapshot = serialize_selection_manifest(manifest)

    assert manifest.work_root == work_root.resolve()
    assert manifest.latest_only is False
    assert snapshot["selections"] == [{"target": "NAG4", "protein_id": "Q7SCE9"}]
    assert snapshot["clustering_pilot"]["label"] == "unit-test-pilot"


def test_load_pilot_protein_overview_tsv_parses_run_modes(tmp_path: Path) -> None:
    overview_path = tmp_path / "overview.tsv"
    overview_path.write_text(
        "\t".join(
            [
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
                "run_domain_only",
                "run_full_length",
                "rationale",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "1",
                "Q7SCE9",
                "AA9",
                "monooxygenase A",
                "Panus similis",
                "1.14.99.54; 1.14.99.56",
                "C1+C4",
                "cellulose",
                "xylo_oligosaccharide",
                "xylohexaose",
                "Nei",
                "None",
                "Kun katalytisk domene",
                "true",
                "false",
                "coverage",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "2",
                "B2ADG1",
                "AA9",
                "monooxygenase H",
                "Podospora anserina",
                "1.14.99.54; 1.14.99.56",
                "C1+C4",
                "cellulose",
                "cello_oligosaccharide",
                "cellohexaose",
                "Ja",
                "PFAM:PF00734",
                "Katalytisk + full-lengde",
                "true",
                "true",
                "full-length representative",
            ]
        )
        + "\n"
    )

    rows = load_pilot_protein_overview_tsv(overview_path)

    assert [row.uniprot_id for row in rows] == ["Q7SCE9", "B2ADG1"]
    assert rows[0].run_domain_only is True
    assert rows[0].run_full_length is False
    assert rows[1].run_domain_only is True
    assert rows[1].run_full_length is True


def test_build_selection_manifest_from_overview_filters_full_length_subset(tmp_path: Path) -> None:
    overview_path = tmp_path / "overview.tsv"
    overview_path.write_text(
        "\t".join(
            [
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
                "run_domain_only",
                "run_full_length",
                "rationale",
            ]
        )
        + "\n"
        + "1\tQ7SCE9\tAA9\tA\tOrg\tEC\tC1\tcellulose\t\t\tNei\tNone\tKun katalytisk domene\ttrue\tfalse\tr1\n"
        + "2\tB2ADG1\tAA9\tB\tOrg\tEC\tC1\tcellulose\t\t\tJa\tPFAM\tKatalytisk + full-lengde\ttrue\ttrue\tr2\n"
    )
    rows = load_pilot_protein_overview_tsv(overview_path)

    domain_manifest = build_selection_manifest_from_overview(
        rows,
        construct_type="domain_only",
        work_root=tmp_path / "work_core",
    )
    full_manifest = build_selection_manifest_from_overview(
        rows,
        construct_type="full_length",
        work_root=tmp_path / "work_full_length",
    )

    assert [selection.protein_id for selection in domain_manifest.selections] == ["Q7SCE9", "B2ADG1"]
    assert all(selection.target is None for selection in domain_manifest.selections)
    assert [selection.protein_id for selection in full_manifest.selections] == ["B2ADG1"]
    assert full_manifest.construct_type == "full_length"


def test_stage_selected_pose_inputs_symlinks_confidence_json_and_latest_link(tmp_path: Path) -> None:
    source_dir = tmp_path / "source"
    source_cif = _touch(
        source_dir / "NAG4" / "af3" / "runs" / "408002" / "Q7SCE9_NAG4" / "seed-1_sample-0" / "Q7SCE9_NAG4_seed-1_sample-0_model.cif"
    )
    source_confidence = _touch(
        source_dir / "NAG4" / "af3" / "runs" / "408002" / "Q7SCE9_NAG4" / "seed-1_sample-0" / "Q7SCE9_NAG4_seed-1_sample-0_model.json",
        content="{}\n",
    )

    payload = [
        {
            "pose_id": "Q7SCE9_NAG4_seed-1_sample-0_model",
            "protein_id": "Q7SCE9",
            "ligand_id": "NAG4",
            "model": "af3",
            "source_run_id": "408002",
            "discovered_run_id": "408002",
            "confidence_json_path": str(source_confidence),
            "run_status": "succeeded",
            "seed": 1,
            "sample": 0,
            "cif_path": str(source_cif),
        }
    ]

    staging = stage_selected_pose_inputs(payload, tmp_path / "staged_work_root")

    assert staging["n_staged_poses"] == 1
    staged_case = staging["staged_cases"][0]
    assert Path(staged_case["staged_cif"]).is_symlink()
    assert Path(staged_case["staged_confidence_json"]).is_symlink()
    latest_link = tmp_path / "staged_work_root" / "NAG4" / "af3" / "latest"
    assert latest_link.is_symlink()
    assert latest_link.resolve() == (
        tmp_path / "staged_work_root" / "NAG4" / "af3" / "runs" / "408002"
    ).resolve()


def test_prepare_real_case_pilot_run_writes_checkpoints_and_config(tmp_path: Path) -> None:
    work_root = _make_work_root(tmp_path)
    manifest_path = tmp_path / "selection.yaml"
    _write_manifest(manifest_path, work_root)

    prepared = prepare_real_case_pilot_run(
        run_dir=tmp_path / "pilot_run",
        manifest_path=manifest_path,
        run_id="pilot-run",
        del_branch="del_a",
        python_executable=Path("/usr/bin/python3"),
        script_path=Path("tests/run_tests_scripts/run_clustering_pilot_real_case.py"),
    )

    assert prepared["ready_to_run"] is True
    assert prepared["n_selected_pairs"] == 1
    assert prepared["n_discovered_poses"] == 2
    assert Path(str(prepared["config_path"])).exists()
    assert Path(str(prepared["checkpoints"]["discovery"])).exists()
    assert Path(str(prepared["checkpoints"]["staging"])).exists()
    assert prepared["step_statuses"]["discovery"] == "written"
    assert prepared["step_statuses"]["staging"] == "written"
    assert "--execute" in str(prepared["next_command"])

    config = yaml.safe_load(Path(str(prepared["config_path"])).read_text())
    assert config["production"]["clustering_pilot"]["enabled"] is True
    assert config["production"]["include_targets"] == ["NAG4"]
    assert config["production"]["include_proteins"] == ["Q7SCE9"]


def test_prepare_real_case_pilot_run_supports_protein_only_selection(tmp_path: Path) -> None:
    work_root = _make_work_root(tmp_path)
    manifest_path = tmp_path / "selection.yaml"
    manifest_path.write_text(
        yaml.safe_dump(
            {
                "work_root": str(work_root),
                "construct_type": "domain_only",
                "latest_only": False,
                "clustering_pilot": {"label": "unit-test-pilot"},
                "selections": [
                    {"protein_id": "Q7SCE9"},
                ],
            },
            sort_keys=False,
        )
    )

    prepared = prepare_real_case_pilot_run(
        run_dir=tmp_path / "pilot_run",
        manifest_path=manifest_path,
        run_id="pilot-run",
        del_branch="del_a",
        python_executable=Path("/usr/bin/python3"),
        script_path=Path("tests/run_tests_scripts/run_clustering_pilot_real_case.py"),
    )

    assert prepared["ready_to_run"] is True
    assert prepared["n_discovered_poses"] == 3
    assert prepared["discovery"]["selection_entries"][0]["discovered_targets"] == ["NAG4", "STA6"]
    config = yaml.safe_load(Path(str(prepared["config_path"])).read_text())
    assert config["production"]["include_targets"] == []
    assert config["production"]["include_proteins"] == ["Q7SCE9"]


def test_prepare_real_case_pilot_run_reuses_discovery_checkpoint(tmp_path: Path) -> None:
    work_root = _make_work_root(tmp_path)
    manifest_path = tmp_path / "selection.yaml"
    _write_manifest(manifest_path, work_root)

    first = prepare_real_case_pilot_run(
        run_dir=tmp_path / "pilot_run",
        manifest_path=manifest_path,
        run_id="pilot-run",
        del_branch="del_a",
        python_executable=Path("/usr/bin/python3"),
        script_path=Path("tests/run_tests_scripts/run_clustering_pilot_real_case.py"),
    )
    discovery_checkpoint = Path(str(first["checkpoints"]["discovery"]))
    discovery_payload = json.loads(discovery_checkpoint.read_text())
    discovery_payload["checkpoint_marker"] = "keep-me"
    discovery_checkpoint.write_text(json.dumps(discovery_payload, indent=2))

    second = prepare_real_case_pilot_run(
        run_dir=tmp_path / "pilot_run",
        manifest_path=manifest_path,
        run_id="pilot-run",
        del_branch="del_a",
        python_executable=Path("/usr/bin/python3"),
        script_path=Path("tests/run_tests_scripts/run_clustering_pilot_real_case.py"),
    )

    reused_payload = json.loads(discovery_checkpoint.read_text())
    assert second["step_statuses"]["discovery"] == "reused"
    assert reused_payload["checkpoint_marker"] == "keep-me"