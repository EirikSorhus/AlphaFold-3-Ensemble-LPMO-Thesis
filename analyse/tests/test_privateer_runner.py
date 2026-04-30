from __future__ import annotations

from pathlib import Path

import pytest

from lpmo_pipeline.qc.privateer_runner import (
    PRIVATEER_DEFAULT_MODE,
    PrivateerBatchInput,
    _aggregate_privateer_results,
    _as_bool,
    _parse_privateer_summary,
    _parse_validation_rows,
    _privateer_subprocess_env,
    build_privateer_invocation,
    get_privateer_version,
    prepare_privateer_input,
    run_privateer,
    run_privateer_batch,
)
from lpmo_pipeline.utils.exceptions import PrivateerCCDError


def _write_minimal_privateer_source(path: Path, residue_name: str = "NAG") -> Path:
    path.write_text(
        "\n".join(
            [
                "data_test",
                "#",
                "loop_",
                "_atom_site.group_PDB",
                "_atom_site.id",
                "_atom_site.type_symbol",
                "_atom_site.label_atom_id",
                "_atom_site.label_alt_id",
                "_atom_site.label_comp_id",
                "_atom_site.label_asym_id",
                "_atom_site.label_entity_id",
                "_atom_site.label_seq_id",
                "_atom_site.pdbx_PDB_ins_code",
                "_atom_site.Cartn_x",
                "_atom_site.Cartn_y",
                "_atom_site.Cartn_z",
                "_atom_site.occupancy",
                "_atom_site.B_iso_or_equiv",
                "_atom_site.auth_seq_id",
                "_atom_site.auth_asym_id",
                "_atom_site.pdbx_PDB_model_num",
                f"HETATM 1 C C1 . {residue_name} C 1 1 ? 0.000 0.000 0.000 1.00 20.00 1 C 1",
                f"HETATM 2 O O5 . {residue_name} C 1 1 ? 1.200 0.000 0.000 1.00 20.00 1 C 1",
                "#",
                "",
            ]
        )
    )
    return path


def test_prepare_privateer_input_writes_validated_mmcif(tmp_path: Path) -> None:
    source = _write_minimal_privateer_source(tmp_path / "input.cif")

    output = prepare_privateer_input(source, tmp_path / "privateer_input.cif")

    assert output.exists()
    assert output.name == "privateer_input.cif"
    assert "data_" in output.read_text()


def test_prepare_privateer_input_rejects_invalid_comp_id(tmp_path: Path) -> None:
    source = _write_minimal_privateer_source(tmp_path / "bad_input.cif", residue_name="XXX")

    with pytest.raises(PrivateerCCDError, match="non-CCD glycan residues"):
        prepare_privateer_input(source, tmp_path / "privateer_input.cif")


def test_parse_privateer_summary_counts() -> None:
    summary = _parse_privateer_summary(
        """
        Wrong anomer: 1
        Wrong configuration: 0
        Unphysical puckering amplitude: 2
        In higher-energy conformations: 3
        Privateer has identified 4 issues, with 2 of 5 sugars affected.
        """
    )

    assert summary.wrong_anomer == 1
    assert summary.unphysical_puckering == 2
    assert summary.high_energy_conformations == 3
    assert summary.issues_detected == 4
    assert summary.affected_sugars == 2
    assert summary.total_sugars_reported == 5


def test_parse_validation_rows() -> None:
    rows = _parse_validation_rows(
        "Q7SC\tNAG-C-1\t0.571\t308.68\t5.97\tbeta-D-aldopyranose\t4c1\t69.10\t(l)\tyes\n"
        "Q7SC\tNAG-C-2\t0.579\t353.75\t4.10\tbeta-D-aldopyranose\t4c1\t76.23\t(l)\tno\n"
    )

    assert len(rows) == 2
    assert rows[0]["sugar_label"] == "NAG-C-1"
    assert rows[1]["ok"] == "no"


def test_build_privateer_invocation_includes_binds_and_optional_inputs(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    from lpmo_pipeline.qc import privateer_runner as mod

    fake_sif = tmp_path / "privateer.sif"
    fake_sif.write_text("fake")
    monkeypatch.setattr(mod, "PRIVATEER_SIF_CANDIDATES", (fake_sif,))

    input_dir = tmp_path / "input dir"
    output_dir = tmp_path / "out dir"
    sf_dir = tmp_path / "density dir"
    input_dir.mkdir()
    output_dir.mkdir()
    sf_dir.mkdir()
    cif_path = _write_minimal_privateer_source(input_dir / "input odd.cif")
    sf_path = sf_dir / "density odd.mtz"
    sf_path.write_text("mtz")

    invocation = build_privateer_invocation(
        cif_path=cif_path,
        output_dir=output_dir,
        structure_factor_mtz=sf_path,
        colin_fo="F,SIGF",
        codein="NAG",
        showgeom=True,
        radiusin=2.5,
    )

    assert invocation.command[:4] == ["apptainer", "run", "--cleanenv", "--bind"]
    assert "-pdbin" in invocation.command
    assert str(cif_path.resolve()) in invocation.command
    assert "-mtzin" in invocation.command
    assert str(sf_path.resolve()) in invocation.command
    assert "-colin-fo" in invocation.command
    assert "F,SIGF" in invocation.command
    assert "-codein" in invocation.command
    assert "NAG" in invocation.command
    assert "-showgeom" in invocation.command
    assert "-radiusin" in invocation.command
    assert str(fake_sif) in invocation.command
    assert any(str(input_dir.resolve()) in bind for bind in invocation.bind_mounts)
    assert any(str(output_dir.resolve()) in bind for bind in invocation.bind_mounts)
    assert any(str(sf_dir.resolve()) in bind for bind in invocation.bind_mounts)


def test_run_privateer_returns_failed_result_on_runner_error(tmp_path: Path) -> None:
    result = run_privateer(tmp_path / "missing.cif", pose_id="pose_1", probe_container=False)

    assert result.all_pass is False
    assert result.total_sugars == 0
    assert result.recognition_rate == 0.0
    assert "does not exist" in result.runner_error


def test_run_privateer_dry_run_returns_invocation_metadata(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    from lpmo_pipeline.qc import privateer_runner as mod

    fake_sif = tmp_path / "privateer.sif"
    fake_sif.write_text("fake")
    monkeypatch.setattr(mod, "PRIVATEER_SIF_CANDIDATES", (fake_sif,))

    prepared = prepare_privateer_input(_write_minimal_privateer_source(tmp_path / "input.cif"))
    result = run_privateer(
        prepared,
        pose_id="pose dry run",
        output_dir=tmp_path / "dry output",
        structure_factor_cif=tmp_path / "factors with spaces.cif",
        dry_run=True,
    )

    assert result.dry_run is True
    assert result.command[:3] == ["apptainer", "run", "--cleanenv"]
    assert "-pdbin" in result.command
    assert "-cifin" in result.command
    assert "-mode" in result.command
    assert PRIVATEER_DEFAULT_MODE in result.command
    assert result.work_dir.endswith("dry output")


def test_privateer_subprocess_env_is_minimal(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("HOME", "/home/tester")
    monkeypatch.setenv("USER", "tester")
    monkeypatch.setenv("LOGNAME", "tester")
    monkeypatch.setenv("LANG", "en_US.UTF-8")
    monkeypatch.setenv("TERM", "screen")
    monkeypatch.setenv("CONDA_PREFIX", "/tmp/conda")

    env = _privateer_subprocess_env()

    assert env == {
        "PATH": "/usr/bin:/bin:/usr/local/bin",
        "HOME": "/home/tester",
        "USER": "tester",
        "LOGNAME": "tester",
        "LANG": "en_US.UTF-8",
        "LC_ALL": "C.UTF-8",
        "TERM": "screen",
    }


def test_run_privateer_parses_ccp4i2_outputs(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    from lpmo_pipeline.qc import privateer_runner as mod

    fake_sif = tmp_path / "privateer.sif"
    fake_sif.write_text("fake")
    monkeypatch.setattr(mod, "PRIVATEER_SIF_CANDIDATES", (fake_sif,))
    monkeypatch.setattr(mod, "_probe_privateer_container", lambda *args, **kwargs: None)

    prepared = prepare_privateer_input(_write_minimal_privateer_source(tmp_path / "input.cif"))
    output_dir = tmp_path / "privateer_output"

    class Proc:
        def __init__(self, stdout: str, stderr: str = "", returncode: int = 0) -> None:
            self.stdout = stdout
            self.stderr = stderr
            self.returncode = returncode

    def _fake_run(cmd, capture_output, text, check, timeout, cwd, env):
        validation_path = Path(cwd) / "validation_data-privateer"
        validation_path.parent.mkdir(parents=True, exist_ok=True)
        validation_path.write_text(
            "Q7SC\tNAG-C-1\t0.571\t308.68\t5.97\tbeta-D-aldopyranose\t4c1\t69.10\t(l)\tyes\n"
        )
        return Proc(
            """
            Wrong anomer: 0
            Wrong configuration: 0
            Unphysical puckering amplitude: 0
            In higher-energy conformations: 0
            Privateer has identified 0 issues, with 0 of 1 sugars affected.
            """
        )

    monkeypatch.setattr(mod.subprocess, "run", _fake_run)

    result = run_privateer(prepared, pose_id="pose_parse", output_dir=output_dir, probe_container=False)

    assert result.pose_id == "pose_parse"
    assert result.returncode == 0
    assert result.total_sugars == 1
    assert result.recognized == 1
    assert result.recognition_rate == 1.0
    assert result.anomer_pass == 1
    assert result.all_pass is True
    assert result.validation_data_path.endswith("validation_data-privateer")


def test_run_privateer_batch_preserves_input_order(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    from lpmo_pipeline.qc import privateer_runner as mod

    calls: list[str] = []

    def _fake_run_privateer(cif_path, pose_id="", **kwargs):
        calls.append(pose_id)
        return mod.PrivateerResult(pose_id=pose_id, all_pass=True)

    monkeypatch.setattr(mod, "run_privateer", _fake_run_privateer)

    results = run_privateer_batch(
        [
            PrivateerBatchInput(cif_path=tmp_path / "a.cif", pose_id="pose_a", dry_run=True),
            PrivateerBatchInput(cif_path=tmp_path / "b.cif", pose_id="pose_b", dry_run=True),
        ],
        max_workers=2,
    )

    assert [result.pose_id for result in results] == ["pose_a", "pose_b"]
    assert sorted(calls) == ["pose_a", "pose_b"]


def test_get_privateer_version_reads_from_sif(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    from lpmo_pipeline.qc import privateer_runner as mod

    fake_sif = tmp_path / "privateer.sif"
    fake_sif.write_text("fake")
    monkeypatch.setattr(mod, "PRIVATEER_SIF_CANDIDATES", (fake_sif,))

    class Proc:
        stdout = "Privateer Version MKV. Copyright 2013-2024"
        stderr = ""
        returncode = 0

    monkeypatch.setattr(mod.subprocess, "run", lambda *args, **kwargs: Proc())

    assert get_privateer_version() == "MKV"


def test_as_bool() -> None:
    assert _as_bool(True) is True
    assert _as_bool("true") is True
    assert _as_bool("YES") is True
    assert _as_bool("1") is True
    assert _as_bool(1) is True
    assert _as_bool("0") is False
    assert _as_bool("no") is False
    assert _as_bool(None) is False


def test_aggregate_privateer_results_requires_all_categories_to_pass(tmp_path: Path) -> None:
    from lpmo_pipeline.qc import privateer_runner as mod

    invocation = build_privateer_invocation(
        cif_path=prepare_privateer_input(_write_minimal_privateer_source(tmp_path / "agg.cif")),
        output_dir=tmp_path / "agg_out",
    )
    summary = _parse_privateer_summary(
        "Privateer has identified 1 issues, with 1 of 1 sugars affected.\nWrong anomer: 0\nWrong configuration: 0\nUnphysical puckering amplitude: 1\nIn higher-energy conformations: 0\n"
    )
    rows = _parse_validation_rows(
        "Q7SC\tNAG-C-1\t0.571\t308.68\t5.97\tbeta-D-aldopyranose\t4c1\t69.10\t(l)\tno\n"
    )
    residues = [mod._row_to_residue_result(rows[0], summary)]
    result = _aggregate_privateer_results(
        pose_id="agg_pose",
        residue_results=residues,
        summary=summary,
        invocation=invocation,
        returncode=0,
    )

    assert result.ring_pucker_pass == 0
    assert result.all_pass is False
