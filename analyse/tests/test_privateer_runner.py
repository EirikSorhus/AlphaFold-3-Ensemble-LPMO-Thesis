from __future__ import annotations

from pathlib import Path

import pytest

from lpmo_pipeline.qc.privateer_runner import (
    _aggregate_privateer_results,
    _parse_privateer_output,
  _privateer_subprocess_env,
    _as_bool,
    run_privateer,
)


def test_parse_privateer_output_top_level_list() -> None:
    raw = """
    [
      {
        "chain": "B",
        "resname": "NAG",
        "resnum": 1,
        "sugar_recognized": true,
        "anomer_ok": true,
        "ring_pucker_ok": true
      },
      {
        "chain": "B",
        "resname": "NAG",
        "resnum": 2,
        "sugar_recognized": false,
        "anomer_ok": false,
        "ring_pucker_ok": false
      }
    ]
    """

    residues = _parse_privateer_output(raw)
    assert len(residues) == 2
    assert residues[0].chain == "B"
    assert residues[0].resname == "NAG"
    assert residues[0].sugar_recognized is True
    assert residues[1].sugar_recognized is False


def test_parse_privateer_output_nested_results_key() -> None:
    raw = """
    {
      "results": [
        {
          "chain_id": "C",
          "residue_name": "BGC",
          "residue_number": 4,
          "is_recognized": "true",
          "anomer_correct": "yes",
          "puckering_ok": "1",
          "linkage_ok": "0"
        }
      ]
    }
    """

    residues = _parse_privateer_output(raw)
    assert len(residues) == 1
    r = residues[0]
    assert r.chain == "C"
    assert r.resname == "BGC"
    assert r.resnum == 4
    assert r.sugar_recognized is True
    assert r.anomer_ok is True
    assert r.ring_pucker_ok is True
    assert r.linkage_ok is False


def test_aggregate_privateer_results_gate_logic() -> None:
    residues = _parse_privateer_output(
        """
        [
          {"chain": "B", "resname": "NAG", "resnum": 1, "sugar_recognized": true, "anomer_ok": true, "ring_pucker_ok": true},
          {"chain": "B", "resname": "NAG", "resnum": 2, "sugar_recognized": true, "anomer_ok": false, "ring_pucker_ok": true}
        ]
        """
    )
    result = _aggregate_privateer_results(residues)

    assert result.total_sugars == 2
    assert result.recognized == 2
    assert result.recognition_rate == 1.0
    assert result.anomer_pass == 1
    assert result.all_pass is False


def test_as_bool() -> None:
    assert _as_bool(True) is True
    assert _as_bool("true") is True
    assert _as_bool("YES") is True
    assert _as_bool("1") is True
    assert _as_bool(1) is True
    assert _as_bool("0") is False
    assert _as_bool("no") is False
    assert _as_bool(None) is False


def test_run_privateer_returns_failed_result_on_runner_error(monkeypatch, tmp_path) -> None:
    from lpmo_pipeline.qc import privateer_runner as mod

    def _boom(_):
        raise RuntimeError("privateer unavailable")

    monkeypatch.setattr(mod, "_run_privateer_cli", _boom)

    cif_path = tmp_path / "input.cif"
    cif_path.write_text("data_test\n")
    result = run_privateer(cif_path, pose_id="pose_1")

    assert result.all_pass is False
    assert result.total_sugars == 0
    assert result.recognition_rate == 0.0
    assert result.runner_error == "privateer unavailable"


def test_run_privateer_cli_requires_sif(monkeypatch, tmp_path) -> None:
    from lpmo_pipeline.qc import privateer_runner as mod

    monkeypatch.setattr(mod, "PRIVATEER_SIF_CANDIDATES", (Path("/does/not/exist/privateer.sif"),))

    cif_path = tmp_path / "input.cif"
    cif_path.write_text("data_test\n")

    with pytest.raises(RuntimeError, match="Privateer SIF not found"):
      mod._run_privateer_cli(cif_path)


def test_run_privateer_cli_uses_apptainer_run(monkeypatch, tmp_path) -> None:
    from lpmo_pipeline.qc import privateer_runner as mod

    class Proc:
        def __init__(self) -> None:
            self.stdout = "[]"

    calls: list[list[str]] = []

    def _fake_run(cmd, capture_output, text, check, timeout, env):
      calls.append(cmd)
      return Proc()

    fake_sif = tmp_path / "privateer.sif"
    fake_sif.write_text("fake")
    cif_path = tmp_path / "input.cif"
    cif_path.write_text("data_test\n")

    monkeypatch.setattr(mod, "PRIVATEER_SIF_CANDIDATES", (fake_sif,))
    monkeypatch.setattr(mod.subprocess, "run", _fake_run)

    raw = mod._run_privateer_cli(cif_path)

    assert raw == "[]"
    assert len(calls) == 1
    assert calls[0][:4] == ["apptainer", "run", "--cleanenv", str(fake_sif)]
    assert "-pdbin" in calls[0]
    assert str(cif_path) in calls[0]


def test_privateer_subprocess_env_is_minimal(monkeypatch) -> None:
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


def test_run_privateer_cli_rejects_non_json_stdout(monkeypatch, tmp_path) -> None:
    from lpmo_pipeline.qc import privateer_runner as mod

    class Proc:
        def __init__(self) -> None:
            self.stdout = "<html>Privateer text output</html>"

    def _fake_run(cmd, capture_output, text, check, timeout, env):
        return Proc()

    fake_sif = tmp_path / "privateer.sif"
    fake_sif.write_text("fake")
    cif_path = tmp_path / "input.cif"
    cif_path.write_text("data_test\n")

    monkeypatch.setattr(mod, "PRIVATEER_SIF_CANDIDATES", (fake_sif,))
    monkeypatch.setattr(mod.subprocess, "run", _fake_run)

    with pytest.raises(RuntimeError, match="did not emit JSON"):
        mod._run_privateer_cli(cif_path)
