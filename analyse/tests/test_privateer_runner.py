from __future__ import annotations

from lpmo_pipeline.qc.privateer_runner import (
    _aggregate_privateer_results,
    _parse_privateer_output,
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
