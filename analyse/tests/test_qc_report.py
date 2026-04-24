from __future__ import annotations

import json
from pathlib import Path

from jsonschema import validate

from lpmo_pipeline.qc.custom_geometry_checks import CuHisMeasurement, GeometryResult
from lpmo_pipeline.qc.posebusters_runner import PoseBustersSingleResult
from lpmo_pipeline.qc.privateer_runner import PrivateerResult, PrivateerResidueResult
from lpmo_pipeline.qc.qc_report import build_qc_report, compute_verdict, write_qc_report


def test_write_qc_report_matches_schema(tmp_path: Path) -> None:
    pb_pass = PoseBustersSingleResult(pose_id="pose_pass", passed=True)
    geom_pass = GeometryResult(
        pose_id="pose_pass",
        cu_found=True,
        cu_his_all_in_range=True,
        cu_his_measurements=[
            CuHisMeasurement(
                his_chain="A",
                his_resnum=1,
                his_atom="NE2",
                cu_chain="E",
                cu_resnum=1,
                distance_angstrom=2.1,
                in_range=True,
            )
        ],
        min_cu_c1=3.2,
        min_cu_c4=4.4,
        passed=True,
    )
    verdict_pass = compute_verdict("pose_pass", pb_pass, None, geom_pass)

    pb_flag = PoseBustersSingleResult(
        pose_id="pose_flag",
        passed=True,
        warnings=["aromatic_ring_flatness"],
    )
    geom_flag = GeometryResult(
        pose_id="pose_flag",
        cu_found=True,
        cu_his_all_in_range=True,
        passed=True,
    )
    verdict_flag = compute_verdict("pose_flag", pb_flag, None, geom_flag)

    priv_fail = PrivateerResult(
        residues=[
            PrivateerResidueResult(
                chain="B",
                resname="XXX",
                resnum=1,
                sugar_recognized=False,
                anomer_ok=False,
                ring_pucker_ok=False,
            )
        ],
        total_sugars=1,
        recognized=0,
        recognition_rate=0.0,
        anomer_pass=0,
        ring_pucker_pass=0,
        linkage_pass=1,
        all_pass=False,
    )
    verdict_fail = compute_verdict("pose_fail", None, priv_fail, None)

    report = build_qc_report("qc_run", [verdict_pass, verdict_flag, verdict_fail])
    out_path = tmp_path / "qc_report.json"
    write_qc_report(report, out_path)

    payload = json.loads(out_path.read_text())
    schema = json.loads(
        (Path(__file__).resolve().parents[1] / "schemas" / "qc_report_schema.json").read_text()
    )
    validate(instance=payload, schema=schema)

    assert payload["total"] == 3
    assert len(payload["poses"]) == 3
    assert {p["overall_status"] for p in payload["poses"]} == {"pass", "soft_flag", "hard_fail"}


def test_privateer_runner_error_becomes_explicit_drop_reason() -> None:
    priv_result = PrivateerResult(
        residues=[],
        total_sugars=0,
        recognized=0,
        recognition_rate=0.0,
        anomer_pass=0,
        ring_pucker_pass=0,
        linkage_pass=0,
        all_pass=False,
        runner_error="privateer unavailable",
    )

    verdict = compute_verdict("pose_err", None, priv_result, None)

    assert verdict.status == "dropped"
    assert "privateer_runner_error:privateer unavailable" in verdict.drop_reasons
    assert all(not reason.startswith("privateer_unrecognized") for reason in verdict.drop_reasons)
    assert verdict.privateer["errors"]