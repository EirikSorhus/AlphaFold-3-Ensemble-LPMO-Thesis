from __future__ import annotations

from lpmo_pipeline.report.build_summary_json import build_summary


def test_build_summary_ignores_missing_crystal_tanimotos() -> None:
    summary = build_summary(
        run_id="test-run",
        mode="production",
        del_variant="a",
        qc_reports=[],
        cluster_results=[],
        geometry_stats=[],
        crystal_reports=[
            {"condition_id": "c1", "best_tanimoto": 0.25},
            {"condition_id": "c2", "best_tanimoto": None},
        ],
    )

    assert summary["crystal_stats"] == {
        "n_comparisons": 2,
        "mean_best_tanimoto": 0.25,
    }