from __future__ import annotations

from pathlib import Path

import yaml

from lpmo_pipeline.report.build_report_html import build_report_html


def test_build_report_html_renders_descriptive_sections_without_placeholders(tmp_path: Path) -> None:
    metrics_csv = tmp_path / "metrics.csv"
    metrics_csv.write_text(
        "\n".join(
            [
                "run_id,protein_id,ligand_id,model,seed,pose_id,cluster_id,qc_status,min_cu_c1,min_cu_c4,pocket_rmsd_vs_crystal,ifp_similarity_crystal,n_ifp_contacts",
                "run-1,P1,NAG4,af3,1,pose-1,0,passed,3.2,4.1,1.4,0.72,5",
                "run-1,P1,NAG4,af3,1,pose-2,1,flagged,,5.2,,0.31,2",
                "",
            ]
        )
    )
    output_path = tmp_path / "report.html"
    summary = {
        "pipeline_version": "2.1",
        "run_id": "run-1",
        "timestamp": "2026-05-22T00:00:00Z",
        "mode": "production",
        "del_variant": "a",
        "dataset_stats": {
            "n_proteins": 1,
            "n_ligands": 1,
            "n_models": 1,
            "n_total_poses": 3,
        },
        "qc_stats": {
            "total_poses_qc": 3,
            "passed": 1,
            "flagged": 1,
            "dropped": 1,
            "pass_rate": 1 / 3,
        },
        "cluster_stats": {
            "mean_n_clusters": 2.0,
            "mean_outlier_rate": 0.1,
            "n_protein_ligand_combinations": 1,
        },
        "geometry_stats": {
            "cu_c1_range": [3.2, 3.2],
            "cu_c4_range": [4.1, 5.2],
            "n_measurements": 2,
        },
        "crystal_stats": {
            "n_comparisons": 2,
            "mean_best_tanimoto": 0.515,
        },
    }

    build_report_html(summary, metrics_csv, output_path)

    html = output_path.read_text()
    assert "Dataset" in html
    assert "QC Summary" in html
    assert "Cluster Landscape" in html
    assert "Geometry" in html
    assert "Crystal Comparison" in html
    assert "33.3%" in html
    assert "pose-1" not in html
    assert "PLACEHOLDER" not in html
    assert "Jinja2" not in html


def test_thresholds_do_not_advertise_inactive_cluster_inclusion_contract() -> None:
    thresholds = yaml.safe_load(Path("configs/thresholds.yaml").read_text())

    assert "cluster_inclusion" not in thresholds
    assert "min_occupancy_for_main_summary" not in thresholds.get("clustering", {})
