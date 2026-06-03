from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import sys

import yaml


SPLIT_SCRIPT = Path(__file__).parent / "run_tests_scripts" / "split_clustering_pilot_selection_manifest.py"
COLLECT_SCRIPT = Path(__file__).parent / "run_tests_scripts" / "collect_clustering_pilot_shards.py"
SRC_PATH = Path(__file__).resolve().parents[1] / "src"


def _load_module(script_path: Path, module_name: str):
    if str(SRC_PATH) not in sys.path:
        sys.path.insert(0, str(SRC_PATH))
    spec = importlib.util.spec_from_file_location(module_name, script_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _run_main(module, argv: list[str]) -> int:
    previous_argv = sys.argv[:]
    sys.argv = argv
    try:
        return int(module.main())
    finally:
        sys.argv = previous_argv


def _write_selection_manifest(path: Path, *, construct_type: str, proteins: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.safe_dump(
            {
                "work_root": str(path.parent / f"work_{construct_type}"),
                "construct_type": construct_type,
                "latest_only": True,
                "selections": [{"protein_id": protein_id} for protein_id in proteins],
            },
            sort_keys=False,
        )
    )


def _write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key in seen:
                continue
            seen.add(key)
            fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _seed_fake_production_output(production_output: Path, *, shard_id: str, construct_type: str) -> None:
    production_output.mkdir(parents=True, exist_ok=True)

    common_rows = [
        {
            "pose_id": f"{construct_type}_{shard_id}_pose_1",
            "protein_id": f"{construct_type}_{shard_id}_protein",
            "ligand_id": "NAG4",
            "cluster_id": "0",
            "qc_status": "pass",
            "min_cu_c1": "4.1",
            "min_cu_c4": "6.2",
        }
    ]

    surfaces_with_rows = {
        "pose_manifest.tsv": [{"pose_id": common_rows[0]["pose_id"], "protein_id": common_rows[0]["protein_id"]}],
        "pose_confidence.tsv": [{"pose_id": common_rows[0]["pose_id"], "ptm": "0.82"}],
        "qc_attrition_table.tsv": [{"protein_id": common_rows[0]["protein_id"], "n_total": "1", "n_pass": "1"}],
        "pose_geometry.tsv": [{"pose_id": common_rows[0]["pose_id"], "min_cu_c1": "4.1", "min_cu_c4": "6.2"}],
        "pose_ifp_table.tsv": [{"pose_id": common_rows[0]["pose_id"], "feature": "Hydrophobic"}],
        "pose_residue_contact_table.tsv": [{"pose_id": common_rows[0]["pose_id"], "residue": "HIS1"}],
        "cluster_assignments.tsv": [{"pose_id": common_rows[0]["pose_id"], "cluster_id": "0"}],
        "medoid_manifest.tsv": [{"cluster_id": "0", "medoid_pose_id": common_rows[0]["pose_id"]}],
        "condition_table.tsv": [{"condition_id": f"{construct_type}_{shard_id}_cond", "protein_id": common_rows[0]["protein_id"]}],
        "cluster_table.tsv": [{"condition_id": f"{construct_type}_{shard_id}_cond", "cluster_id": "0"}],
        "protein_summary_table.tsv": [{"protein_id": common_rows[0]["protein_id"], "n_conditions": "1"}],
        "protein_condition_residue_scores.tsv": [{"protein_id": common_rows[0]["protein_id"], "residue": "HIS1", "score": "0.9"}],
        "protein_residue_regio_delta.tsv": [{"protein_id": common_rows[0]["protein_id"], "residue": "HIS1", "delta": "0.1"}],
        "condition_patch_summary.tsv": [{"condition_id": f"{construct_type}_{shard_id}_cond", "patch": "A"}],
        "protein_patch_summary.tsv": [{"protein_id": common_rows[0]["protein_id"], "patch": "A"}],
        "crystal_anchor_table.tsv": [{"pose_id": common_rows[0]["pose_id"], "pdb_id": "5ACI", "tanimoto": "0.42"}],
        "crystal_geometry_table.tsv": [{"pose_id": common_rows[0]["pose_id"], "crystal_min_cu_c1": "4.0"}],
        "crystal_ifp_diagnostic_summary.tsv": [{"condition_id": f"{construct_type}_{shard_id}_cond", "n_nonempty": "1"}],
        "metrics.csv": common_rows,
    }

    for filename, rows in surfaces_with_rows.items():
        _write_tsv(production_output / filename, rows)

    summary_json = production_output / "summary.json"
    summary_json.write_text(
        json.dumps(
            {
                "run_id": f"{construct_type}_{shard_id}",
                "dataset_stats": {"n_total_poses": 1},
                "cluster_stats": {"n_protein_ligand_combinations": 1},
                "geometry_stats": {"cu_c1_range": [4.1, 4.1]},
                "crystal_stats": {"n_comparisons": 1, "mean_best_tanimoto": 0.42},
            },
            indent=2,
        )
    )
    report_html = production_output / "report.html"
    report_html.write_text("<html><body>synthetic report</body></html>\n")

    analysis_summary = {
        "n_discovered": 1,
        "n_prepared": 1,
        "n_analyzed": 1,
        "pose_manifest_tsv": str(production_output / "pose_manifest.tsv"),
        "cluster_assignments_tsv": str(production_output / "cluster_assignments.tsv"),
        "medoid_manifest_tsv": str(production_output / "medoid_manifest.tsv"),
        "crystal_anchor_tsv": str(production_output / "crystal_anchor_table.tsv"),
        "metrics_csv": str(production_output / "metrics.csv"),
        "summary_json": str(summary_json),
        "report_html": str(report_html),
    }
    (production_output / "analysis_core_summary.json").write_text(json.dumps(analysis_summary, indent=2))


def test_staged_pipeline_synthetic_end_to_end_smoke(tmp_path: Path, monkeypatch) -> None:
    split_module = _load_module(SPLIT_SCRIPT, "split_selection_manifest_smoke")
    collect_module = _load_module(COLLECT_SCRIPT, "collect_shards_smoke")

    run_root = tmp_path / "staged_run"
    manifests_dir = run_root / "manifests"
    domain_manifest = manifests_dir / "clustering_pilot_domain_only_selection.yaml"
    full_manifest = manifests_dir / "clustering_pilot_full_length_selection.yaml"
    _write_selection_manifest(domain_manifest, construct_type="domain_only", proteins=["P001", "P002"])
    _write_selection_manifest(full_manifest, construct_type="full_length", proteins=["P101", "P102"])

    for construct_type, manifest_path in (
        ("domain_only", domain_manifest),
        ("full_length", full_manifest),
    ):
        shard_dir = run_root / "manifests" / "shards" / construct_type
        exit_code = _run_main(
            split_module,
            [
                str(SPLIT_SCRIPT),
                "--input-manifest",
                str(manifest_path),
                "--output-dir",
                str(shard_dir),
                "--proteins-per-shard",
                "2",
                "--balance-by",
                "protein_count",
                "--output-root-template",
                str(run_root / "{construct_type}_shards" / "{shard_id}" / "production_output"),
            ],
        )
        assert exit_code == 0

        with (shard_dir / "shard_index.tsv").open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                shard_id = str(row["shard_id"])
                production_output = Path(str(row["output_root"]))
                _seed_fake_production_output(production_output, shard_id=shard_id, construct_type=construct_type)
                shard_run_dir = run_root / f"{construct_type}_shards" / shard_id
                shard_run_dir.mkdir(parents=True, exist_ok=True)
                (shard_run_dir / "clustering_pilot_real_case_summary.json").write_text(
                    json.dumps(
                        {
                            "run_step_status": "executed",
                            "summary_step_status": "written",
                            "cli_exit_code": 0,
                            "n_discovered_poses": 1,
                        },
                        indent=2,
                    )
                )

    protein_metadata = tmp_path / "protein_metadata.tsv"
    _write_tsv(
        protein_metadata,
        [{"Protein_ID": "P001", "Family": "AA9", "regioselectivity": "C1"}],
    )
    core_fasta = tmp_path / "core.fasta"
    core_fasta.write_text(">P001\nMSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGV\n")

    def _fake_predictive(**kwargs):
        root = Path(kwargs["output_dir"]) / "10_predictive"
        root.mkdir(parents=True, exist_ok=True)
        summary_path = root / "predictive_summary.json"
        summary_path.write_text(json.dumps({"status": "synthetic"}, indent=2))
        return SimpleNamespace(
            summary_path=summary_path,
            modeling_table_paths={},
            metrics_paths={},
            predictions_paths={},
        )

    def _fake_cbm(**kwargs):
        root = Path(kwargs["output_dir"]) / "15_cbm_paired_analysis"
        root.mkdir(parents=True, exist_ok=True)
        summary_path = root / "cbm_paired_analysis_summary.json"
        summary_path.write_text(json.dumps({"status": "synthetic"}, indent=2))
        return SimpleNamespace(summary_path=summary_path, table_paths={})

    def _fake_family(**kwargs):
        root = Path(kwargs["output_dir"]) / "08_family_residue_enrichment"
        root.mkdir(parents=True, exist_ok=True)
        summary_path = root / "family_enrichment_summary.json"
        aligned = root / "family_aligned_residue_table.tsv"
        enrichment = root / "family_residue_enrichment.tsv"
        substrate_enrichment = root / "family_substrate_residue_enrichment.tsv"
        wrong_ligand_enrichment = root / "family_wrong_ligand_residue_enrichment.tsv"
        aligned.write_text("family\tpos\nAA9\t1\n")
        enrichment.write_text("family\tresidue\nAA9\tHIS1\n")
        substrate_enrichment.write_text("family\tsubstrate\tresidue\nAA9\tchitin\tHIS1\n")
        wrong_ligand_enrichment.write_text("family\tactive_substrate\tresidue\nAA9\tchitin\tHIS1\n")
        summary_path.write_text(json.dumps({"status": "synthetic"}, indent=2))
        return SimpleNamespace(
            summary_path=summary_path,
            processed_families={"AA9": {"n_rows": 1}},
            skipped_families={},
            family_aligned_residue_table_path=aligned,
            family_residue_enrichment_path=enrichment,
            family_substrate_residue_enrichment_path=substrate_enrichment,
            family_wrong_ligand_residue_enrichment_path=wrong_ligand_enrichment,
        )

    monkeypatch.setattr(collect_module, "run_predictive_postprocess", _fake_predictive)
    monkeypatch.setattr(collect_module, "run_cbm_paired_analysis", _fake_cbm)
    monkeypatch.setattr(collect_module, "run_family_enrichment_postprocess", _fake_family)

    output_json = run_root / "summaries" / "staged_array_summary.json"
    output_tsv = run_root / "summaries" / "staged_array_shards.tsv"
    merge_output_root = run_root / "merged_production_output"
    exit_code = _run_main(
        collect_module,
        [
            str(COLLECT_SCRIPT),
            "--run-root",
            str(run_root),
            "--output-json",
            str(output_json),
            "--output-tsv",
            str(output_tsv),
            "--merge-output-root",
            str(merge_output_root),
            "--run-global-postprocess",
            "--protein-metadata",
            str(protein_metadata),
            "--core-fasta",
            str(core_fasta),
        ],
    )

    assert exit_code == 0
    summary = json.loads(output_json.read_text())
    assert summary["all_shards_complete"] is True
    assert summary["n_shards"] == 2
    assert summary["n_failed_or_incomplete_shards"] == 0

    assert (merge_output_root / "pose_manifest.tsv").exists()
    assert (merge_output_root / "cluster_assignments.tsv").exists()
    assert (merge_output_root / "medoid_manifest.tsv").exists()
    assert (merge_output_root / "crystal_anchor_table.tsv").exists()
    assert (merge_output_root / "summary.json").exists()
    assert (merge_output_root / "report.html").exists()

    with (merge_output_root / "cluster_assignments.tsv").open(newline="") as handle:
        merged_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(merged_rows) == 2
    assert {row["_construct_type"] for row in merged_rows} == {"domain_only", "full_length"}

    postprocess = summary["postprocess"]
    assert "summary_path" in postprocess["predictive"]
    assert "summary_path" in postprocess["cbm_paired"]
    assert "summary_path" in postprocess["family_enrichment"]
    assert Path(postprocess["predictive"]["summary_path"]).exists()
    assert Path(postprocess["cbm_paired"]["summary_path"]).exists()
    assert Path(postprocess["family_enrichment"]["summary_path"]).exists()
