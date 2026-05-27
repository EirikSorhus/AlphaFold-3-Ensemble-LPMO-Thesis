from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from lpmo_pipeline.analysis.cluster_signatures import write_cluster_table_tsv
from lpmo_pipeline.analysis.cbm_comparison import (
    build_cbm_construct_condition_summary_rows,
    build_cbm_paired_comparison_rows,
    build_cbm_primary_metric_summary_rows,
    build_cbm_representative_example_rows,
    build_cbm_secondary_descriptive_summary_rows,
    build_cbm_stratified_summary_rows,
    run_cbm_paired_analysis,
    write_cbm_construct_condition_summary,
    write_cbm_paired_comparison_table,
    write_cbm_primary_metric_summary,
)
from lpmo_pipeline.analysis.condition_summary import (
    build_condition_table_rows,
    read_tsv,
    write_condition_table,
)


def _resolve_production_outputs_for_test(path: Path) -> list[Path]:
    if (path / "qc_attrition_table.tsv").exists() or (path / "cluster_signatures.json").exists():
        return [path]

    shard_outputs = sorted(
        child / "production_output"
        for child in path.iterdir()
        if child.is_dir() and child.name.startswith("shard_") and (child / "production_output").exists()
    )
    return shard_outputs


def _load_cluster_rows_for_test(production_output: Path) -> list[dict[str, str]]:
    cluster_table_path = production_output / "cluster_table.tsv"
    if cluster_table_path.exists():
        return read_tsv(cluster_table_path)

    cluster_signatures_path = production_output / "cluster_signatures.json"
    if not cluster_signatures_path.exists():
        raise FileNotFoundError(
            f"Missing both cluster_table.tsv and cluster_signatures.json under {production_output}"
        )

    import json

    payload = json.loads(cluster_signatures_path.read_text())
    if isinstance(payload, dict):
        clusters = payload.get("clusters") or []
    else:
        clusters = payload
    return list(clusters)


def _build_summary_tables_for_test(production_root: Path, output_dir: Path) -> tuple[Path, Path]:
    existing_condition_path = output_dir / "generated_condition_table.tsv"
    existing_cluster_path = output_dir / "generated_cluster_table.tsv"
    if existing_condition_path.exists() and existing_cluster_path.exists():
        return existing_condition_path, existing_cluster_path

    qc_attrition_rows: list[dict[str, str]] = []
    condition_cluster_summary_rows: list[dict[str, str]] = []
    condition_convergence_summary_rows: list[dict[str, str]] = []
    condition_patch_summary_rows: list[dict[str, str]] = []
    pose_confidence_rows: list[dict[str, str]] = []
    cluster_rows: list[dict[str, str]] = []

    for production_output in _resolve_production_outputs_for_test(production_root):
        try:
            cluster_rows.extend(_load_cluster_rows_for_test(production_output))
        except FileNotFoundError:
            continue
        qc_attrition_rows.extend(read_tsv(production_output / "qc_attrition_table.tsv"))
        condition_cluster_summary_rows.extend(read_tsv(production_output / "condition_cluster_summary.tsv"))
        condition_convergence_summary_rows.extend(read_tsv(production_output / "condition_convergence_summary.tsv"))
        condition_patch_summary_rows.extend(read_tsv(production_output / "condition_patch_summary.tsv"))
        pose_confidence_rows.extend(read_tsv(production_output / "pose_confidence.tsv"))

    if not qc_attrition_rows or not cluster_rows:
        raise FileNotFoundError(f"Could not build summary tables from {production_root}")

    cluster_table_path = output_dir / "generated_cluster_table.tsv"
    condition_table_path = output_dir / "generated_condition_table.tsv"
    write_cluster_table_tsv(cluster_rows, cluster_table_path)
    condition_rows = build_condition_table_rows(
        qc_attrition_rows=qc_attrition_rows,
        condition_cluster_summary_rows=condition_cluster_summary_rows,
        condition_convergence_summary_rows=condition_convergence_summary_rows,
        condition_patch_summary_rows=condition_patch_summary_rows,
        pose_confidence_rows=pose_confidence_rows,
        cluster_table_rows=cluster_rows,
    )
    write_condition_table(condition_rows, condition_table_path)
    return condition_table_path, cluster_table_path


def _write_merged_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_cbm_condition_summary_and_pairs_are_condition_level(tmp_path: Path) -> None:
    condition_rows = [
        {
            "condition_id": "P1__domain_only__chitin_DP4",
            "protein_id": "P1",
            "construct_type": "domain_only",
            "substrate_class": "chitin",
            "dp": "4",
            "n_generated": "10",
            "n_stage1_pass": "8",
            "n_ifp_success": "6",
            "n_clusters": "1",
            "noise_fraction": "0.1",
            "top_cluster_occupancy": "0.7",
            "cluster_entropy": "0.2",
            "occupancy_weighted_c1_plausible_fraction": "0.8",
            "occupancy_weighted_c4_plausible_fraction": "0.1",
            "catalytic_surface_contact_fraction": "0.5",
            "aromatic_contact_fraction": "0.25",
            "polar_contact_fraction": "0.3",
            "active_site_ifp_weighted_vector": "[0.2, 0.0, 0.5]",
        },
        {
            "condition_id": "P1__full_length__chitin_DP4",
            "protein_id": "P1",
            "construct_type": "full_length",
            "substrate_class": "chitin",
            "dp": "4",
            "n_generated": "10",
            "n_stage1_pass": "7",
            "n_ifp_success": "7",
            "n_clusters": "2",
            "noise_fraction": "0.2",
            "top_cluster_occupancy": "0.6",
            "cluster_entropy": "0.4",
            "occupancy_weighted_c1_plausible_fraction": "0.2",
            "occupancy_weighted_c4_plausible_fraction": "0.6",
            "catalytic_surface_contact_fraction": "0.4",
            "aromatic_contact_fraction": "0.35",
            "polar_contact_fraction": "0.2",
            "active_site_ifp_weighted_vector": "[0.0, 0.3, 0.5]",
        },
        {
            "condition_id": "P2__domain_only__chitin_DP4",
            "protein_id": "P2",
            "construct_type": "domain_only",
            "substrate_class": "chitin",
            "dp": "4",
            "n_generated": "4",
            "n_stage1_pass": "4",
            "n_ifp_success": "0",
            "n_clusters": "0",
        },
    ]
    cluster_rows = [
        {
            "condition_id": "P1__domain_only__chitin_DP4",
            "cluster_id": "0",
            "cluster_type": "C1_compatible",
            "occupancy": "0.7",
            "c1_geometry_computable_fraction": "1.0",
            "c4_geometry_computable_fraction": "1.0",
        },
        {
            "condition_id": "P1__full_length__chitin_DP4",
            "cluster_id": "0",
            "cluster_type": "C4_compatible",
            "occupancy": "0.6",
            "c1_geometry_computable_fraction": "1.0",
            "c4_geometry_computable_fraction": "1.0",
        },
        {
            "condition_id": "P1__full_length__chitin_DP4",
            "cluster_id": "1",
            "cluster_type": "uncertain",
            "occupancy": "0.4",
            "c1_geometry_computable_fraction": "1.0",
            "c4_geometry_computable_fraction": "1.0",
        },
    ]
    cluster_residue_rows = [
        {
            "condition_id": "P1__full_length__chitin_DP4",
            "cluster_id": "0",
            "contact_frequency": "1.0",
            "is_core_region": "True",
            "is_non_core_region": "False",
        },
        {
            "condition_id": "P1__full_length__chitin_DP4",
            "cluster_id": "0",
            "contact_frequency": "1.0",
            "is_core_region": "False",
            "is_non_core_region": "True",
        },
        {
            "condition_id": "P1__full_length__chitin_DP4",
            "cluster_id": "1",
            "contact_frequency": "1.0",
            "is_core_region": "True",
            "is_non_core_region": "False",
        },
    ]

    construct_rows = build_cbm_construct_condition_summary_rows(
        condition_rows,
        cluster_rows=cluster_rows,
        cluster_residue_rows=cluster_residue_rows,
        protein_metadata_rows=[{"protein_id": "P1", "family": "AA9", "cbm_type": "CBM1"}],
    )
    assert len(construct_rows) == 3
    domain = next(row for row in construct_rows if row["construct_type"] == "domain_only" and row["protein_id"] == "P1")
    full = next(row for row in construct_rows if row["construct_type"] == "full_length")
    assert domain["qc_pass_fraction"] == pytest.approx(0.8)
    assert domain["ifp_success_fraction"] == pytest.approx(0.75)
    assert domain["C1_compatible_fraction"] == pytest.approx(1.0)
    assert full["C4_compatible_fraction"] == pytest.approx(0.6)
    assert full["non_core_ligand_contact_fraction"] == pytest.approx(0.6)
    assert full["bridge_fraction"] == pytest.approx(0.6)
    assert full["non_core_recruitment_score"] == pytest.approx(0.6)

    pair_rows = build_cbm_paired_comparison_rows(construct_rows)
    assert len(pair_rows) == 1
    pair = pair_rows[0]
    assert pair["protein_id"] == "P1"
    assert pair["domain_only_condition_id"] == "P1__domain_only__chitin_DP4"
    assert pair["full_length_condition_id"] == "P1__full_length__chitin_DP4"
    assert pair["delta_qc_pass_fraction"] == pytest.approx(-0.1)
    assert pair["delta_C4_minus_C1_geometry_bias"] == pytest.approx(1.6)
    assert pair["non_core_ligand_contact_fraction_full_length"] == pytest.approx(0.6)
    assert pair["bridge_fraction_full_length"] == pytest.approx(0.6)
    assert pair["catalytic_domain_ifp_jaccard_distance"] == pytest.approx(2 / 3)

    construct_path = tmp_path / "cbm_construct_condition_summary.tsv"
    pair_path = tmp_path / "cbm_paired_comparison_table.tsv"
    write_cbm_construct_condition_summary(construct_rows, construct_path)
    write_cbm_paired_comparison_table(pair_rows, pair_path)
    with pair_path.open(newline="") as handle:
        written_pairs = list(csv.DictReader(handle, delimiter="\t"))
    assert written_pairs[0]["protein_id"] == "P1"
    assert written_pairs[0]["catalytic_domain_ifp_jaccard_distance"] != ""
    assert construct_path.exists()


def test_cbm_paired_summaries_follow_sample_size_rules_and_figure_contract(tmp_path: Path) -> None:
    paired_rows = []
    for index in range(12):
        paired_rows.append(
            {
                "protein_id": f"P{index:02d}",
                "family_label": "AA9",
                "cbm_type": "CBM1" if index < 6 else "CBM2",
                "substrate_class": "chitin" if index % 2 == 0 else "cellulose",
                "dp": "4" if index < 8 else "6",
                "domain_only_condition_id": f"P{index:02d}__domain_only",
                "full_length_condition_id": f"P{index:02d}__full_length",
                "bridge_fraction_full_length": 0.15 + index * 0.01,
                "non_core_ligand_contact_fraction_full_length": 0.25 + index * 0.01,
                "non_core_recruitment_score_full_length": 0.25 + index * 0.01,
                "catalytic_domain_ifp_jaccard_distance": 0.1 + index * 0.02,
                "delta_C4_minus_C1_geometry_bias": 0.1 + index * 0.03,
                "delta_qc_pass_fraction": -0.05 if index < 6 else 0.05,
                "delta_cluster_entropy": 0.2 + index * 0.01,
                "delta_C1_compatible_fraction": -0.1,
                "delta_C4_compatible_fraction": 0.1,
                "delta_geometry_plausible_fraction": 0.02,
                "delta_ifp_success_fraction": 0.03,
                "delta_noise_fraction": -0.04,
                "delta_top_cluster_occupancy": 0.08,
                "delta_n_clusters": 1,
                "delta_catalytic_surface_contact_fraction": -0.03,
                "delta_aromatic_contact_fraction": 0.04,
                "delta_polar_contact_fraction": -0.02,
            }
        )

    primary_rows = build_cbm_primary_metric_summary_rows(paired_rows, random_state=7)
    secondary_rows = build_cbm_secondary_descriptive_summary_rows(paired_rows)
    substrate_rows = build_cbm_stratified_summary_rows(paired_rows, stratum_name="substrate_class")
    example_rows = build_cbm_representative_example_rows(paired_rows, top_n=1)

    bridge_summary = next(row for row in primary_rows if row["metric_name"] == "bridge_fraction")
    geometry_summary = next(
        row for row in primary_rows if row["metric_name"] == "delta_C4_minus_C1_geometry_bias"
    )
    assert bridge_summary["sample_size_label"] == "adequate_for_simple_nonparametric_paired_summary"
    assert bridge_summary["statistical_test"] == "descriptive_full_length_only"
    assert bridge_summary["wilcoxon_p_value"] == ""
    assert geometry_summary["statistical_test"] == "paired_wilcoxon_signed_rank"
    assert geometry_summary["wilcoxon_p_value"] != ""
    assert geometry_summary["bh_fdr_p_value"] != ""
    assert geometry_summary["bootstrap_ci_low"] != ""
    assert len(secondary_rows) > 5
    assert {row["stratum_value"] for row in substrate_rows} == {"cellulose", "chitin"}
    assert example_rows[0]["example_type"] == "full_length_medoid_with_ligand_bridging_catalytic_domain_and_cbm"

    primary_path = tmp_path / "cbm_primary_metric_summary.tsv"
    write_cbm_primary_metric_summary(primary_rows, primary_path)
    with primary_path.open(newline="") as handle:
        written_primary = list(csv.DictReader(handle, delimiter="\t"))
    assert written_primary[0]["sample_size_label"]
    assert "interpretation_scope" not in written_primary[0]


def test_run_cbm_paired_analysis_writes_complete_side_analysis(tmp_path: Path) -> None:
    condition_table = tmp_path / "condition_table.tsv"
    cluster_table = tmp_path / "cluster_table.tsv"
    cluster_residue_table = tmp_path / "cluster_residue_signature.tsv"
    metadata_table = tmp_path / "protein_metadata.tsv"
    output_dir = tmp_path / "results"

    condition_table.write_text(
        "\t".join(
            [
                "condition_id",
                "protein_id",
                "construct_type",
                "substrate_class",
                "dp",
                "n_generated",
                "n_stage1_pass",
                "n_ifp_success",
                "n_clusters",
                "noise_fraction",
                "top_cluster_occupancy",
                "cluster_entropy",
                "occupancy_weighted_c1_plausible_fraction",
                "occupancy_weighted_c4_plausible_fraction",
                "active_site_ifp_weighted_vector",
            ]
        )
        + "\n"
        + "P1__domain\tP1\tdomain_only\tchitin\t4\t10\t8\t8\t1\t0.1\t0.7\t0.2\t0.8\t0.1\t[0.2, 0.0]\n"
        + "P1__full\tP1\tfull_length\tchitin\t4\t10\t7\t7\t1\t0.2\t0.6\t0.4\t0.2\t0.7\t[0.0, 0.3]\n"
    )
    cluster_table.write_text(
        "condition_id\tcluster_id\tcluster_type\toccupancy\tc1_geometry_computable_fraction\tc4_geometry_computable_fraction\n"
        "P1__domain\t0\tC1_compatible\t0.7\t1.0\t1.0\n"
        "P1__full\t0\tC4_compatible\t0.6\t1.0\t1.0\n"
    )
    cluster_residue_table.write_text(
        "condition_id\tcluster_id\tcontact_frequency\tis_core_region\tis_non_core_region\n"
        "P1__full\t0\t1.0\tTrue\tFalse\n"
        "P1__full\t0\t1.0\tFalse\tTrue\n"
    )
    metadata_table.write_text("protein_id\tfamily\tcbm_type\nP1\tAA9\tCBM1\n")

    result = run_cbm_paired_analysis(
        condition_table_path=condition_table,
        cluster_table_path=cluster_table,
        cluster_residue_signature_table_path=cluster_residue_table,
        protein_metadata_path=metadata_table,
        output_dir=output_dir,
        random_state=5,
    )

    assert result.summary_path.exists()
    expected_tables = {
        "construct_condition_summary",
        "paired_comparison",
        "primary_metric_summary",
        "secondary_descriptive_summary",
        "stratified_summary_by_substrate",
        "stratified_summary_by_dp",
        "stratified_summary_by_cbm_type",
        "representative_examples",
    }
    assert set(result.table_paths) == expected_tables
    for path in result.table_paths.values():
        assert path.exists()
    assert not (result.output_dir / "cbm_figure_manifest.tsv").exists()
    summary = json.loads(result.summary_path.read_text())
    assert summary["figures_to_generate"]
    assert summary["n_full_length_conditions_with_bridge_flag"] == 1


def test_run_cbm_paired_analysis_on_staged_real_summary_tables(tmp_path: Path) -> None:
    project_root = Path(__file__).resolve().parents[1]
    domain_root = project_root / "tests/tests_results/clustering_pilot_staged/domain_only_shards"
    full_root = project_root / "tests/tests_results/clustering_pilot_staged/full_length_shards"
    metadata_path = project_root / "input_data/metadata_final_ec_fixed.tsv"

    if not domain_root.exists() or not full_root.exists() or not metadata_path.exists():
        pytest.skip("staged real CBM validation inputs are not available in this workspace")

    domain_condition_path, domain_cluster_path = _build_summary_tables_for_test(
        domain_root,
        tmp_path / "domain_summary_validation",
    )
    full_condition_path, full_cluster_path = _build_summary_tables_for_test(
        full_root,
        tmp_path / "full_summary_validation",
    )

    domain_condition_rows = read_tsv(domain_condition_path)
    full_condition_rows = read_tsv(full_condition_path)
    overlap_keys = {
        (row["protein_id"], row["substrate_class"], row["dp"])
        for row in domain_condition_rows
    } & {
        (row["protein_id"], row["substrate_class"], row["dp"])
        for row in full_condition_rows
    }
    assert overlap_keys

    merged_condition_path = tmp_path / "merged_condition_table.tsv"
    merged_cluster_path = tmp_path / "merged_cluster_table.tsv"
    _write_merged_tsv(merged_condition_path, domain_condition_rows + full_condition_rows)
    _write_merged_tsv(
        merged_cluster_path,
        read_tsv(domain_cluster_path) + read_tsv(full_cluster_path),
    )

    result = run_cbm_paired_analysis(
        condition_table_path=merged_condition_path,
        cluster_table_path=merged_cluster_path,
        protein_metadata_path=metadata_path,
        output_dir=tmp_path,
        random_state=11,
    )

    paired_rows = list(csv.DictReader(result.table_paths["paired_comparison"].open(), delimiter="\t"))
    primary_rows = list(csv.DictReader(result.table_paths["primary_metric_summary"].open(), delimiter="\t"))

    assert len(paired_rows) == len(overlap_keys)
    assert all(row["domain_only_condition_id"] for row in paired_rows)
    assert all(row["full_length_condition_id"] for row in paired_rows)
    assert all(row["family_label"] for row in paired_rows)
    assert all(row["cbm_type"] for row in paired_rows)
    assert len(
        {
            (row["protein_id"], row["substrate_class"], row["dp"])
            for row in paired_rows
        }
    ) == len(paired_rows)

    delta_qc_summary = next(row for row in primary_rows if row["metric_name"] == "delta_qc_pass_fraction")
    bridge_summary = next(row for row in primary_rows if row["metric_name"] == "bridge_fraction")
    assert int(delta_qc_summary["n_pairs"]) == len(paired_rows)
    assert delta_qc_summary["statistical_test"] in {"sign_test", "paired_wilcoxon_signed_rank"}
    assert bridge_summary["statistical_test"] == "descriptive_full_length_only"
    assert result.summary_path.exists()
