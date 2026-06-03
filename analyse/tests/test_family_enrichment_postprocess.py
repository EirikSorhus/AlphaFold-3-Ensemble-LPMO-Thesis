from __future__ import annotations

import csv
import json
from pathlib import Path

from lpmo_pipeline.analysis.family_enrichment_postprocess import (
    run_family_enrichment_postprocess,
)


def _write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_run_family_enrichment_postprocess_builds_outputs_from_precomputed_alignment(tmp_path: Path) -> None:
    scores_path = tmp_path / "protein_condition_residue_scores.tsv"
    delta_path = tmp_path / "protein_residue_regio_delta.tsv"
    metadata_path = tmp_path / "protein_metadata.tsv"
    core_fasta_path = tmp_path / "core.fasta"
    alignment_dir = tmp_path / "alignments"
    output_dir = tmp_path / "results"

    _write_tsv(
        scores_path,
        [
            {
                "protein_id": "P1",
                "condition_id": "P1__domain_only__chitin_DP4",
                "residue_chain": "A",
                "residue_number": "2",
                "residue_name": "ASN",
                "residue_label": "ASN2.A",
                "residue_contact_score": "0.8",
            },
            {
                "protein_id": "P2",
                "condition_id": "P2__domain_only__chitin_DP4",
                "residue_chain": "A",
                "residue_number": "2",
                "residue_name": "ASP",
                "residue_label": "ASP2.A",
                "residue_contact_score": "0.4",
            },
            {
                "protein_id": "P3",
                "condition_id": "P3__full_length__chitin_DP4",
                "residue_chain": "A",
                "residue_number": "2",
                "residue_name": "ASN",
                "residue_label": "ASN2.A",
                "residue_contact_score": "0.9",
            },
        ],
        [
            "protein_id",
            "condition_id",
            "residue_chain",
            "residue_number",
            "residue_name",
            "residue_label",
            "residue_contact_score",
        ],
    )
    _write_tsv(
        delta_path,
        [
            {
                "protein_id": "P1",
                "residue_chain": "A",
                "residue_number": "2",
                "residue_name": "ASN",
                "residue_label": "ASN2.A",
                "mean_c1_weighted_residue_score": "0.6",
                "mean_c4_weighted_residue_score": "0.1",
                "c1_minus_c4_weighted_delta": "0.5",
            },
            {
                "protein_id": "P2",
                "residue_chain": "A",
                "residue_number": "2",
                "residue_name": "ASP",
                "residue_label": "ASP2.A",
                "mean_c1_weighted_residue_score": "0.2",
                "mean_c4_weighted_residue_score": "0.3",
                "c1_minus_c4_weighted_delta": "-0.1",
            },
            {
                "protein_id": "P3",
                "residue_chain": "A",
                "residue_number": "2",
                "residue_name": "ASN",
                "residue_label": "ASN2.A",
                "mean_c1_weighted_residue_score": "0.9",
                "mean_c4_weighted_residue_score": "0.0",
                "c1_minus_c4_weighted_delta": "0.9",
            },
        ],
        [
            "protein_id",
            "residue_chain",
            "residue_number",
            "residue_name",
            "residue_label",
            "mean_c1_weighted_residue_score",
            "mean_c4_weighted_residue_score",
            "c1_minus_c4_weighted_delta",
        ],
    )
    _write_tsv(
        metadata_path,
        [
            {
                "UniProt_ID": "P1",
                "CAZy_family": "AA10",
                "Cat_Seq_Group": "AA10_G1",
            },
            {
                "UniProt_ID": "P2",
                "CAZy_family": "AA10",
                "Cat_Seq_Group": "AA10_G2",
            },
            {
                "UniProt_ID": "P3",
                "CAZy_family": "AA10",
                "Cat_Seq_Group": "AA10_G3",
            },
        ],
        ["UniProt_ID", "CAZy_family", "Cat_Seq_Group"],
    )
    core_fasta_path.write_text(
        ">UniProtIDs|P1|Example organism|AA10 enzyme\n"
        "MNAS\n"
        ">UniProtIDs|P2|Example organism|AA10 enzyme\n"
        "MDAS\n"
        ">UniProtIDs|P3|Example organism|AA10 enzyme\n"
        "MQAS\n"
    )
    alignment_dir.mkdir(parents=True, exist_ok=True)
    (alignment_dir / "AA10.aligned.fasta").write_text(
        ">AA10_G1\n"
        "MNAS\n"
        ">AA10_G2\n"
        "MDAS\n"
        ">AA10_G3\n"
        "MQAS\n"
    )

    result = run_family_enrichment_postprocess(
        protein_condition_residue_scores_path=scores_path,
        protein_residue_regio_delta_path=delta_path,
        protein_metadata_path=metadata_path,
        core_fasta_path=core_fasta_path,
        output_dir=output_dir,
        alignment_dir=alignment_dir,
    )

    aligned_rows = _read_tsv(result.family_aligned_residue_table_path)
    enrichment_rows = _read_tsv(result.family_residue_enrichment_path)
    substrate_enrichment_rows = _read_tsv(result.family_substrate_residue_enrichment_path)
    wrong_ligand_rows = _read_tsv(result.family_wrong_ligand_residue_enrichment_path)
    manifest_rows = _read_tsv(result.alignment_manifest_path)
    summary = json.loads(result.summary_path.read_text())

    assert len(aligned_rows) == 2
    assert {row["protein_id"] for row in aligned_rows} == {"P1", "P2"}
    assert all(row["construct_type"] == "domain_only" for row in aligned_rows)
    assert len(enrichment_rows) == 1
    assert enrichment_rows[0]["family_label"] == "AA10"
    assert enrichment_rows[0]["n_proteins_observed"] == "2"
    assert enrichment_rows[0]["n_family_aggregation_units_observed"] == "2"
    assert {row["target_substrate"] for row in substrate_enrichment_rows} == {
        "cellulose",
        "chitin",
    }
    assert {row["active_substrate"] for row in wrong_ligand_rows} == {"cellulose", "chitin"}
    assert all(row["family_label"] == "AA10" for row in substrate_enrichment_rows)
    assert len(manifest_rows) == 2
    assert summary["processed_families"]["AA10"]["alignment_source"] == "precomputed"
    assert summary["n_family_substrate_residue_enrichment_rows"] == len(substrate_enrichment_rows)
    assert summary["n_family_wrong_ligand_residue_enrichment_rows"] == len(wrong_ligand_rows)
    assert "AA9" in summary["skipped_families"]


def test_run_family_enrichment_postprocess_skips_when_alignment_is_unavailable(tmp_path: Path) -> None:
    scores_path = tmp_path / "protein_condition_residue_scores.tsv"
    delta_path = tmp_path / "protein_residue_regio_delta.tsv"
    metadata_path = tmp_path / "protein_metadata.tsv"
    core_fasta_path = tmp_path / "core.fasta"
    output_dir = tmp_path / "results"

    _write_tsv(
        scores_path,
        [
            {
                "protein_id": "P1",
                "condition_id": "P1__domain_only__chitin_DP4",
                "residue_chain": "A",
                "residue_number": "2",
                "residue_name": "ASN",
                "residue_label": "ASN2.A",
                "residue_contact_score": "0.8",
            },
            {
                "protein_id": "P2",
                "condition_id": "P2__domain_only__chitin_DP4",
                "residue_chain": "A",
                "residue_number": "2",
                "residue_name": "ASP",
                "residue_label": "ASP2.A",
                "residue_contact_score": "0.4",
            },
        ],
        [
            "protein_id",
            "condition_id",
            "residue_chain",
            "residue_number",
            "residue_name",
            "residue_label",
            "residue_contact_score",
        ],
    )
    _write_tsv(
        delta_path,
        [
            {
                "protein_id": "P1",
                "residue_chain": "A",
                "residue_number": "2",
                "residue_name": "ASN",
                "residue_label": "ASN2.A",
                "mean_c1_weighted_residue_score": "0.6",
                "mean_c4_weighted_residue_score": "0.1",
                "c1_minus_c4_weighted_delta": "0.5",
            },
            {
                "protein_id": "P2",
                "residue_chain": "A",
                "residue_number": "2",
                "residue_name": "ASP",
                "residue_label": "ASP2.A",
                "mean_c1_weighted_residue_score": "0.2",
                "mean_c4_weighted_residue_score": "0.3",
                "c1_minus_c4_weighted_delta": "-0.1",
            },
        ],
        [
            "protein_id",
            "residue_chain",
            "residue_number",
            "residue_name",
            "residue_label",
            "mean_c1_weighted_residue_score",
            "mean_c4_weighted_residue_score",
            "c1_minus_c4_weighted_delta",
        ],
    )
    _write_tsv(
        metadata_path,
        [
            {
                "UniProt_ID": "P1",
                "CAZy_family": "AA10",
                "Cat_Seq_Group": "AA10_G1",
            },
            {
                "UniProt_ID": "P2",
                "CAZy_family": "AA10",
                "Cat_Seq_Group": "AA10_G2",
            },
        ],
        ["UniProt_ID", "CAZy_family", "Cat_Seq_Group"],
    )
    core_fasta_path.write_text(
        ">UniProtIDs|P1|Example organism|AA10 enzyme\n"
        "MNAS\n"
        ">UniProtIDs|P2|Example organism|AA10 enzyme\n"
        "MDAS\n"
    )

    result = run_family_enrichment_postprocess(
        protein_condition_residue_scores_path=scores_path,
        protein_residue_regio_delta_path=delta_path,
        protein_metadata_path=metadata_path,
        core_fasta_path=core_fasta_path,
        output_dir=output_dir,
        mafft_executable="definitely-not-a-real-mafft",
    )

    summary = json.loads(result.summary_path.read_text())
    assert summary["skipped_families"]["AA10"]["reason"] == "mafft_not_available"
    assert _read_tsv(result.family_aligned_residue_table_path) == []
    assert _read_tsv(result.family_residue_enrichment_path) == []
    assert _read_tsv(result.family_substrate_residue_enrichment_path) == []
    assert _read_tsv(result.family_wrong_ligand_residue_enrichment_path) == []
