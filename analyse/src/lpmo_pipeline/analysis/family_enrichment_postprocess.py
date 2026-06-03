"""Optional within-family residue enrichment postprocess over Stage 16b outputs."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from lpmo_pipeline.analysis.family_alignment import (
    TARGET_FAMILY_LABELS,
    FamilySequenceGroup,
    FamilySequenceMember,
    build_alignment_column_lookup,
    build_core_sequence_index,
    collapse_family_sequence_members,
    construct_type_from_condition_id,
    protein_id_from_metadata_row,
    read_alignment_sequences,
    resolve_family_sequence_members,
    write_family_alignment_input_fasta,
)
from lpmo_pipeline.analysis.family_residue_enrichment import (
    compute_family_residue_enrichment_outputs,
    write_family_aligned_residue_table,
    write_family_residue_enrichment,
    write_family_substrate_residue_enrichment,
    write_family_wrong_ligand_residue_enrichment,
)


FAMILY_ALIGNMENT_MANIFEST_COLUMNS = [
    "family_label",
    "family_aggregation_id",
    "n_proteins",
    "protein_ids",
    "matched_accessions",
    "alignment_source",
    "alignment_fasta_path",
    "sequence_length",
]


@dataclass
class FamilyEnrichmentPostprocessResult:
    output_dir: Path
    summary_path: Path
    family_aligned_residue_table_path: Path
    family_residue_enrichment_path: Path
    family_substrate_residue_enrichment_path: Path
    family_wrong_ligand_residue_enrichment_path: Path
    alignment_manifest_path: Path
    processed_families: dict[str, dict[str, Any]] = field(default_factory=dict)
    skipped_families: dict[str, dict[str, Any]] = field(default_factory=dict)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in fieldnames})


def _as_int(value: Any) -> int:
    if value in {None, "", "None"}:
        return 0
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def _residue_label_from_row(row: dict[str, Any]) -> str:
    label = str(row.get("residue_label", "")).strip()
    if label:
        return label
    residue_name = str(row.get("residue_name", "")).strip()
    residue_number = _as_int(row.get("residue_number"))
    residue_chain = str(row.get("residue_chain", "")).strip()
    if not residue_name or residue_number <= 0:
        return ""
    return f"{residue_name}{residue_number}.{residue_chain}"


def _condition_rows_for_construct(
    rows: list[dict[str, str]],
    *,
    construct_type: str,
) -> tuple[list[dict[str, str]], set[str]]:
    filtered_rows: list[dict[str, str]] = []
    protein_ids: set[str] = set()
    for row in rows:
        protein_id = str(row.get("protein_id", "")).strip()
        if not protein_id:
            continue
        if construct_type_from_condition_id(row.get("condition_id", "")) != construct_type:
            continue
        filtered_rows.append(row)
        protein_ids.add(protein_id)
    return filtered_rows, protein_ids


def _index_members_by_protein(
    members: list[FamilySequenceMember],
) -> dict[str, FamilySequenceMember]:
    return {member.protein_id: member for member in members}


def _build_residue_catalog(
    protein_condition_residue_score_rows: list[dict[str, Any]],
    protein_residue_regio_delta_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    catalog: dict[tuple[str, str], dict[str, Any]] = {}
    for row in protein_residue_regio_delta_rows + protein_condition_residue_score_rows:
        protein_id = str(row.get("protein_id", "")).strip()
        residue_label = _residue_label_from_row(row)
        if not protein_id or not residue_label:
            continue
        entry = catalog.setdefault(
            (protein_id, residue_label),
            {
                "protein_id": protein_id,
                "residue_label": residue_label,
                "residue_chain": "",
                "residue_number": "",
                "residue_name": "",
            },
        )
        for key in ("residue_chain", "residue_number", "residue_name"):
            value = row.get(key, "")
            if value not in {None, "", "None"} and not entry.get(key):
                entry[key] = value
    return [catalog[key] for key in sorted(catalog)]


def _build_residue_alignment_rows(
    *,
    family_members: list[FamilySequenceMember],
    aligned_sequences: dict[str, str],
    protein_condition_residue_score_rows: list[dict[str, Any]],
    protein_residue_regio_delta_rows: list[dict[str, Any]],
    construct_type: str,
) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    member_by_protein = _index_members_by_protein(family_members)
    alignment_lookup = build_alignment_column_lookup(aligned_sequences)
    residue_catalog = _build_residue_catalog(
        protein_condition_residue_score_rows,
        protein_residue_regio_delta_rows,
    )

    residue_alignment_rows: list[dict[str, Any]] = []
    unresolved_rows: list[dict[str, str]] = []
    for residue_row in residue_catalog:
        protein_id = str(residue_row.get("protein_id", "")).strip()
        residue_label = str(residue_row.get("residue_label", "")).strip()
        member = member_by_protein.get(protein_id)
        if member is None:
            unresolved_rows.append(
                {
                    "protein_id": protein_id,
                    "residue_label": residue_label,
                    "reason": "missing_family_sequence_member",
                }
            )
            continue

        catalytic_core_residue_index = _as_int(residue_row.get("residue_number"))
        if catalytic_core_residue_index <= 0:
            unresolved_rows.append(
                {
                    "protein_id": protein_id,
                    "residue_label": residue_label,
                    "reason": "invalid_residue_number",
                }
            )
            continue

        group_alignment_lookup = alignment_lookup.get(member.family_aggregation_id)
        if group_alignment_lookup is None:
            unresolved_rows.append(
                {
                    "protein_id": protein_id,
                    "residue_label": residue_label,
                    "reason": "missing_alignment_lookup",
                }
            )
            continue

        alignment_column = group_alignment_lookup.get(catalytic_core_residue_index)
        if alignment_column is None:
            unresolved_rows.append(
                {
                    "protein_id": protein_id,
                    "residue_label": residue_label,
                    "reason": "residue_outside_alignment",
                }
            )
            continue

        alignment_residue = aligned_sequences[member.family_aggregation_id][alignment_column - 1]
        residue_alignment_rows.append(
            {
                "family_label": member.family_label,
                "family_aggregation_id": member.family_aggregation_id,
                "protein_id": protein_id,
                "construct_type": construct_type,
                "alignment_column": alignment_column,
                "alignment_residue": alignment_residue,
                "catalytic_core_residue_index": catalytic_core_residue_index,
                "residue_chain": residue_row.get("residue_chain", ""),
                "residue_number": catalytic_core_residue_index,
                "residue_name": residue_row.get("residue_name", ""),
                "residue_label": residue_label,
            }
        )

    residue_alignment_rows.sort(
        key=lambda row: (
            str(row.get("family_label", "")),
            str(row.get("alignment_column", "")),
            str(row.get("family_aggregation_id", "")),
            str(row.get("protein_id", "")),
            str(row.get("residue_label", "")),
        )
    )
    unresolved_rows.sort(key=lambda row: (row.get("protein_id", ""), row.get("residue_label", "")))
    return residue_alignment_rows, unresolved_rows


def _validate_alignment_sequences(
    groups: list[FamilySequenceGroup],
    aligned_sequences: dict[str, str],
) -> str:
    if not aligned_sequences:
        return "empty_alignment"

    alignment_lengths = {len(sequence) for sequence in aligned_sequences.values()}
    if len(alignment_lengths) > 1:
        return "alignment_length_mismatch"

    for group in groups:
        aligned_sequence = aligned_sequences.get(group.family_aggregation_id)
        if aligned_sequence is None:
            return f"missing_alignment_record:{group.family_aggregation_id}"
        ungapped_sequence = aligned_sequence.replace("-", "").replace(".", "")
        if ungapped_sequence != group.sequence:
            return f"alignment_sequence_mismatch:{group.family_aggregation_id}"
    return ""


def _resolve_alignment_path(
    *,
    family_label: str,
    groups: list[FamilySequenceGroup],
    alignments_dir: Path,
    precomputed_alignment_dir: Path | None,
    mafft_executable: str,
) -> tuple[Path | None, str, str]:
    if precomputed_alignment_dir is not None:
        precomputed_path = precomputed_alignment_dir / f"{family_label}.aligned.fasta"
        if precomputed_path.exists():
            return precomputed_path, "precomputed", ""

    resolved_mafft = shutil.which(mafft_executable)
    if resolved_mafft is None:
        return None, "", "mafft_not_available"

    input_path = alignments_dir / f"{family_label}.input.fasta"
    output_path = alignments_dir / f"{family_label}.aligned.fasta"
    write_family_alignment_input_fasta(groups, input_path)

    result = subprocess.run(
        [resolved_mafft, "--localpair", "--maxiterate", "1000", str(input_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        return None, "", "mafft_failed"
    output_path.write_text(result.stdout)
    return output_path, "mafft_linsi", ""


def run_family_enrichment_postprocess(
    *,
    protein_condition_residue_scores_path: Path,
    protein_residue_regio_delta_path: Path,
    protein_metadata_path: Path,
    core_fasta_path: Path,
    output_dir: Path,
    alignment_dir: Path | None = None,
    families: tuple[str, ...] = TARGET_FAMILY_LABELS,
    substrate_classes: tuple[str, ...] = ("cellulose", "chitin"),
    construct_type: str = "domain_only",
    min_family_aggregation_units: int = 2,
    mafft_executable: str = "mafft",
) -> FamilyEnrichmentPostprocessResult:
    protein_condition_rows = read_tsv(protein_condition_residue_scores_path)
    protein_delta_rows = read_tsv(protein_residue_regio_delta_path)
    protein_metadata_rows = read_tsv(protein_metadata_path)

    family_root = output_dir / "08_family_residue_enrichment"
    alignments_root = family_root / "alignments"
    family_root.mkdir(parents=True, exist_ok=True)
    alignments_root.mkdir(parents=True, exist_ok=True)

    family_aligned_residue_table_path = family_root / "family_aligned_residue_table.tsv"
    family_residue_enrichment_path = family_root / "family_residue_enrichment.tsv"
    family_substrate_residue_enrichment_path = family_root / "family_substrate_residue_enrichment.tsv"
    family_wrong_ligand_residue_enrichment_path = family_root / "family_wrong_ligand_residue_enrichment.tsv"
    alignment_manifest_path = family_root / "family_alignment_manifest.tsv"
    summary_path = family_root / "family_enrichment_summary.json"

    filtered_condition_rows, construct_protein_ids = _condition_rows_for_construct(
        protein_condition_rows,
        construct_type=construct_type,
    )
    filtered_delta_rows = [
        row
        for row in protein_delta_rows
        if str(row.get("protein_id", "")).strip() in construct_protein_ids
    ]
    metadata_rows = [
        row
        for row in protein_metadata_rows
        if protein_id_from_metadata_row(row) in construct_protein_ids
    ]

    core_sequence_index = build_core_sequence_index(core_fasta_path)
    resolution = resolve_family_sequence_members(
        metadata_rows,
        core_sequence_index,
        allowed_families=families,
    )

    members_by_family: dict[str, list[FamilySequenceMember]] = {}
    for member in resolution.members:
        members_by_family.setdefault(member.family_label, []).append(member)

    all_family_aligned_rows: list[dict[str, Any]] = []
    all_family_enrichment_rows: list[dict[str, Any]] = []
    all_family_substrate_enrichment_rows: list[dict[str, Any]] = []
    all_family_wrong_ligand_enrichment_rows: list[dict[str, Any]] = []
    alignment_manifest_rows: list[dict[str, Any]] = []
    processed_families: dict[str, dict[str, Any]] = {}
    skipped_families: dict[str, dict[str, Any]] = {}

    for family_label in families:
        family_members = members_by_family.get(family_label, [])
        if not family_members:
            skipped_families[family_label] = {"reason": "no_family_members"}
            continue

        family_groups = collapse_family_sequence_members(family_members)
        if len(family_groups) < min_family_aggregation_units:
            skipped_families[family_label] = {
                "reason": "insufficient_family_aggregation_units",
                "n_family_aggregation_units": len(family_groups),
            }
            continue

        alignment_path, alignment_source, skip_reason = _resolve_alignment_path(
            family_label=family_label,
            groups=family_groups,
            alignments_dir=alignments_root,
            precomputed_alignment_dir=alignment_dir,
            mafft_executable=mafft_executable,
        )
        if alignment_path is None:
            skipped_families[family_label] = {"reason": skip_reason}
            continue

        aligned_sequences = read_alignment_sequences(alignment_path)
        validation_error = _validate_alignment_sequences(family_groups, aligned_sequences)
        if validation_error:
            skipped_families[family_label] = {
                "reason": validation_error,
                "alignment_fasta_path": str(alignment_path),
            }
            continue

        family_protein_ids = {member.protein_id for member in family_members}
        family_condition_rows = [
            row for row in filtered_condition_rows if str(row.get("protein_id", "")).strip() in family_protein_ids
        ]
        family_delta_rows = [
            row for row in filtered_delta_rows if str(row.get("protein_id", "")).strip() in family_protein_ids
        ]
        residue_alignment_rows, unresolved_residue_rows = _build_residue_alignment_rows(
            family_members=family_members,
            aligned_sequences=aligned_sequences,
            protein_condition_residue_score_rows=family_condition_rows,
            protein_residue_regio_delta_rows=family_delta_rows,
            construct_type=construct_type,
        )
        if not residue_alignment_rows:
            skipped_families[family_label] = {
                "reason": "no_residue_alignment_rows",
                "alignment_fasta_path": str(alignment_path),
            }
            continue

        family_outputs = compute_family_residue_enrichment_outputs(
            protein_condition_residue_score_rows=family_condition_rows,
            protein_residue_regio_delta_rows=family_delta_rows,
            residue_alignment_rows=residue_alignment_rows,
            protein_metadata_rows=metadata_rows,
            substrate_classes=substrate_classes,
        )
        all_family_aligned_rows.extend(family_outputs.family_aligned_residue_table)
        all_family_enrichment_rows.extend(family_outputs.family_residue_enrichment)
        all_family_substrate_enrichment_rows.extend(
            family_outputs.family_substrate_residue_enrichment
        )
        all_family_wrong_ligand_enrichment_rows.extend(
            family_outputs.family_wrong_ligand_residue_enrichment
        )

        for group in family_groups:
            alignment_manifest_rows.append(
                {
                    "family_label": family_label,
                    "family_aggregation_id": group.family_aggregation_id,
                    "n_proteins": len(group.protein_ids),
                    "protein_ids": ";".join(group.protein_ids),
                    "matched_accessions": ";".join(group.matched_accessions),
                    "alignment_source": alignment_source,
                    "alignment_fasta_path": str(alignment_path),
                    "sequence_length": len(group.sequence),
                }
            )

        processed_families[family_label] = {
            "alignment_source": alignment_source,
            "alignment_fasta_path": str(alignment_path),
            "n_family_members": len(family_members),
            "n_family_aggregation_units": len(family_groups),
            "n_residue_alignment_rows": len(residue_alignment_rows),
            "n_unresolved_residue_rows": len(unresolved_residue_rows),
            "n_family_aligned_residue_rows": len(family_outputs.family_aligned_residue_table),
            "n_family_enrichment_rows": len(family_outputs.family_residue_enrichment),
            "n_family_substrate_enrichment_rows": len(
                family_outputs.family_substrate_residue_enrichment
            ),
            "n_family_wrong_ligand_residue_enrichment_rows": len(
                family_outputs.family_wrong_ligand_residue_enrichment
            ),
        }

    write_family_aligned_residue_table(all_family_aligned_rows, family_aligned_residue_table_path)
    write_family_residue_enrichment(all_family_enrichment_rows, family_residue_enrichment_path)
    write_family_substrate_residue_enrichment(
        all_family_substrate_enrichment_rows,
        family_substrate_residue_enrichment_path,
    )
    write_family_wrong_ligand_residue_enrichment(
        all_family_wrong_ligand_enrichment_rows,
        family_wrong_ligand_residue_enrichment_path,
    )
    _write_tsv(alignment_manifest_path, alignment_manifest_rows, FAMILY_ALIGNMENT_MANIFEST_COLUMNS)

    summary_data = {
        "protein_condition_residue_scores": str(protein_condition_residue_scores_path),
        "protein_residue_regio_delta": str(protein_residue_regio_delta_path),
        "protein_metadata": str(protein_metadata_path),
        "core_fasta": str(core_fasta_path),
        "alignment_dir": str(alignment_dir) if alignment_dir is not None else None,
        "families": list(families),
        "substrate_classes": list(substrate_classes),
        "construct_type_filter": construct_type,
        "min_family_aggregation_units": min_family_aggregation_units,
        "unresolved_core_sequences": resolution.unresolved_rows,
        "processed_families": processed_families,
        "skipped_families": skipped_families,
        "family_aligned_residue_table": str(family_aligned_residue_table_path),
        "family_residue_enrichment": str(family_residue_enrichment_path),
        "family_substrate_residue_enrichment": str(family_substrate_residue_enrichment_path),
        "family_wrong_ligand_residue_enrichment": str(family_wrong_ligand_residue_enrichment_path),
        "family_alignment_manifest": str(alignment_manifest_path),
        "n_family_aligned_residue_rows": len(all_family_aligned_rows),
        "n_family_residue_enrichment_rows": len(all_family_enrichment_rows),
        "n_family_substrate_residue_enrichment_rows": len(all_family_substrate_enrichment_rows),
        "n_family_wrong_ligand_residue_enrichment_rows": len(all_family_wrong_ligand_enrichment_rows),
    }
    summary_path.write_text(json.dumps(summary_data, indent=2))

    return FamilyEnrichmentPostprocessResult(
        output_dir=family_root,
        summary_path=summary_path,
        family_aligned_residue_table_path=family_aligned_residue_table_path,
        family_residue_enrichment_path=family_residue_enrichment_path,
        family_substrate_residue_enrichment_path=family_substrate_residue_enrichment_path,
        family_wrong_ligand_residue_enrichment_path=family_wrong_ligand_residue_enrichment_path,
        alignment_manifest_path=alignment_manifest_path,
        processed_families=processed_families,
        skipped_families=skipped_families,
    )
