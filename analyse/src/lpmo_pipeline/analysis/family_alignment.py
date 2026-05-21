"""Sequence and alignment helpers for optional AA9/AA10 family enrichment."""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any


TARGET_FAMILY_LABELS = ("AA9", "AA10")
_GAP_CHARACTERS = {"-", "."}


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    description: str
    sequence: str
    accessions: tuple[str, ...]


@dataclass(frozen=True)
class FamilySequenceMember:
    protein_id: str
    family_label: str
    family_aggregation_id: str
    sequence: str
    matched_accession: str
    sequence_header: str


@dataclass(frozen=True)
class FamilySequenceGroup:
    family_label: str
    family_aggregation_id: str
    sequence: str
    protein_ids: tuple[str, ...]
    matched_accessions: tuple[str, ...]


@dataclass(frozen=True)
class FamilySequenceResolution:
    members: list[FamilySequenceMember]
    unresolved_rows: list[dict[str, str]]


def _unique_tokens(values: list[str]) -> tuple[str, ...]:
    seen: set[str] = set()
    unique: list[str] = []
    for value in values:
        token = value.strip()
        if not token or token in seen:
            continue
        seen.add(token)
        unique.append(token)
    return tuple(unique)


def split_accession_values(value: Any) -> tuple[str, ...]:
    if value in {None, "", "None"}:
        return ()
    return _unique_tokens(re.split(r"[;,\s]+", str(value).strip()))


def parse_accessions_from_header(header: str) -> tuple[str, ...]:
    leading_field = header.strip().split(maxsplit=1)[0]
    if leading_field.startswith("UniProtIDs|"):
        fields = leading_field.split("|", 2)
        if len(fields) >= 2:
            return split_accession_values(fields[1])
    if "|" in leading_field:
        fields = leading_field.split("|")
        for field in fields[1:]:
            accessions = split_accession_values(field)
            if accessions:
                return accessions
    return split_accession_values(leading_field)


def parse_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    description = ""
    sequence_chunks: list[str] = []

    def _flush() -> None:
        if not description:
            return
        sequence = "".join(sequence_chunks).strip().upper()
        identifier = description.split(maxsplit=1)[0]
        records.append(
            FastaRecord(
                identifier=identifier,
                description=description,
                sequence=sequence,
                accessions=parse_accessions_from_header(description),
            )
        )

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                _flush()
                description = line[1:].strip()
                sequence_chunks = []
                continue
            sequence_chunks.append(line)

    _flush()
    return records


def build_core_sequence_index(core_fasta_path: Path) -> dict[str, FastaRecord]:
    index: dict[str, FastaRecord] = {}
    for record in parse_fasta(core_fasta_path):
        for accession in record.accessions:
            existing = index.get(accession)
            if existing is not None and existing.sequence != record.sequence:
                raise ValueError(
                    f"Conflicting catalytic-core sequences for accession '{accession}'"
                )
            index[accession] = record
    return index


def protein_id_from_metadata_row(row: dict[str, Any]) -> str:
    for key in ("protein_id", "UniProt_ID", "uniprot_id"):
        value = str(row.get(key, "")).strip()
        if value:
            return value
    return ""


def family_label_from_metadata_row(row: dict[str, Any]) -> str:
    for key in ("family_label", "family", "CAZy_family", "cazy_family"):
        value = str(row.get(key, "")).strip()
        if value:
            return value
    return ""


def family_aggregation_id_from_metadata_row(row: dict[str, Any]) -> str:
    for key in (
        "family_aggregation_id",
        "Cat_Seq_Group",
        "cat_seq_group",
        "Sequence_Group",
        "sequence_group",
    ):
        value = str(row.get(key, "")).strip()
        if value:
            return value
    return protein_id_from_metadata_row(row)


def accession_candidates_from_metadata_row(row: dict[str, Any]) -> tuple[str, ...]:
    values: list[str] = []
    for key in ("protein_id", "UniProt_ID", "uniprot_id", "CoAccessions", "Cat_CoAccessions"):
        values.extend(split_accession_values(row.get(key, "")))
    return _unique_tokens(values)


def resolve_family_sequence_members(
    metadata_rows: list[dict[str, Any]],
    core_sequence_index: dict[str, FastaRecord],
    *,
    allowed_families: tuple[str, ...] = TARGET_FAMILY_LABELS,
) -> FamilySequenceResolution:
    allowed = {family.upper() for family in allowed_families}
    family_rank = {family.upper(): index for index, family in enumerate(allowed_families)}
    members: list[FamilySequenceMember] = []
    unresolved_rows: list[dict[str, str]] = []

    for row in metadata_rows:
        protein_id = protein_id_from_metadata_row(row)
        family_label = family_label_from_metadata_row(row)
        if not protein_id or family_label.upper() not in allowed:
            continue

        matched_accession = ""
        matched_record: FastaRecord | None = None
        for accession in accession_candidates_from_metadata_row(row):
            matched_record = core_sequence_index.get(accession)
            if matched_record is not None:
                matched_accession = accession
                break

        if matched_record is None:
            unresolved_rows.append(
                {
                    "protein_id": protein_id,
                    "family_label": family_label,
                    "reason": "missing_core_sequence",
                }
            )
            continue

        members.append(
            FamilySequenceMember(
                protein_id=protein_id,
                family_label=family_label,
                family_aggregation_id=family_aggregation_id_from_metadata_row(row),
                sequence=matched_record.sequence,
                matched_accession=matched_accession,
                sequence_header=matched_record.description,
            )
        )

    members.sort(
        key=lambda member: (
            family_rank.get(member.family_label.upper(), len(family_rank)),
            member.family_aggregation_id,
            member.protein_id,
        )
    )
    unresolved_rows.sort(
        key=lambda row: (
            family_rank.get(str(row.get("family_label", "")).upper(), len(family_rank)),
            row.get("protein_id", ""),
        )
    )
    return FamilySequenceResolution(members=members, unresolved_rows=unresolved_rows)


def collapse_family_sequence_members(
    members: list[FamilySequenceMember],
) -> list[FamilySequenceGroup]:
    grouped: dict[tuple[str, str], list[FamilySequenceMember]] = defaultdict(list)
    for member in members:
        grouped[(member.family_label, member.family_aggregation_id)].append(member)

    groups: list[FamilySequenceGroup] = []
    for family_label, family_aggregation_id in sorted(grouped):
        grouped_members = grouped[(family_label, family_aggregation_id)]
        sequences = {member.sequence for member in grouped_members}
        if len(sequences) != 1:
            raise ValueError(
                "Inconsistent catalytic-core sequences for "
                f"{family_label}:{family_aggregation_id}"
            )
        groups.append(
            FamilySequenceGroup(
                family_label=family_label,
                family_aggregation_id=family_aggregation_id,
                sequence=grouped_members[0].sequence,
                protein_ids=tuple(sorted({member.protein_id for member in grouped_members})),
                matched_accessions=tuple(
                    sorted({member.matched_accession for member in grouped_members if member.matched_accession})
                ),
            )
        )
    return groups


def write_family_alignment_input_fasta(
    groups: list[FamilySequenceGroup],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        for group in groups:
            protein_ids = ",".join(group.protein_ids)
            handle.write(
                f">{group.family_aggregation_id} family={group.family_label} proteins={protein_ids}\n"
            )
            handle.write(f"{group.sequence}\n")


def read_alignment_sequences(aligned_fasta_path: Path) -> dict[str, str]:
    return {
        record.identifier: record.sequence
        for record in parse_fasta(aligned_fasta_path)
    }


def build_alignment_column_index(gapped_sequence: str) -> dict[int, int]:
    mapping: dict[int, int] = {}
    residue_index = 0
    for alignment_column, residue in enumerate(gapped_sequence, start=1):
        if residue in _GAP_CHARACTERS:
            continue
        residue_index += 1
        mapping[residue_index] = alignment_column
    return mapping


def build_alignment_column_lookup(
    aligned_sequences: dict[str, str],
) -> dict[str, dict[int, int]]:
    return {
        identifier: build_alignment_column_index(sequence)
        for identifier, sequence in aligned_sequences.items()
    }


def construct_type_from_condition_id(condition_id: str) -> str:
    parts = str(condition_id).split("__")
    if len(parts) < 2:
        return ""
    return parts[1].strip()