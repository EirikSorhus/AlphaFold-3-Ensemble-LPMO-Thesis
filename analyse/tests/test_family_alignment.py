from __future__ import annotations

import pytest

from lpmo_pipeline.analysis.family_alignment import (
    FamilySequenceMember,
    build_alignment_column_index,
    build_core_sequence_index,
    collapse_family_sequence_members,
    construct_type_from_condition_id,
    parse_fasta,
    parse_accessions_from_header,
    resolve_family_sequence_members,
)


def test_parse_accessions_from_header_supports_uniprotids_headers() -> None:
    header = "UniProtIDs|P9WEM3;O83009|Serratia marcescens|AA10 family enzyme"

    assert parse_accessions_from_header(header) == ("P9WEM3", "O83009")


def test_parse_fasta_reads_records_and_builds_accession_index(tmp_path) -> None:
    fasta_path = tmp_path / "core.fasta"
    fasta_path.write_text(
        ">UniProtIDs|P9WEM3;O83009|Serratia marcescens|AA10 family enzyme\n"
        "ACDEFG\n"
        ">UniProtIDs|A0A0S2GKZ1|Example organism|AA9 family enzyme\n"
        "FGHIKL\n"
    )

    records = parse_fasta(fasta_path)
    index = build_core_sequence_index(fasta_path)

    assert [record.sequence for record in records] == ["ACDEFG", "FGHIKL"]
    assert index["P9WEM3"].sequence == "ACDEFG"
    assert index["O83009"].sequence == "ACDEFG"
    assert index["A0A0S2GKZ1"].sequence == "FGHIKL"


def test_resolve_family_sequence_members_uses_coaccessions_and_filters_to_aa9_aa10(tmp_path) -> None:
    fasta_path = tmp_path / "core.fasta"
    fasta_path.write_text(
        ">UniProtIDs|P9WEM3;O83009|Serratia marcescens|AA10 family enzyme\n"
        "ACDEFG\n"
        ">UniProtIDs|AA9SEQ|Example fungus|AA9 family enzyme\n"
        "FGHIKL\n"
    )
    core_index = build_core_sequence_index(fasta_path)

    resolution = resolve_family_sequence_members(
        [
            {
                "UniProt_ID": "O83009",
                "CAZy_family": "AA10",
                "Cat_Seq_Group": "AA10_CORE_1",
            },
            {
                "UniProt_ID": "AA9SEQ",
                "CAZy_family": "AA9",
                "Cat_Seq_Group": "AA9_CORE_1",
            },
            {
                "UniProt_ID": "SKIPME",
                "CAZy_family": "AA11",
                "Cat_Seq_Group": "AA11_CORE_1",
            },
            {
                "UniProt_ID": "MISSING",
                "CAZy_family": "AA10",
                "Cat_Seq_Group": "AA10_CORE_2",
                "CoAccessions": "MISSING_ALT",
            },
        ],
        core_index,
    )

    assert resolution.members == [
        FamilySequenceMember(
            protein_id="AA9SEQ",
            family_label="AA9",
            family_aggregation_id="AA9_CORE_1",
            sequence="FGHIKL",
            matched_accession="AA9SEQ",
            sequence_header="UniProtIDs|AA9SEQ|Example fungus|AA9 family enzyme",
        ),
        FamilySequenceMember(
            protein_id="O83009",
            family_label="AA10",
            family_aggregation_id="AA10_CORE_1",
            sequence="ACDEFG",
            matched_accession="O83009",
            sequence_header="UniProtIDs|P9WEM3;O83009|Serratia marcescens|AA10 family enzyme",
        ),
    ]
    assert resolution.unresolved_rows == [
        {
            "protein_id": "MISSING",
            "family_label": "AA10",
            "reason": "missing_core_sequence",
        }
    ]


def test_collapse_family_sequence_members_groups_duplicate_catalytic_cores() -> None:
    members = [
        FamilySequenceMember(
            protein_id="P1",
            family_label="AA10",
            family_aggregation_id="AA10_CORE_1",
            sequence="ACDEFG",
            matched_accession="P1",
            sequence_header="header1",
        ),
        FamilySequenceMember(
            protein_id="P2",
            family_label="AA10",
            family_aggregation_id="AA10_CORE_1",
            sequence="ACDEFG",
            matched_accession="P2",
            sequence_header="header2",
        ),
    ]

    groups = collapse_family_sequence_members(members)

    assert len(groups) == 1
    assert groups[0].family_label == "AA10"
    assert groups[0].family_aggregation_id == "AA10_CORE_1"
    assert groups[0].protein_ids == ("P1", "P2")
    assert groups[0].matched_accessions == ("P1", "P2")


def test_collapse_family_sequence_members_rejects_inconsistent_group_sequences() -> None:
    members = [
        FamilySequenceMember(
            protein_id="P1",
            family_label="AA10",
            family_aggregation_id="AA10_CORE_1",
            sequence="ACDEFG",
            matched_accession="P1",
            sequence_header="header1",
        ),
        FamilySequenceMember(
            protein_id="P2",
            family_label="AA10",
            family_aggregation_id="AA10_CORE_1",
            sequence="ACDFFG",
            matched_accession="P2",
            sequence_header="header2",
        ),
    ]

    with pytest.raises(ValueError, match="Inconsistent catalytic-core sequences"):
        collapse_family_sequence_members(members)


def test_build_alignment_column_index_maps_ungapped_positions_to_alignment_columns() -> None:
    mapping = build_alignment_column_index("A-CD.EF")

    assert mapping == {1: 1, 2: 3, 3: 4, 4: 6, 5: 7}


def test_construct_type_from_condition_id_reads_second_token() -> None:
    assert construct_type_from_condition_id("P1__domain_only__chitin_DP4") == "domain_only"
    assert construct_type_from_condition_id("bad_condition_id") == ""