"""FASTA parsing and protein manifest generation."""

from __future__ import annotations

import csv
import hashlib
import re
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterator


@dataclass
class ProteinRecord:
    """A protein record from FASTA with parsed metadata."""

    protein_id: str  # Primary UniProt ID (first before ;)
    uniprot_ids_all: str  # All UniProt IDs (semicolon-separated)
    organism: str  # Organism from header
    annotation: str  # Protein annotation/name
    sequence: str  # Amino acid sequence
    sequence_sha256: str  # SHA256 hash of sequence
    sequence_length: int  # Length of sequence
    fasta_header_raw: str  # Original FASTA header

    def to_dict(self) -> dict:
        """Convert to dictionary for CSV export."""
        return asdict(self)


def _parse_header(header: str) -> dict:
    """Parse FASTA header with UniProt format.

    Supports two formats:
    1. Pipe-separated (preferred):
       >UniProtIDs|ID1;ID2|Organism name|Annotation text
       Example: >UniProtIDs|Q1K4Q1|Neurospora crassa|AA9 family LPMO

    2. Underscore-separated (legacy):
       >UniProtIDs_ID1;ID2_Organism_name_Annotation_text
       Example: >UniProtIDs_Q1K4Q1;Q873G1_Neurospora_crassa_AA9_family_LPMO

    Returns:
        Dictionary with parsed fields
    """
    # Remove leading > if present
    header = header.lstrip(">").strip()

    # Try pipe-separated format first (preferred)
    # Pattern: UniProtIDs|{IDs}|{organism}|{annotation}
    if header.startswith("UniProtIDs|"):
        parts = header.split("|", 3)  # Split into max 4 parts
        uniprot_ids = parts[1] if len(parts) > 1 else ""
        organism = parts[2] if len(parts) > 2 else ""
        annotation = parts[3] if len(parts) > 3 else ""

        # Get primary ID (first before semicolon)
        primary_id = uniprot_ids.split(";")[0] if uniprot_ids else ""

        return {
            "protein_id": primary_id,
            "uniprot_ids_all": uniprot_ids,
            "organism": organism,
            "annotation": annotation,
        }

    # Try underscore-separated format (legacy)
    # Pattern: UniProtIDs_{IDs}_{organism}_{annotation}
    parts = header.split("_", 3)  # Split into max 4 parts

    if len(parts) >= 2 and parts[0] == "UniProtIDs":
        uniprot_ids = parts[1] if len(parts) > 1 else ""
        organism = parts[2] if len(parts) > 2 else ""
        annotation = parts[3] if len(parts) > 3 else ""

        # Get primary ID (first before semicolon)
        primary_id = uniprot_ids.split(";")[0] if uniprot_ids else ""

        return {
            "protein_id": primary_id,
            "uniprot_ids_all": uniprot_ids,
            "organism": organism.replace("_", " "),
            "annotation": annotation.replace("_", " "),
        }

    # Fallback: try standard UniProt header format
    # >sp|P12345|PROT_HUMAN Description OS=Homo sapiens
    sp_match = re.match(r"^(?:sp|tr)\|([A-Z0-9]+)\|", header)
    if sp_match:
        return {
            "protein_id": sp_match.group(1),
            "uniprot_ids_all": sp_match.group(1),
            "organism": "",
            "annotation": header,
        }

    # Last resort: use first word as ID
    first_word = header.split()[0] if header else "unknown"
    return {
        "protein_id": first_word,
        "uniprot_ids_all": first_word,
        "organism": "",
        "annotation": header,
    }


def parse_fasta(fasta_path: Path) -> Iterator[ProteinRecord]:
    """Parse a FASTA file and yield ProteinRecords.

    Args:
        fasta_path: Path to FASTA file

    Yields:
        ProteinRecord for each sequence in the file
    """
    current_header = None
    current_sequence_parts: list[str] = []

    def make_record(header: str, seq_parts: list[str]) -> ProteinRecord:
        sequence = "".join(seq_parts).upper()
        parsed = _parse_header(header)
        return ProteinRecord(
            protein_id=parsed["protein_id"],
            uniprot_ids_all=parsed["uniprot_ids_all"],
            organism=parsed["organism"],
            annotation=parsed["annotation"],
            sequence=sequence,
            sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
            sequence_length=len(sequence),
            fasta_header_raw=header,
        )

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # Save previous record if exists
                if current_header is not None:
                    yield make_record(current_header, current_sequence_parts)

                current_header = line[1:]  # Remove >
                current_sequence_parts = []
            else:
                current_sequence_parts.append(line)

        # Don't forget the last record
        if current_header is not None:
            yield make_record(current_header, current_sequence_parts)


def write_proteins_manifest(
    proteins: list[ProteinRecord], output_path: Path
) -> None:
    """Write proteins manifest to CSV.

    Args:
        proteins: List of ProteinRecord objects
        output_path: Path to output CSV file
    """
    if not proteins:
        raise ValueError("No proteins to write")

    fieldnames = list(proteins[0].to_dict().keys())

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for protein in proteins:
            writer.writerow(protein.to_dict())


def load_proteins_manifest(manifest_path: Path) -> list[ProteinRecord]:
    """Load proteins from manifest CSV.

    Args:
        manifest_path: Path to manifest CSV file

    Returns:
        List of ProteinRecord objects
    """
    proteins = []
    with open(manifest_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            proteins.append(
                ProteinRecord(
                    protein_id=row["protein_id"],
                    uniprot_ids_all=row["uniprot_ids_all"],
                    organism=row["organism"],
                    annotation=row["annotation"],
                    sequence=row["sequence"],
                    sequence_sha256=row["sequence_sha256"],
                    sequence_length=int(row["sequence_length"]),
                    fasta_header_raw=row["fasta_header_raw"],
                )
            )
    return proteins
