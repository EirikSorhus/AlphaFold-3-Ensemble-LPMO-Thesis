"""MSA file matching and manifest generation.

Supports:
- Regular a3m files in directories
- a3m files inside squashfs images (mounted or overlayed)
"""

from __future__ import annotations

import csv
import gzip
import logging
import os
import re
import subprocess
import tempfile
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterator

from .proteins import ProteinRecord

logger = logging.getLogger(__name__)


@dataclass
class MSARecord:
    """An MSA file record matched to a protein."""

    protein_id: str  # Protein ID this MSA belongs to
    msa_path: str  # Full path to a3m file (or path inside sqsh)
    msa_filename: str = ""  # Just the filename for matching
    msa_source: str = "mmseqs"  # Source of MSA (mmseqs, colabfold, etc.)
    matched_id: str = ""  # Which UniProt ID matched
    match_method: str = "header_parse"  # How the match was made
    sqsh_file: str = ""  # Path to squashfs if MSA is inside sqsh
    num_sequences: int = 0  # Number of sequences in MSA (optional)

    def to_dict(self) -> dict:
        """Convert to dictionary for CSV export."""
        return asdict(self)


def list_squashfs_contents(sqsh_path: Path) -> list[str]:
    """List files inside a squashfs image.

    Args:
        sqsh_path: Path to .sqsh or .sqfs file

    Returns:
        List of filenames (just the filename, not full path) ending in .a3m
    """
    try:
        result = subprocess.run(
            ["unsquashfs", "-l", str(sqsh_path)],
            capture_output=True,
            text=True,
            timeout=60,
        )
        if result.returncode != 0:
            raise RuntimeError(f"unsquashfs failed: {result.stderr}")

        files = []
        for line in result.stdout.splitlines():
            line = line.strip()
            # Lines look like: squashfs-root/filename.a3m
            if line.endswith(".a3m") or line.endswith(".a3m.gz"):
                # Extract just filename
                filename = os.path.basename(line)
                files.append(filename)
        return files
    except FileNotFoundError:
        raise RuntimeError("unsquashfs not found - install squashfs-tools")
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"Timeout listing {sqsh_path}")


def _parse_a3m_header(a3m_path: Path) -> set[str]:
    """Extract UniProt IDs from the first line (query) of an a3m file.

    The query header often contains UniProt IDs in various formats.

    Args:
        a3m_path: Path to a3m file (can be gzipped)

    Returns:
        Set of potential UniProt IDs found in header
    """
    ids = set()

    try:
        # Handle gzipped files
        if a3m_path.suffix == ".gz":
            opener = gzip.open
            mode = "rt"
        else:
            opener = open
            mode = "r"

        with opener(a3m_path, mode) as f:
            first_line = f.readline().strip()

        if not first_line.startswith(">"):
            return ids

        header = first_line[1:]  # Remove >

        # Extract UniProt-like IDs (6-10 alphanumeric, starts with letter)
        # Pattern matches: A0A223GEC9, Q1K4Q1, P12345
        uniprot_pattern = r"\b([A-Z][A-Z0-9]{5,9})\b"
        matches = re.findall(uniprot_pattern, header, re.IGNORECASE)
        ids.update(m.upper() for m in matches)

    except Exception:
        pass

    return ids


def _extract_ids_from_filename(filename: str) -> set[str]:
    """Extract UniProt IDs from filename.

    Handles patterns like:
    - UniProtIDs_Q1K4Q1_Q873G1_Organism_...a3m
    - Q1K4Q1.a3m
    - some_prefix_Q1K4Q1_suffix.a3m

    Args:
        filename: Filename (without path)

    Returns:
        Set of potential UniProt IDs
    """
    ids = set()

    # Remove extension
    stem = filename.replace(".a3m.gz", "").replace(".a3m", "")

    # Split by common delimiters and check each part
    # This handles underscore-separated filenames like UniProtIDs_B6EQJ6_Organism
    parts = re.split(r"[_\-\s\|]+", stem)
    
    # UniProt ID pattern: starts with letter, 6-10 alphanumeric chars
    uniprot_pattern = r"^[A-Z][A-Z0-9]{5,9}$"
    for part in parts:
        if re.match(uniprot_pattern, part, re.IGNORECASE):
            # Exclude common non-ID words
            if part.upper() not in {"UNIPROTIDS", "UNKNOWN", "STRAIN"}:
                ids.add(part.upper())

    return ids


def match_msa(
    proteins: list[ProteinRecord],
    msa_dir: Path | None = None,
    msa_sqsh_files: list[Path] | None = None,
    strict: bool = False,
) -> tuple[list[MSARecord], list[str]]:
    """Match MSA files to proteins.

    Strategy:
    1. Build index of all a3m files with their potential UniProt IDs
    2. For each protein, find matching MSA by ID

    Supports both:
    - Regular directories with a3m files
    - Squashfs images containing a3m files

    Args:
        proteins: List of ProteinRecord objects
        msa_dir: Directory containing a3m files (optional)
        msa_sqsh_files: List of squashfs files containing a3m files (optional)
        strict: If True, raise error on ambiguous matches

    Returns:
        Tuple of (matched MSA records, list of protein IDs without MSA)
    """
    # Build MSA index: {uniprot_id: [(path_or_filename, sqsh_file)]}
    # sqsh_file is None for regular files, Path for squashfs files
    msa_index: dict[str, list[tuple[str, Path | None]]] = {}

    # Index regular a3m files from directory
    if msa_dir and msa_dir.is_dir():
        a3m_files = list(msa_dir.glob("**/*.a3m")) + list(msa_dir.glob("**/*.a3m.gz"))

        for a3m_path in a3m_files:
            # Get IDs from both header and filename
            ids_from_header = _parse_a3m_header(a3m_path)
            ids_from_filename = _extract_ids_from_filename(a3m_path.name)

            all_ids = ids_from_header | ids_from_filename

            for uid in all_ids:
                if uid not in msa_index:
                    msa_index[uid] = []
                msa_index[uid].append((str(a3m_path), None))

    # Index a3m files from squashfs images
    if msa_sqsh_files:
        for sqsh_path in msa_sqsh_files:
            if not sqsh_path.exists():
                continue

            try:
                filenames = list_squashfs_contents(sqsh_path)
                for filename in filenames:
                    ids_from_filename = _extract_ids_from_filename(filename)

                    for uid in ids_from_filename:
                        if uid not in msa_index:
                            msa_index[uid] = []
                        msa_index[uid].append((filename, sqsh_path))
            except RuntimeError as e:
                # Log warning but continue
                print(f"Warning: Could not list {sqsh_path}: {e}")

    # Match proteins to MSAs
    matched: list[MSARecord] = []
    missing: list[str] = []

    for protein in proteins:
        # Get all possible IDs for this protein
        candidate_ids = [protein.protein_id.upper()]
        if protein.uniprot_ids_all:
            candidate_ids.extend(
                uid.strip().upper() for uid in protein.uniprot_ids_all.split(";")
            )

        # Find matching MSA - collect (path_or_filename, sqsh_file) tuples
        found_entries: list[tuple[str, Path | None]] = []
        matched_id = ""

        for uid in candidate_ids:
            if uid in msa_index:
                found_entries.extend(msa_index[uid])
                if not matched_id:
                    matched_id = uid

        if len(found_entries) == 0:
            missing.append(protein.protein_id)
            continue

        # Deduplicate
        found_entries = list(set(found_entries))

        if len(found_entries) > 1:
            if strict:
                raise ValueError(
                    f"Ambiguous MSA match for {protein.protein_id}: "
                    f"found {len(found_entries)} files"
                )
            # Use first match (sorted for determinism)
            found_entries.sort(key=lambda x: (x[0], str(x[1] or "")))
            logger.warning(
                "Multiple MSA matches for protein %s: found %d, "
                "using first match (%s)",
                protein.protein_id,
                len(found_entries),
                found_entries[0][0],
            )

        path_or_filename, sqsh_file = found_entries[0]

        # Determine the path to use:
        # - For regular files: full path
        # - For squashfs: just the filename (will be mounted at runtime)
        if sqsh_file:
            # MSA is inside squashfs - store filename and sqsh path
            matched.append(
                MSARecord(
                    protein_id=protein.protein_id,
                    msa_path=path_or_filename,  # Just filename
                    msa_filename=path_or_filename,
                    msa_source="mmseqs",
                    matched_id=matched_id,
                    match_method="filename_in_sqsh",
                    sqsh_file=str(sqsh_file),
                )
            )
        else:
            # Regular file
            matched.append(
                MSARecord(
                    protein_id=protein.protein_id,
                    msa_path=path_or_filename,
                    msa_filename=Path(path_or_filename).name,
                    msa_source="mmseqs",
                    matched_id=matched_id,
                    match_method="header_and_filename",
                    sqsh_file="",
                )
            )

    return matched, missing


def write_msa_manifest(msa_records: list[MSARecord], output_path: Path) -> None:
    """Write MSA manifest to CSV.

    Args:
        msa_records: List of MSARecord objects
        output_path: Path to output CSV file
    """
    if not msa_records:
        # Write empty file with headers
        fieldnames = list(MSARecord.__dataclass_fields__.keys())
        with open(output_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
        return

    fieldnames = list(msa_records[0].to_dict().keys())

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for record in msa_records:
            writer.writerow(record.to_dict())


def load_msa_manifest(manifest_path: Path) -> list[MSARecord]:
    """Load MSA records from manifest CSV.

    Args:
        manifest_path: Path to manifest CSV file

    Returns:
        List of MSARecord objects
    """
    records = []
    with open(manifest_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            records.append(
                MSARecord(
                    protein_id=row["protein_id"],
                    msa_path=row["msa_path"],
                    msa_filename=row.get("msa_filename", ""),
                    msa_source=row.get("msa_source", "mmseqs"),
                    matched_id=row.get("matched_id", ""),
                    match_method=row.get("match_method", ""),
                    sqsh_file=row.get("sqsh_file", ""),
                    num_sequences=int(row.get("num_sequences", 0)),
                )
            )
    return records
