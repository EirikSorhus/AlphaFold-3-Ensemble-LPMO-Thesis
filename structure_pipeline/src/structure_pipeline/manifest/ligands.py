"""Ligand scanning and manifest generation."""

from __future__ import annotations

import csv
import hashlib
import logging
import re
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterator

logger = logging.getLogger(__name__)


@dataclass
class LigandRecord:
    """A ligand record from CIF file or CCD code list."""

    ligand_id: str  # Generated ID (L001, L002, ...)
    category: str  # Subdirectory name (amylose, cellulose, chitin)
    ccd_code: str  # CCD code extracted from filename or file
    filename: str  # Original filename (empty for CCD-list mode)
    cif_path: str  # Full path to CIF file (empty for CCD-list mode)
    cif_sha256: str  # SHA256 hash of file contents (empty for CCD-list mode)

    def to_dict(self) -> dict:
        """Convert to dictionary for CSV export."""
        return asdict(self)


def _extract_ccd_code(cif_path: Path) -> str:
    """Extract CCD code from CIF file or filename.

    Tries:
    1. Parse _chem_comp.id from CIF file
    2. Use filename stem (without extension)

    Args:
        cif_path: Path to CIF file

    Returns:
        CCD code string
    """
    # Try to extract from CIF file content
    try:
        with open(cif_path) as f:
            content = f.read(4096)  # Read first 4KB

        # Look for _chem_comp.id or data_XXX
        match = re.search(r"_chem_comp\.id\s+(\S+)", content)
        if match:
            return match.group(1).strip("'\"")

        match = re.search(r"^data_(\S+)", content, re.MULTILINE)
        if match:
            return match.group(1)
    except Exception:
        pass

    # Fallback to filename
    return cif_path.stem.upper()


def _compute_file_hash(file_path: Path) -> str:
    """Compute SHA256 hash of file contents."""
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def scan_ligands(ligands_dir: Path) -> Iterator[LigandRecord]:
    """Scan ligands directory for CIF files.

    Supports multiple structures:
    1. ligands/{category}/*.cif
    2. ligands/{category}/cif/*.cif  (with cif subdirectory)
    3. ligands/*.cif (no category)

    Args:
        ligands_dir: Path to ligands directory

    Yields:
        LigandRecord for each CIF file found
    """
    if not ligands_dir.is_dir():
        raise FileNotFoundError(f"Ligands directory not found: {ligands_dir}")

    # Collect all CIF files with their categories
    cif_files: list[tuple[str, Path]] = []

    for category_dir in sorted(ligands_dir.iterdir()):
        if not category_dir.is_dir():
            continue

        # Skip non-category directories like boltz_ccd_lib
        if category_dir.name.startswith(".") or category_dir.name == "boltz_ccd_lib":
            continue

        category = category_dir.name

        # Check for CIF files directly in category dir
        for cif_file in sorted(category_dir.glob("*.cif")):
            cif_files.append((category, cif_file))

        # Check for CIF files in cif/ subdirectory
        cif_subdir = category_dir / "cif"
        if cif_subdir.is_dir():
            for cif_file in sorted(cif_subdir.glob("*.cif")):
                cif_files.append((category, cif_file))

    # Also check for CIF files directly in ligands_dir (no category)
    for cif_file in sorted(ligands_dir.glob("*.cif")):
        cif_files.append(("uncategorized", cif_file))

    # Generate sequential IDs and yield records
    for idx, (category, cif_path) in enumerate(cif_files, start=1):
        ligand_id = f"L{idx:03d}"
        ccd_code = _extract_ccd_code(cif_path)
        file_hash = _compute_file_hash(cif_path)

        yield LigandRecord(
            ligand_id=ligand_id,
            category=category,
            ccd_code=ccd_code,
            filename=cif_path.name,
            cif_path=str(cif_path.resolve()),  # Use absolute path
            cif_sha256=file_hash,
        )


def write_ligands_manifest(ligands: list[LigandRecord], output_path: Path) -> None:
    """Write ligands manifest to CSV.

    Args:
        ligands: List of LigandRecord objects
        output_path: Path to output CSV file
    """
    if not ligands:
        raise ValueError("No ligands to write")

    fieldnames = list(ligands[0].to_dict().keys())

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for ligand in ligands:
            writer.writerow(ligand.to_dict())


def load_ligands_manifest(manifest_path: Path) -> list[LigandRecord]:
    """Load ligands from manifest CSV.

    Args:
        manifest_path: Path to manifest CSV file

    Returns:
        List of LigandRecord objects
    """
    ligands = []
    with open(manifest_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ligands.append(
                LigandRecord(
                    ligand_id=row["ligand_id"],
                    category=row["category"],
                    ccd_code=row["ccd_code"],
                    filename=row["filename"],
                    cif_path=row["cif_path"],
                    cif_sha256=row["cif_sha256"],
                )
            )
    return ligands


def load_ccd_list(ccd_list_path: Path) -> list[LigandRecord]:
    """Load ligands from a plain-text CCD code list.

    Each non-empty, non-comment line is treated as a single CCD code.
    No CIF files are associated – runners will use standard CCD lookups.

    Args:
        ccd_list_path: Path to ``.txt`` file with one CCD code per line.

    Returns:
        List of LigandRecord objects with empty ``cif_path``.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If no valid CCD codes are found.
    """
    if not ccd_list_path.is_file():
        raise FileNotFoundError(f"CCD list file not found: {ccd_list_path}")

    codes: list[str] = []
    with open(ccd_list_path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            codes.append(stripped)

    if not codes:
        raise ValueError(f"No CCD codes found in {ccd_list_path}")

    seen: set[str] = set()
    records: list[LigandRecord] = []
    for idx, code in enumerate(codes, start=1):
        if code in seen:
            logger.warning("Duplicate CCD code '%s' in %s – skipping", code, ccd_list_path)
            continue
        seen.add(code)
        records.append(
            LigandRecord(
                ligand_id=f"L{idx:03d}",
                category="ccd_list",
                ccd_code=code,
                filename="",
                cif_path="",
                cif_sha256="",
            )
        )

    logger.info("Loaded %d CCD codes from %s", len(records), ccd_list_path)
    return records
