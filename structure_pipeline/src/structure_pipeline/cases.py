"""Case generation - expand proteins × ligands × models into job cases."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass, asdict
from enum import Enum
from pathlib import Path
from typing import Iterator

from .manifest.proteins import ProteinRecord
from .manifest.ligands import LigandRecord
from .manifest.msa import MSARecord


class ModelType(str, Enum):
    """Supported structure prediction models."""

    AF3 = "af3"
    BOLTZ = "boltz"
    RF3 = "rf3"


class CaseStatus(str, Enum):
    """Status of a prediction case."""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    SKIPPED = "skipped"  # e.g., missing MSA for boltz/rf3


@dataclass
class Case:
    """A single prediction case: one protein + one ligand + one model."""

    case_id: str  # Unique hash-based ID
    protein_id: str
    ligand_id: str
    ligand_ccd_code: str
    model: str  # af3, boltz, rf3
    msa_source: str  # native (af3), mmseqs (boltz/rf3), or none
    msa_path: str  # Path to MSA file (empty for AF3 native)
    status: str = CaseStatus.PENDING.value
    error_message: str = ""

    def to_dict(self) -> dict:
        """Convert to dictionary for CSV export."""
        return asdict(self)

    @property
    def output_dir_name(self) -> str:
        """Generate output directory name for this case.

        Format: {protein_id}_{ligand_ccd_code}, e.g. B6EQJ6_STA6
        """
        return f"{self.protein_id}_{self.ligand_ccd_code}"


def _generate_case_id(
    protein_id: str,
    ligand_id: str,
    model: str,
    msa_source: str,
) -> str:
    """Generate a deterministic case ID based on inputs.

    Args:
        protein_id: Protein identifier
        ligand_id: Ligand identifier
        model: Model name (af3, boltz, rf3)
        msa_source: MSA source (native, mmseqs)

    Returns:
        Short hash-based case ID
    """
    key = f"{protein_id}:{ligand_id}:{model}:{msa_source}"
    return hashlib.sha256(key.encode()).hexdigest()[:12]


def generate_cases(
    proteins: list[ProteinRecord],
    ligands: list[LigandRecord],
    msa_records: list[MSARecord],
    models: list[str],
    missing_msa_proteins: list[str] | None = None,
) -> Iterator[Case]:
    """Generate prediction cases for all combinations.

    Args:
        proteins: List of protein records
        ligands: List of ligand records
        msa_records: List of MSA records (for boltz/rf3)
        models: List of models to run (af3, boltz, rf3)
        missing_msa_proteins: List of protein IDs without MSA

    Yields:
        Case objects for each valid combination
    """
    if missing_msa_proteins is None:
        missing_msa_proteins = []

    # Build MSA lookup
    msa_by_protein = {r.protein_id: r for r in msa_records}

    for protein in proteins:
        for ligand in ligands:
            for model in models:
                model_lower = model.lower()

                # Determine MSA source and path
                if model_lower == "af3":
                    # AF3 uses native MSA generation
                    msa_source = "native"
                    msa_path = ""
                else:
                    # Boltz and RF3 use pre-computed MSA
                    if protein.protein_id in missing_msa_proteins:
                        # Skip this protein for non-AF3 models
                        yield Case(
                            case_id=_generate_case_id(
                                protein.protein_id,
                                ligand.ligand_id,
                                model_lower,
                                "mmseqs",
                            ),
                            protein_id=protein.protein_id,
                            ligand_id=ligand.ligand_id,
                            ligand_ccd_code=ligand.ccd_code,
                            model=model_lower,
                            msa_source="mmseqs",
                            msa_path="",
                            status=CaseStatus.SKIPPED.value,
                            error_message="Missing MSA",
                        )
                        continue

                    msa_record = msa_by_protein.get(protein.protein_id)
                    msa_source = "mmseqs"
                    msa_path = msa_record.msa_path if msa_record else ""

                case_id = _generate_case_id(
                    protein.protein_id,
                    ligand.ligand_id,
                    model_lower,
                    msa_source,
                )

                yield Case(
                    case_id=case_id,
                    protein_id=protein.protein_id,
                    ligand_id=ligand.ligand_id,
                    ligand_ccd_code=ligand.ccd_code,
                    model=model_lower,
                    msa_source=msa_source,
                    msa_path=msa_path,
                    status=CaseStatus.PENDING.value,
                )


def write_cases_manifest(cases: list[Case], output_path: Path) -> None:
    """Write cases manifest to CSV.

    Args:
        cases: List of Case objects
        output_path: Path to output CSV file
    """
    if not cases:
        raise ValueError("No cases to write")

    fieldnames = list(cases[0].to_dict().keys())

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for case in cases:
            writer.writerow(case.to_dict())


def load_cases_manifest(manifest_path: Path) -> list[Case]:
    """Load cases from manifest CSV.

    Args:
        manifest_path: Path to manifest CSV file

    Returns:
        List of Case objects
    """
    cases = []
    with open(manifest_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cases.append(
                Case(
                    case_id=row["case_id"],
                    protein_id=row["protein_id"],
                    ligand_id=row["ligand_id"],
                    ligand_ccd_code=row["ligand_ccd_code"],
                    model=row["model"],
                    msa_source=row["msa_source"],
                    msa_path=row["msa_path"],
                    status=row.get("status", CaseStatus.PENDING.value),
                    error_message=row.get("error_message", ""),
                )
            )
    return cases


def group_cases_by_protein_model(cases: list[Case]) -> dict[tuple[str, str], list[Case]]:
    """Group cases by (protein_id, model) for batch processing.

    Args:
        cases: List of Case objects

    Returns:
        Dictionary mapping (protein_id, model) to list of cases
    """
    grouped: dict[tuple[str, str], list[Case]] = {}

    for case in cases:
        if case.status == CaseStatus.SKIPPED.value:
            continue

        key = (case.protein_id, case.model)
        if key not in grouped:
            grouped[key] = []
        grouped[key].append(case)

    return grouped


def group_cases_by_ligand_model(cases: list[Case]) -> dict[tuple[str, str], list[Case]]:
    """Group cases by (ligand_ccd_code, model) for batch processing.

    One SLURM job per (ligand, model): batches all proteins together.

    Args:
        cases: List of Case objects

    Returns:
        Dictionary mapping (ligand_ccd_code, model) to list of cases
    """
    grouped: dict[tuple[str, str], list[Case]] = {}

    for case in cases:
        if case.status == CaseStatus.SKIPPED.value:
            continue

        key = (case.ligand_ccd_code, case.model)
        if key not in grouped:
            grouped[key] = []
        grouped[key].append(case)

    return grouped
