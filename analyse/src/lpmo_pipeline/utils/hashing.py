# src/lpmo_pipeline/utils/hashing.py
"""
Responsibility: Deterministic hashing of configs, inputs, and artifacts.
Input:  File paths or dicts
Output: SHA-256 hex digests (str)
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


def hash_file(path: Path) -> str:
    """Return SHA-256 hex digest of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def hash_dict(d: dict[str, Any]) -> str:
    """Deterministic SHA-256 of a dict (sorted-key JSON serialization)."""
    blob = json.dumps(d, sort_keys=True, default=str).encode()
    return hashlib.sha256(blob).hexdigest()


def hash_config(config_path: Path) -> str:
    """Hash a YAML/JSON config file.

    Reads as bytes to be encoding-agnostic.
    """
    return hash_file(config_path)


def compute_input_checksums(
    config_path: Path,
    protein_paths: list[Path],
    ligand_paths: list[Path],
) -> dict[str, str]:
    """Compute checksums for all pipeline inputs.

    Returns:
        {"config_yaml": "sha256:...", "protein_sequences": "sha256:...", ...}
    """
    checksums: dict[str, str] = {}
    checksums["config_yaml"] = f"sha256:{hash_file(config_path)}"

    # Protein: hash all fasta/pdb files concatenated (sorted by name for determinism)
    h_prot = hashlib.sha256()
    for p in sorted(protein_paths):
        h_prot.update(p.read_bytes())
    checksums["protein_sequences"] = f"sha256:{h_prot.hexdigest()}"

    # Ligands
    h_lig = hashlib.sha256()
    for p in sorted(ligand_paths):
        h_lig.update(p.read_bytes())
    checksums["ligand_inputs"] = f"sha256:{h_lig.hexdigest()}"

    return checksums
