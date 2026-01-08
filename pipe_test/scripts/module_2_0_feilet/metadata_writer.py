"""
Metadata Writer – Write enriched CSV and FASTA outputs.

Skriver ren, strukturert CSV med 14 definerte kolonner og FASTA-filer.
"""

from __future__ import annotations

import csv
from typing import List, Dict, Any, Optional
from pathlib import Path


# ============================================================================
# COLUMN DEFINITIONS
# ============================================================================

CSV_COLUMNS = [
    # Identitet / opprinnelse
    "seq_id",
    "family",
    "source_db",
    "taxonomy_id",
    "kingdom",
    "organism_name",
    "reviewed",
    "is_fragment",
    
    # Sekvens
    "sequence_length",
    
    # 3D / struktur
    "has_experimental_structure",
    "pdb_ids",
    "best_pdb_resolution",
    "has_alphafold_model",
    "alphafold_accession",
]


# ============================================================================
# FORMATTING HELPERS
# ============================================================================

def _format_value(val: Any) -> str:
    """Format value for CSV output."""
    
    if val is None:
        return "NA"
    
    if isinstance(val, bool):
        return "Y" if val else "N"
    
    if isinstance(val, (list, tuple)):
        return ";".join(str(v) for v in val if v)
    
    if isinstance(val, float):
        return f"{val:.2f}"
    
    return str(val)


# ============================================================================
# WRITE CSV
# ============================================================================

def write_enriched_csv(
    rows: List[Dict[str, Any]],
    output_path: Path,
) -> int:
    """
    Write enriched metadata to CSV.
    
    Args:
        rows: List of enriched sequence dicts
        output_path: Output CSV path
        
    Returns:
        Number of rows written
    """
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    n_written = 0
    
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        
        for row in rows:
            # Ensure all columns exist
            csv_row = {}
            for col in CSV_COLUMNS:
                val = row.get(col)
                csv_row[col] = _format_value(val)
            
            writer.writerow(csv_row)
            n_written += 1
    
    return n_written


# ============================================================================
# WRITE FASTA
# ============================================================================

def write_sequences_fasta(
    rows: List[Dict[str, Any]],
    output_path: Path,
) -> int:
    """
    Write sequences to FASTA file.
    
    Args:
        rows: List of enriched sequence dicts (must have "seq_id", "family", "sequence")
        output_path: Output FASTA path
        
    Returns:
        Number of sequences written
    """
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    n_written = 0
    
    with open(output_path, "w", encoding="utf-8") as f:
        for row in rows:
            seq_id = row.get("seq_id")
            family = row.get("family")
            seq = row.get("sequence")
            
            if not seq_id or not seq:
                continue
            
            # Header: >seq_id family
            header = f">{seq_id}"
            if family:
                header += f" {family}"
            
            f.write(header + "\n")
            
            # Sequence (70 chars per line)
            for i in range(0, len(seq), 70):
                f.write(seq[i:i+70] + "\n")
            
            n_written += 1
    
    return n_written


# ============================================================================
# VALIDATE
# ============================================================================

def validate_rows(rows: List[Dict[str, Any]]) -> List[str]:
    """
    Validate rows for common issues.
    
    Returns:
        List of warning messages
    """
    
    warnings = []
    
    for idx, row in enumerate(rows, start=1):
        seq_id = row.get("seq_id")
        seq = row.get("sequence")
        seq_len = row.get("sequence_length")
        
        if not seq_id:
            warnings.append(f"Row {idx}: missing seq_id")
        
        if not seq:
            warnings.append(f"Row {idx} ({seq_id}): missing sequence")
        
        if seq and seq_len:
            if len(seq) != seq_len:
                warnings.append(
                    f"Row {idx} ({seq_id}): sequence_length mismatch "
                    f"(declared={seq_len}, actual={len(seq)})"
                )
    
    return warnings


if __name__ == "__main__":
    print("Testing metadata_writer module...")
    
    # Example row
    test_rows = [
        {
            "seq_id": "P12345",
            "family": "AA13",
            "source_db": "UniProt",
            "taxonomy_id": 5596,
            "kingdom": "Eukaryota",
            "organism_name": "Alternaria alternata",
            "reviewed": True,
            "is_fragment": False,
            "sequence_length": 347,
            "has_experimental_structure": True,
            "pdb_ids": ["1ABC", "2DEF"],
            "best_pdb_resolution": 1.85,
            "has_alphafold_model": True,
            "alphafold_accession": "AF-P12345-F1",
            "sequence": "MRAISALAGAGAG" * 26 + "MR",
        }
    ]
    
    # Test CSV write
    print("\n1. Writing test CSV...")
    csv_path = Path("/tmp/test_metadata.csv")
    n = write_enriched_csv(test_rows, csv_path)
    print(f"   Wrote {n} rows to {csv_path}")
    
    # Test FASTA write
    print("\n2. Writing test FASTA...")
    fasta_path = Path("/tmp/test_sequences.fasta")
    n = write_sequences_fasta(test_rows, fasta_path)
    print(f"   Wrote {n} sequences to {fasta_path}")
    
    # Test validation
    print("\n3. Validation...")
    warnings = validate_rows(test_rows)
    if warnings:
        for w in warnings:
            print(f"   ⚠ {w}")
    else:
        print("   ✓ All validations passed")
