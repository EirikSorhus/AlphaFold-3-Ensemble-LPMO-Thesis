# scripts/module_2/merge_sequences_improved.py
"""
FASTA file merging and deduplication module.

Merges multiple FASTA files and removes duplicate sequences.
"""

from typing import List, Tuple, Optional
from pathlib import Path
from Bio import SeqIO

from scripts.module_2.config_improved import (
    DATA_RAW_DIR,
    DATA_SEQUENCES_DIR,
    ensure_data_directory,
    get_merged_output_filename,
    get_family_output_filename
)


def merge_fastas(
    input_files: List[Path],
    output_file: Path,
    deduplicate_by: str = "sequence"
) -> Tuple[int, int]:
    """
    Merge multiple FASTA files and optionally deduplicate.
    
    Args:
        input_files: List of input FASTA file paths
        output_file: Output FASTA file path
        deduplicate_by: Deduplication strategy:
            - "sequence": Only sequence must match (different IDs with same sequence = duplicate)
            - "id": Only ID must match (same ID with different sequences = duplicate, keeps first)
            - "both": Both ID and sequence must match
            - "none": No deduplication
            
    Returns:
        tuple: (total records before deduplication, unique records after deduplication)
    """
    ensure_data_directory()
    
    print(f"[MERGE] Merging {len(input_files)} FASTA files")
    print(f"[MERGE] Deduplication strategy: {deduplicate_by}")
    
    records = []
    files_successfully_read = 0
    
    # Read all input files
    for filepath in input_files:
        if not filepath.exists():
            print(f"[MERGE] WARNING: File does not exist, skipping: {filepath}")
            continue
        
        try:
            file_records = list(SeqIO.parse(filepath, "fasta"))
            
            if not file_records:
                print(f"[MERGE] WARNING: No sequences found in {filepath.name}")
            else:
                records.extend(file_records)
                files_successfully_read += 1
                print(f"[MERGE] Read {len(file_records)} sequences from {filepath.name}")
                
        except Exception as e:
            print(f"[MERGE] ERROR: Failed to read {filepath}: {e}")
    
    if files_successfully_read == 0:
        print(f"[MERGE] ERROR: No valid input files could be read")
        return 0, 0
    
    total_before = len(records)
    print(f"[MERGE] Total sequences before deduplication: {total_before}")
    
    # Deduplicate based on strategy
    if deduplicate_by == "none":
        unique_records = records
    else:
        seen = {}
        
        for rec in records:
            # Determine key based on strategy
            if deduplicate_by == "sequence":
                key = str(rec.seq)
            elif deduplicate_by == "id":
                key = rec.id
            elif deduplicate_by == "both":
                key = (rec.id, str(rec.seq))
            else:
                raise ValueError(f"Invalid deduplicate_by value: {deduplicate_by}")
            
            # Keep first occurrence
            if key not in seen:
                seen[key] = rec
        
        unique_records = list(seen.values())
    
    duplicates_removed = total_before - len(unique_records)
    print(f"[MERGE] Sequences after deduplication: {len(unique_records)}")
    print(f"[MERGE] Duplicates removed: {duplicates_removed}")
    
    # Write merged file
    try:
        SeqIO.write(unique_records, output_file, "fasta")
        print(f"[MERGE] SUCCESS: Wrote merged FASTA to {output_file}")
    except Exception as e:
        print(f"[MERGE] ERROR: Failed to write output file: {e}")
        return total_before, 0
    
    return total_before, len(unique_records)


def merge_family_files(
    families: List[str],
    run_id: str,
    sources: List[str] = ['uniprot', 'cazy'],
    output_dir: Optional[Path] = None,
    deduplicate_by: str = "sequence"
) -> Tuple[int, int]:
    """
    Merge all family-specific FASTA files into one merged file.
    
    Args:
        families: List of family names to merge
        run_id: Run ID used for the input files
        sources: List of sources to include (e.g., ['uniprot', 'cazy'])
        output_dir: Output directory. If None, uses DATA_RAW_DIR
        deduplicate_by: Deduplication strategy
        
    Returns:
        tuple: (total records, unique records)
    """
    if output_dir is None:
        output_dir = DATA_SEQUENCES_DIR
    
    print(f"[MERGE] Merging files from {len(families)} families and {len(sources)} sources")
    print(f"[MERGE] Families: {', '.join(families)}")
    print(f"[MERGE] Sources: {', '.join(sources)}")
    print()
    
    # Collect all input files
    input_files = []
    for family in families:
        for source in sources:
            filename = get_family_output_filename(family, run_id, source)
            filepath = output_dir / filename
            
            if filepath.exists():
                input_files.append(filepath)
                print(f"[MERGE] Will include: {filename}")
            else:
                print(f"[MERGE] Skipping (not found): {filename}")
    
    if not input_files:
        print("[MERGE] ERROR: No input files found to merge")
        return 0, 0
    
    print(f"\n[MERGE] Total files to merge: {len(input_files)}")
    print()
    
    # Generate output filename
    output_filename = get_merged_output_filename(run_id)
    output_file = output_dir / output_filename
    
    # Merge
    total, unique = merge_fastas(input_files, output_file, deduplicate_by)
    
    return total, unique


if __name__ == "__main__":
    print("=== Merge Sequences Module ===")
    from scripts.module_2.config_improved import DEFAULT_CAZY_FAMILIES, get_run_id
    
    # Test with a run_id (you would normally get this from actual run)
    run_id = "20241202_120000"  # Example run ID
    
    total, unique = merge_family_files(
        families=DEFAULT_CAZY_FAMILIES,
        run_id=run_id
    )
    
    if unique > 0:
        print(f"\n✓ Successfully merged: {total} → {unique} unique sequences")
    else:
        print(f"\n✗ Merge failed or no sequences found")
