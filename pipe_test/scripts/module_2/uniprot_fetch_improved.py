# scripts/module_2/uniprot_fetch_improved.py
"""
UniProt sequence fetching module.

Fetches LPMO sequences from UniProt using the bioservices library.
"""

from typing import Optional, List, Dict
from pathlib import Path
from bioservices import UniProt

from scripts.module_2.config_improved import (
    DATA_RAW_DIR,
    DATA_SEQUENCES_DIR,
    ensure_data_directory,
    build_uniprot_query,
    get_family_output_filename
)


def fetch_uniprot_lpmos(
    query: str,
    output_file: Path
) -> int:
    """
    Fetch LPMO sequences from UniProt based on query.
    
    Args:
        query: UniProt search query
        output_file: Output FASTA file path
        
    Returns:
        int: Number of sequences fetched, or 0 if failed
        
    Raises:
        Exception: If API call fails critically
    """
    # Ensure output directory exists
    ensure_data_directory()
    
    print(f"[UniProt] Starting fetch with query: {query}")
    
    try:
        u = UniProt(verbose=False)
        fasta_text = u.search(query, frmt="fasta")
    except Exception as e:
        print(f"[UniProt] ERROR: API call failed: {e}")
        return 0
    
    # Validate response
    if not fasta_text or not isinstance(fasta_text, str):
        print(f"[UniProt] WARNING: Search returned no data (type: {type(fasta_text)})")
        return 0
    
    if not fasta_text.strip():
        print(f"[UniProt] WARNING: Search returned empty string")
        return 0
    
    # Count sequences (lines starting with '>')
    seq_count = fasta_text.count('\n>')
    if fasta_text.startswith('>'):
        seq_count += 1
    
    if seq_count == 0:
        print(f"[UniProt] WARNING: No sequences found in response")
        return 0
    
    # Write to file
    try:
        with open(output_file, "w", encoding="utf-8") as f:
            f.write(fasta_text)
        print(f"[UniProt] SUCCESS: Saved {seq_count} sequences to {output_file}")
    except IOError as e:
        print(f"[UniProt] ERROR: Failed to write file: {e}")
        return 0
    
    return seq_count


def fetch_uniprot_by_families(
    families: List[str],
    run_id: str,
    output_dir: Optional[Path] = None
) -> Dict[str, int]:
    """
    Fetch UniProt sequences for multiple families, saving each to separate files.
    
    Args:
        families: List of CAZy family names (e.g., ['AA9', 'AA10'])
        run_id: Run ID for filename generation
        output_dir: Output directory. If None, uses DATA_RAW_DIR
        
    Returns:
        dict: Dictionary mapping family name to sequence count
    """
    if output_dir is None:
        output_dir = DATA_SEQUENCES_DIR
    
    results = {}
    
    print(f"[UniProt] Fetching sequences for {len(families)} families")
    print(f"[UniProt] Families: {', '.join(families)}")
    print()
    
    for i, family in enumerate(families, 1):
        print(f"[UniProt] Processing family {i}/{len(families)}: {family}")
        print("-" * 60)
        
        # Build query for single family
        query = build_uniprot_query([family])
        
        # Generate output filename
        filename = get_family_output_filename(family, run_id, 'uniprot')
        output_file = output_dir / filename
        print(f"[UniProt] Output: {output_file.name}")
        
        # Fetch sequences
        count = fetch_uniprot_lpmos(query, output_file)
        results[family] = count
        
        print()
    
    # Summary
    print("=" * 60)
    print("[UniProt] Summary:")
    total = 0
    for family, count in results.items():
        print(f"  {family}: {count} sequences")
        total += count
    print(f"  TOTAL: {total} sequences across {len(families)} families")
    print("=" * 60)
    
    return results


if __name__ == "__main__":
    print("=== UniProt Fetch Module ===")
    from scripts.module_2.config_improved import DEFAULT_CAZY_FAMILIES, get_run_id
    
    # Test with default families
    run_id = get_run_id()
    results = fetch_uniprot_by_families(DEFAULT_CAZY_FAMILIES, run_id)
    
    if any(count > 0 for count in results.values()):
        print(f"\n✓ Successfully fetched sequences")
    else:
        print(f"\n✗ Failed to fetch sequences")
