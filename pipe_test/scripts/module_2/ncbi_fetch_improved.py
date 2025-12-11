# scripts/module_2/ncbi_fetch_improved.py
"""
NCBI sequence fetching module.

Fetches protein sequences from NCBI Protein database using Biopython's Entrez.
"""

import time
from typing import List, Tuple, Optional, Dict
from pathlib import Path
from Bio import Entrez

from scripts.module_2.config_improved import (
    get_ncbi_email,
    DATA_RAW_DIR,
    DATA_SEQUENCES_DIR,
    NCBI_BATCH_SIZE,
    NCBI_REQUESTS_PER_SECOND,
    ensure_data_directory,
    get_family_output_filename
)
from scripts.module_2.cazy_fetch_improved import (
    fetch_all_cazy_genbank_ids,
    fetch_cazy_genbank_ids_by_family,
    fetch_cazy_data_by_family
)

# Global store for taxonomy info from CAZy (used by run_module2_improved.py)
_cazy_taxonomy_cache: Dict[str, Tuple[str, str]] = {}


def fetch_ncbi_fasta(
    genbank_ids: List[str],
    output_file: Path,
    ncbi_email: Optional[str] = None,
    batch_size: int = NCBI_BATCH_SIZE
) -> Tuple[int, int]:
    """
    Fetch FASTA sequences from NCBI Protein database.
    
    Args:
        genbank_ids: List of GenBank accession IDs to fetch
        output_file: Output FASTA file path
        ncbi_email: Email address for NCBI API. If None, uses get_ncbi_email()
        batch_size: Number of sequences to fetch per API call
        
    Returns:
        tuple: (number of successful sequences, number of failed sequences)
        
    Raises:
        ValueError: If NCBI email is not properly configured
    """
    
    # Ensure output directory exists
    ensure_data_directory()
    
    # Get and validate NCBI email
    try:
        email = get_ncbi_email(ncbi_email)
        Entrez.email = email
        print(f"[NCBI] Using email: {email}")
    except ValueError as e:
        print(f"[NCBI] ERROR: {e}")
        raise
    
    if not genbank_ids:
        print("[NCBI] WARNING: No GenBank IDs provided")
        return 0, 0
    
    print(f"[NCBI] Fetching {len(genbank_ids)} sequences in batches of {batch_size}")
    
    success_count = 0
    fail_count = 0
    
    # Calculate delay between requests for rate limiting
    # NCBI allows max 3 requests/second without API key
    delay = 1.0 / NCBI_REQUESTS_PER_SECOND
    
    try:
        with open(output_file, "w", encoding="utf-8") as outfile:
            # Process in batches
            total_batches = (len(genbank_ids) + batch_size - 1) // batch_size
            
            for batch_num in range(total_batches):
                start_idx = batch_num * batch_size
                end_idx = min(start_idx + batch_size, len(genbank_ids))
                batch = genbank_ids[start_idx:end_idx]
                
                batch_ids_str = ",".join(batch)
                
                print(f"[NCBI] Batch {batch_num + 1}/{total_batches}: "
                      f"Fetching {len(batch)} sequences...")
                
                try:
                    handle = Entrez.efetch(
                        db="protein",
                        id=batch_ids_str,
                        rettype="fasta",
                        retmode="text"
                    )
                    fasta_text = handle.read()
                    handle.close()
                    
                    if fasta_text and fasta_text.strip():
                        outfile.write(fasta_text)
                        
                        # Count sequences in this batch
                        batch_seq_count = fasta_text.count('\n>')
                        if fasta_text.startswith('>'):
                            batch_seq_count += 1
                        
                        success_count += batch_seq_count
                        
                        # If we got fewer sequences than requested, some failed
                        if batch_seq_count < len(batch):
                            fail_count += len(batch) - batch_seq_count
                    else:
                        print(f"[NCBI] WARNING: Empty response for batch {batch_num + 1}")
                        fail_count += len(batch)
                        
                except Exception as e:
                    print(f"[NCBI] ERROR: Batch {batch_num + 1} failed: {e}")
                    fail_count += len(batch)
                
                # Rate limiting: wait between requests (except after last batch)
                if batch_num < total_batches - 1:
                    time.sleep(delay)
        
    except IOError as e:
        print(f"[NCBI] ERROR: Failed to write output file: {e}")
        return success_count, fail_count + len(genbank_ids) - success_count
    
    print(f"[NCBI] SUCCESS: Saved {success_count} sequences to {output_file}")
    
    if fail_count > 0:
        print(f"[NCBI] WARNING: {fail_count} sequences could not be fetched")
    
    return success_count, fail_count


def fetch_ncbi_by_families(
    families: List[str],
    run_id: str,
    ncbi_email: Optional[str] = None,
    output_dir: Optional[Path] = None
) -> Dict[str, Tuple[int, int]]:
    """
    Fetch NCBI sequences for multiple families, saving each to separate files.
    Also extracts and caches taxonomy information from CAZy.
    
    Args:
        families: List of CAZy family names
        run_id: Run ID for filename generation
        ncbi_email: Email address for NCBI API
        output_dir: Output directory. If None, uses DATA_RAW_DIR
        
    Returns:
        dict: Dictionary mapping family name to (success_count, fail_count) tuple
    """
    if output_dir is None:
        output_dir = DATA_SEQUENCES_DIR
    
    print(f"[NCBI] Fetching sequences for {len(families)} families")
    print(f"[NCBI] Families: {', '.join(families)}")
    print()
    
    # Step 1: Fetch GenBank IDs per family from CAZy (with taxonomy)
    print("[NCBI] Step 1: Fetching GenBank IDs and taxonomy from CAZy per family...")
    family_ids, all_taxonomy = fetch_cazy_data_by_family(families)
    
    # Cache taxonomy for later use in metadata generation
    global _cazy_taxonomy_cache
    _cazy_taxonomy_cache.update(all_taxonomy)
    print(f"[NCBI] Cached {len(_cazy_taxonomy_cache)} taxonomy entries from CAZy")
    print()
    
    # Step 2: Fetch sequences per family from NCBI
    print("[NCBI] Step 2: Fetching sequences from NCBI per family...")
    results = {}
    
    for i, family in enumerate(families, 1):
        genbank_ids = family_ids.get(family, [])
        
        if not genbank_ids:
            print(f"[NCBI] Family {family}: No GenBank IDs to fetch")
            results[family] = (0, 0)
            continue
        
        print(f"[NCBI] Processing family {i}/{len(families)}: {family}")
        print(f"[NCBI] Fetching {len(genbank_ids)} sequences from CAZy/NCBI...")
        print("-" * 60)
        
        # Generate output filename
        filename = get_family_output_filename(family, run_id, 'ncbi')
        output_file = output_dir / filename
        print(f"[NCBI] Output: {output_file.name}")
        
        # Fetch sequences
        success, failed = fetch_ncbi_fasta(
            genbank_ids,
            output_file,
            ncbi_email=ncbi_email
        )
        
        results[family] = (success, failed)
        print()
    
    # Summary
    print("=" * 60)
    print("[NCBI] Summary:")
    total_success = 0
    total_failed = 0
    for family, (success, failed) in results.items():
        print(f"  {family}: {success} sequences ({failed} failed)")
        total_success += success
        total_failed += failed
    print(f"  TOTAL: {total_success} sequences ({total_failed} failed)")
    print("=" * 60)
    
    return results


def get_cazy_taxonomy_cache() -> Dict[str, Tuple[str, str]]:
    """
    Get the cached taxonomy information from CAZy (populated by fetch_ncbi_by_families).
    
    Returns:
        dict: Mapping of protein_id → (kingdom, organism) tuples
        
    Note:
        This cache is populated during fetch_ncbi_by_families() execution.
        For each CAZy family, taxonomy info is extracted and stored here.
        Keys follow the naming convention:
        - "organism" column should contain Kingdom (Eukaryota, Bacteria, etc.)
        - "species" column should contain Organism/Species (Alternaria alternata, etc.)
    """
    return _cazy_taxonomy_cache


if __name__ == "__main__":
    print("=== NCBI Fetch Module ===")
    from scripts.module_2.config_improved import DEFAULT_CAZY_FAMILIES, get_run_id
    
    try:
        run_id = get_run_id()
        results = fetch_ncbi_by_families(DEFAULT_CAZY_FAMILIES, run_id)
        
        if any(success > 0 for success, _ in results.values()):
            print(f"\n✓ Successfully fetched sequences")
        else:
            print(f"\n✗ Failed to fetch sequences")
            
    except ValueError as e:
        print(f"\n✗ Configuration error: {e}")
