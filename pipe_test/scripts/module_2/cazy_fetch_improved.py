# scripts/module_2/cazy_fetch_improved.py
"""
CAZy web scraping module.

Scrapes GenBank accession IDs from CAZy family pages.
"""

import time
from typing import List, Optional, Dict, Tuple
import requests
from bs4 import BeautifulSoup

from scripts.module_2.config_improved import (
    REQUEST_TIMEOUT
)
from scripts.module_2.taxonomy_info import parse_cazy_taxonomy_from_file


BASE_URL = "https://www.cazy.org"


def fetch_cazy_data_for_family(family: str, timeout: int = REQUEST_TIMEOUT) -> Tuple[List[str], Dict[str, Tuple[str, str]]]:
    """
    Fetch GenBank accession IDs and taxonomy info for a single CAZy family.
    
    Args:
        family: CAZy family name (e.g., 'AA9')
        timeout: HTTP request timeout in seconds
        
    Returns:
        tuple: (list of GenBank IDs, dict mapping protein_id → (kingdom, organism))
        
    Note:
        CAZy now provides direct text files at /IMG/cazy_data/{family}.txt
        These contain tab-separated data:
        - Column 4 (index 3): Protein_ID
        - Column 5 (index 4): Source
        - Column 2 (index 1): Kingdom
        - Column 3 (index 2): Organism/Species
    """
    # CAZy now provides direct .txt files with all data
    txt_url = f"{BASE_URL}/IMG/cazy_data/{family}.txt"
    candidate_urls = [txt_url]
    response = None
    
    for url in candidate_urls:
        print(f"[CAZy] Attempting: {url}")
        try:
            response = requests.get(url, timeout=timeout, headers={"User-Agent": "Mozilla/5.0 (LPMO Pipeline)"})
            response.raise_for_status()
            print(f"[CAZy] SUCCESS: Got response {response.status_code} from {url}")
            break  # success
        except requests.exceptions.Timeout:
            print(f"[CAZy] ERROR: Request timeout for {family} (>{timeout}s) at {url}")
            return [], {}
        except requests.exceptions.HTTPError as e:
            print(f"[CAZy] WARNING: HTTP {e.response.status_code} for {url}")
            response = None
        except requests.exceptions.RequestException as e:
            print(f"[CAZy] ERROR: Request failed for {url}: {e}")
            return [], {}

    if response is None:
        print(f"[CAZy] ERROR: All URL patterns failed for {family}")
        return [], {}
    
    print(f"[CAZy] Data file length: {len(response.text)} bytes")
    
    # Parse tab-separated text file
    # Format: Family\tKingdom\tOrganism\tProtein_ID\tSource
    genbank_ids = []
    taxonomy_map = {}  # protein_id → (kingdom, organism)
    
    try:
        lines = response.text.strip().split('\n')
        print(f"[CAZy] Parsing {len(lines)} lines from data file...")
        
        for i, line in enumerate(lines):
            # Skip header line or comments
            if i == 0 or line.startswith('#') or not line.strip():
                continue
            
            parts = line.split('\t')
            if len(parts) >= 5:
                protein_id = parts[3].strip()
                source = parts[4].strip().lower()
                kingdom = parts[1].strip() if len(parts) > 1 else ""
                organism = parts[2].strip() if len(parts) > 2 else ""
                
                # Only take NCBI entries (skip JGI, Uniprot if we want pure GenBank)
                # NCBI GenBank/Protein IDs typically: letters + numbers (e.g., RYN72100.1)
                # or pure numbers representing GI numbers
                if source == 'ncbi' and protein_id:
                    if protein_id not in genbank_ids:
                        genbank_ids.append(protein_id)
                        taxonomy_map[protein_id] = (kingdom, organism)
                # Also accept entries that look like GenBank accessions even if source unclear
                elif protein_id and '.' in protein_id and not source == 'jgi':
                    # Pattern like ABC12345.1 (letters followed by digits and version)
                    if protein_id not in genbank_ids:
                        genbank_ids.append(protein_id)
                        taxonomy_map[protein_id] = (kingdom, organism)
        
        print(f"[CAZy] Extracted {len(genbank_ids)} unique GenBank/NCBI IDs from data file")
        if genbank_ids:
            print(f"[CAZy] Example IDs: {', '.join(genbank_ids[:5])}")
            if genbank_ids[0] in taxonomy_map:
                kingdom, organism = taxonomy_map[genbank_ids[0]]
                print(f"[CAZy] Example taxonomy: {genbank_ids[0]} → Kingdom: {kingdom}, Organism: {organism}")
        
    except Exception as e:
        print(f"[CAZy] ERROR: Failed to parse data file for {family}: {e}")
        return [], {}
    
    print(f"[CAZy] {family}: Found {len(genbank_ids)} GenBank IDs with taxonomy info")
    return genbank_ids, taxonomy_map


def fetch_cazy_genbank_ids_for_family(family: str, timeout: int = REQUEST_TIMEOUT) -> List[str]:
    """
    Fetch GenBank accession IDs for a single CAZy family.
    
    DEPRECATED: Use fetch_cazy_data_for_family() instead to get taxonomy too.
    
    Args:
        family: CAZy family name (e.g., 'AA9')
        timeout: HTTP request timeout in seconds
        
    Returns:
        list: List of GenBank accession IDs found in the CAZy data file
    """
    ids, _ = fetch_cazy_data_for_family(family, timeout)
    return ids


def fetch_cazy_data_by_family(
    families: List[str],
    delay_between_requests: float = 1.0
) -> Tuple[Dict[str, List[str]], Dict[str, Tuple[str, str]]]:
    """
    Fetch GenBank IDs and taxonomy from CAZy families, keeping results per family.
    
    Args:
        families: List of CAZy family names
        delay_between_requests: Seconds to wait between requests
        
    Returns:
        tuple: (
            dict mapping family name to list of GenBank IDs,
            dict mapping protein_id to (kingdom, organism) tuples
        )
    """
    print(f"[CAZy] Fetching GenBank IDs and taxonomy per family from {len(families)} families")
    print(f"[CAZy] Families: {', '.join(families)}")
    print()
    
    results = {}
    all_taxonomy = {}
    
    for i, family in enumerate(families):
        print(f"[CAZy] Processing family {i+1}/{len(families)}: {family}")
        ids, taxonomy = fetch_cazy_data_for_family(family)
        results[family] = ids
        all_taxonomy.update(taxonomy)
        
        # Rate limiting
        if i < len(families) - 1:
            time.sleep(delay_between_requests)
        print()
    
    # Summary
    print("=" * 60)
    print("[CAZy] Summary:")
    total = 0
    for family, ids in results.items():
        print(f"  {family}: {len(ids)} GenBank IDs")
        total += len(ids)
    print(f"  TOTAL: {total} GenBank IDs across {len(families)} families")
    print(f"  Taxonomy entries: {len(all_taxonomy)}")
    print("=" * 60)
    
    return results, all_taxonomy


def fetch_all_cazy_genbank_ids(
    families: List[str],
    delay_between_requests: float = 1.0
) -> List[str]:
    """
    Fetch all GenBank IDs from multiple CAZy families.
    
    Args:
        families: List of CAZy family names
        delay_between_requests: Seconds to wait between requests (be nice to CAZy servers)
        
    Returns:
        list: Sorted list of unique GenBank accession IDs
    """
    print(f"[CAZy] Fetching GenBank IDs from {len(families)} families: {families}")
    
    all_ids = []
    
    for i, family in enumerate(families):
        ids = fetch_cazy_genbank_ids_for_family(family)
        all_ids.extend(ids)
        
        # Rate limiting: wait between requests (except after last one)
        if i < len(families) - 1:
            time.sleep(delay_between_requests)
    
    # Remove duplicates and sort
    unique_ids = sorted(set(all_ids))
    
    duplicates_removed = len(all_ids) - len(unique_ids)
    print(f"[CAZy] Total unique GenBank IDs: {len(unique_ids)}")
    if duplicates_removed > 0:
        print(f"[CAZy] Removed {duplicates_removed} duplicate IDs")
    
    return unique_ids


def fetch_cazy_genbank_ids_by_family(
    families: List[str],
    delay_between_requests: float = 1.0
) -> Dict[str, List[str]]:
    """
    Fetch GenBank IDs from CAZy families, keeping results per family.
    
    DEPRECATED: Use fetch_cazy_data_by_family() instead to get taxonomy too.
    
    Args:
        families: List of CAZy family names
        delay_between_requests: Seconds to wait between requests
        
    Returns:
        dict: Dictionary mapping family name to list of GenBank IDs
    """
    results, _ = fetch_cazy_data_by_family(families, delay_between_requests)
    return results


if __name__ == "__main__":
    print("=== CAZy Fetch Module ===")
    from scripts.module_2.config_improved import DEFAULT_CAZY_FAMILIES
    
    # Test per-family fetching
    results = fetch_cazy_genbank_ids_by_family(DEFAULT_CAZY_FAMILIES)
    
    if results and any(len(ids) > 0 for ids in results.values()):
        print(f"\n✓ Successfully fetched GenBank IDs")
        for family, ids in results.items():
            if ids:
                print(f"  {family}: First 5 IDs: {ids[:5]}")
    else:
        print(f"\n✗ Failed to fetch GenBank IDs")
