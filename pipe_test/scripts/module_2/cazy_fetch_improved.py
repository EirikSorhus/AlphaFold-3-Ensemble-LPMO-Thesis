# scripts/module_2/cazy_fetch_improved.py
"""
CAZy web scraping module.

Scrapes GenBank accession IDs from CAZy family pages.
"""

import time
from typing import List, Optional, Dict
import requests
from bs4 import BeautifulSoup

from scripts.module_2.config_improved import (
    REQUEST_TIMEOUT
)


BASE_URL = "https://www.cazy.org"


def fetch_cazy_genbank_ids_for_family(family: str, timeout: int = REQUEST_TIMEOUT) -> List[str]:
    """
    Fetch GenBank accession IDs for a single CAZy family.
    
    Args:
        family: CAZy family name (e.g., 'AA9')
        timeout: HTTP request timeout in seconds
        
    Returns:
        list: List of GenBank accession IDs found on the page
        
    Note:
        This function uses web scraping and may break if CAZy changes their HTML structure.
        Consider checking CAZy for API access or alternative data sources.
    """
    url = f"{BASE_URL}/{family}_family.html"
    print(f"[CAZy] Fetching: {url}")
    
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
    except requests.exceptions.Timeout:
        print(f"[CAZy] ERROR: Request timeout for {family} (>{timeout}s)")
        return []
    except requests.exceptions.HTTPError as e:
        print(f"[CAZy] ERROR: HTTP error for {family}: {e}")
        return []
    except requests.exceptions.RequestException as e:
        print(f"[CAZy] ERROR: Request failed for {family}: {e}")
        return []
    
    # Parse HTML
    try:
        soup = BeautifulSoup(response.text, "html.parser")
    except Exception as e:
        print(f"[CAZy] ERROR: Failed to parse HTML for {family}: {e}")
        return []
    
    genbank_ids = []
    
    # Find all links pointing to NCBI Protein database
    for link in soup.find_all("a", href=True):
        href = link["href"]
        if "ncbi.nlm.nih.gov/protein" in href:
            acc = link.text.strip()
            if acc and acc not in genbank_ids:  # Avoid immediate duplicates
                genbank_ids.append(acc)
    
    print(f"[CAZy] {family}: Found {len(genbank_ids)} GenBank IDs")
    return genbank_ids


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
    
    Args:
        families: List of CAZy family names
        delay_between_requests: Seconds to wait between requests
        
    Returns:
        dict: Dictionary mapping family name to list of GenBank IDs
    """
    print(f"[CAZy] Fetching GenBank IDs per family from {len(families)} families")
    print(f"[CAZy] Families: {', '.join(families)}")
    print()
    
    results = {}
    
    for i, family in enumerate(families):
        print(f"[CAZy] Processing family {i+1}/{len(families)}: {family}")
        ids = fetch_cazy_genbank_ids_for_family(family)
        results[family] = ids
        
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
    print("=" * 60)
    
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
