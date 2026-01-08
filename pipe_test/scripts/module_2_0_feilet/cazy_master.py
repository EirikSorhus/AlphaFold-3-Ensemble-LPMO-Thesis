"""
CAZy Master data retrieval – Fetch family table from CAZy .txt files.

Henter alle proteiner for en familie (f.eks. AA13) fra CAZy sine
tab-separerte .txt-filer og returnerer strukturert master-liste.
"""

from __future__ import annotations

import re
from typing import List, Dict, Optional, Any
import requests

# ============================================================================
# CONSTANTS
# ============================================================================

CAZY_BASE_URL = "https://www.cazy.org"
USER_AGENT = "LPMO-Pipeline/1.0 (eirik.sorhus@nmbu.no)"
REQUEST_TIMEOUT = 30


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def get_cazy_family_table(family: str) -> List[Dict[str, Any]]:
    """
    Fetch and parse CAZy family table from .txt file.
    
    Args:
        family: CAZy family name (e.g., "AA13", "GH5")
        
    Returns:
        List of dicts, each representing one protein entry:
        {
            "family": "AA13",
            "taxonomy_id": "5596",
            "organism_name": "Alternaria alternata",
            "protein_name": "laccase",
            "genbank_acc": "RYN72100.1",
            "uniprot_acc": "P12345",
            "pdb_acc": "1ABC;2DEF",
            "source_db": "ncbi"
        }
    """
    
    rows = []
    
    # Try .txt file
    txt_url = f"{CAZY_BASE_URL}/IMG/cazy_data/{family}.txt"
    rows = _parse_cazy_txt(txt_url, family)
    
    if rows:
        return rows
    
    # Fallback: Try HTML (legacy)
    html_url = f"{CAZY_BASE_URL}/{family}.html"
    rows = _parse_cazy_html(html_url, family)
    
    return rows


# ============================================================================
# PARSE TXT FILE (preferred)
# ============================================================================

def _parse_cazy_txt(url: str, family: str) -> List[Dict[str, Any]]:
    """
    Parse CAZy .txt file (tab-separated).
    
    Format:
      Family | Kingdom | Organism | Protein_ID | Source
      AA13   | Eukaryota | Alternaria alternata | RYN72100.1 | ncbi
    """
    
    rows = []
    
    try:
        resp = requests.get(url, headers={"User-Agent": USER_AGENT}, timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
    except Exception:
        return rows
    
    lines = resp.text.strip().split("\n")
    
    for line_num, line in enumerate(lines, start=1):
        # Skip empty / comment lines
        if not line.strip() or line.startswith("#"):
            continue
        
        # Skip header line (first line typically)
        if line_num == 1:
            continue
        
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        
        family_val = parts[0].strip()
        kingdom = parts[1].strip()
        organism = parts[2].strip()
        protein_id = parts[3].strip()
        source_db = parts[4].strip().lower()
        
        # Skip invalid entries
        if not protein_id:
            continue
        
        # Filter: only NCBI or accessions with dots (GenBank-like)
        if source_db == "jgi":
            continue
        
        if source_db != "ncbi" and "." not in protein_id:
            continue
        
        row = {
            "family": family_val or family,
            "taxonomy_id": None,  # Will be filled later
            "kingdom": kingdom or None,
            "organism_name": organism or None,
            "protein_name": None,  # Not in .txt
            "genbank_acc": protein_id if source_db == "ncbi" else None,
            "uniprot_acc": protein_id if source_db == "uniprot" else None,
            "pdb_acc": None,  # Not in .txt (would need separate lookup)
            "source_db": source_db,
        }
        
        rows.append(row)
    
    return rows


# ============================================================================
# PARSE HTML FILE (fallback / legacy)
# ============================================================================

def _parse_cazy_html(url: str, family: str) -> List[Dict[str, Any]]:
    """
    Parse CAZy HTML family page (fallback if .txt unavailable).
    
    Note: This is a fallback; .txt is preferred.
    HTML structure may change, so this is not reliable.
    """
    
    rows = []
    
    try:
        from bs4 import BeautifulSoup
    except ImportError:
        return rows
    
    try:
        resp = requests.get(url, headers={"User-Agent": USER_AGENT}, timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
    except Exception:
        return rows
    
    try:
        soup = BeautifulSoup(resp.text, "html.parser")
    except Exception:
        return rows
    
    # Find main data table
    for table in soup.find_all("table"):
        for tr in table.find_all("tr"):
            tds = [td.get_text(" ", strip=True) for td in tr.find_all(["td"])]
            
            if len(tds) < 6:
                continue
            
            # Expected columns: Family, Tax_id, Organism, Protein, GenBank, UniProt, PDB, ...
            # This is approximate; structure varies
            
            if not re.match(r"^\d+$", tds[1]):
                continue  # Skip non-numeric taxid
            
            row = {
                "family": tds[0] or family,
                "taxonomy_id": int(tds[1]) if tds[1].isdigit() else None,
                "kingdom": None,
                "organism_name": tds[2] or None,
                "protein_name": tds[3] or None,
                "genbank_acc": tds[4] if tds[4] and tds[4] != "-" else None,
                "uniprot_acc": tds[5] if len(tds) > 5 and tds[5] and tds[5] != "-" else None,
                "pdb_acc": tds[6] if len(tds) > 6 and tds[6] and tds[6] != "-" else None,
                "source_db": "cazy_html",
            }
            
            rows.append(row)
    
    return rows


if __name__ == "__main__":
    print("Testing CAZy master retrieval...")
    
    print("\n1. Fetch AA13 from CAZy:")
    rows = get_cazy_family_table("AA13")
    print(f"   Found {len(rows)} entries")
    if rows:
        print(f"   First entry: {rows[0]}")
