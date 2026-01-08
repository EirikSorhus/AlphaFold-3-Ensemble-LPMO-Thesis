"""
Enricher module – Data enrichment fra UniProt, NCBI, RCSB PDB, NCBI Taxonomy.

Dette modulet håndterer henting av metadata fra multiple kilder med
robust error-handling, timeout, retry-logikk og safe defaults.
"""

from __future__ import annotations

import re
import time
import json
import xml.etree.ElementTree as ET
from typing import Dict, Optional, List, Tuple, Any
from pathlib import Path
import requests

# ============================================================================
# CONSTANTS & CONFIG
# ============================================================================

REQUEST_TIMEOUT = 30  # seconds
MAX_RETRIES = 3
RETRY_DELAY = 1  # seconds

USER_AGENT = "LPMO-Pipeline/1.0 (eirik.sorhus@nmbu.no)"

# API URLs
UNIPROT_JSON_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"
NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCBI_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_TAXONOMY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
RCSB_PDB_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb}"


# ============================================================================
# HELPER: ROBUST HTTP REQUEST
# ============================================================================

def _http_get(url: str, params: Optional[Dict] = None, max_retries: int = MAX_RETRIES) -> Optional[requests.Response]:
    """
    Robust HTTP GET with retry, timeout, error-safe handling.
    
    Args:
        url: URL to fetch
        params: Query parameters
        max_retries: Number of retries on timeout/5xx
        
    Returns:
        Response object or None on persistent error
    """
    for attempt in range(1, max_retries + 1):
        try:
            resp = requests.get(
                url,
                params=params,
                headers={"User-Agent": USER_AGENT},
                timeout=REQUEST_TIMEOUT
            )
            if resp.status_code == 200:
                return resp
            elif resp.status_code == 404:
                return None  # Not found (safe failure)
            elif 500 <= resp.status_code < 600:
                # Server error; retry
                if attempt < max_retries:
                    time.sleep(RETRY_DELAY * attempt)
                    continue
                return None
            else:
                return None
        except requests.exceptions.Timeout:
            if attempt < max_retries:
                time.sleep(RETRY_DELAY * attempt)
                continue
            return None
        except requests.exceptions.RequestException:
            return None
    return None


# ============================================================================
# 1. UNIPROT ENRICHMENT
# ============================================================================

def enrich_from_uniprot(uniprot_acc: str) -> Dict[str, Any]:
    """
    Enrich from UniProt JSON API.
    
    Returns dict with:
      - reviewed (bool)
      - is_fragment (bool)
      - sequence_length (int)
      - sequence (str)
      - taxonomy_id (int|None)
      - organism_name (str|None)
      - pdb_ids (List[str])
      - has_alphafold_model (bool)
      - alphafold_accession (str|None)
      - source_db (str)
    """
    
    result = {
        "source_db": "UniProt",
        "reviewed": None,
        "is_fragment": None,
        "sequence_length": None,
        "sequence": None,
        "taxonomy_id": None,
        "organism_name": None,
        "pdb_ids": [],
        "has_alphafold_model": False,
        "alphafold_accession": None,
    }
    
    url = UNIPROT_JSON_URL.format(acc=uniprot_acc)
    resp = _http_get(url)
    
    if not resp:
        return result
    
    try:
        data = resp.json()
    except Exception:
        return result
    
    # reviewed: check entryType
    entry_type = data.get("entryType", "")
    result["reviewed"] = "reviewed" in entry_type.lower()
    
    # sequence info
    seq_obj = data.get("sequence", {})
    result["is_fragment"] = bool(seq_obj.get("fragment"))
    result["sequence_length"] = seq_obj.get("length")
    result["sequence"] = seq_obj.get("value")
    
    # organism / taxonomy
    org = data.get("organism", {})
    result["taxonomy_id"] = org.get("taxonId")
    result["organism_name"] = org.get("scientificName")
    
    # cross-references: PDB + AlphaFold
    for xref in (data.get("uniProtKBCrossReferences") or []):
        db = xref.get("database")
        xid = xref.get("id")
        if db == "PDB" and xid:
            result["pdb_ids"].append(xid)
        elif db == "AlphaFoldDB" and xid:
            result["has_alphafold_model"] = True
            result["alphafold_accession"] = xid
    
    result["pdb_ids"] = sorted(set(result["pdb_ids"]))
    
    return result


# ============================================================================
# 2. NCBI PROTEIN ENRICHMENT
# ============================================================================

def enrich_from_ncbi_protein(genbank_acc: str, ncbi_email: str = "eirik.sorhus@nmbu.no") -> Dict[str, Any]:
    """
    Enrich from NCBI Protein database (Entrez).
    
    Returns dict with:
      - source_db (str)
      - sequence_length (int|None)
      - sequence (str|None)
      - organism_name (str|None)
      - taxonomy_id (int|None)
    """
    
    result = {
        "source_db": "NCBI",
        "sequence_length": None,
        "sequence": None,
        "organism_name": None,
        "taxonomy_id": None,
    }
    
    # Fetch FASTA from Entrez
    params = {
        "db": "protein",
        "id": genbank_acc,
        "rettype": "fasta",
        "retmode": "text",
        "email": ncbi_email,
    }
    
    resp = _http_get(NCBI_EFETCH_URL, params=params)
    
    if not resp:
        return result
    
    try:
        lines = resp.text.strip().split("\n")
        if not lines:
            return result
        
        # Parse FASTA
        header = lines[0]  # >gi|...|ref|NP_000001.3|...
        seq = "".join(lines[1:])
        
        result["sequence"] = seq
        result["sequence_length"] = len(seq)
        
        # Try to extract organism from header or fetch via taxonomy
        # Format often: >gi|123|ref|NP_000001.3| Homo sapiens ...
        parts = header.split("|")
        if len(parts) > 4:
            desc = "|".join(parts[4:])
            # Extract organism (first part before [ or space)
            match = re.search(r"^([A-Z][a-z]+(?:\s+[a-z]+)*)", desc)
            if match:
                result["organism_name"] = match.group(1)
    except Exception:
        pass
    
    return result


# ============================================================================
# 3. NCBI TAXONOMY – KINGDOM
# ============================================================================

def get_kingdom_from_taxid(taxid: int | str) -> Optional[str]:
    """
    Fetch kingdom/superkingdom from NCBI Taxonomy.
    
    Args:
        taxid: NCBI taxonomy ID
        
    Returns:
        Kingdom name (Eukaryota, Bacteria, Archaea, Viruses, etc.) or None
    """
    
    if not taxid:
        return None
    
    params = {
        "db": "taxonomy",
        "id": str(taxid),
        "retmode": "xml",
    }
    
    resp = _http_get(NCBI_TAXONOMY_URL, params=params)
    
    if not resp:
        return None
    
    try:
        root = ET.fromstring(resp.text)
        
        # Find LineageEx
        lineage_ex = root.find(".//LineageEx")
        if lineage_ex is None:
            return None
        
        # Look for Rank=superkingdom or Rank=kingdom
        for taxon in lineage_ex.findall("Taxon"):
            rank = taxon.findtext("Rank")
            name = taxon.findtext("ScientificName")
            if rank in ("superkingdom", "kingdom") and name:
                return name
        
    except Exception:
        pass
    
    return None


# ============================================================================
# 4. RCSB PDB – RESOLUTION
# ============================================================================

def get_pdb_resolution(pdb_ids: List[str]) -> Optional[float]:
    """
    Fetch best (minimum) resolution from RCSB PDB.
    
    Args:
        pdb_ids: List of PDB IDs (e.g., ["1ABC", "2DEF"])
        
    Returns:
        Best resolution in Ångström or None
    """
    
    if not pdb_ids:
        return None
    
    resolutions = []
    
    for pdb_id in pdb_ids:
        url = RCSB_PDB_URL.format(pdb=pdb_id)
        resp = _http_get(url)
        
        if not resp:
            continue
        
        try:
            data = resp.json()
            rcsb_info = data.get("rcsb_entry_info", {})
            
            # Try multiple resolution fields
            res_combined = rcsb_info.get("resolution_combined")
            if res_combined:
                if isinstance(res_combined, list):
                    resolutions.extend([r for r in res_combined if isinstance(r, (int, float))])
                elif isinstance(res_combined, (int, float)):
                    resolutions.append(res_combined)
            
            # Fallback: try individual resolution field
            res_single = rcsb_info.get("resolution")
            if res_single and isinstance(res_single, (int, float)):
                resolutions.append(res_single)
        
        except Exception:
            pass
    
    return min(resolutions) if resolutions else None


# ============================================================================
# 5. ORCHESTRATOR – COMBINE SOURCES
# ============================================================================

def enrich_sequence_metadata(
    uniprot_acc: Optional[str] = None,
    genbank_acc: Optional[str] = None,
    taxonomy_id: Optional[int | str] = None,
    organism_name: Optional[str] = None,
    pdb_acc_from_cazy: Optional[str] = None,
    ncbi_email: str = "eirik.sorhus@nmbu.no",
) -> Dict[str, Any]:
    """
    Orchestrate enrichment from multiple sources.
    
    Priority:
      1. If UniProt: enrich_from_uniprot()
      2. Else if GenBank: enrich_from_ncbi_protein()
      3. Fill kingdom via NCBI Taxonomy (taxid)
      4. Fill PDB resolution via RCSB
    
    Args:
        uniprot_acc: UniProt accession (e.g., "P12345")
        genbank_acc: GenBank accession (e.g., "NP_000001.3")
        taxonomy_id: NCBI taxonomy ID for kingdom lookup
        organism_name: Organism name (fallback)
        pdb_acc_from_cazy: PDB ID(s) from CAZy (semicolon-separated)
        ncbi_email: Email for NCBI API
        
    Returns:
        Dict with all enriched fields
    """
    
    result = {
        "source_db": None,
        "reviewed": None,
        "is_fragment": None,
        "sequence_length": None,
        "sequence": None,
        "taxonomy_id": None,
        "kingdom": None,
        "organism_name": None,
        "pdb_ids": [],
        "has_experimental_structure": False,
        "best_pdb_resolution": None,
        "has_alphafold_model": False,
        "alphafold_accession": None,
    }
    
    # Priority 1: UniProt
    if uniprot_acc:
        up_data = enrich_from_uniprot(uniprot_acc)
        result.update(up_data)
        
        # Merge PDB IDs from CAZy if provided
        pdb_from_cazy = []
        if pdb_acc_from_cazy:
            pdb_from_cazy = [p.strip() for p in re.split(r"[;,\s]+", pdb_acc_from_cazy) if p.strip()]
        
        result["pdb_ids"] = sorted(set(result.get("pdb_ids", []) + pdb_from_cazy))
    
    # Priority 2: GenBank (if no UniProt)
    elif genbank_acc:
        nb_data = enrich_from_ncbi_protein(genbank_acc, ncbi_email=ncbi_email)
        result.update(nb_data)
        
        # Add PDB from CAZy
        pdb_from_cazy = []
        if pdb_acc_from_cazy:
            pdb_from_cazy = [p.strip() for p in re.split(r"[;,\s]+", pdb_acc_from_cazy) if p.strip()]
        result["pdb_ids"] = pdb_from_cazy
    
    # Fill kingdom from taxid
    if taxonomy_id and not result.get("kingdom"):
        kingdom = get_kingdom_from_taxid(taxonomy_id)
        result["kingdom"] = kingdom
    
    # Fallback organism name
    if not result.get("organism_name") and organism_name:
        result["organism_name"] = organism_name
    
    # Fallback taxonomy_id
    if not result.get("taxonomy_id") and taxonomy_id:
        result["taxonomy_id"] = taxonomy_id
    
    # PDB resolution
    if result.get("pdb_ids"):
        result["has_experimental_structure"] = True
        result["best_pdb_resolution"] = get_pdb_resolution(result["pdb_ids"])
    
    return result


if __name__ == "__main__":
    # Test example
    print("Testing enricher module...")
    
    # Test UniProt
    print("\n1. UniProt enrichment (P12345):")
    up = enrich_from_uniprot("P12345")
    print(f"   reviewed={up.get('reviewed')}, length={up.get('sequence_length')}")
    
    # Test kingdom
    print("\n2. Kingdom from taxid (9606 = Homo sapiens):")
    kingdom = get_kingdom_from_taxid(9606)
    print(f"   kingdom={kingdom}")
