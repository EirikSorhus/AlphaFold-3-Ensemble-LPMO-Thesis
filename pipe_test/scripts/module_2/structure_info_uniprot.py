import requests
from typing import Optional, Dict, Any

def fetch_uniprot_structure_info(uniprot_id: str) -> Dict[str, Any]:
    """
    Fetch structure info for a UniProt ID from UniProt JSON API.
    Returns dict with structure fields (see pipeline spec).
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        resp = requests.get(url, timeout=20)
        resp.raise_for_status()
        data = resp.json()
    except Exception:
        return {
            "has_experimental_structure": False,
            "pdb_ids": [],
            "best_method": None,
            "best_resolution": None,
            "best_coverage": None,
            "has_alphafold_model": False,
            "alphafold_confidence_summary": None
        }
    # Parse experimental structures (PDB)
    pdb_ids = []
    best_method = None
    best_resolution = None
    best_coverage = None
    if "uniProtKBCrossReferences" in data:
        pdb_entries = [x for x in data["uniProtKBCrossReferences"] if x["database"] == "PDB"]
        for entry in pdb_entries:
            pdb_ids.append(entry["id"])
            # Method, resolution, coverage
            props = {p["key"]: p["value"] for p in entry.get("properties", [])}
            if not best_method and "Method" in props:
                best_method = props["Method"]
            if not best_resolution and "Resolution" in props:
                try:
                    best_resolution = float(props["Resolution"].split()[0])
                except Exception:
                    pass
            if not best_coverage and "Chains" in props:
                # Example: "A=30-250"
                try:
                    chain_range = props["Chains"].split("=")[1]
                    start, end = map(int, chain_range.split("-")[:2])
                    best_coverage = end - start + 1
                except Exception:
                    pass
    # AlphaFold
    has_af = False
    af_conf = None
    if "uniProtKBCrossReferences" in data:
        af_entries = [x for x in data["uniProtKBCrossReferences"] if x["database"] == "AlphaFold DB"]
        if af_entries:
            has_af = True
            # Try to extract confidence summary if present
            af_conf = af_entries[0].get("properties", [{}])[0].get("value")
    return {
        "has_experimental_structure": bool(pdb_ids),
        "pdb_ids": pdb_ids,
        "best_method": best_method,
        "best_resolution": best_resolution,
        "best_coverage": best_coverage,
        "has_alphafold_model": has_af,
        "alphafold_confidence_summary": af_conf
    }
