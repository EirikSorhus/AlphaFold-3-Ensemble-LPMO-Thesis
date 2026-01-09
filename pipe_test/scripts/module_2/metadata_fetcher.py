# metadata_fetcher.py
import requests
import time
import csv
from typing import Dict, List


def fetch_uniprot_metadata_batch(uniprot_ids: List[str], contact_email=None, batch_size: int = 100) -> List[Dict[str, str]]:
    """Fetch UniProt metadata in batches and normalise field names using csv.DictReader."""
    headers = {
        "User-Agent": f"BioPipeline/1.0 ({contact_email})" if contact_email else "BioPipeline/1.0"
    }
    base_url = "https://rest.uniprot.org/uniprotkb/search"

    all_rows: List[Dict[str, str]] = []
    for i in range(0, len(uniprot_ids), batch_size):
        chunk = uniprot_ids[i:i + batch_size]
        if not chunk:
            continue
        query = " OR ".join([f"accession:{uid}" for uid in chunk])
        params = {
            "query": query,
            "fields": "accession,protein_name,gene_names,organism_name,ec,reviewed,xref_pdb,xref_pfam,xref_interpro",
            "format": "tsv",
            "size": len(chunk),
        }

        try:
            r = requests.get(base_url, headers=headers, params=params, timeout=60)
            if not r.ok:
                print(f"[ERROR] UniProt batch {i//batch_size + 1} failed: {r.status_code}")
                continue

            lines = r.text.strip().split("\n")
            if not lines or len(lines) < 2:
                continue

            reader = csv.DictReader(lines, delimiter="\t")
            for raw in reader:
                all_rows.append({
                    "UniProt_ID": raw.get("Entry", raw.get("accession", "")),
                    "Protein_Name": raw.get("Protein names", raw.get("protein_name", "")),
                    "Gene_Name": raw.get("Gene Names", raw.get("gene_names", "")),
                    "Organism": raw.get("Organism", raw.get("organism_name", "")),
                    "EC_Number": raw.get("EC number", raw.get("ec", "")),
                    "Reviewed": raw.get("Reviewed", raw.get("reviewed", "")),
                    "PDB_IDs": raw.get("Cross-reference (PDB)", raw.get("xref_pdb", "")),
                    "Pfam_IDs": raw.get("Cross-reference (Pfam)", raw.get("xref_pfam", "")),
                    "InterPro_IDs": raw.get("Cross-reference (InterPro)", raw.get("xref_interpro", "")),
                    "Match_Status": "matched_uniprot",
                    "Source_DB": "UniProt",
                })
            time.sleep(0.5)
        except Exception as exc:
            print(f"[ERROR] UniProt batch {i//batch_size + 1} exception: {exc}")
            continue

    return all_rows
