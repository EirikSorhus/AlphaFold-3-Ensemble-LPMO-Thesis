##### FILE: scripts/module_2/metadata_from_fasta.py #####
import argparse
import hashlib
import json
import csv
import requests
import re
import time
from pathlib import Path
from datetime import datetime
from xml.etree import ElementTree

# === Configuration ===
BASE_DIR = Path(".")
METADATA_DIR = BASE_DIR / "data" / "metadata"
RUN_DIR = BASE_DIR / "data" / "run"
SEQUENCES_DIR = BASE_DIR / "data" / "sequences"
LOGS_DIR = BASE_DIR / "logs"

for d in [METADATA_DIR, RUN_DIR, SEQUENCES_DIR, LOGS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

def timestamp():
    return datetime.now().strftime("%Y%m%d_%H%M%S")

# === I/O Functions ===
def read_fasta(filepath):
    with open(filepath) as f:
        seqs = {}
        current_header = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_header = line[1:]
                seqs[current_header] = ""
            elif current_header:
                seqs[current_header] += line
        return seqs

def write_json(data, path):
    with open(path, "w") as f:
        json.dump(data, f, indent=2)

def write_tsv(rows, fieldnames, out_path):
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

def seq_hash(sequence):
    return hashlib.md5(sequence.encode()).hexdigest()

# === Helper: Extract ID from Header ===
def extract_id_from_header(header):
    """
    Attempts to parse UniProt ID from header using Regex.
    Format example: >sp|P12345|ID_NAME
    Returns ID string (e.g., P12345) or None.
    """
    match = re.search(r"\|([A-Z0-9]+)\|", header)
    if match:
        return match.group(1)
    return None

# === UniProt Metadata Lookup ===
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{}?format=json"

def fetch_uniprot_metadata(uniprot_id):
    """Fetches detailed metadata for a known UniProt ID."""
    print(f"   [API] Fetching metadata for: {uniprot_id}...")
    try:
        r = requests.get(UNIPROT_ENTRY_URL.format(uniprot_id), timeout=10)
        if not r.ok:
            print(f"   [API] Error: {r.status_code}")
            return {}
        data = r.json()

        # 1. Protein Name
        protein_name = data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
        if not protein_name:
             if "submissionNames" in data.get("proteinDescription", {}):
                 protein_name = data["proteinDescription"]["submissionNames"][0]["fullName"]["value"]
        
        # 2. Gene Names
        gene_names = ",".join([g["value"] for g in data.get("genes", []) if g.get("value")])
        
        # 3. Organism
        organism = data.get("organism", {}).get("scientificName", "")
        
        # 4. EC Number
        ec_numbers = [x["value"] for x in data.get("proteinDescription", {}).get("recommendedName", {}).get("ecNumbers", [])]
        
        # 5. Reviewed Status (Swiss-Prot vs TrEMBL)
        reviewed = data.get("entryType", "") == "Swiss-Prot"

        # 6. External IDs (PDB, Pfam, InterPro)
        pdb_ids = []
        pfam_ids = []
        interpro_ids = []

        for xref in data.get("uniProtKBCrossReferences", []):
            if xref["database"] == "PDB":
                pdb_ids.append(xref["id"])
            elif xref["database"] == "Pfam":
                pfam_ids.append(xref["id"])
            elif xref["database"] == "InterPro":
                interpro_ids.append(xref["id"])

        return {
            "Protein_Name": protein_name,
            "Gene_Name": gene_names,
            "Organism": organism,
            "EC_Number": ";".join(ec_numbers) if ec_numbers else "",
            "Reviewed": str(reviewed),
            "PDB_IDs": ";".join(pdb_ids),
            "Pfam_IDs": ";".join(pfam_ids),
            "InterPro_IDs": ";".join(interpro_ids)
        }
    except Exception as e:
        print(f"   [ERROR] Failed to fetch metadata for {uniprot_id}: {e}")
        return {}

# === BLAST Fallback (NCBI API) ===
NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

def run_blast_search(sequence):
    """
    Submits sequence to NCBI BLAST (blastp) against 'swissprot'.
    Returns ID only if 100% Identity and covers full query length.
    """
    print("   [BLAST] ID missing. Running BLAST search against 'swissprot'...")
    
    query_len = len(sequence)
    
    # 1. Submit BLAST job
    # We use 'swissprot' because 'nr' returns GenBank IDs (WP_/XP_) which do not work with UniProt API.
    params = {
        "CMD": "Put",
        "PROGRAM": "blastp",
        "DATABASE": "swissprot", 
        "QUERY": sequence
    }
    
    try:
        response = requests.post(NCBI_BLAST_URL, data=params)
        if not response.ok:
            return None
        
        rid_match = re.search(r"RID = (.*)", response.text)
        if not rid_match:
            print("   [BLAST] Could not retrieve RID.")
            return None
        rid = rid_match.group(1).strip()
        
        # 2. Wait for results
        print(f"   [BLAST] Job submitted. RID: {rid}. Waiting...", end="", flush=True)
        while True:
            time.sleep(5) 
            check_params = {"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid}
            check_r = requests.get(NCBI_BLAST_URL, params=check_params)
            if "Status=WAITING" in check_r.text:
                print(".", end="", flush=True)
                continue
            elif "Status=FAILED" in check_r.text or "Status=UNKNOWN" in check_r.text:
                print(" Failed.")
                return None
            elif "Status=READY" in check_r.text:
                print(" Ready!")
                break
        
        # 3. Retrieve and parse results
        result_params = {"CMD": "Get", "FORMAT_TYPE": "XML", "RID": rid}
        result_r = requests.get(NCBI_BLAST_URL, params=result_params)
        
        root = ElementTree.fromstring(result_r.content)
        
        # Iterate through hits to find an EXACT match
        for hit in root.findall(".//Hit"):
            # Check the first HSP (High-scoring Segment Pair) of the hit
            hsp = hit.find("Hit_hsps/Hsp")
            if hsp is None: 
                continue

            identity = int(hsp.find("Hsp_identity").text)
            align_len = int(hsp.find("Hsp_align-len").text)
            
            # Exact Match Logic:
            # Identity must equal alignment length (100% identity)
            # Alignment length must equal Query length (100% coverage)
            if identity == align_len and align_len == query_len:
                accession = hit.find("Hit_accession").text
                # Clean up accession if needed (sometimes NCBI returns 'sp|ID|NAME', we just want ID)
                # Usually Hit_accession in XML is just the ID (e.g. "P12345")
                print(f"   [BLAST] Found EXACT UniProt match: {accession}")
                return accession
        
        print("   [BLAST] No exact UniProt match (100% Identity + 100% Coverage) found.")
        return None
            
    except Exception as e:
        print(f"\n   [ERROR] BLAST failed: {e}")
        return None

# === CLI ===
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-fasta", required=True, help="Path to FASTA file")
    args = parser.parse_args()

    run_id = timestamp()
    fasta_seqs = read_fasta(args.input_fasta)
    
    metadata_rows = []
    unmatched_count = 0
    
    total_seqs = len(fasta_seqs)
    print(f"--- Starting processing of {total_seqs} sequences ---")

    for idx, (header, seq) in enumerate(fasta_seqs.items(), 1):
        print(f"Processing {idx}/{total_seqs}: {header[:40]}...")
        
        sid = seq_hash(seq)
        
        # 1. Try to get ID from Header (Regex)
        uniprot_id = extract_id_from_header(header)
        source_method = "header_extract"
        
        # 2. If no ID, try BLAST (NCBI SwissProt)
        if not uniprot_id:
            uniprot_id = run_blast_search(seq)
            source_method = "blast_search" if uniprot_id else "failed"

        # 3. Construct Row
        row = {
            "Input_Header": header,
            "Seq_Hash": sid,
            "Match_Status": "matched" if uniprot_id else "not_matched",
            "UniProt_ID": uniprot_id or "NA",
            "Source_DB": source_method
        }

        # 4. Fetch Metadata (For ALL valid IDs found)
        if uniprot_id:
            meta = fetch_uniprot_metadata(uniprot_id)
            if meta:
                row.update(meta)
            else:
                row["Match_Status"] = "id_found_metadata_failed"
        else:
            unmatched_count += 1

        metadata_rows.append(row)
        # Sleep to respect API rate limits
        time.sleep(0.5)

    # === Write Output ===
    out_path = METADATA_DIR / f"metadata_expanded_{run_id}.tsv"
    
    fieldnames = [
        "Input_Header", "Seq_Hash", "Match_Status", "UniProt_ID", "Source_DB",
        "Protein_Name", "Gene_Name", "Organism", "EC_Number", 
        "Reviewed", "PDB_IDs", "Pfam_IDs", "InterPro_IDs"
    ]
    
    # Ensure all rows have all fields
    for r in metadata_rows:
        for f in fieldnames:
            if f not in r:
                r[f] = ""

    write_tsv(metadata_rows, fieldnames, out_path)

    run_meta = {
        "run_id": run_id,
        "input_file": str(Path(args.input_fasta).resolve()),
        "total": len(metadata_rows),
        "matched_uniprot": len(metadata_rows) - unmatched_count,
        "unmatched": unmatched_count,
        "output_file": str(out_path.resolve())
    }

    run_path = RUN_DIR / f"run_metadata_{run_id}.json"
    write_json(run_meta, run_path)

    print(f"\n[✓] Metadata saved: {out_path}")
    print(f"[✓] Run metadata saved: {run_path}")