import requests
from typing import List, Dict, Optional

def map_ncbi_to_uniprot(genbank_ids: List[str]) -> Dict[str, Optional[str]]:
    """
    Map a list of GenBank protein IDs to UniProt IDs using UniProt ID Mapping API.
    Returns a dict: {genbank_id: uniprot_id or None}
    """
    # Step 1: Submit mapping job
    url = "https://rest.uniprot.org/idmapping/run"
    payload = {
        "from": "RefSeq_Protein",
        "to": "UniProtKB",
        "ids": ",".join(genbank_ids)
    }
    resp = requests.post(url, data=payload)
    resp.raise_for_status()
    job_id = resp.json()["jobId"]

    # Step 2: Poll for job completion
    import time
    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    while True:
        status_resp = requests.get(status_url)
        status_resp.raise_for_status()
        status = status_resp.json()
        if status.get("jobStatus") == "FINISHED":
            break
        time.sleep(1)

    # Step 3: Download results
    results_url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/{job_id}"
    results_resp = requests.get(results_url)
    results_resp.raise_for_status()
    results = results_resp.json()
    mapping = {}
    for r in results.get("results", []):
        from_id = r["from"]
        to_id = r["to"]
        mapping[from_id] = to_id
    # Fill in None for unmapped
    for gid in genbank_ids:
        if gid not in mapping:
            mapping[gid] = None
    return mapping
