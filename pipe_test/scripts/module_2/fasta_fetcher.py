import requests
import time
from pathlib import Path
from typing import Dict, List

UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/search"


def fetch_uniprot_fasta_bulk(uniprot_ids: List[str], output_path: Path, chunk_size: int = 100, contact_email=None) -> Dict[str, object]:
    """Fetch FASTA sequences from UniProt in batches and write to disk."""
    headers = {"User-Agent": f"BioPipeline/1.0 ({contact_email})" if contact_email else "BioPipeline/1.0"}
    output_path = Path(output_path)
    failed: List[str] = []

    with output_path.open("w") as handle:
        for i in range(0, len(uniprot_ids), chunk_size):
            chunk = uniprot_ids[i:i + chunk_size]
            query = " OR ".join(f"accession:{uid}" for uid in chunk)
            params = {"query": query, "format": "fasta", "size": len(chunk)}
            try:
                r = requests.get(UNIPROT_FASTA_URL, params=params, headers=headers, timeout=60)
                r.raise_for_status()
                handle.write(r.text)
            except Exception as exc:
                failed.extend(chunk)
                print(f"[!] UniProt FASTA fetch failed for chunk {i // chunk_size + 1}: {exc}")
            time.sleep(0.5)

    return {"requested": len(uniprot_ids), "written": len(uniprot_ids) - len(failed), "failed_ids": failed}


def fetch_ncbi_fasta_bulk(ncbi_ids: List[str], output_path: Path, contact_email=None, delay: float = 0.35, mode: str = "w") -> Dict[str, object]:
    """Fetch FASTA sequences from NCBI efetch endpoint one-by-one with throttling."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    output_path = Path(output_path)
    failed: List[str] = []

    with output_path.open(mode) as handle:
        for acc in ncbi_ids:
            params = {"db": "protein", "id": acc, "rettype": "fasta", "retmode": "text"}
            if contact_email:
                params["email"] = contact_email
            try:
                r = requests.get(base_url, params=params, timeout=30)
                if r.ok and r.text.startswith(">"):
                    handle.write(r.text if r.text.endswith("\n") else r.text + "\n")
                else:
                    failed.append(acc)
            except Exception as exc:
                failed.append(acc)
                print(f"[!] NCBI FASTA fetch failed for {acc}: {exc}")
            time.sleep(max(delay, 0.3))

    return {"requested": len(ncbi_ids), "written": len(ncbi_ids) - len(failed), "failed_ids": failed}
