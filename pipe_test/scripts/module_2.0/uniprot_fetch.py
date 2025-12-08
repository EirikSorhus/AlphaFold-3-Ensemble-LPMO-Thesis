# scripts/module2/uniprot_fetch.py

# Fetch LPMO sequences from UniProt
import os
from bioservices import UniProt
from .config import UNIPROT_QUERY, UNIPROT_FASTA, DATA_RAW_DIR

def fetch_uniprot_lpmos():
    os.makedirs(DATA_RAW_DIR, exist_ok=True)

    u = UniProt()
    # 'fasta' gir oss sekvensene direkte
    print(f"Kjører UniProt-søk: {UNIPROT_QUERY}")
    fasta_text = u.search(UNIPROT_QUERY, format="fasta")

    with open(UNIPROT_FASTA, "w") as f:
        f.write(fasta_text)

    print(f"[UniProt] Lagret sekvenser til: {UNIPROT_FASTA}")

if __name__ == "__main__":
    fetch_uniprot_lpmos()
