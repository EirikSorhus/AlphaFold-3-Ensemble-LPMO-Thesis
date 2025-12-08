# scripts/module2/ncbi_fetch.py

# Fetch sequences from NCBI given GenBank IDs
from time import sleep
from typing import List
from Bio import Entrez
from .config import NCBI_EMAIL, CAZY_FASTA, DATA_RAW_DIR
from .cazy_fetch import fetch_all_cazy_genbank_ids
import os

def fetch_ncbi_fasta(genbank_ids: List[str], out_fasta: str):
    os.makedirs(DATA_RAW_DIR, exist_ok=True)
    Entrez.email = NCBI_EMAIL

    print(f"[NCBI] Henter {len(genbank_ids)} sekvenser...")
    with open(out_fasta, "w") as out:
        for i, acc in enumerate(genbank_ids, start=1):
            try:
                handle = Entrez.efetch(
                    db="protein",
                    id=acc,
                    rettype="fasta",
                    retmode="text"
                )
                seq_text = handle.read()
                handle.close()

                if seq_text.strip():
                    out.write(seq_text)
                else:
                    print(f"[NCBI] Tom respons for {acc}")
            except Exception as e:
                print(f"[NCBI] Feil for {acc}: {e}")

            # litt pause for å være snill mot NCBI
            if i % 10 == 0:
                sleep(0.5)

    print(f"[NCBI] Lagret CAZy/NCBI-sekvenser til: {out_fasta}")

def fetch_cazy_lpmos_to_fasta():
    ids = fetch_all_cazy_genbank_ids()
    fetch_ncbi_fasta(ids, CAZY_FASTA)

if __name__ == "__main__":
    fetch_cazy_lpmos_to_fasta()
