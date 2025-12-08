# scripts/module2/run_module2.py

# Run the full pipeline for Module 2: Fetching LPMO sequences
from .uniprot_fetch import fetch_uniprot_lpmos
from .ncbi_fetch import fetch_cazy_lpmos_to_fasta
from .merge_sequences import merge_uniprot_and_cazy

def run_module2_pipeline():
    print("=== Modul 2: Hente LPMO-sekvenser ===")
    print("Steg 1: UniProt")
    fetch_uniprot_lpmos()

    print("\nSteg 2: CAZy + NCBI")
    fetch_cazy_lpmos_to_fasta()

    print("\nSteg 3: Merge UniProt + CAZy")
    merge_uniprot_and_cazy()

    print("\n=== Modul 2 ferdig ===")

if __name__ == "__main__":
    run_module2_pipeline()
