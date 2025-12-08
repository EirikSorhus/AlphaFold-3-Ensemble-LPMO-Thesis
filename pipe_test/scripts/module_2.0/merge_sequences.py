# scripts/module2/merge_sequences.py

# Merge FASTA files from UniProt and CAZy, removing duplicates
import os
from Bio import SeqIO
from .config import UNIPROT_FASTA, CAZY_FASTA, MERGED_FASTA, DATA_RAW_DIR

def merge_fastas(input_fastas, out_fasta):
    os.makedirs(DATA_RAW_DIR, exist_ok=True)

    records = []
    for f in input_fastas:
        if not os.path.exists(f):
            print(f"[MERGE] Advarsel: {f} finnes ikke, hopper over.")
            continue
        print(f"[MERGE] Leser: {f}")
        for rec in SeqIO.parse(f, "fasta"):
            records.append(rec)

    print(f"[MERGE] Totalt antall sekvenser før enkel deduplisering: {len(records)}")

    # Enkel deduplisering basert på (id, sequence)
    unique = {}
    for rec in records:
        key = (str(rec.seq), rec.id)
        if key not in unique:
            unique[key] = rec

    unique_records = list(unique.values())
    print(f"[MERGE] Antall sekvenser etter enkel deduplisering: {len(unique_records)}")

    SeqIO.write(unique_records, out_fasta, "fasta")
    print(f"[MERGE] Skrev sammenslått FASTA til: {out_fasta}")

def merge_uniprot_and_cazy():
    merge_fastas([UNIPROT_FASTA, CAZY_FASTA], MERGED_FASTA)

if __name__ == "__main__":
    merge_uniprot_and_cazy()
