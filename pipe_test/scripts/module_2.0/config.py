# scripts/module2/config.py

# CAZy-familier for LPMOs
CAZY_FAMILIES = ["AA9", "AA10", "AA11", "AA13"]

# UniProt-query for LPMOs (kan justeres underveis)
UNIPROT_QUERY = '(family:"AA9" OR family:"AA10" OR family:"AA11" OR family:"AA13")'

# NCBI krever en e-postadresse
NCBI_EMAIL = "din.epost@institusjon.no"

# Paths (relativt til project_root)
DATA_RAW_DIR = "data_raw"

UNIPROT_FASTA = f"{DATA_RAW_DIR}/uniprot_LPMO_raw.fasta"
CAZY_FASTA = f"{DATA_RAW_DIR}/cazy_LPMO_raw.fasta"
MERGED_FASTA = f"{DATA_RAW_DIR}/LPMO_all_raw.fasta"