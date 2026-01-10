import sys
import os
import logging

# Ensure we can import the local module
sys.path.append(os.getcwd())

from scripts.module_2_new.uniprot_client import UniProtClient

logging.basicConfig(level=logging.INFO)
# Mute connection logs
logging.getLogger("urllib3").setLevel(logging.WARNING)

client = UniProtClient()

# Data from the error log
organism = "Aureococcus anophagefferens clone 1984"
ids = ["61551", "66375"]

def test_query(name, q):
    print(f"\n--- {name} ---")
    print(f"Query: {q}")
    res = client.search_by_query(q)
    print(f"Hit count: {len(res)}")

# 1. Exact reproduction of failure
id_part_1 = " OR ".join([f'"{jid}"' for jid in ids])
q1 = f'organism:"{organism}" AND ({id_part_1})'
test_query("Original Failing Query", q1)

# 2. Try 'organism_id' if we can find it? No we have name. Use 'taxonomy_name'? 
# UniProt field is 'organism_name' or 'organism'.
# Try replacing 'organism' with 'organism_name'
q2 = f'organism_name:"{organism}" AND ({id_part_1})'
test_query("Using organism_name field", q2)

# 3. Maybe the organism name has issues? Try searching just the ID.
q3 = f'"{ids[0]}"'
test_query("Search ID only", q3)

# 4. Try removing quotes from ID
id_part_4 = " OR ".join(ids)
q4 = f'organism:"{organism}" AND ({id_part_4})'
test_query("No quotes on IDs", q4)
