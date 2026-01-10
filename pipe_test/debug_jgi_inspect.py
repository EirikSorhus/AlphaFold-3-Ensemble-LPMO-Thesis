import sys
import os
import logging
import json

sys.path.append(os.getcwd())
from scripts.module_2_new.uniprot_client import UniProtClient

logging.basicConfig(level=logging.INFO)
logging.getLogger("urllib3").setLevel(logging.WARNING)

client = UniProtClient()

# Test simplified organism name
full_org = "Aureococcus anophagefferens clone 1984"
simple_org = " ".join(full_org.split()[:2])
ids = ["61551", "66375"] # 61551 is verified to be in Aureococcus anophagefferens

print(f"Testing simplified organism: '{simple_org}'")

# Test 1: organism_name:"Genus species" AND (ID OR ID)
id_part = " OR ".join([f'"{jid}"' for jid in ids])
q1 = f'organism_name:"{simple_org}" AND ({id_part})'
print(f"Query 1: {q1}")
res1 = client.search_by_query(q1)
print(f"Result 1 count: {len(res1)}")

# Test 2: organism:"Genus species" AND (ID OR ID) (using scientific name field)
q2 = f'organism:"{simple_org}" AND ({id_part})'
print(f"Query 2: {q2}")
res2 = client.search_by_query(q2)
print(f"Result 2 count: {len(res2)}")

