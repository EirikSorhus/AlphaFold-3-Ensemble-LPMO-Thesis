# scripts/module2/cazy_fetch.py

# Fetch genebank ids from CAZy
import requests
from bs4 import BeautifulSoup
from .config import CAZY_FAMILIES

BASE_URL = "https://www.cazy.org"

def fetch_cazy_genbank_ids_for_family(family: str):
    """
    Henter GenBank-accessions for én CAZy-familie (f.eks. 'AA9').
    NB: CAZy kan endre HTML-struktur – dette er et praktisk utgangspunkt.
    """
    url = f"{BASE_URL}/{family}_family.html"
    print(f"[CAZy] Henter: {url}")
    r = requests.get(url)
    r.raise_for_status()

    soup = BeautifulSoup(r.text, "html.parser")
    genbank_ids = []

    # Mange av GenBank-linkene peker mot NCBI Protein
    for link in soup.find_all("a", href=True):
        href = link["href"]
        if "ncbi.nlm.nih.gov/protein" in href:
            acc = link.text.strip()
            if acc:
                genbank_ids.append(acc)

    print(f"[CAZy] {family}: fant {len(genbank_ids)} GenBank-IDer")
    return genbank_ids

def fetch_all_cazy_genbank_ids():
    all_ids = []
    for fam in CAZY_FAMILIES:
        ids = fetch_cazy_genbank_ids_for_family(fam)
        all_ids.extend(ids)
    # unike IDer
    all_ids = sorted(set(all_ids))
    print(f"[CAZy] Totalt unike GenBank-IDer: {len(all_ids)}")
    return all_ids

if __name__ == "__main__":
    ids = fetch_all_cazy_genbank_ids()
    # bare for test
    print(ids[:20])
