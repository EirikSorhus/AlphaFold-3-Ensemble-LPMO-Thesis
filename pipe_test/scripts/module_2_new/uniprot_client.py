import requests
import time
import logging
import json
import os
from tqdm import tqdm

logger = logging.getLogger(__name__)

class UniProtClient:
    def __init__(self, cache_path="discovery_cache.json", batch_size=50):
        self.base_url = "https://rest.uniprot.org/uniprotkb/search"
        self.cache_path = cache_path
        self.batch_size = batch_size
        self.cache = self._load_cache()

    def _load_cache(self):
        if os.path.exists(self.cache_path):
            with open(self.cache_path, 'r') as f:
                try:
                    return json.load(f)
                except json.JSONDecodeError:
                    return {}
        return {}

    def _save_cache(self):
        with open(self.cache_path, 'w') as f:
            json.dump(self.cache, f)

    def fetch_batch(self, ids):
        """
        Hovedmetode for å hente metadata. Bruker bisection hvis batchen feiler.
        """
        # Filtrer ut ID-er som allerede er i cache
        to_fetch = [i for i in ids if i not in self.cache]
        results = [self.cache[i] for i in ids if i in self.cache]

        if not to_fetch:
            return results

        fetched_data = self._request_with_bisection(to_fetch)
        
        # Oppdater cache med nye funn
        for entry in fetched_data:
            acc = entry.get('primaryAccession')
            if acc:
                self.cache[acc] = entry
        
        self._save_cache()
        return results + fetched_data

    def _request_with_bisection(self, id_list):
        """
        Rekursiv bisection-logikk: Hvis en batch feiler, splitt den i to 
        helt til vi isolerer den korrupte ID-en.
        """
        if not id_list:
            return []

        # Bygg query: accession:ID1 OR accession:ID2...
        query = " OR ".join([f"accession:{i}" for i in id_list])
        params = {
            "query": query,
            "fields": "accession,sequence,organism_name,protein_name,xref_interpro,ft_signal,ft_domain,ft_region,ft_binding,ft_site,ft_act_site,ft_zn_fing,ft_motif",
            "format": "json",
            "size": len(id_list)
        }

        try:
            response = requests.get(self.base_url, params=params, timeout=30)
            
            if response.status_code == 200:
                results = response.json().get('results', [])
                return results
            
            elif response.status_code in [400, 500, 503]:
                # Hvis vi har bare én ID og den feiler, er det en "bad ID"
                if len(id_list) == 1:
                    logger.error(f"Permanent feil for ID {id_list[0]} (HTTP {response.status_code})")
                    return []
                
                # Bisection: Del listen i to og prøv hver del for seg
                logger.warning(f"Batch feilet (HTTP {response.status_code}). Splitter {len(id_list)} ID-er...")
                mid = len(id_list) // 2
                return self._request_with_bisection(id_list[:mid]) + \
                       self._request_with_bisection(id_list[mid:])
            
            else:
                response.raise_for_status()

        except Exception as e:
            logger.error(f"Nettverksfeil: {e}. Prøver exponential backoff...")
            time.sleep(5) # Enkel backoff
            return self._request_with_bisection(id_list)

    def search_by_sequence(self, sequence):
        """
        Søker i UniProt etter eksakt sekvens-match (for ukjente FASTA headere).
        """
        # Merk: UniProt API for sekvens-søk kan kreve POST eller spesifikt endpoint
        # Her bruker vi en forenklet versjon av deres krysreferanse-søk
        params = {
            "query": f"sequence:{sequence}",
            "format": "json"
        }
        try:
            res = requests.get(self.base_url, params=params, timeout=30)
            if res.status_code == 200:
                results = res.json().get('results', [])
                return results[0] if results else None
        except Exception:
            return None
