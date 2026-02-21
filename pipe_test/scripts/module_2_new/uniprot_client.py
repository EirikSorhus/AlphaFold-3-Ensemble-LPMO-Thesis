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
        self.accession_url = "https://rest.uniprot.org/uniprotkb/accessions"
        self.cache_path = cache_path
        self.batch_size = batch_size
        self.cache = self._load_cache()
        self.obsolete_mapping_cache = {}

    def _load_cache(self):
        if os.path.exists(self.cache_path):
            with open(self.cache_path, 'r') as f:
                try:
                    return json.load(f)
                except json.JSONDecodeError:
                    return {}
        return {}

    def _save_cache(self):
        import tempfile
        try:
            # Write to a temp file in the same directory as the cache file,
            # then atomically rename to avoid corruption on crash.
            cache_dir = os.path.dirname(os.path.abspath(self.cache_path))
            with tempfile.NamedTemporaryFile(mode='w', dir=cache_dir, delete=False, suffix='.tmp') as tmp:
                json.dump(self.cache, tmp)
                tmp_path = tmp.name
            os.replace(tmp_path, self.cache_path)
        except Exception as e:
            logger.warning(f"Failed to save UniProt cache: {e}")

    def _lookup_primary_via_search(self, acc):
        """Fallback: try search endpoint to retrieve a matching entry for an accession."""
        try:
            params = {
                "query": f"accession:{acc}",
                "fields": "accession",
                "format": "json",
                "size": 1,
            }
            response = requests.get(self.base_url, params=params, timeout=15)
            if response.status_code == 200:
                results = response.json().get('results', [])
                if results:
                    return results[0].get('primaryAccession')
            else:
                logger.debug(f"Search fallback failed for {acc} (HTTP {response.status_code})")
        except Exception as e:
            logger.debug(f"Search fallback error for {acc}: {e}")
        return None

    def resolve_obsolete_ids(self, ids):
        """
        Resolves obsolete/secondary UniProt accessions to their current primary accessions using the
        UniProt accessions endpoint. Returns a dict mapping: {original_id: primary_id or None}.
        """
        mapping = {}

        # Use cached mappings first
        pending = []
        for acc in ids:
            if acc in self.obsolete_mapping_cache:
                mapping[acc] = self.obsolete_mapping_cache[acc]
            else:
                pending.append(acc)

        if pending:
            try:
                params = {"accessions": ",".join(pending)}
                response = requests.get(self.accession_url, params=params, timeout=20)
                if response.status_code == 200:
                    data = response.json()
                    for item in data.get("results", []):
                        old = item.get("from")
                        new = item.get("to")
                        if old:
                            mapping[old] = new
                            self.obsolete_mapping_cache[old] = new
                    # Anything not returned is unmapped
                    for acc in pending:
                        if acc not in mapping:
                            mapping[acc] = None
                            self.obsolete_mapping_cache[acc] = None
                else:
                    logger.warning(f"Accessions endpoint failed (HTTP {response.status_code}); falling back to search for {len(pending)} IDs")
            except Exception as e:
                logger.warning(f"Accessions endpoint error: {e}; falling back to search for {len(pending)} IDs")

        # Fallback search for any still unresolved
        for acc in ids:
            if mapping.get(acc) is None:
                primary = self._lookup_primary_via_search(acc)
                mapping[acc] = primary
                self.obsolete_mapping_cache[acc] = primary
                if primary and primary != acc:
                    logger.info(f"Obsolete ID {acc} -> {primary}")

        return mapping

    def fetch_batch(self, ids):
        """
        Hovedmetode for å hente metadata. Bruker bisection hvis batchen feiler.
        """
        # First resolve any obsolete IDs to their primary accessions
        obsolete_mapping = self.resolve_obsolete_ids(ids)
        
        # Map original IDs to resolved (or keep original if not obsolete)
        resolved_ids = []
        original_to_resolved = {}
        for orig_id in ids:
            resolved = obsolete_mapping.get(orig_id, orig_id)
            if resolved:
                resolved_ids.append(resolved)
                original_to_resolved[orig_id] = resolved
            else:
                # Could not resolve, skip
                logger.warning(f"Could not resolve obsolete ID: {orig_id}")
        
        # Filtrer ut ID-er som allerede er i cache
        to_fetch = [i for i in resolved_ids if i not in self.cache]
        results = [self.cache[i] for i in resolved_ids if i in self.cache]

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

        # Fix: Only use 'accession:' for valid UniProt Accessions (alphanumeric, no dots).
        # GenBank IDs (e.g. AGE49160.1) cause HTTP 400 if used with accession field.
        query_parts = []
        for i in id_list:
             # Basic heuristic: UniProt IDs are alphanumeric and 6-10 chars. GenBank usually has dots or differs.
             if '.' in i or not i.isalnum():
                 query_parts.append(f'"{i}"')
             else:
                 query_parts.append(f"accession:{i}")

        query = " OR ".join(query_parts)
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

    def search_by_query(self, query):
        """
        Executes a direct search query (e.g. for JGI grouped searches).
        """
        params = {
            "query": query,
            "fields": "accession,sequence,organism_name,protein_name,xref_interpro,ft_signal,ft_domain,ft_region,ft_binding,ft_site,ft_act_site,ft_zn_fing,ft_motif",
            "format": "json"
        }
        
        try:
            response = requests.get(self.base_url, params=params, timeout=30)
            if response.status_code == 200:
                results = response.json().get('results', [])
                
                # Cache results found via query
                for entry in results:
                    acc = entry.get('primaryAccession')
                    if acc:
                        self.cache[acc] = entry
                self._save_cache()
                
                return results
            else:
                logger.warning(f"Query failed: {query} (HTTP {response.status_code})")
                return []
        except Exception as e:
            logger.error(f"Network error in search_by_query: {e}")
            return []

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
