import requests
import time
import logging
import json

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class InterProClient:
    """
    Fetches domain-type entries from InterPro7 member databases with source prioritization.
    Deduplicates overlapping domains and adjusts LPMO domain starts to first Histidine.
    """
    
    BASE_URL = "https://www.ebi.ac.uk/interpro/api"
    # Source prioritization: Pfam > CDD > SMART > PROSITE profiles
    MEMBER_DATABASES = ["pfam", "cdd", "smart", "profile"]
    SOURCE_PRIORITY = {"pfam": 1, "cdd": 2, "smart": 3, "profile": 4, "ncbifam": 5}

    def __init__(self, max_retries=3, delay=1.0, cache_path="interpro_cache.json"):
        self.max_retries = max_retries
        self.delay = delay
        self.cache_path = cache_path
        self.cache = self._load_cache()
    
    def _load_cache(self):
        """Load InterPro cache from JSON file."""
        import os
        if os.path.exists(self.cache_path):
            try:
                with open(self.cache_path, 'r') as f:
                    return json.load(f)
            except Exception as e:
                logger.warning(f"Could not load InterPro cache: {e}")
                return {}
        return {}
    
    def _save_cache(self):
        """Save InterPro cache to JSON file."""
        try:
            with open(self.cache_path, 'w') as f:
                json.dump(self.cache, f)
        except Exception as e:
            logger.warning(f"Could not save InterPro cache: {e}")

    def fetch_domains(self, uniprot_acc, sequence=None, signal_end=None):
        """
        Fetches deduplicated domain-type entries with source prioritization.
        Adjusts LPMO domain starts to first H after signal peptide.
        Uses cache to avoid redundant API calls.
        
        Args:
            uniprot_acc: UniProt accession
            sequence: Protein sequence (for H-adjustment)
            signal_end: Signal peptide cleavage site (0-indexed position)
        """
        # Check cache first
        cache_key = f"{uniprot_acc}_domains"
        if cache_key in self.cache:
            logger.debug(f"Using cached domains for {uniprot_acc}")
            return self.cache[cache_key]
        
        # Fetch from all sources
        all_domains = []
        for db in self.MEMBER_DATABASES:
            domains = self._fetch_from_source(uniprot_acc, db)
            all_domains.extend(domains)
            if domains:
                logger.debug(f"Fetched {len(domains)} domains from {db} for {uniprot_acc}")
        
        if not all_domains:
            return []
        
        logger.debug(f"Total domains fetched for {uniprot_acc}: {len(all_domains)}")
        for d in all_domains:
            logger.debug(f"  - {d['source']}:{d['model']} ({d['start']}-{d['end']}) {d['name']}")
        
        # Deduplicate overlapping domains by priority
        deduped = self._deduplicate_domains(all_domains)
        
        logger.info(f"After dedup for {uniprot_acc}: {len(deduped)} domains kept")
        for d in deduped:
            logger.info(f"  - {d['source']}:{d['model']} ({d['start']}-{d['end']}) {d['name']}")
        
        # Adjust LPMO domain starts to first H
        if sequence and signal_end is not None:
            deduped = self._adjust_lpmo_starts(deduped, sequence, signal_end)
        
        # Cache the result
        self.cache[cache_key] = deduped
        self._save_cache()
        
        return deduped

    def fetch_signal_peptide(self, uniprot_acc):
        """
        Fetches SignalP predictions from InterPro extra_features as fallback.
        Returns dict with 'start', 'end', 'source' or None.
        Uses cache to avoid redundant API calls.
        """
        cache_key = f"{uniprot_acc}_signal"
        if cache_key in self.cache:
            logger.debug(f"Using cached signal peptide for {uniprot_acc}")
            return self.cache[cache_key]
        
        url = f"{self.BASE_URL}/protein/uniprot/{uniprot_acc}"
        params = {"extra_features": "true"}
        headers = {"Accept": "application/json"}
        
        attempt = 0
        while attempt < self.max_retries:
            try:
                response = requests.get(url, headers=headers, params=params, timeout=10)
                
                if response.status_code in [404, 204]:
                    self.cache[cache_key] = None
                    self._save_cache()
                    return None
                
                response.raise_for_status()
                data = response.json()
                
                # Look for SignalP or Phobius
                extra_features = data.get("extra_features", {})
                for source_key in ["signalp_euk", "signalp_gram_positive", "signalp_gram_negative", "phobius"]:
                    if source_key in extra_features:
                        locations = extra_features[source_key].get("locations", [])
                        if locations:
                            loc = locations[0]
                            fragments = loc.get("fragments", [])
                            if fragments:
                                frag = fragments[0]
                                result = {
                                    "start": frag.get("start"),
                                    "end": frag.get("end"),
                                    "source": source_key
                                }
                                self.cache[cache_key] = result
                                self._save_cache()
                                return result
                
                self.cache[cache_key] = None
                self._save_cache()
                return None

            except (json.JSONDecodeError, requests.exceptions.RequestException) as e:
                attempt += 1
                if attempt >= self.max_retries:
                    logger.error(f"Failed to fetch SignalP for {uniprot_acc}: {e}")
                    self.cache[cache_key] = None
                    self._save_cache()
                    return None
                time.sleep(self.delay * attempt)
        
        return None

    def _fetch_from_source(self, uniprot_acc, source_db):
        """Fetches domain-type entries from a specific member database."""
        url = f"{self.BASE_URL}/entry/{source_db}/protein/uniprot/{uniprot_acc}"
        params = {"type": "domain", "page_size": 200}
        headers = {"Accept": "application/json"}
        
        attempt = 0
        while attempt < self.max_retries:
            try:
                response = requests.get(url, headers=headers, params=params, timeout=10)
                
                if response.status_code in [404, 204]:
                    return []
                
                response.raise_for_status()
                data = response.json()
                
                return self._parse_response(data, uniprot_acc, source_db)

            except (json.JSONDecodeError, requests.exceptions.RequestException) as e:
                attempt += 1
                if attempt >= self.max_retries:
                    logger.error(f"Failed to fetch {source_db} for {uniprot_acc}: {e}")
                    return []
                time.sleep(self.delay * attempt)
        
        return []

    def _parse_response(self, data, target_acc, source_db):
        """Parses InterPro7 JSON for domain entries."""
        parsed_entries = []
        results = data.get("results", [])
        
        for res in results:
            metadata = res.get("metadata", {})
            entry_acc = metadata.get("accession")
            entry_name = metadata.get("name", "Unknown")
            entry_type = metadata.get("type", "")
            source = metadata.get("source_database", source_db)
            integrated_ipr = metadata.get("integrated")
            
            # Filter: only domain type
            if entry_type.lower() != "domain":
                continue
            
            proteins = res.get("proteins", [])
            for prot in proteins:
                if prot.get("accession", "").lower() == target_acc.lower():
                    locations = prot.get("entry_protein_locations", [])
                    for loc in locations:
                        model = loc.get("model", entry_acc)
                        representative = loc.get("representative", False)
                        score = loc.get("score")
                        
                        fragments = loc.get("fragments", [])
                        for frag in fragments:
                            parsed_entries.append({
                                "source": source.lower(),
                                "model": model,
                                "entry_id": entry_acc,
                                "integrated_ipr": integrated_ipr,
                                "name": entry_name,
                                "type": entry_type,
                                "start": frag.get("start"),
                                "end": frag.get("end"),
                                "representative": representative,
                                "score": score,
                            })
        
        return parsed_entries

    def _deduplicate_domains(self, domains):
        """
        Deduplicates overlapping domains by source priority.
        For overlapping pairs, always keeps the highest-priority (lowest priority number).
        """
        if not domains:
            return []
        
        # Sort ONLY by priority (Pfam=1, CDD=2, SMART=3, etc.)
        domains_sorted = sorted(
            domains,
            key=lambda d: self.SOURCE_PRIORITY.get(d["source"], 99)
        )
        
        logger.debug(f"Dedup: Sorted by priority: {[(d['source'], d['model'], d['start'], d['end']) for d in domains_sorted]}")
        
        kept = []
        for domain in domains_sorted:
            # Check if this domain overlaps with any already-kept domain
            overlaps_with_kept = False
            for kept_domain in kept:
                if self._domains_overlap(domain, kept_domain):
                    # Since we sorted by priority, current domain has equal or lower priority
                    # Skip it (keep the higher-priority one already in kept)
                    logger.debug(f"Dedup: Dropping {domain['source']}:{domain['model']} ({domain['start']}-{domain['end']}) - overlaps with kept {kept_domain['source']}:{kept_domain['model']}")
                    overlaps_with_kept = True
                    break
            
            if not overlaps_with_kept:
                logger.debug(f"Dedup: Keeping {domain['source']}:{domain['model']} ({domain['start']}-{domain['end']})")
                kept.append(domain)
        
        return kept

    def _domains_overlap(self, d1, d2, threshold=0.5):
        """
        Check if two domains overlap significantly.
        Returns True if overlap > threshold * length of shorter domain.
        """
        start1, end1 = d1["start"], d1["end"]
        start2, end2 = d2["start"], d2["end"]
        
        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)
        
        if overlap_start >= overlap_end:
            return False
        
        overlap_len = overlap_end - overlap_start
        min_len = min(end1 - start1, end2 - start2)
        
        return (overlap_len / min_len) > threshold

    def _adjust_lpmo_starts(self, domains, sequence, signal_end):
        """
        Adjust LPMO domain starts to first Histidine after signal peptide.
        Only modifies domains identified as LPMO-related.
        """
        lpmo_keywords = ["lpmo", "aa9", "aa10", "aa11", "aa13", "aa14", "aa15", 
                         "aa16", "aa17", "polysaccharide monooxygenase", 
                         "cbm33", "gh61"]
        
        adjusted = []
        for domain in domains:
            is_lpmo = any(kw in domain["name"].lower() for kw in lpmo_keywords)
            
            if is_lpmo:
                domain_start = domain["start"] - 1  # Convert to 0-indexed
                
                # Check if starts on H
                if domain_start < len(sequence) and sequence[domain_start].upper() == 'H':
                    adjusted.append(domain)
                    continue
                
                # Search backward for first H (not before signal_end)
                first_h_pos = None
                for pos in range(domain_start - 1, signal_end, -1):
                    if pos < len(sequence) and sequence[pos].upper() == 'H':
                        first_h_pos = pos
                        break
                
                if first_h_pos is not None:
                    domain_copy = domain.copy()
                    domain_copy["start"] = first_h_pos + 1  # Convert back to 1-indexed
                    domain_copy["adjusted"] = True
                    adjusted.append(domain_copy)
                    logger.debug(f"Adjusted LPMO domain '{domain['name']}' start from {domain['start']} to {first_h_pos + 1}")
                else:
                    # No H found, keep original
                    adjusted.append(domain)
            else:
                adjusted.append(domain)
        
        return adjusted


if __name__ == "__main__":
    # Test with example
    client = InterProClient()
    print("Testing with O09185...")
    domains = client.fetch_domains("O09185")
    print(f"Found {len(domains)} deduplicated domains:")
    for d in domains:
        print(f"  {d['source']:8} {d['model']:12} {d['name']:40} ({d['start']}-{d['end']})")
