import requests
import re
import time
import logging
from xml.etree import ElementTree

logger = logging.getLogger(__name__)

NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

def run_blast_search(sequence):
    """
    Submits sequence to NCBI BLAST (blastp) against 'swissprot'.
    Returns ID only if 100% Identity and covers full query length.
    """
    logger.info("   [BLAST] ID missing. Running BLAST search against 'swissprot'...")
    
    query_len = len(sequence)
    
    # 1. Submit BLAST job
    # We use 'swissprot' because it's curated and maps directly to UniProt IDs
    params = {
        "CMD": "Put",
        "PROGRAM": "blastp",
        "DATABASE": "swissprot", 
        "QUERY": sequence
    }
    
    try:
        response = requests.post(NCBI_BLAST_URL, data=params, timeout=30)
        if not response.ok:
            logger.error(f"[BLAST] Submission failed: {response.status_code}")
            return None
        
        rid_match = re.search(r"RID = (.*)", response.text)
        if not rid_match:
            logger.error("   [BLAST] Could not retrieve RID.")
            return None
        rid = rid_match.group(1).strip()
        
        # 2. Wait for results
        logger.info(f"   [BLAST] Job submitted. RID: {rid}. Waiting...")
        max_polls = 60  # Max 60 polls = 10 minutes total wait time
        poll_count = 0
        while poll_count < max_polls:
            poll_count += 1
            time.sleep(10)  # Polling interval
            check_params = {"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid}
            try:
                check_r = requests.get(NCBI_BLAST_URL, params=check_params, timeout=30)
                if "Status=WAITING" in check_r.text:
                    continue
                elif "Status=FAILED" in check_r.text or "Status=UNKNOWN" in check_r.text:
                    logger.error("   [BLAST] Job Failed or Unknown status.")
                    return None
                elif "Status=READY" in check_r.text:
                    break
                else:
                    logger.warning(f"   [BLAST] Unexpected status in response. Continuing...")
                    continue
            except requests.exceptions.Timeout:
                # If polling times out, just retry
                continue
        
        if poll_count >= max_polls:
            logger.error(f"   [BLAST] Timeout: Job did not complete after {max_polls} polls (10 min)")
            return None
        
        # 3. Retrieve and parse results
        result_params = {"CMD": "Get", "FORMAT_TYPE": "XML", "RID": rid}
        result_r = requests.get(NCBI_BLAST_URL, params=result_params, timeout=30)
        
        try:
            root = ElementTree.fromstring(result_r.content)
        except ElementTree.ParseError:
            logger.error("[BLAST] Failed to parse XML response.")
            return None
        
        # Iterate through hits to find an EXACT match
        for hit in root.findall(".//Hit"):
            # Check the first HSP (High-scoring Segment Pair) of the hit
            hsp = hit.find("Hit_hsps/Hsp")
            if hsp is None: 
                continue

            identity = int(hsp.find("Hsp_identity").text)
            align_len = int(hsp.find("Hsp_align-len").text)
            
            # Exact Match Logic:
            # Identity must equal alignment length (100% identity)
            # Alignment length must equal Query length (100% coverage)
            # Alignment length should also match identity (implied)
            # Note: Sometimes align_len includes gaps, but for 100% identity gaps should be 0.
            
            if identity == align_len and align_len == query_len:
                accession = hit.find("Hit_accession").text
                logger.info(f"   [BLAST] Found EXACT UniProt match: {accession}")
                return accession
        
        logger.info("   [BLAST] No exact UniProt match (100% Identity + 100% Coverage) found.")
        return None
            
    except Exception as e:
        logger.error(f"   [ERROR] BLAST failed: {e}")
        return None
