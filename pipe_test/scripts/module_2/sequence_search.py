import re
import time
from typing import Optional
import requests
from xml.etree import ElementTree

NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def run_blast_search(
    sequence: str,
    poll_seconds: float = 5.0,
    database: str = "swissprot",
    max_seconds: Optional[float] = None,
    max_polls: Optional[int] = None,
) -> Optional[str]:
    """Submit sequence to NCBI BLAST (blastp) against SwissProt and return a UniProt accession
    only if an exact match is found (100% identity and 100% coverage).

    By default, no hard wall-clock cap is enforced; set max_seconds or max_polls to limit runtime.
    """
    params = {
        "CMD": "Put",
        "PROGRAM": "blastp",
        "DATABASE": database,
        "QUERY": sequence,
    }

    try:
        response = requests.post(NCBI_BLAST_URL, data=params, timeout=30)
        if not response.ok:
            return None
        rid_match = re.search(r"RID = (.*)", response.text)
        if not rid_match:
            return None
        rid = rid_match.group(1).strip()

        polls = 0
        while True:
            time.sleep(max(poll_seconds, 1.0))
            polls += 1
            if max_polls is not None and polls > max_polls:
                return None
            # If a wall-clock limit is provided, enforce it
            # (We check elapsed by storing the start time lazily to avoid overhead when not needed.)
            #
            # Note: start_time is captured only if max_seconds is set to keep the common path fast.
            #
            # pragma: no cover
            if max_seconds is not None:
                if polls == 1:
                    start_time = time.time()
                else:
                    if time.time() - start_time > max_seconds:
                        return None
            check_params = {"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid}
            check_r = requests.get(NCBI_BLAST_URL, params=check_params, timeout=20)
            text = check_r.text
            if "Status=WAITING" in text:
                continue
            if "Status=FAILED" in text or "Status=UNKNOWN" in text:
                return None
            if "Status=READY" in text:
                break

        result_params = {"CMD": "Get", "FORMAT_TYPE": "XML", "RID": rid}
        result_r = requests.get(NCBI_BLAST_URL, params=result_params, timeout=60)
        if not result_r.ok:
            return None

        root = ElementTree.fromstring(result_r.content)
        query_len = None
        qlen_node = root.find(".//Iteration_query-len")
        if qlen_node is not None:
            try:
                query_len = int(qlen_node.text)
            except Exception:
                query_len = None
        if query_len is None:
            query_len = len(sequence)

        for hit in root.findall(".//Hit"):
            hsp = hit.find("Hit_hsps/Hsp")
            if hsp is None:
                continue
            try:
                identity = int(hsp.find("Hsp_identity").text)
                align_len = int(hsp.find("Hsp_align-len").text)
            except Exception:
                continue
            if identity == align_len and align_len == query_len:
                accession = hit.findtext("Hit_accession")
                return accession
        return None
    except Exception:
        return None
