#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple
import time
import requests

UNIPROT_REST = "https://rest.uniprot.org"


@dataclass
class MappingStats:
    submitted: int
    mapped: int
    unmapped: int
    job_id: str


class UniProtIdMapper:
    """
    UniProt ID mapping (programmatic) – kjører:
      1) POST /idmapping/run
      2) GET  /idmapping/status/{jobId} (poll)
      3) GET  /idmapping/results/{jobId}?format=tsv

    NB: 'from' database-koder kan variere; driveren prøver flere.
    """

    def __init__(
        self,
        user_agent: str,
        rate_limit_s: float = 0.25,
        timeout_s: int = 60,
        max_poll_s: int = 600,
    ) -> None:
        self.session = requests.Session()
        self.headers = {"User-Agent": user_agent}
        self.rate_limit_s = rate_limit_s
        self.timeout_s = timeout_s
        self.max_poll_s = max_poll_s

    def _sleep(self) -> None:
        time.sleep(self.rate_limit_s)

    def submit_job(self, from_db: str, to_db: str, ids: List[str]) -> str:
        url = f"{UNIPROT_REST}/idmapping/run"
        data = {"from": from_db, "to": to_db, "ids": ",".join(ids)}
        r = self.session.post(url, data=data, headers=self.headers, timeout=self.timeout_s)
        r.raise_for_status()
        self._sleep()
        js = r.json()
        job_id = js.get("jobId")
        if not job_id:
            raise RuntimeError(f"Unexpected submit response (no jobId): {js}")
        return job_id

    def wait_for_job(self, job_id: str) -> None:
        url = f"{UNIPROT_REST}/idmapping/status/{job_id}"
        t0 = time.time()
        while True:
            r = self.session.get(url, headers=self.headers, timeout=self.timeout_s)
            r.raise_for_status()
            self._sleep()
            js = r.json()
            status = js.get("jobStatus")

            if status in (None, "FINISHED"):
                return
            if status == "FAILED":
                raise RuntimeError(f"ID mapping job failed: {js}")

            if time.time() - t0 > self.max_poll_s:
                raise TimeoutError(f"Timed out waiting for UniProt mapping job {job_id}")

            time.sleep(1.0)

    def fetch_results_tsv(self, job_id: str) -> str:
        url = f"{UNIPROT_REST}/idmapping/results/{job_id}"
        params = {"format": "tsv"}
        r = self.session.get(url, params=params, headers=self.headers, timeout=self.timeout_s)
        r.raise_for_status()
        self._sleep()
        return r.text

    @staticmethod
    def parse_results_tsv(tsv_text: str) -> Dict[str, str]:
        lines = [ln for ln in tsv_text.splitlines() if ln.strip()]
        if not lines:
            return {}

        header = lines[0].split("\t")
        try:
            i_from = header.index("From")
            i_to = header.index("To")
        except ValueError:
            raise RuntimeError(f"Unexpected ID mapping TSV header: {header}")

        mapping: Dict[str, str] = {}
        for ln in lines[1:]:
            cols = ln.split("\t")
            if len(cols) <= max(i_from, i_to):
                continue
            f = cols[i_from].strip()
            t = cols[i_to].strip()
            if f and t and f not in mapping:
                mapping[f] = t
        return mapping

    def map_ids(
        self,
        ids: List[str],
        from_db: str,
        to_db: str = "UniProtKB",
        batch_size: int = 500,
    ) -> Tuple[Dict[str, str], MappingStats]:
        """
        Batch mapping for store ID-lister.
        Returnerer:
          mapping: from_id -> uniprot_id (første mapping hvis flere)
          stats: summerte tellinger
        """
        full_map: Dict[str, str] = {}
        submitted = 0
        last_job = None

        for start in range(0, len(ids), batch_size):
            chunk = ids[start : start + batch_size]
            submitted += len(chunk)

            job_id = self.submit_job(from_db=from_db, to_db=to_db, ids=chunk)
            last_job = job_id
            self.wait_for_job(job_id)
            tsv = self.fetch_results_tsv(job_id)
            m = self.parse_results_tsv(tsv)
            full_map.update(m)

        mapped = len(full_map)
        unmapped = submitted - mapped
        stats = MappingStats(submitted=submitted, mapped=mapped, unmapped=unmapped, job_id=last_job or "NONE")
        return full_map, stats
