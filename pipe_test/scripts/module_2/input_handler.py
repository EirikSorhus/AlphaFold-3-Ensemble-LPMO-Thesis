from pathlib import Path
import time
from typing import Dict, Iterable, List, Optional, Tuple
import requests

from cazy_parser import (
    parse_cazy_txt_lines,
    filter_jgi_seeds,
    filter_ncbi_protein_seeds,
    unique_seed_ids,
)
from ncbi_fallback import fetch_from_ncbi_with_metadata
from paths_config import PipelinePaths
from utils_io import read_fasta


def _download_cazy_txt(family: str, contact_email: Optional[str]) -> Path:
    PipelinePaths.cazy_dir.mkdir(parents=True, exist_ok=True)
    txt_path = PipelinePaths.cazy_dir / f"{family}.txt"
    if txt_path.exists():
        return txt_path

    url = f"https://www.cazy.org/IMG/cazy_data/{family}.txt"
    headers = {"User-Agent": f"BioPipeline/1.0 ({contact_email})" if contact_email else "BioPipeline/1.0"}
    r = requests.get(url, headers=headers, timeout=60)
    r.raise_for_status()
    txt_path.write_text(r.text)
    return txt_path


def bulk_uniprot_search(id_org_list: List[Tuple[str, Optional[str]]], contact_email: Optional[str], chunk_size: int = 100) -> Dict[str, str]:
    """Resolve JGI IDs to UniProt accessions in batches."""
    found: Dict[str, str] = {}
    headers = {"User-Agent": f"BioPipeline/1.0 ({contact_email})" if contact_email else "BioPipeline/1.0"}

    for i in range(0, len(id_org_list), chunk_size):
        chunk = id_org_list[i:i + chunk_size]
        query_parts = []
        for sid, org in chunk:
            if org:
                query_parts.append(f'(protein:"{sid}" AND organism:"{org}")')
            else:
                query_parts.append(f'all:"{sid}"')
        query = " OR ".join(query_parts)

        try:
            r = requests.get(
                "https://rest.uniprot.org/uniprotkb/search",
                params={"query": query, "fields": "accession", "format": "json", "size": 500},
                headers=headers,
                timeout=30,
            )
            if r.ok:
                for entry in r.json().get("results", []):
                    acc = entry.get("primaryAccession")
                    if not acc:
                        continue
                    entry_str = str(entry)
                    for sid, _ in chunk:
                        if sid in entry_str:
                            found.setdefault(sid, acc)
            time.sleep(0.5)
        except Exception as exc:  # pragma: no cover - best-effort mapping
            print(f"[!] UniProt batch search failed: {exc}")
    return found


def bulk_ncbi_to_uniprot_map(ncbi_ids: Iterable[str], contact_email: Optional[str]) -> Dict[str, str]:
    """Attempt a lightweight mapping of NCBI protein IDs to UniProt accessions."""
    ncbi_list = list(ncbi_ids)
    if not ncbi_list:
        return {}

    mapped: Dict[str, str] = {}
    headers = {"User-Agent": f"BioPipeline/1.0 ({contact_email})" if contact_email else "BioPipeline/1.0"}
    for i in range(0, len(ncbi_list), 100):
        chunk = ncbi_list[i:i + 100]
        query = " OR ".join(f'database:(type:refseq {cid})' for cid in chunk)
        try:
            r = requests.get(
                "https://rest.uniprot.org/uniprotkb/search",
                params={"query": query, "fields": "accession", "format": "json", "size": 500},
                headers=headers,
                timeout=30,
            )
            if r.ok:
                for entry in r.json().get("results", []):
                    acc = entry.get("primaryAccession")
                    entry_str = str(entry)
                    for cid in chunk:
                        if cid in entry_str:
                            mapped.setdefault(cid, acc)
            time.sleep(0.3)
        except Exception as exc:
            print(f"[!] NCBI->UniProt search failed: {exc}")
    return mapped


def handle_input(
    families: Optional[List[str]] = None,
    id_file: Optional[str] = None,
    fasta_path: Optional[str] = None,
    allow_ncbi: bool = False,
    max_seq_search: int = 10,
    contact_email: Optional[str] = None,
    run_id: Optional[str] = None,
):
    ids = set()
    ncbi_metadata_dict: Dict[str, Dict[str, str]] = {}
    failed: List[Tuple[str, str, str]] = []  # (ID, Organism, Reason)
    source_log = {
        "ncbi_mapped": 0,
        "jgi_mapped": 0,
        "ncbi_fallback": 0,
        "input_file": 0,
        "fasta_matched": 0,
        "fasta_unmatched": 0,
        "failed": 0,
    }

    # === CAZy-family input ===
    if families:
        seeds: List = []
        for fam in families:
            print(f"[*] Processing family: {fam}")
            txt_path = _download_cazy_txt(fam, contact_email)
            seeds.extend(parse_cazy_txt_lines(txt_path.read_text().splitlines(), fam))

        jgi_seeds = filter_jgi_seeds(seeds)
        ncbi_seeds = filter_ncbi_protein_seeds(seeds)

        # Map JGI -> UniProt in batches
        jgi_map = bulk_uniprot_search([(s.seed_id, s.organism) for s in jgi_seeds], contact_email)
        for seed in jgi_seeds:
            mapped = jgi_map.get(seed.seed_id)
            if mapped:
                ids.add(mapped)
                source_log["jgi_mapped"] += 1
            elif allow_ncbi:
                seq, meta = fetch_from_ncbi_with_metadata(seed.seed_id, contact_email=contact_email)
                if seq and meta:
                    meta["Source_DB"] = "NCBI_Fallback"
                    ncbi_id = f"NCBI_{meta['UniProt_ID']}"
                    ncbi_metadata_dict[ncbi_id] = meta
                    ids.add(ncbi_id)
                    source_log["ncbi_fallback"] += 1
                else:
                    failed.append((seed.seed_id, seed.organism or "", "ncbi_fallback_failed"))
                    source_log["failed"] += 1
            else:
                failed.append((seed.seed_id, seed.organism or "", "uniprot_not_found"))
                source_log["failed"] += 1

        # Map NCBI protein IDs -> UniProt
        ncbi_ids = unique_seed_ids(ncbi_seeds)
        ncbi_map = bulk_ncbi_to_uniprot_map(ncbi_ids, contact_email)
        for cid in ncbi_ids:
            mapped = ncbi_map.get(cid)
            if mapped:
                ids.add(mapped)
                source_log["ncbi_mapped"] += 1
            elif allow_ncbi:
                seq, meta = fetch_from_ncbi_with_metadata(cid, contact_email=contact_email)
                if seq and meta:
                    ncbi_id = f"NCBI_{meta['UniProt_ID']}"
                    ncbi_metadata_dict[ncbi_id] = meta
                    ids.add(ncbi_id)
                    source_log["ncbi_fallback"] += 1
                else:
                    failed.append((cid, "", "ncbi_fallback_failed"))
                    source_log["failed"] += 1
            else:
                failed.append((cid, "", "uniprot_not_found"))
                source_log["failed"] += 1

    # === UniProt ID file ===
    elif id_file:
        for line in Path(id_file).read_text().splitlines():
            val = line.strip()
            if val:
                ids.add(val)
                source_log["input_file"] += 1

    # === Fasta input ===
    elif fasta_path:
        seqs = read_fasta(Path(fasta_path))
        for i, (header, seq) in enumerate(seqs.items()):
            if i >= max_seq_search:
                print(f"[!] Max sequence search limit ({max_seq_search}) reached.")
                break
            uid = None
            if "|" in header:
                parts = header.split("|")
                if len(parts) >= 2 and 2 <= len(parts[1]) <= 12:
                    uid = parts[1]
            if uid:
                ids.add(uid)
                source_log["fasta_matched"] += 1
            elif allow_ncbi:
                seq_text, meta = fetch_from_ncbi_with_metadata(header.strip(), contact_email=contact_email)
                if seq_text and meta:
                    ncbi_id = f"NCBI_{meta['UniProt_ID']}"
                    ncbi_metadata_dict[ncbi_id] = meta
                    ids.add(ncbi_id)
                    source_log["ncbi_fallback"] += 1
                else:
                    source_log["fasta_unmatched"] += 1
                    failed.append((header, "", "fasta_unmatched"))
            else:
                source_log["fasta_unmatched"] += 1
                failed.append((header, "", "fasta_unmatched"))

    return sorted(list(ids)), ncbi_metadata_dict, source_log, failed
