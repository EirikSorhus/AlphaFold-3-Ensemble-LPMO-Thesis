# FILE: cazy_parser.py

from dataclasses import dataclass
import re
from typing import Iterable, List, Optional, Set


CAZY_SOURCE_TOKENS = {"ncbi", "jgi", "uniprot", "genbank", "embl", "ddb", "refseq"}
CAZY_KINGDOM_TOKENS = {"Bacteria", "Archaea", "Eukaryota", "Viruses"}
_ncbi_protein_like = re.compile(r"^(?:[A-Z]{3}\d+\.[0-9]+|[A-Z]{2}_[0-9]+\.[0-9]+|[A-Z]{2,}[0-9]+\.[0-9]+)$")


@dataclass(frozen=True)
class CazySeed:
    seed_id: str
    seed_source: str
    family: str
    organism: Optional[str]
    kingdom: Optional[str]
    raw_line: str


def parse_cazy_txt_lines(lines: Iterable[str], family: str) -> List[CazySeed]:
    """Parse CAZy .txt content and normalise into CazySeed objects.

    The CAZy dumps are mostly tab-separated but can occasionally be space
    separated. We try both formats and keep whatever metadata we can recover.
    """
    seeds: List[CazySeed] = []
    fam = family.strip()
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue

        if "\t" in line:
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            fam_col = cols[0].strip()
            kingdom_col = cols[1].strip() if len(cols) > 1 else None
            organism = cols[2].strip() if len(cols) > 2 else None
            seed_id = cols[-2].strip()
            src = cols[-1].strip().lower()
            if fam_col != fam or not seed_id or not src:
                continue
            if kingdom_col not in CAZY_KINGDOM_TOKENS:
                kingdom_col = None
            seeds.append(CazySeed(seed_id, src, fam, organism, kingdom_col, line))
            continue

        # Fallback for space-separated legacy lines
        toks = line.split()
        if len(toks) < 4:
            continue
        for i in range(1, len(toks) - 2):
            src = toks[i].lower()
            if src not in CAZY_SOURCE_TOKENS:
                continue
            if toks[i + 1] != fam:
                continue
            kingdom = toks[i + 2] if toks[i + 2] in CAZY_KINGDOM_TOKENS else None
            seed_id = toks[i - 1]
            seeds.append(CazySeed(seed_id, src, fam, None, kingdom, line))
            break
    return seeds


def filter_ncbi_protein_seeds(seeds: Iterable[CazySeed]) -> List[CazySeed]:
    return [s for s in seeds if s.seed_source == "ncbi" and _ncbi_protein_like.match(s.seed_id)]


def filter_jgi_seeds(seeds: Iterable[CazySeed]) -> List[CazySeed]:
    return [s for s in seeds if s.seed_source == "jgi"]


def unique_seed_ids(seeds: Iterable[CazySeed]) -> List[str]:
    seen: Set[str] = set()
    ids: List[str] = []
    for s in seeds:
        if s.seed_id not in seen:
            seen.add(s.seed_id)
            ids.append(s.seed_id)
    return ids


# Note: Single-ID JGI lookups removed; pipeline uses batched UniProt searches.
