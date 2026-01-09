#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Set
import re

# Kilder vi kan støte på i CAZy txt. (Vi bruker først og fremst ncbi/jgi her.)
CAZY_SOURCE_TOKENS = {"ncbi", "jgi", "uniprot", "genbank", "embl", "ddb", "refseq"}

# CAZy i disse filene bruker typisk superkingdom (Bacteria/Archaea/Eukaryota/Viruses).
CAZY_KINGDOM_TOKENS = {"Bacteria", "Archaea", "Eukaryota", "Viruses"}

# NCBI/INSDC protein accessions med versjon (mange varianter finnes; vi matcher konservativt).
# Eksempler: WDK13118.1, CAQ16217.1, WQF89147.1, WP_012345678.1
_ncbi_protein_like = re.compile(
    r"^(?:[A-Z]{3}\d+\.\d+|[A-Z]{2}_\d+\.\d+|[A-Z]{2,}\d+\.\d+)$"
)


@dataclass(frozen=True)
class CazySeed:
    seed_id: str              # f.eks. WDK13118.1 eller 428069
    seed_source: str          # "ncbi" / "jgi" / ...
    family: str               # f.eks. AA9
    kingdom: Optional[str]    # f.eks. Eukaryota (kan være None)
    raw_line: str             # for debugging/trace


def parse_cazy_txt_lines(lines: Iterable[str], family: str) -> List[CazySeed]:
    """
    CAZy TXT finnes i minst to formater:

      A) Tab-separert (som du viste):
         family  kingdom  organism_name  taxid?  seed_id  source_db
         AA9     Eukaryota  ...          428069  WDK13118.1  ncbi
         (seed_id og source_db ligger ofte i de to siste kolonnene)

      B) Whitespace-format (eldre/komprimert):
         <seed_id> <source> <family> <kingdom> ...

    Denne parseren prøver A først per linje, og faller tilbake til B.
    """
    seeds: List[CazySeed] = []
    fam = family.strip()

    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue

        # --- Format A: tab-separert ---
        if "\t" in line:
            cols = line.split("\t")
            # Vi gjør det defensivt: minst 5 kolonner, og antar seed_id + source er de to siste.
            if len(cols) >= 5:
                fam_col = cols[0].strip()
                kingdom_col = cols[1].strip() if len(cols) > 1 else None
                seed_id = cols[-2].strip()
                src = cols[-1].strip().lower()

                if fam_col == fam and seed_id and src:
                    # Normaliser kingdom hvis det er “rart”
                    if kingdom_col not in CAZY_KINGDOM_TOKENS:
                        kingdom_col = None

                    seeds.append(
                        CazySeed(
                            seed_id=seed_id,
                            seed_source=src,
                            family=fam,
                            kingdom=kingdom_col,
                            raw_line=line,
                        )
                    )
            continue

        # --- Format B: whitespace fallback ---
        toks = line.split()
        if len(toks) < 4:
            continue

        for i in range(1, len(toks) - 2):
            src = toks[i].lower()
            if src not in CAZY_SOURCE_TOKENS:
                continue
            if toks[i + 1] != fam:
                continue

            kingdom = toks[i + 2]
            if kingdom not in CAZY_KINGDOM_TOKENS:
                kingdom = None

            seed_id = toks[i - 1]
            seeds.append(
                CazySeed(
                    seed_id=seed_id,
                    seed_source=src,
                    family=fam,
                    kingdom=kingdom,
                    raw_line=line,
                )
            )

    return seeds


def filter_ncbi_protein_seeds(seeds: Iterable[CazySeed]) -> List[CazySeed]:
    """
    Beholder bare seeds der:
      - source == "ncbi"
      - seed_id ser ut som NCBI/INSDC protein accession med versjon (.1)
    """
    out: List[CazySeed] = []
    for s in seeds:
        if s.seed_source != "ncbi":
            continue
        if _ncbi_protein_like.match(s.seed_id):
            out.append(s)
    return out


def unique_seed_ids(seeds: Iterable[CazySeed]) -> List[str]:
    seen: Set[str] = set()
    ids: List[str] = []
    for s in seeds:
        if s.seed_id not in seen:
            seen.add(s.seed_id)
            ids.append(s.seed_id)
    return ids
