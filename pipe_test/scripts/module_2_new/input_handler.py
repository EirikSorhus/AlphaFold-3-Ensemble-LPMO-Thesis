import csv
import re
from typing import Dict, List, Set, Tuple

class CAZyHandler:
    def __init__(self):
        # Regex to remove version numbers like " v1.0", " v2.0", " v1.1" etc.
        self.version_pattern = re.compile(r"\s+v\d+(\.\d+)*.*")

    def process_cazy_file(self, file_path):
        """
        Parses CAZy export file and returns data sorted for UniProt lookup.
        Expected columns: Family | Kingdom | Organism | Accession | Source
        """
        ncbi_ids = []
        jgi_groups = {} # { "Cleaned Organism": [ID1, ID2, ...] }

        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue

                # Indices based on: Family, Kingdom, Organism, Accession, Source
                organism = parts[2].strip()
                prot_id = parts[3].strip()
                source = parts[4].strip().lower()

                if source in ["ncbi", "genbank"]:
                    ncbi_ids.append(prot_id)
                
                elif source == "jgi":
                    # Clean organism name: Remove versions, then take first 2 words (Genus species)
                    clean_org = self.version_pattern.sub("", organism).strip()
                    # Heuristic: Use first two words to match UniProt 'organism_name'
                    # e.g. "Aureococcus anophagefferens clone 1984" -> "Aureococcus anophagefferens"
                    # e.g. "Bacillariophyceae sp. MOSAICH1_1" -> "Bacillariophyceae sp."
                    words = clean_org.split()
                    if len(words) >= 2:
                        clean_org = f"{words[0]} {words[1]}"
                    else:
                        clean_org = words[0] # Fallback for single word names

                    if clean_org not in jgi_groups:
                        jgi_groups[clean_org] = []
                    jgi_groups[clean_org].append(prot_id)

        return ncbi_ids, jgi_groups

    def generate_jgi_queries(self, jgi_groups):
        """
        Generates a list of UniProt search queries based on JGI groupings.
        Format: organism:"Name" AND (all:"ID1" OR all:"ID2")
        """
        queries = []
        for organism, ids in jgi_groups.items():
            # UniProt has a query length limit, split if too many IDs per organism
            chunk_size = 50 
            for i in range(0, len(ids), chunk_size):
                chunk = ids[i:i + chunk_size]
                # Removing "all:" prefix as it causes HTTP 400. Using quoted ID for exact phrase match in any field.
                id_part = " OR ".join([f'"{jid}"' for jid in chunk])
                # Use organism_name field instead of organism (which causes HTTP 400)
                query = f'organism_name:"{organism}" AND ({id_part})'
                queries.append(query)
        
        return queries

def parse_fasta_header(header):
    """
    Attempts to parse UniProt ID from header using Regex.
    Works for >sp|P12345|... and >tr|...
    """
    match = re.search(r"\|([A-Z0-9]+)\|", header)
    if match:
        return match.group(1)
    return None

def parse_fasta_file(file_path):
    """
    Generator that yields (header, sequence).
    """
    current_header = None
    current_seq = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            
            if line.startswith('>'):
                if current_header:
                    yield current_header, "".join(current_seq)
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_header:
            yield current_header, "".join(current_seq)


SEMICOLON_REQUIRED_HEADERS = {
    'protein name': 'Protein Name',
    'ec#': 'EC#',
    'reference': 'Reference',
    'organism': 'Organism',
    'genbank': 'GenBank',
    'uniprot': 'Uniprot',
    'pdb/3d': 'PDB/3D',
}


def _normalize_uniprot_cell(cell: str) -> List[str]:
    """Extract UniProt accessions; tolerant to spaces/commas/semicolons."""
    if not cell:
        return []
    tokens = re.findall(r'[A-Z0-9]{6,10}', cell.upper())
    return tokens


def _normalize_genbank_cell(cell: str) -> List[str]:
    """Extract GenBank-like accessions, keeping version suffixes."""
    if not cell:
        return []
    toks: List[str] = []
    for raw in re.split(r'[;,\s]+', cell.strip()):
        tok = re.sub(r'[^A-Za-z0-9_.]', '', raw)
        if tok:
            toks.append(tok)
    return toks


class CharacterizedHandler:
    """Parse semicolon-delimited characterized CSV (CAZy layout)."""

    def __init__(self, file_path: str):
        self.file_path = file_path

    def process_characterized_file(self) -> Tuple[List[str], List[str], List[Dict]]:
        """
        Strictly parses a semicolon-delimited CSV with the CAZy characterized layout.
        Returns (uniprot_ids, genbank_ids, rows_meta).
        rows_meta keeps per-row context for later failure reporting and grouping.
        """

        # Quick delimiter sanity check: require semicolon, reject comma as delimiter
        with open(self.file_path, 'r', encoding='utf-8') as fh:
            first_line = fh.readline()
            if ';' not in first_line or ',' in first_line:
                raise ValueError("Characterized CSV must be semicolon-delimited (;) and must not contain commas.")

        with open(self.file_path, newline='', encoding='utf-8') as handle:
            reader = csv.reader(handle, delimiter=';')
            try:
                header = next(reader)
            except StopIteration:
                raise ValueError("Characterized CSV is empty; expected header row.")

            normalized = [h.strip().lower() for h in header]
            header_map = {name: idx for idx, name in enumerate(normalized)}

            missing = [col for col in SEMICOLON_REQUIRED_HEADERS if col not in header_map]
            if missing:
                raise ValueError(
                    "Characterized CSV header missing columns or wrong delimiter (expected ';'): "
                    + ", ".join(sorted(missing))
                )

            def _get(col_name: str, row: List[str]) -> str:
                idx = header_map[col_name]
                if idx < len(row):
                    return row[idx].strip()
                return ""

            uniprot_ids: List[str] = []
            genbank_ids: List[str] = []
            seen: Set[str] = set()
            rows_meta: List[Dict] = []

            for line_no, row in enumerate(reader, start=2):  # account for header line
                if not row or all(not cell.strip() for cell in row):
                    continue

                # Pad short rows if needed
                if len(row) < len(header):
                    row = row + [""] * (len(header) - len(row))

                uni_cell = _get('uniprot', row)
                gb_cell = _get('genbank', row)

                uni_tokens = _normalize_uniprot_cell(uni_cell)
                gb_tokens = _normalize_genbank_cell(gb_cell)

                for tok in uni_tokens:
                    if tok not in seen:
                        uniprot_ids.append(tok)
                        seen.add(tok)
                for tok in gb_tokens:
                    if tok not in seen:
                        genbank_ids.append(tok)
                        seen.add(tok)

                rows_meta.append({
                    "line_no": line_no,
                    "raw": ";".join(row),
                    "protein": _get('protein name', row),
                    "organism": _get('organism', row),
                    "ec": _get('ec#', row),
                    "reference": _get('reference', row),
                    "pdb": _get('pdb/3d', row),
                    "uniprot_ids": uni_tokens,
                    "genbank_ids": gb_tokens,
                })

        return uniprot_ids, genbank_ids, rows_meta

