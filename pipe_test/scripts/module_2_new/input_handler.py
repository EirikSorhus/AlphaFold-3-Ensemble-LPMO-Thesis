import re

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

