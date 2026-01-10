import re

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

def format_jgi_query(organism, ids):
    """
    Removes version from organism and builds query string.
    """
    clean_org = re.sub(r" v\d+\.\d+.*", "", organism).strip()
    id_part = " OR ".join([f'all:"{jid}"' for jid in ids])
    return f'organism:"{clean_org}" AND ({id_part})'
