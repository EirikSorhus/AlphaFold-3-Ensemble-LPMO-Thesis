import re
import logging

# Oppsett av logging
logger = logging.getLogger(__name__)

# Utvidet liste med mønstre for å fange opp alle typer bindingsmoduler og LPMO-familier
CBM_PATTERNS = [
    r"CBM\d+",                          # Standard CAZy (CBM1, CBM2, etc.)
    r"Chitin[- ]binding",               # Beskrivende navn
    r"Cellulose[- ]binding",            # Beskrivende navn
    r"Carbohydrate[- ]binding",         # Generell betegnelse
    r"Glycan[- ]binding",               # Spesifikk for sukkerkjeder
    r"WSC[- ]domain",                   # Cellevegg-interaksjon
    r"Lectin",                          # Sukkerbindende
    r"LysM",                            # Kitin-bindende (vanlig i bakterier)
    r"X8[- ]domain",                    # Assosiert med CAZymes
    r"CBD",                             # Legacy: Cellulose Binding Domain
    r"F5/8 type C domain",              # InterPro alias for visse CBMs
    r"Dockerin",                        # Cellulosom-relatert
    r"Fibronectin type III",            # Spacer/binding domene
    r"Chitinase_insertion_domain"       # Eldre annoteringer
]

LPMO_FAMILIES = [
    r"AA9", r"AA10", r"AA11", r"AA13", r"AA14", 
    r"AA15", r"AA16", r"AA17", r"LPMO", 
    r"Lytic polysaccharide monooxygenase"
]

def parse_uniprot_features(json_entry):
    """
    Parser UniProt JSON-data for å trekke ut LPMO-spesifikke features.
    """
    acc = json_entry.get('primaryAccession', 'Unknown')
    sequence = json_entry.get('sequence', {}).get('value', '')
    features = json_entry.get('features', [])
    
    # --- 1. SIGNAL PEPTIDE (Max-End logikk) ---
    signal_ends = [
        f['location']['end']['value'] 
        for f in features if f.get('type') == 'Signal'
    ]
    # Vi velger det lengste signalpeptidet for å sikre ren N-terminal
    sig_end = max(signal_ends) if signal_ends else 0
    
    # --- 2. H1-VALIDERING (LPMO-spesifikk) ---
    # Sjekker om aminosyren rett etter signalpeptidet er Histidin (H)
    h1_verified = False
    found_aa = None
    if sig_end > 0 and len(sequence) > sig_end:
        found_aa = sequence[sig_end] # Indexing 0-based: sig_end er posisjonen til AA (sig_end+1)
        if found_aa.upper() == 'H':
            h1_verified = True
    
    # --- 3. DOMENE OG CBM MAPPING ---
    core_data = {"type": "Unknown", "start": None, "end": None}
    found_cbms = []
    
    for f in features:
        f_type = f.get('type')
        # Hent beskrivelse fra 'description', 'note' eller 'comment'
        desc = ""
        if 'description' in f:
            desc = f['description']
        elif 'note' in f: # Noen ganger ligger navnet i notes
            desc = str(f['note'])
            
        start = f['location']['start']['value']
        end = f['location']['end']['value']
        
        # Sjekk for LPMO kjerne (skal starte rett etter signalpeptidet)
        # Vi tillater et slingringsmonn på +/- 1 aminosyre
        if abs(start - (sig_end + 1)) <= 1:
            if any(re.search(p, desc, re.I) for p in LPMO_FAMILIES):
                core_data = {"type": desc, "start": start, "end": end}
        
        # Sjekk for CBMs ved bruk av den utvidede mønsterlisten
        if any(re.search(p, desc, re.I) for p in CBM_PATTERNS):
            found_cbms.append(f"{desc} ({start}-{end})")

    # --- 4. INTERPRO IDS ---
    interpro_ids = [
        db.get('id') for db in json_entry.get('uniProtKBCrossReferences', [])
        if db.get('database') == 'InterPro'
    ]

    return {
        "UniProt_ID": acc,
        "Signal_End": sig_end,
        "H1_Verified": h1_verified,
        "H1_AminoAcid": found_aa,
        "LPMO_Core_Type": core_data["type"],
        "LPMO_Core_Start": core_data["start"],
        "LPMO_Core_End": core_data["end"],
        "Binding_Modules": "; ".join(found_cbms) if found_cbms else "None",
        "InterPro_IDs": "; ".join(interpro_ids) if interpro_ids else "None"
    }

if __name__ == "__main__":
    # Enkel test-logikk
    print("Feature Parser modul lastet med", len(CBM_PATTERNS), "CBM mønstre.")
