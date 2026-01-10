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
    # Modern CAZy names
    r"AA9", r"AA10", r"AA11", r"AA13", 
    r"AA14", r"AA15", r"AA16", r"AA17",
    
    # Legacy names (Critical)
    r"GH61", r"Glycoside[ -]hydrolase[ -]family[ -]61",
    r"CBM33", r"Carbohydrate[ -]binding[ -]module[ -]33",
    r"GH18", r"Glycoside[ -]hydrolase[ -]family[ -]18", # AA15 association
    
    # InterPro IDs (Highly specific, can appear in detailed description)
    r"IPR005123", r"IPR004302", r"IPR027003", r"IPR026998",
    r"IPR031548", r"IPR031154", r"IPR043232", r"IPR044679",
    
    # Generic
    r"LPMO", r"Lytic[ -]polysaccharide[ -]monooxygenase",
    r"Fusolin", r"Spindle" # AA15 aliases
]

def parse_uniprot_features(json_entry, interpro_entries=None):
    """
    Parser UniProt JSON-data med mulighet for å forbedre domenegrenser ved hjelp av InterPro-data.
    
    Args:
        json_entry (dict): Rå JSON fra UniProt API.
        interpro_entries (list, optional): Liste med domener fra InterProClient.
    """
    acc = json_entry.get('primaryAccession', 'Unknown')
    sequence = json_entry.get('sequence', {}).get('value', '')
    features = json_entry.get('features', [])
    
    # --- 1. SIGNAL PEPTIDE (Max-End logikk) ---
    signal_ends = [
        f['location']['end']['value'] 
        for f in features if f.get('type') == 'Signal'
    ]
    sig_end = max(signal_ends) if signal_ends else 0
    
    # --- 2. H1-VALIDERING (LPMO-spesifikk) ---
    h1_verified = False
    found_aa = None
    if sig_end > 0 and len(sequence) > sig_end:
        found_aa = sequence[sig_end] 
        if found_aa.upper() == 'H':
            h1_verified = True
            
    # --- 3. DOMENE OG CBM MAPPING (Hybrid Approach) ---
    core_data = {"type": "Unknown", "start": None, "end": None}
    found_cbms = []
    
    # A. Parse UniProt internal features (Legacy/Fallback)
    for f in features:
        desc = f.get('description') or str(f.get('note', ''))
        start = f['location']['start']['value']
        end = f['location']['end']['value']
        
        # CBM Detection
        if any(re.search(p, desc, re.I) for p in CBM_PATTERNS):
            found_cbms.append(f"{desc} ({start}-{end})")
            
        # LPMO Core (UniProt internal) - only use if we don't prefer InterPro later
        is_lpmo = any(re.search(p, desc, re.I) for p in LPMO_FAMILIES)
        if is_lpmo:
             current_is_better = core_data["start"] is not None and abs(core_data["start"] - (sig_end + 1)) <= 1
             new_is_n_term = abs(start - (sig_end + 1)) <= 1
             if not current_is_better or new_is_n_term:
                 core_data = {"type": desc, "start": start, "end": end}

    # B. Apply InterPro Data (Priority)
    if interpro_entries:
        # Filter for relevant domains (Families or Domains, ignore minimal repeats if better exist)
        # We look for something that contains "LPMO", "Polysaccharide monooxygenase", "Auxiliary Activity" or specific ID
        # Or simply the largest domain that starts near N-term.
        
        best_ipr = None
        
        # 1. Look for specific LPMO text match in name
        lpmo_candidates = [
            e for e in interpro_entries 
            if any(re.search(p, e['name'], re.I) for p in LPMO_FAMILIES) or 
               any(re.search(p, e['id'], re.I) for p in LPMO_FAMILIES) # Check ID too
        ]
        
        if lpmo_candidates:
            # Pick one that starts earliest (closest to signal) but is after signal start
            # Sort by start position
            lpmo_candidates.sort(key=lambda x: x['start'])
            best_ipr = lpmo_candidates[0]
        
        if best_ipr:
            # LOGIC REFINEMENT:
            # If H1 is verified, we FORCE the start to be exactly at H1 (SigEnd + 1).
            # Because InterPro often includes the Signal Peptide region (e.g. starts at 1, 11) 
            # or starts inside the domain (e.g. 25).
            
            ipr_start = best_ipr['start']
            ipr_end = best_ipr['end']
            
            final_start = ipr_start
            
            if h1_verified:
                final_start = sig_end + 1
            else:
                # If no H1, but InterPro starts BEFORE signal end, we must trim it
                if ipr_start <= sig_end:
                    final_start = sig_end + 1
            
            core_data = {
                "type": f"{best_ipr['name']} ({best_ipr['id']})",
                "start": final_start,
                "end": ipr_end
            }

    # C. Fallback: Implicit
    if core_data["start"] is None and h1_verified:
        core_data["start"] = sig_end + 1
        core_data["type"] = "Inferred (H1)"
        if core_data["end"] is None:
            core_data["end"] = len(sequence)

    # --- 4. INTERPRO IDS LIST ---
    # Merge both UniProt refs and fetched refs
    uni_refs = [
        db.get('id') for db in json_entry.get('uniProtKBCrossReferences', [])
        if db.get('database') == 'InterPro'
    ]
    if interpro_entries:
        uni_refs.extend([e['id'] for e in interpro_entries])
    
    unique_refs = sorted(list(set(uni_refs)))

    return {
        "UniProt_ID": acc,
        "Signal_End": int(sig_end),
        "H1_Verified": h1_verified,
        "H1_AminoAcid": found_aa,
        "LPMO_Core_Type": core_data["type"],
        "LPMO_Core_Start": int(core_data["start"]) if core_data["start"] is not None else None,
        "LPMO_Core_End": int(core_data["end"]) if core_data["end"] is not None else None,
        "Binding_Modules": "; ".join(found_cbms) if found_cbms else "None",
        "InterPro_IDs": "; ".join(unique_refs) if unique_refs else "None"
    }


if __name__ == "__main__":
    # Enkel test-logikk
    print("Feature Parser modul lastet med", len(CBM_PATTERNS), "CBM mønstre.")
