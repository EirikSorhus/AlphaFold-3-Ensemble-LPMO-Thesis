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

def parse_uniprot_features(json_entry, interpro_domains=None):
    """
    Parser UniProt JSON-data med forbedret domenegrenser fra InterPro domain-type entries.
    
    Args:
        json_entry (dict): Rå JSON fra UniProt API.
        interpro_domains (list, optional): Liste med domener fra InterProClient.fetch_domains()
                                           (already deduplicated and H-adjusted)
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
            
    # --- 3. DOMENE OG CBM MAPPING (InterPro Priority) ---
    core_data = {"type": "Unknown", "start": None, "end": None}
    found_cbms = []
    
    # A. Use InterPro domains if available (already deduplicated and H-adjusted)
    if interpro_domains:
        # Separate LPMO core and CBM domains
        lpmo_candidates = []
        cbm_candidates = []
        
        for domain in interpro_domains:
            domain_name = domain.get("name", "")
            
            # Check if CBM
            is_cbm = any(re.search(p, domain_name, re.I) for p in CBM_PATTERNS)
            if is_cbm:
                cbm_candidates.append(domain)
                continue
            
            # Check if LPMO
            is_lpmo = any(re.search(p, domain_name, re.I) for p in LPMO_FAMILIES)
            if is_lpmo:
                lpmo_candidates.append(domain)
        
        # Pick LPMO core: prefer one that starts closest to signal end
        if lpmo_candidates:
            lpmo_candidates.sort(key=lambda d: abs(d["start"] - (sig_end + 1)))
            best_lpmo = lpmo_candidates[0]
            source_label = f"{best_lpmo['source'].upper()}:{best_lpmo['model']}"
            ipr_label = f" (IPR:{best_lpmo['integrated_ipr']})" if best_lpmo.get('integrated_ipr') else ""
            core_data = {
                "type": f"{best_lpmo['name']} [{source_label}]{ipr_label}",
                "start": best_lpmo["start"],
                "end": best_lpmo["end"]
            }
        
        # Collect CBMs
        for cbm in cbm_candidates:
            source_label = f"{cbm['source'].upper()}:{cbm['model']}"
            cbm_desc = f"{cbm['name']} [{source_label}] ({cbm['start']}-{cbm['end']})"
            found_cbms.append(cbm_desc)
    
    # B. Fallback: Parse UniProt internal features if InterPro didn't provide domains
    if core_data["start"] is None:
        for f in features:
            desc = f.get('description') or str(f.get('note', ''))
            start = f['location']['start']['value']
            end = f['location']['end']['value']
            
            # CBM Detection
            if any(re.search(p, desc, re.I) for p in CBM_PATTERNS):
                found_cbms.append(f"{desc} ({start}-{end})")
                
            # LPMO Core (UniProt internal)
            is_lpmo = any(re.search(p, desc, re.I) for p in LPMO_FAMILIES)
            if is_lpmo:
                current_is_better = core_data["start"] is not None and abs(core_data["start"] - (sig_end + 1)) <= 1
                new_is_n_term = abs(start - (sig_end + 1)) <= 1
                if not current_is_better or new_is_n_term:
                    core_data = {"type": desc, "start": start, "end": end}

    # C. Fallback: Implicit domain if H1 verified but no domain found
    if core_data["start"] is None and h1_verified:
        core_data["start"] = sig_end + 1
        core_data["type"] = "Inferred (H1)"
        if core_data["end"] is None:
            core_data["end"] = len(sequence)

    # --- 4. RE-CHECK H1 FOR CASES WHERE sig_end=0 BUT LPMO DOMAIN FOUND ---
    if sig_end == 0 and core_data["start"] is not None and not h1_verified:
        lpmo_start_0based = core_data["start"] - 1  # Convert 1-based to 0-based
        if 0 <= lpmo_start_0based < len(sequence):
            found_aa = sequence[lpmo_start_0based]
            if found_aa.upper() == 'H':
                h1_verified = True

    # --- 5. INTERPRO IDS LIST ---
    # Merge UniProt refs and InterPro domain IDs
    uni_refs = [
        db.get('id') for db in json_entry.get('uniProtKBCrossReferences', [])
        if db.get('database') == 'InterPro'
    ]
    if interpro_domains:
        for domain in interpro_domains:
            if domain.get("integrated_ipr"):
                uni_refs.append(domain["integrated_ipr"])
            if domain.get("entry_id"):
                uni_refs.append(domain["entry_id"])
    
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
