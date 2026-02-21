import argparse
import logging
import os
import sys
import pandas as pd
import json
import requests
import re
from datetime import datetime
from tqdm import tqdm
import hashlib

# Import local modules - robust check for package vs script execution
try:
    from .input_handler import parse_fasta_file, parse_fasta_header, CAZyHandler, CharacterizedHandler
    from .uniprot_client import UniProtClient
    from .interpro_client import InterProClient
    from .feature_parser import parse_uniprot_features
    from .blast_client import run_blast_search
except ImportError:
    from input_handler import parse_fasta_file, parse_fasta_header, CAZyHandler, CharacterizedHandler
    from uniprot_client import UniProtClient
    from interpro_client import InterProClient
    from feature_parser import parse_uniprot_features
    from blast_client import run_blast_search

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Constants for CAZy download
CAZY_BASE_URL = "https://www.cazy.org"
USER_AGENT = "LPMO-Pipeline/1.0 (eirik.sorhus@nmbu.no)"
REQUEST_TIMEOUT = 30

def main():
    parser = argparse.ArgumentParser(description="Module 2: Metadata Enrichment")
    parser.add_argument(
        "--mode",
        choices=['fasta', 'cazy', 'list', 'characterized'],
        required=True,
        help=(
            "Input mode: fasta (FASTA headers with UniProt IDs), cazy (CAZy export), "
            "list (plain ID list), characterized (semikolonskilt CSV med samme oppsett som CAZy characterized: "
            "Protein Name;EC#;Reference;Organism;GenBank;Uniprot;PDB/3D)."
        ),
    )
    parser.add_argument("--input", required=True, help="Input file path")
    parser.add_argument("--output_dir", default="data", help="Base output directory (default: data)")
    parser.add_argument("--allow-ncbi-fallback", action="store_true", help="Enable NCBI fallback (not currently active)")
    parser.add_argument("--allow-sequence-search", action="store_true", help="Allow BLAST sequence search for unknown FASTA headers")
    parser.add_argument("--max-sequence-searches", type=int, default=20, help="Maximum number of sequences to search via BLAST (default: 20)")
    parser.add_argument(
        "--cazy-family",
        default=None,
        help=(
            "CAZy family name written to the CAZy_family metadata column (e.g. AA9). "
            "Required for modes fasta, list, and characterized. "
            "Optional for cazy mode: if omitted, the family name is auto-detected from --input."
        ),
    )
    
    args = parser.parse_args()

    # --cazy-family is required for all modes except cazy (where it is auto-detected)
    if args.mode != 'cazy' and args.cazy_family is None:
        parser.error(f"--cazy-family is required when --mode is '{args.mode}'")
    
    # Setup Output Directories with Metadata
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base_dir = args.output_dir
    
    metadata_dir = os.path.join(base_dir, "metadata")
    sequences_dir = os.path.join(base_dir, "sequences")
    run_dir = os.path.join(base_dir, "run")
    
    for d in [metadata_dir, sequences_dir, run_dir]:
        os.makedirs(d, exist_ok=True)
        
    out_tsv = os.path.join(metadata_dir, f"metadata_expanded_{timestamp}.tsv")
    out_fasta = os.path.join(sequences_dir, f"all_sequences_{timestamp}.fasta")
    out_failed = os.path.join(run_dir, f"failed_ids_{timestamp}.txt")
    out_run_json = os.path.join(run_dir, f"run_metadata_{timestamp}.json")

    # Use absolute cache paths resolved relative to this script file,
    # so cache is reused regardless of the working directory (e.g. sbatch vs interactive).
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    client = UniProtClient(cache_path=os.path.join(_script_dir, "discovery_cache.json"))
    ip_client = InterProClient(cache_path=os.path.join(_script_dir, "interpro_cache.json"))

    successful_entries = []
    failed_ids = []
    raw_fasta_entries = []
    start_time = datetime.now()
    
    logger.info(f"Starting Module 2 in {args.mode} mode. Output will be timestamped: {timestamp}")
    
    # --- MODE: CAZy ---
    if args.mode == 'cazy':
        # Resolve one or more families from --input (comma-separated family names or a single file path)
        input_tokens = [t.strip() for t in args.input.split(',')]
        families_to_process = []  # list of (local_file_path, family_name)

        for token in input_tokens:
            if os.path.exists(token):
                # Existing file: derive family name from filename
                fname = os.path.splitext(os.path.basename(token))[0].upper()
                families_to_process.append((token, fname))
            else:
                family_match = re.match(r'^((AA|GH|GT|PL|CE|CBM)\d+)(\.txt)?$', token, re.I)
                if family_match:
                    family = family_match.group(1).upper()
                    raw_dir = os.path.join(args.output_dir, "cazy_raw")
                    os.makedirs(raw_dir, exist_ok=True)
                    target_file = os.path.join(raw_dir, f"{family}.txt")
                    if os.path.exists(target_file):
                        logger.info(f"Found existing CAZy file for {family} at {target_file}")
                    else:
                        url = f"{CAZY_BASE_URL}/IMG/cazy_data/{family}.txt"
                        logger.info(f"Downloading CAZy data for {family} from {url}...")
                        try:
                            resp = requests.get(url, headers={"User-Agent": USER_AGENT}, timeout=REQUEST_TIMEOUT)
                            resp.raise_for_status()
                            if not resp.text.strip():
                                logger.error(f"Downloaded file for {family} is empty.")
                                sys.exit(1)
                            with open(target_file, 'w', encoding='utf-8') as f:
                                f.write(resp.text)
                            logger.info(f"Saved to {target_file}")
                        except Exception as e:
                            logger.error(f"Failed to download CAZy data for {family}: {e}")
                            sys.exit(1)
                    families_to_process.append((target_file, family))
                else:
                    logger.error(f"Input '{token}' is not a valid file path or CAZy family ID.")
                    sys.exit(1)

        for cazy_file, current_family in families_to_process:
            # --cazy-family overrides auto-detected family label for all rows if provided
            cazy_family_label = args.cazy_family if args.cazy_family else current_family

            logger.info(f"Reading CAZy input file: {cazy_file} (CAZy_family label: {cazy_family_label})")
            cazy_handler = CAZyHandler()
            ncbi_ids, jgi_groups = cazy_handler.process_cazy_file(cazy_file)

            logger.info(f"Parsers found: {len(ncbi_ids)} NCBI IDs and {len(jgi_groups)} JGI organism groups.")

            # 1. Process NCBI IDs (Batch fetch)
            if ncbi_ids:
                chunk_size = client.batch_size
                for i in tqdm(range(0, len(ncbi_ids), chunk_size), desc=f"Fetching NCBI-linked IDs ({current_family})"):
                    batch = ncbi_ids[i:i+chunk_size]
                    results = client.fetch_batch(batch)

                    for data in results:
                        acc = data.get('primaryAccession')
                        seq = data.get('sequence', {}).get('value', '')
                        # Extract signal peptide end from UniProt features
                        signal_ends = [f['location']['end']['value'] for f in data.get('features', []) if f.get('type') == 'Signal']
                        sig_end = max(signal_ends) if signal_ends else 0
                        # Fetch InterPro domains with sequence and signal_end for H-adjustment
                        ipr_domains = ip_client.fetch_domains(acc, sequence=seq, signal_end=sig_end)

                        features = parse_uniprot_features(data, ipr_domains)
                        uid = features['UniProt_ID']
                        features['Protein_Name'] = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')
                        features['Organism'] = data.get('organism', {}).get('scientificName', 'Unknown')
                        ecs = [db.get('id') for db in data.get('uniProtKBCrossReferences', []) if db.get('database') == 'EC']
                        features['EC_Number'] = ";".join(ecs)
                        features['Match_Status'] = "Success_CAZy_NCBI"
                        features['CAZy_family'] = cazy_family_label
                        successful_entries.append(features)

                        seq = data.get('sequence', {}).get('value', '')
                        new_header = f">UniProtID|{uid}|{features['Organism']}|{features['Protein_Name']}"
                        raw_fasta_entries.append((new_header, seq))

            # 2. Process JGI Groups (Query)
            if jgi_groups:
                queries = cazy_handler.generate_jgi_queries(jgi_groups)
                logger.info(f"Generated {len(queries)} grouped queries for JGI entries.")

                for query in tqdm(queries, desc=f"Searching JGI Groups ({current_family})"):
                    results = client.search_by_query(query)

                    if not results:
                        failed_ids.append(f"JGI_Query_Failed: {query}")

                    for data in results:
                        acc = data.get('primaryAccession')
                        seq = data.get('sequence', {}).get('value', '')
                        # Extract signal peptide end
                        signal_ends = [f['location']['end']['value'] for f in data.get('features', []) if f.get('type') == 'Signal']
                        sig_end = max(signal_ends) if signal_ends else 0
                        ipr_domains = ip_client.fetch_domains(acc, sequence=seq, signal_end=sig_end)

                        features = parse_uniprot_features(data, ipr_domains)
                        uid = features['UniProt_ID']
                        features['Protein_Name'] = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')
                        features['Organism'] = data.get('organism', {}).get('scientificName', 'Unknown')
                        ecs = [db.get('id') for db in data.get('uniProtKBCrossReferences', []) if db.get('database') == 'EC']
                        features['EC_Number'] = ";".join(ecs)
                        features['Match_Status'] = "Success_CAZy_JGI"
                        features['CAZy_family'] = cazy_family_label
                        successful_entries.append(features)

                        new_header = f">UniProtID|{uid}|{features['Organism']}|{features['Protein_Name']}"
                        raw_fasta_entries.append((new_header, seq))

    # --- MODE: CHARACTERIZED CSV ---
    elif args.mode == 'characterized':
        try:
            handler = CharacterizedHandler(args.input)
            uni_ids, genbank_ids, rows_meta = handler.process_characterized_file()
        except ValueError as exc:
            logger.error(str(exc))
            sys.exit(1)

        # Split rows into those with UniProt IDs and those without
        rows_with_uniprot = [r for r in rows_meta if r['uniprot_ids']]
        rows_without_uniprot = [r for r in rows_meta if not r['uniprot_ids'] and r['genbank_ids']]
        rows_no_ids = [r for r in rows_meta if not r['uniprot_ids'] and not r['genbank_ids']]

        # Ordered unique UniProt IDs from rows that have UniProt
        uni_like_ordered = []
        seen_uni = set()
        for r in rows_with_uniprot:
            for uid in r['uniprot_ids']:
                if uid not in seen_uni:
                    uni_like_ordered.append(uid)
                    seen_uni.add(uid)

        logger.info(
            f"Characterized CSV parsed: {len(uni_like_ordered)} UniProt IDs (rows with UniProt), "
            f"{len(rows_without_uniprot)} rows without UniProt but with GenBank/RefSeq, {len(rows_no_ids)} rows without IDs."
        )

        features_by_uid = {}
        seq_groups = {}  # seq string -> list of ids in order
        id_to_sequence = {}  # Map UniProt ID -> amino acid sequence string
        failed_id_set = set()
        success_id_set = set()

        def choose_primary(id_list):
            for cid in id_list:
                if re.fullmatch(r"[A-Z0-9]{6,10}", cid):
                    return cid
            return id_list[0] if id_list else None

        # Helper: fetch UniProt entries for UniProt IDs
        def fetch_uniprot_ids(id_list):
            chunk_size = client.batch_size
            for i in tqdm(range(0, len(id_list), chunk_size), desc="Fetching UniProt IDs"):
                batch = id_list[i:i+chunk_size]
                results = client.fetch_batch(batch)
                for data in results:
                    acc = data.get('primaryAccession')
                    if not acc:
                        continue
                    seq = data.get('sequence', {}).get('value', '')
                    if not seq:
                        failed_id_set.add(acc)
                        continue

                    # Extract signal peptide end
                    signal_ends = [f['location']['end']['value'] for f in data.get('features', []) if f.get('type') == 'Signal']
                    sig_end = max(signal_ends) if signal_ends else 0
                    ipr_domains = ip_client.fetch_domains(acc, sequence=seq, signal_end=sig_end)

                    features = parse_uniprot_features(data, ipr_domains)
                    uid = features['UniProt_ID']
                    features['Protein_Name'] = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')
                    features['Organism'] = data.get('organism', {}).get('scientificName', 'Unknown')
                    ecs = [db.get('id') for db in data.get('uniProtKBCrossReferences', []) if db.get('database') == 'EC']
                    features['EC_Number'] = ";".join(ecs)
                    features['Match_Status'] = "Success_Characterized"
                    features_by_uid[uid] = {**features, "_seq": seq, "CAZy_family": args.cazy_family}
                    success_id_set.add(uid)

                    if seq not in seq_groups:
                        seq_groups[seq] = []
                    seq_groups[seq].append(uid)
                    id_to_sequence[uid] = seq

        # Helper: map GenBank/RefSeq IDs to UniProt via crossref query
        def map_non_uniprot_to_uniprot(id_token):
            token = id_token.strip()
            # RefSeq proteins: NP_, XP_, YP_, WP_, ZP_
            is_refseq = bool(re.match(r"^(NP_|XP_|YP_|WP_|ZP_)", token, re.IGNORECASE))
            candidates = []
            if is_refseq:
                candidates.append(f"xref:RefSeq:{token}")
                if '.' in token:
                    candidates.append(f"xref:RefSeq:{token.split('.')[0]}")
            else:
                # GenBank protein_id via EMBL-CDS cross-ref
                candidates.append(f"xref:EMBL-CDS:{token}")
                if '.' in token:
                    candidates.append(f"xref:EMBL-CDS:{token.split('.')[0]}")

            for q in candidates:
                results = client.search_by_query(q)
                if results:
                    return results[0]
            return None

        # Helper: NCBI FASTA fallback
        def fetch_ncbi_fasta(genbank_acc):
            try:
                params = {
                    "db": "protein",
                    "id": genbank_acc,
                    "rettype": "fasta",
                    "retmode": "text",
                }
                resp = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=params, timeout=30)
                if not resp.ok:
                    return None, None, None
                lines = [ln.strip() for ln in resp.text.splitlines() if ln.strip()]
                if not lines or not lines[0].startswith('>'):
                    return None, None, None
                header = lines[0]
                seq = "".join(lines[1:])
                org = None
                m = re.search(r"\[(.*?)\]\s*$", header)
                if m:
                    org = m.group(1)
                prot = header.lstrip('>').split(' ', 1)[1] if ' ' in header else genbank_acc
                return seq, org or "Unknown", prot or "Unknown"
            except Exception:
                return None, None, None

        # Fetch UniProt IDs from rows that already have UniProt
        if uni_like_ordered:
            fetch_uniprot_ids(uni_like_ordered)
        
        # Fallback: Try InterPro SignalP for entries with sig_end=0
        for uid, feat_data in list(features_by_uid.items()):
            if feat_data.get('Signal_End') == 0:
                ipr_signal = ip_client.fetch_signal_peptide(uid)
                if ipr_signal and ipr_signal.get('end') is not None:
                    # Re-parse features with corrected signal_end
                    seq = feat_data.get('_seq', '')
                    sig_end_corrected = ipr_signal['end']
                    
                    # Validate sig_end_corrected is within sequence bounds
                    if not isinstance(sig_end_corrected, int) or sig_end_corrected < 0:
                        logger.warning(f"Invalid signal_end from InterPro for {uid}: {sig_end_corrected}")
                        continue
                    
                    # Fetch domains again with correct signal_end
                    try:
                        ipr_domains = ip_client.fetch_domains(uid, sequence=seq, signal_end=sig_end_corrected)
                    except Exception:
                        ipr_domains = []
                    # Find original UniProt entry to re-parse
                    # We don't have it cached, so just update signal_end manually
                    feat_data['Signal_End'] = sig_end_corrected
                    # Re-check H1 with corrected signal end
                    if len(seq) > sig_end_corrected:
                        found_aa = seq[sig_end_corrected]
                        feat_data['H1_AminoAcid'] = found_aa
                        # ALWAYS update H1_Verified to match the actual amino acid (fixes consistency bug)
                        feat_data['H1_Verified'] = (found_aa.upper() == 'H')
                    logger.info(f"Applied InterPro SignalP fallback for {uid}: sig_end={sig_end_corrected}")

        mapped_non_uni_ids = set()

        # Map rows without UniProt: try GenBank/RefSeq -> UniProt; optional NCBI fallback
        for row in tqdm(rows_without_uniprot, desc="Mapping GenBank/RefSeq IDs"):
            mapped_any = False
            for nid in row['genbank_ids']:
                mapped = map_non_uniprot_to_uniprot(nid)
                if mapped:
                    acc = mapped.get('primaryAccession')
                    seq = mapped.get('sequence', {}).get('value', '')
                    if not acc or not seq:
                        failed_id_set.add(nid)
                        continue
                    
                    # Extract signal peptide end
                    signal_ends = [f['location']['end']['value'] for f in mapped.get('features', []) if f.get('type') == 'Signal']
                    sig_end = max(signal_ends) if signal_ends else 0
                    
                    try:
                        ipr_domains = ip_client.fetch_domains(acc, sequence=seq, signal_end=sig_end)
                    except Exception:
                        ipr_domains = []

                    feats = parse_uniprot_features(mapped, ipr_domains)
                    uid = feats['UniProt_ID']
                    feats['Protein_Name'] = mapped.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')
                    feats['Organism'] = mapped.get('organism', {}).get('scientificName', 'Unknown')
                    ecs = [db.get('id') for db in mapped.get('uniProtKBCrossReferences', []) if db.get('database') == 'EC']
                    feats['EC_Number'] = ";".join(ecs)
                    feats['Match_Status'] = "Success_Characterized_Mapped"
                    features_by_uid[uid] = {**feats, "_seq": seq, "CAZy_family": args.cazy_family}
                    success_id_set.add(uid)
                    mapped_non_uni_ids.add(nid)
                    mapped_any = True
                    if seq not in seq_groups:
                        seq_groups[seq] = []
                    id_to_sequence[uid] = seq
                    seq_groups[seq].append(uid)
                else:
                    if args.allow_ncbi_fallback:
                        seq, org, prot = fetch_ncbi_fasta(nid)
                        if seq:
                            new_header = f">UniProtIDs|{nid}|{org}|{prot}"
                            raw_fasta_entries.append((new_header, seq))
                            md = {
                                'UniProt_ID': nid,
                                'InterPro_IDs': '',
                                'Protein_Name': prot,
                                'Organism': org,
                                'EC_Number': '',
                                'Signal_End': None,
                                'LPMO_Core_Start': None,
                                'LPMO_Core_End': None,
                                'Binding_Modules': None,
                                'H1_Verified': None,
                                'H1_AminoAcid': None,
                                'Match_Status': 'NCBI_Fallback',
                                'Sequence_Group': nid,
                                'CoAccessions': '',
                                'CAZy_family': args.cazy_family,
                            }
                            successful_entries.append(md)
                            success_id_set.add(nid)
                            mapped_non_uni_ids.add(nid)
                            mapped_any = True
                        else:
                            failed_id_set.add(nid)
                    else:
                        failed_id_set.add(nid)

            if not mapped_any:
                failed_ids.append(f"No valid IDs (line {row['line_no']}): {row['raw']}")

        # Missing IDs (not returned) + empty sequence cases
        # Build metadata rows with Sequence_Group and CoAccessions; dedupe FASTA
        for seq, ids_in_group in seq_groups.items():
            if not ids_in_group:
                continue
            primary = choose_primary(ids_in_group)
            if not primary:
                primary = ids_in_group[0]
            co_ids = [cid for cid in ids_in_group if cid != primary]

            # FASTA header uses primary + coaccessions
            primary_feat = features_by_uid.get(primary)
            org = primary_feat.get('Organism', 'Unknown') if primary_feat else 'Unknown'
            pname = primary_feat.get('Protein_Name', 'Unknown') if primary_feat else 'Unknown'
            header_ids = [primary] + co_ids
            new_header = f">UniProtIDs|{';'.join(header_ids)}|{org}|{pname}"
            raw_fasta_entries.append((new_header, seq))

            # Metadata rows
            for cid in ids_in_group:
                feats = features_by_uid.get(cid)
                if not feats:
                    continue
                feats = feats.copy()
                feats.pop('_seq', None)
                feats['Sequence_Group'] = primary
                feats['CoAccessions'] = ";".join([i for i in ids_in_group if i != cid])
                successful_entries.append(feats)

        for mid in sorted(failed_id_set):
            failed_ids.append(f"{mid} (Not found or unmapped)")

        # Detect CSV rows that split into multiple sequences
        row_splits = []  # List of {"line_no": int, "uniprot_ids": list, "sequences_produced": int}
        if 'rows_meta' in locals():
            for row in rows_meta:
                row_uniprot_ids = row.get('uniprot_ids', [])
                if not row_uniprot_ids:
                    # Skip rows without UniProt IDs (GenBank-only rows)
                    continue
                
                # Collect unique sequences for all UniProt IDs in this row
                sequences_in_row = set()
                valid_ids_in_row = []
                for uid in row_uniprot_ids:
                    if uid in id_to_sequence:
                        sequences_in_row.add(id_to_sequence[uid])
                        valid_ids_in_row.append(uid)
                
                # If more than 1 unique sequence, this row was split
                if len(sequences_in_row) > 1:
                    row_splits.append({
                        "line_no": row['line_no'],
                        "uniprot_ids": valid_ids_in_row,
                        "sequences_produced": len(sequences_in_row)
                    })

    # --- MODE: FASTA ---
    elif args.mode == 'fasta':
        entries = []
        unknown_headers = []
        
        logger.info(f"Reading FASTA file: {args.input}")
        for header, seq in parse_fasta_file(args.input):
            uid = parse_fasta_header(header)
            if uid:
                entries.append((uid, header, seq))
            else:
                unknown_headers.append((header, seq))
        
        # Unique IDs to fetch
        known_ids = list(set([e[0] for e in entries]))
        logger.info(f"Found {len(known_ids)} unique UniProt IDs in headers.")
        
        # Batch Fetch
        metadata_results = []
        # Chunking handled inside fetched batch? No, UniProtClient.fetch_batch handles list, but better to chunk here if list is huge?
        # The client code says "fetch_batch(self, ids)". And it has "batch_size" in init.
        # But fetch_batch implementation doesn't loop over batches!
        # It says "if len(id_list) > 1: ... bisection".
        # It takes a list of IDs.
        # So we should feed it in chunks.
        
        chunk_size = client.batch_size
        for i in tqdm(range(0, len(known_ids), chunk_size), desc="Fetching Metadata"):
            batch = known_ids[i:i+chunk_size]
            results = client.fetch_batch(batch)
            metadata_results.extend(results)
            
        metadata_map = {res['primaryAccession']: res for res in metadata_results}
        
        # Process Known
        for uid, header, seq in entries:
            if uid in metadata_map:
                data = metadata_map[uid]
                
                # Extract signal peptide
                signal_ends = [f['location']['end']['value'] for f in data.get('features', []) if f.get('type') == 'Signal']
                sig_end = max(signal_ends) if signal_ends else 0
                
                # Fetch InterPro
                ipr_domains = ip_client.fetch_domains(uid, sequence=seq, signal_end=sig_end)
                
                features = parse_uniprot_features(data, ipr_domains)
                
                # Add extra fields standardisation
                features['Protein_Name'] = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')
                features['Organism'] = data.get('organism', {}).get('scientificName', 'Unknown')
                ecs = [db.get('id') for db in data.get('uniProtKBCrossReferences', []) if db.get('database') == 'EC']
                features['EC_Number'] = ";".join(ecs)
                features['Match_Status'] = "Success"
                features['CAZy_family'] = args.cazy_family

                successful_entries.append(features)

                # Normalized Header
                new_header = f">UniProtID|{uid}|{features['Organism']}|{features['Protein_Name']}"
                raw_fasta_entries.append((new_header, seq))
            else:
                failed_ids.append(f"{uid} (Not found in UniProt)")
        
        # Process Unknowns via BLAST (if enabled)
        if unknown_headers:
            if args.allow_sequence_search:
                count_ok = 0
                count_total = len(unknown_headers)
                limit = args.max_sequence_searches
                
                logger.info(f"Processing unknown headers via BLAST sequence match (Limit: {limit})...")
                
                # Iterate only up to the limit
                to_process = unknown_headers[:limit]
                skipped = unknown_headers[limit:]
                
                if skipped:
                   logger.warning(f"Skipping {len(skipped)} sequences due to --max-sequence-searches limit ({limit}).")
                   for h, s in skipped:
                       failed_ids.append(f"SkippedSequence: {h} (Max Limit Reached)")

                for header, seq in tqdm(to_process, desc="Sequence Matching"):
                    # Use BLAST fallback for sequence matching
                    match_id = run_blast_search(seq)
                    
                    if match_id:
                        # Fetch metadata for the newly found ID
                        results = client.fetch_batch([match_id])
                        
                        if results:
                            data = results[0]
                            acc = data.get('primaryAccession')
                            seq = data.get('sequence', {}).get('value', '')
                            # Extract signal peptide
                            signal_ends = [f['location']['end']['value'] for f in data.get('features', []) if f.get('type') == 'Signal']
                            sig_end = max(signal_ends) if signal_ends else 0
                            
                            # Precision fetch
                            ipr_domains = ip_client.fetch_domains(acc, sequence=seq, signal_end=sig_end)
                            
                            features = parse_uniprot_features(data, ipr_domains)
                            uid = features['UniProt_ID']
                            
                            features['Protein_Name'] = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')
                            features['Organism'] = data.get('organism', {}).get('scientificName', 'Unknown')
                            ecs = [db.get('id') for db in data.get('uniProtKBCrossReferences', []) if db.get('database') == 'EC']
                            features['EC_Number'] = ";".join(ecs)
                            features['Match_Status'] = "Success_SeqMatch"
                            features['CAZy_family'] = args.cazy_family

                            successful_entries.append(features)

                            # Add to fasta output with new header
                            new_header = f">UniProtID|{uid}|{features['Organism']}|{features['Protein_Name']}"
                            raw_fasta_entries.append((new_header, seq))
                            count_ok += 1
                        else:
                            failed_ids.append(f"UnknownHeader: {header} (BLAST ID {match_id} Fetch Failed)")

                    else:
                        failed_ids.append(f"UnknownHeader: {header} (No Seq Match)")
                        
            else:
                logger.info("Skipping unknown headers (enable --allow-sequence-search to perform BLAST)")
                for h, s in unknown_headers:
                    failed_ids.append(f"UnknownHeader: {h} (Sequence Search Disabled)")

    # --- MODE: LIST ---
    elif args.mode == 'list':
        with open(args.input, 'r') as f:
            ids = [line.strip() for line in f if line.strip()]
        
        logger.info(f"Loaded {len(ids)} IDs from list.")
        chunk_size = client.batch_size
        metadata_results = []
        for i in tqdm(range(0, len(ids), chunk_size), desc="Fetching Metadata"):
            batch = ids[i:i+chunk_size]
            results = client.fetch_batch(batch)
            metadata_results.extend(results)
            
        for data in metadata_results:
            acc = data.get('primaryAccession')
            seq = data.get('sequence', {}).get('value', '')
            # Extract signal peptide
            signal_ends = [f['location']['end']['value'] for f in data.get('features', []) if f.get('type') == 'Signal']
            sig_end = max(signal_ends) if signal_ends else 0
            
            ipr_domains = ip_client.fetch_domains(acc, sequence=seq, signal_end=sig_end)
            
            features = parse_uniprot_features(data, ipr_domains)
            uid = features['UniProt_ID']
            features['Protein_Name'] = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')
            features['Organism'] = data.get('organism', {}).get('scientificName', 'Unknown')
            ecs = [db.get('id') for db in data.get('uniProtKBCrossReferences', []) if db.get('database') == 'EC']
            features['EC_Number'] = ";".join(ecs)
            features['Match_Status'] = "Success"
            features['CAZy_family'] = args.cazy_family
            successful_entries.append(features)

            new_header = f">UniProtID|{uid}|{features['Organism']}|{features['Protein_Name']}"
            raw_fasta_entries.append((new_header, seq))

    # --- OUTPUT ---
    if successful_entries:
        df = pd.DataFrame(successful_entries)
        # Reorder columns
        desired_cols = ['UniProt_ID', 'CAZy_family', 'InterPro_IDs', 'Protein_Name', 'Organism', 'EC_Number',
                        'Signal_End', 'Transmembrane_Regions', 'LPMO_Core_Start', 'LPMO_Core_End', 'Domain_Provenance', 'Binding_Modules',
                        'H1_Verified', 'H1_AminoAcid', 'Match_Status', 'Sequence_Group', 'CoAccessions']
        
        # Ensure all cols exist
        for c in desired_cols:
            if c not in df.columns:
                df[c] = None
        
        df = df[desired_cols]

        # Force integer types for coordinate columns (handling NaNs via Int64)
        for int_col in ['Signal_End', 'LPMO_Core_Start', 'LPMO_Core_End']:
            if int_col in df.columns:
                df[int_col] = df[int_col].astype('Int64')

        df.to_csv(out_tsv, sep='\t', index=False)
        logger.info(f"Wrote metadata to {out_tsv}")
    
    if raw_fasta_entries:
        with open(out_fasta, 'w') as f:
            for header, seq in raw_fasta_entries:
                f.write(f"{header}\n{seq}\n")
        logger.info(f"Wrote sequences to {out_fasta}")
        
    if failed_ids:
        with open(out_failed, 'w') as f:
            for fail in failed_ids:
                f.write(f"{fail}\n")
        logger.info(f"Wrote failed IDs to {out_failed}")
        
    # Write JSON Run Stats
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    
    # Calculate processed count based on mode
    total_processed = 0
    if args.mode == 'fasta':
        total_processed = len(entries) + len(unknown_headers)
    elif args.mode == 'cazy':
        # Safely calculate total from available variables
        c_ncbi = len(ncbi_ids) if 'ncbi_ids' in locals() else 0
        c_jgi = sum(len(v) for v in jgi_groups.values()) if 'jgi_groups' in locals() else 0
        total_processed = c_ncbi + c_jgi
    elif args.mode == 'characterized':
        # For characterized mode, total_processed should be the number of ROWS in the CSV, not unique IDs
        total_processed = len(rows_meta) if 'rows_meta' in locals() else 0
    elif args.mode == 'list':
         total_processed = len(ids) if 'ids' in locals() else 0

    # Derive success/failed counts more accurately for characterized mode
    if args.mode == 'characterized':
        success_ids_count = len(success_id_set) if 'success_id_set' in locals() else len(successful_entries)
        rows_no_ids_count = len(rows_no_ids) if 'rows_no_ids' in locals() else 0
        failed_ids_count = len(failed_id_set) if 'failed_id_set' in locals() else len(failed_ids)
        scount = success_ids_count
        fcount = failed_ids_count + rows_no_ids_count
        
        # Calculate FASTA-related statistics
        fasta_sequences_count = len(seq_groups) if 'seq_groups' in locals() else 0
        row_splits_count = len(row_splits) if 'row_splits' in locals() else 0
        # Count unique CSV rows that produced at least one sequence
        rows_with_sequences = set()
        if 'rows_meta' in locals() and 'id_to_sequence' in locals():
            for row in rows_meta:
                for uid in row.get('uniprot_ids', []):
                    if uid in id_to_sequence:
                        rows_with_sequences.add(row['line_no'])
                        break
        rows_with_sequences_count = len(rows_with_sequences)
        additional_sequences_from_splits = fasta_sequences_count - rows_with_sequences_count
    else:
        scount = len(successful_entries)
        fcount = len(failed_ids)

    run_stats = {
        "run_id": timestamp,
        "input_file": args.input,
        "mode": args.mode,
        "total_processed": total_processed,
        "success_count": scount,
        "failed_count": fcount,
        "duration_seconds": duration,
        "output_files": {
            "metadata": out_tsv,
            "sequences": out_fasta,
            "failed_ids": out_failed
        }
    }
    
    # Add FASTA and split tracking info for characterized mode
    if args.mode == 'characterized':
        run_stats["fasta_sequences_count"] = fasta_sequences_count
        run_stats["row_splits_count"] = row_splits_count
        run_stats["additional_sequences_from_splits"] = additional_sequences_from_splits
        run_stats["split_rows"] = row_splits
    
    with open(out_run_json, 'w') as f:
        json.dump(run_stats, f, indent=2)
    logger.info(f"Wrote run stats to {out_run_json}")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("\n⚠️  Pipeline interrupted by user (Ctrl+C)")
        logger.warning("   Partial results may have been written to output files.")
        sys.exit(130)  # Standard exit code for SIGINT
    except Exception as e:
        logger.error("\n❌ CRITICAL ERROR: Pipeline crashed unexpectedly")
        logger.error(f"   Error type: {type(e).__name__}")
        logger.error(f"   Error message: {e}")
        logger.error("\n   Traceback:", exc_info=True)
        logger.error("\n   ℹ️  Partial results may have been written to output files.")
        logger.error("   Check 'data/metadata/' and 'data/sequences/' directories.")
        sys.exit(1)
