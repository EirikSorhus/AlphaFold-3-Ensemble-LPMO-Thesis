import argparse
import logging
import os
import sys
import pandas as pd
import json
from datetime import datetime
from tqdm import tqdm

# Import local modules - robust check for package vs script execution
try:
    from .input_handler import parse_fasta_file, parse_fasta_header, format_jgi_query
    from .uniprot_client import UniProtClient
    from .feature_parser import parse_uniprot_features
    from .blast_client import run_blast_search
except ImportError:
    from input_handler import parse_fasta_file, parse_fasta_header, format_jgi_query
    from uniprot_client import UniProtClient
    from feature_parser import parse_uniprot_features
    from blast_client import run_blast_search

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(description="Module 2: Metadata Enrichment")
    parser.add_argument("--mode", choices=['fasta', 'cazy', 'list'], required=True, help="Input mode")
    parser.add_argument("--input", required=True, help="Input file path")
    parser.add_argument("--output_dir", default="data", help="Base output directory (default: data)")
    parser.add_argument("--allow-ncbi-fallback", action="store_true", help="Enable NCBI fallback (not currently active)")
    parser.add_argument("--allow-sequence-search", action="store_true", help="Allow BLAST sequence search for unknown FASTA headers")
    parser.add_argument("--max-sequence-searches", type=int, default=20, help="Maximum number of sequences to search via BLAST (default: 20)")
    
    args = parser.parse_args()
    
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

    client = UniProtClient()
    
    successful_entries = []
    failed_ids = []
    raw_fasta_entries = []
    start_time = datetime.now()
    
    logger.info(f"Starting Module 2 in {args.mode} mode. Output will be timestamped: {timestamp}")
    
    # --- MODE: FASTA ---
    if args.mode == 'fasta':
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
                features = parse_uniprot_features(data)
                
                # Add extra fields standardisation
                features['Protein_Name'] = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')
                features['Organism'] = data.get('organism', {}).get('scientificName', 'Unknown')
                ecs = [db.get('id') for db in data.get('uniProtKBCrossReferences', []) if db.get('database') == 'EC']
                features['EC_Number'] = ";".join(ecs)
                features['Match_Status'] = "Success"
                
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
                            features = parse_uniprot_features(data)
                            uid = features['UniProt_ID']
                            
                            features['Protein_Name'] = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')
                            features['Organism'] = data.get('organism', {}).get('scientificName', 'Unknown')
                            ecs = [db.get('id') for db in data.get('uniProtKBCrossReferences', []) if db.get('database') == 'EC']
                            features['EC_Number'] = ";".join(ecs)
                            features['Match_Status'] = "Success_SeqMatch"
                            
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
             features = parse_uniprot_features(data)
             uid = features['UniProt_ID']
             features['Protein_Name'] = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')
             features['Organism'] = data.get('organism', {}).get('scientificName', 'Unknown')
             ecs = [db.get('id') for db in data.get('uniProtKBCrossReferences', []) if db.get('database') == 'EC']
             features['EC_Number'] = ";".join(ecs)
             features['Match_Status'] = "Success"
             successful_entries.append(features)
             
             seq = data.get('sequence', {}).get('value', '')
             new_header = f">UniProtID|{uid}|{features['Organism']}|{features['Protein_Name']}"
             raw_fasta_entries.append((new_header, seq))
             
    # --- OUTPUT ---
    if successful_entries:
        df = pd.DataFrame(successful_entries)
        # Reorder columns
        desired_cols = ['UniProt_ID', 'InterPro_IDs', 'Protein_Name', 'Organism', 'EC_Number', 
                        'Signal_End', 'LPMO_Core_Start', 'LPMO_Core_End', 'Binding_Modules', 
                        'H1_Verified', 'H1_AminoAcid', 'Match_Status']
        
        # Ensure all cols exist
        for c in desired_cols:
            if c not in df.columns:
                df[c] = None
        
        df = df[desired_cols]
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
    run_stats = {
        "run_id": timestamp,
        "input_file": args.input,
        "mode": args.mode,
        "total_processed": len(entries) + len(unknown_headers) if args.mode == 'fasta' else len(ids),
        "success_count": len(successful_entries),
        "failed_count": len(failed_ids),
        "duration_seconds": duration,
        "output_files": {
            "metadata": out_tsv,
            "sequences": out_fasta,
            "failed_ids": out_failed
        }
    }
    with open(out_run_json, 'w') as f:
        json.dump(run_stats, f, indent=2)
    logger.info(f"Wrote run stats to {out_run_json}")

if __name__ == "__main__":
    main()
