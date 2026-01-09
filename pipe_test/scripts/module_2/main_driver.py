import argparse
from pathlib import Path
from datetime import datetime

from input_handler import handle_input
from metadata_fetcher import fetch_uniprot_metadata_batch
from fasta_fetcher import fetch_uniprot_fasta_bulk, fetch_ncbi_fasta_bulk
from paths_config import PipelinePaths
from utils_io import write_json, write_tsv


def timestamp():
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def main():
    parser = argparse.ArgumentParser(description="Main driver for protein metadata pipeline")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--families", nargs="*", help="List of CAZy families to include")
    group.add_argument("--uniprot-id-file", help="File with UniProt IDs (one per line)")
    group.add_argument("--input-fasta", help="Input FASTA file for sequence lookup")
    parser.add_argument("--project-dir", default=".", help="Project root folder for outputs")
    parser.add_argument("--contact-email", help="Email for API courtesy header")
    parser.add_argument("--allow-ncbi-fallback", action="store_true", help="Use NCBI for unmatched proteins")
    parser.add_argument("--batch-size", type=int, default=100, help="Batch size for UniProt queries")
    parser.add_argument("--max-seq-search", type=int, default=10, help="Max FASTA sequence lookups")
    args = parser.parse_args()

    PipelinePaths.configure_base(Path(args.project_dir))
    PipelinePaths.make_all_dirs()

    run_id = timestamp()

    print("[*] Parsing input sources...")
    all_ids, ncbi_metadata_dict, source_log, failed_items = handle_input(
        families=args.families,
        id_file=args.uniprot_id_file,
        fasta_path=args.input_fasta,
        allow_ncbi=args.allow_ncbi_fallback,
        contact_email=args.contact_email,
        max_seq_search=args.max_seq_search,
        run_id=run_id,
    )

    print(f"[*] Resolved {len(all_ids)} total input IDs")
    for k, v in source_log.items():
        print(f"   - {k}: {v}")

    uniprot_ids = [uid for uid in all_ids if not uid.startswith("NCBI_")]
    ncbi_ids = [uid.replace("NCBI_", "") for uid in all_ids if uid.startswith("NCBI_")]

    # Metadata
    metadata_rows = fetch_uniprot_metadata_batch(uniprot_ids, contact_email=args.contact_email, batch_size=args.batch_size)
    metadata_rows.extend(ncbi_metadata_dict.values())

    fieldnames = [
        "UniProt_ID",
        "Protein_Name",
        "Gene_Name",
        "Organism",
        "EC_Number",
        "Reviewed",
        "PDB_IDs",
        "Pfam_IDs",
        "InterPro_IDs",
        "Match_Status",
        "Source_DB",
    ]

    for row in metadata_rows:
        for f in fieldnames:
            row.setdefault(f, "")

    meta_path = PipelinePaths.metadata_dir / f"metadata_expanded_{run_id}.tsv"
    write_tsv(metadata_rows, fieldnames, meta_path)

    # FASTA
    print("[*] Fetching FASTA sequences...")
    fasta_path = PipelinePaths.sequences_dir / f"all_sequences_{run_id}.fasta"
    fasta_stats = {}
    if uniprot_ids:
        fasta_stats["uniprot"] = fetch_uniprot_fasta_bulk(uniprot_ids, fasta_path, contact_email=args.contact_email)
    if ncbi_ids:
        mode = "a" if uniprot_ids else "w"
        fasta_stats["ncbi"] = fetch_ncbi_fasta_bulk(ncbi_ids, fasta_path, contact_email=args.contact_email, mode=mode)

    # Collect failed ID info
    failed_rows = list(failed_items)  # [(ID, Organism, Reason)] from input handling
    # Add FASTA failures with generic reason
    for fid in fasta_stats.get("uniprot", {}).get("failed_ids", []):
        failed_rows.append((fid, "", "fasta_fetch_failed"))
    for fid in fasta_stats.get("ncbi", {}).get("failed_ids", []):
        failed_rows.append((f"NCBI_{fid}", "", "fasta_fetch_failed"))

    failed_path = None
    if failed_rows:
        failed_path = PipelinePaths.run_dir / f"failed_ids_{run_id}.txt"
        with failed_path.open("w", encoding="utf-8") as f:
            f.write("ID\tOrganism\tReason\n")
            for idv, org, reason in failed_rows:
                f.write(f"{idv}\t{org}\t{reason}\n")

    run_json = {
        "run_id": run_id,
        "timestamp": run_id,
        "input_counts": {
            "total_ids": len(all_ids),
            "uniprot": len(uniprot_ids),
            "ncbi": len(ncbi_ids),
        },
        "source_breakdown": source_log,
        "output_files": {
            "metadata_tsv": str(meta_path),
            "fasta_all": str(fasta_path),
            "failed_ids_txt": str(failed_path) if failed_path else None,
        },
        "fasta_stats": fasta_stats,
    }

    run_path = PipelinePaths.run_dir / f"run_metadata_{run_id}.json"
    write_json(run_json, run_path)

    print("[✓] Done. Output written:")
    print(f"   - Metadata: {meta_path}")
    print(f"   - FASTA (combined): {fasta_path}")
    if failed_path:
        print(f"   - Failed IDs: {failed_path}")
    print(f"   - Run JSON: {run_path}")


if __name__ == "__main__":
    main()
