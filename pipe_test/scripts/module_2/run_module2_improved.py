# scripts/module_2/run_module2_improved.py
"""
Main pipeline orchestrator for Module 2: LPMO sequence fetching.

This script runs the complete pipeline:
1. Fetch sequences from UniProt
2. Fetch GenBank IDs from CAZy and sequences from NCBI
3. Merge and deduplicate sequences

Usage:
    python -m scripts.module_2.run_module2_improved [--email your.email@domain.com]
    
    Or set environment variable:
    export NCBI_EMAIL="your.email@domain.com"
    python -m scripts.module_2.run_module2_improved
"""

import sys
import argparse
from pathlib import Path

from scripts.module_2.config_improved import (
    check_dependencies,
    validate_config,
    ensure_data_directory,
    DEFAULT_CAZY_FAMILIES,
    DATA_RAW_DIR,
    DATA_SEQUENCES_DIR,
    DATA_METADATA_DIR,
    get_run_id,
    build_uniprot_query,
    save_metadata,
    get_merged_output_filename,
    get_metadata_filename,
    get_family_output_filename
)
from scripts.module_2.uniprot_fetch_improved import fetch_uniprot_by_families
from scripts.module_2.ncbi_fetch_improved import fetch_ncbi_by_families
from scripts.module_2.merge_sequences_improved import merge_family_files


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Module 2: Fetch LPMO sequences from UniProt, CAZy, and NCBI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Using environment variable NCBI_EMAIL with default families
  export NCBI_EMAIL="your.email@institution.no"
  python -m scripts.module_2.run_module2_improved
  
  # Specify families and email
  python -m scripts.module_2.run_module2_improved --families AA9 AA10 --email your.email@institution.no
  
  # Only fetch AA13 family
  python -m scripts.module_2.run_module2_improved --families AA13 --email your@email.com
  
  # Skip certain steps
  python -m scripts.module_2.run_module2_improved --skip-uniprot --families AA9 AA10
        """
    )
    
    parser.add_argument(
        "--families",
        type=str,
        nargs="+",
        default=None,
        help=f"CAZy families to fetch (default: {' '.join(DEFAULT_CAZY_FAMILIES)}). "
             "Example: --families AA9 AA10 AA11"
    )
    
    parser.add_argument(
        "--email",
        type=str,
        help="Email address for NCBI API (required by NCBI)"
    )
    
    parser.add_argument(
        "--skip-uniprot",
        action="store_true",
        help="Skip UniProt fetching step"
    )
    
    parser.add_argument(
        "--skip-cazy",
        action="store_true",
        help="Skip CAZy/NCBI fetching step"
    )
    
    parser.add_argument(
        "--skip-merge",
        action="store_true",
        help="Skip merging step"
    )
    
    parser.add_argument(
        "--deduplicate-by",
        type=str,
        choices=["sequence", "id", "both", "none"],
        default="sequence",
        help="Deduplication strategy for merge step (default: sequence)"
    )
    
    return parser.parse_args()


def run_module2_pipeline(
    families=None,
    ncbi_email=None,
    skip_uniprot=False,
    skip_cazy=False,
    skip_merge=False,
    deduplicate_by="sequence"
) -> int:
    """
    Run the complete Module 2 pipeline with per-family output files.
    
    Args:
        families: List of CAZy families to fetch. If None, uses DEFAULT_CAZY_FAMILIES
        ncbi_email: Email for NCBI API (required if skip_cazy=False)
        skip_uniprot: Skip UniProt fetching
        skip_cazy: Skip CAZy/NCBI fetching
        skip_merge: Skip merge step
        deduplicate_by: Deduplication strategy for merging
        
    Returns:
        int: Exit code (0 for success, 1 for failure)
    """
    # Use default families if not specified
    if families is None:
        families = DEFAULT_CAZY_FAMILIES
    
    # Generate run ID for this pipeline execution
    run_id = get_run_id()
    
    print("=" * 70)
    print("MODULE 2: LPMO SEQUENCE FETCHING PIPELINE")
    print("=" * 70)
    print(f"Run ID: {run_id}")
    print(f"Families: {', '.join(families)}")
    print(f"Output directories:")
    print(f"  - Sequences: {DATA_SEQUENCES_DIR}")
    print(f"  - Metadata: {DATA_METADATA_DIR}")
    print()
    
    # Check dependencies
    print("Pre-flight checks:")
    print("-" * 70)
    
    if not check_dependencies():
        print("\n✗ Missing dependencies. Please install required packages.")
        return 1
    
    if not validate_config(families):
        print("\n✗ Configuration validation failed.")
        return 1
    
    ensure_data_directory()
    print()
    
    # Build UniProt query for metadata
    uniprot_query = build_uniprot_query(families)
    
    # Track overall success and results
    all_steps_successful = True
    metadata_info = {
        "uniprot_results": {},
        "cazy_results": {},
        "merge_total": 0,
        "merge_unique": 0
    }
    
    # ========================================================================
    # STEP 1: UniProt
    # ========================================================================
    if not skip_uniprot:
        print("STEP 1: UniProt Sequence Fetching (per family)")
        print("-" * 70)
        
        try:
            uniprot_results = fetch_uniprot_by_families(families, run_id)
            metadata_info["uniprot_results"] = uniprot_results
            
            total_seqs = sum(uniprot_results.values())
            if total_seqs == 0:
                print("[PIPELINE] WARNING: UniProt returned 0 sequences total")
                all_steps_successful = False
            else:
                print(f"[PIPELINE] ✓ UniProt: {total_seqs} sequences across {len(families)} families")
                
        except Exception as e:
            print(f"[PIPELINE] ✗ ERROR in UniProt step: {e}")
            import traceback
            traceback.print_exc()
            all_steps_successful = False
        
        print()
    else:
        print("STEP 1: UniProt Sequence Fetching [SKIPPED]")
        print("-" * 70)
        print()
    
    # ========================================================================
    # STEP 2: CAZy + NCBI
    # ========================================================================
    if not skip_cazy:
        print("STEP 2: CAZy + NCBI Sequence Fetching (per family)")
        print("-" * 70)
        
        try:
            cazy_results = fetch_ncbi_by_families(families, run_id, ncbi_email=ncbi_email)
            metadata_info["cazy_results"] = {
                fam: {"success": succ, "failed": fail}
                for fam, (succ, fail) in cazy_results.items()
            }
            
            total_success = sum(succ for succ, _ in cazy_results.values())
            total_failed = sum(fail for _, fail in cazy_results.values())
            
            if total_success == 0:
                print("[PIPELINE] WARNING: CAZy/NCBI returned 0 sequences")
                all_steps_successful = False
            else:
                print(f"[PIPELINE] ✓ CAZy/NCBI: {total_success} sequences ({total_failed} failed)")
                    
        except ValueError as e:
            print(f"[PIPELINE] ✗ ERROR: {e}")
            print("[PIPELINE] Please set NCBI_EMAIL environment variable or use --email")
            return 1
        except Exception as e:
            print(f"[PIPELINE] ✗ ERROR in CAZy/NCBI step: {e}")
            import traceback
            traceback.print_exc()
            all_steps_successful = False
        
        print()
    else:
        print("STEP 2: CAZy + NCBI Sequence Fetching [SKIPPED]")
        print("-" * 70)
        print()
    
    # ========================================================================
    # STEP 3: Merge
    # ========================================================================
    if not skip_merge:
        print("STEP 3: Merge and Deduplicate Sequences")
        print("-" * 70)
        
        # Determine which sources to include
        sources = []
        if not skip_uniprot:
            sources.append('uniprot')
        if not skip_cazy:
            sources.append('ncbi')
        
        if not sources:
            print("[PIPELINE] WARNING: No sources to merge (all steps skipped)")
            all_steps_successful = False
        else:
            print(f"[PIPELINE] Merging files from sources: {', '.join(sources)}")
            
            try:
                total, unique = merge_family_files(
                    families,
                    run_id,
                    sources=sources,
                    deduplicate_by=deduplicate_by
                )
                
                metadata_info["merge_total"] = total
                metadata_info["merge_unique"] = unique
                
                if unique == 0:
                    print("[PIPELINE] WARNING: Merge resulted in 0 sequences")
                    all_steps_successful = False
                else:
                    print(f"[PIPELINE] ✓ Merge: {total} → {unique} unique sequences")
                    
            except Exception as e:
                print(f"[PIPELINE] ✗ ERROR in merge step: {e}")
                import traceback
                traceback.print_exc()
                all_steps_successful = False
        
        print()
    else:
        print("STEP 3: Merge and Deduplicate Sequences [SKIPPED]")
        print("-" * 70)
        print()
    
    # ========================================================================
    # STEP 4: Save Metadata (MOVED TO AFTER VALIDATION - see below)
    # ========================================================================
    # Metadata will be saved after all pipeline steps are complete
    # to avoid partial/false success records
    print("STEP 4: Save Metadata [DEFERRED - will save after validation]")
    print("-" * 70)
    print()

    # ========================================================================
    # STEP 5: Add family names to FASTA headers and generate sequence metadata
    # ========================================================================
    print("STEP 5: Annotate sequences with family names and 3D structure info")
    print("-" * 70)
    
    try:
        from scripts.module_2.structure_info_uniprot import fetch_uniprot_structure_info
        from scripts.module_2.structure_info import map_ncbi_to_uniprot
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        import json
        import csv
        
        # Process each family file and add family name to headers
        all_sequences = []
        family_seq_count = {}
        
        for family in families:
            family_file = DATA_SEQUENCES_DIR / get_family_output_filename(family, run_id)
            if family_file.exists():
                records = list(SeqIO.parse(str(family_file), "fasta"))
                family_seq_count[family] = len(records)
                
                # Add family name to description
                for rec in records:
                    rec.description = f"{rec.id} {family}"
                    all_sequences.append((rec, family))
                
                # Rewrite file with updated headers
                SeqIO.write(records, family_file, "fasta")
                print(f"[PIPELINE] ✓ Updated {family_file.name} with family names ({len(records)} seqs)")
            else:
                print(f"[PIPELINE] ⚠ File not found: {family_file.name}")
                family_seq_count[family] = 0
        
        # Create merged file with all sequences
        merged_file = DATA_SEQUENCES_DIR / get_merged_output_filename(run_id)
        if all_sequences:
            merged_records = [rec for rec, _ in all_sequences]
            SeqIO.write(merged_records, merged_file, "fasta")
            print(f"[PIPELINE] ✓ Created {merged_file.name} ({len(merged_records)} total sequences)")
        
        # Generate CSV metadata with 3D structure info (include run_id in filename)
        csv_metadata_path = DATA_METADATA_DIR / f"m2_sequence_3d_metadata_{run_id}.csv"
        print(f"[PIPELINE] Fetching 3D structure info for {len(all_sequences)} sequences...")
        
        csv_rows = []
        for rec, family in all_sequences:
            seq_id = rec.id
            
            # Try to get UniProt ID for structure lookup
            uniprot_id = None
            if seq_id.startswith("sp|") or seq_id.startswith("tr|"):
                uniprot_id = seq_id.split("|")[1]
            elif len(seq_id) == 6 and seq_id[0].isalpha():
                uniprot_id = seq_id
            
            # Fetch 3D structure info
            if uniprot_id:
                try:
                    struct_info = fetch_uniprot_structure_info(uniprot_id)
                except Exception as e:
                    print(f"[PIPELINE] WARNING: Could not fetch structure for {uniprot_id}: {e}")
                    struct_info = {
                        "has_experimental_structure": False,
                        "pdb_ids": [],
                        "best_method": None,
                        "best_resolution": None,
                        "best_coverage": None,
                        "has_alphafold_model": False,
                        "alphafold_confidence_summary": None
                    }
            else:
                struct_info = {
                    "has_experimental_structure": False,
                    "pdb_ids": [],
                    "best_method": None,
                    "best_resolution": None,
                    "best_coverage": None,
                    "has_alphafold_model": False,
                    "alphafold_confidence_summary": None
                }
            
            # Build CSV row
            csv_rows.append({
                "protein_id": seq_id,
                "family": family,
                "sequence_length": len(rec.seq),
                "has_experimental_structure": struct_info["has_experimental_structure"],
                "pdb_ids": ",".join(struct_info["pdb_ids"]) if struct_info["pdb_ids"] else "",
                "best_pdb_resolution": struct_info["best_resolution"] if struct_info["best_resolution"] else "",
                "has_alphafold_model": struct_info["has_alphafold_model"],
                "alphafold_accession": struct_info["alphafold_confidence_summary"] if struct_info["alphafold_confidence_summary"] else ""
            })
        
        # Write CSV
        if csv_rows:
            with open(csv_metadata_path, "w", newline="", encoding="utf-8") as f:
                fieldnames = ["protein_id", "family", "sequence_length", "has_experimental_structure", 
                             "pdb_ids", "best_pdb_resolution", "has_alphafold_model", "alphafold_accession"]
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(csv_rows)
            print(f"[PIPELINE] ✓ Sequence metadata saved: {csv_metadata_path.name}")
        
        # Update run metadata with family counts
        metadata_path = DATA_METADATA_DIR / get_metadata_filename(run_id)
        if metadata_path.exists():
            with open(metadata_path, "r", encoding="utf-8") as f:
                meta = json.load(f)
            meta["n_sequences_per_family"] = family_seq_count
            meta["source_databases"] = {
                "uniprot": "REST API",
                "cazy": "Web scraping",
                "ncbi": "Entrez API",
                "structure": ["PDB", "AlphaFold DB"]
            }
            with open(metadata_path, "w", encoding="utf-8") as f:
                json.dump(meta, f, indent=2, ensure_ascii=False)
        
    except Exception as e:
        print(f"[PIPELINE] ✗ ERROR in annotation step: {e}")
        import traceback
        traceback.print_exc()
        all_steps_successful = False
    
    print()
    
    # ========================================================================
    # STEP 6: Validation checks
    # ========================================================================
    print("STEP 6: Validation Checks")
    print("-" * 70)
    
    validation_passed = True
    
    try:
        from Bio import SeqIO
        import json
        import csv
        
        # Check 1: All IDs in metadata CSV exist in lpmo_all_raw.fasta
        merged_file = DATA_SEQUENCES_DIR / get_merged_output_filename(run_id)
        csv_metadata_path = DATA_METADATA_DIR / f"m2_sequence_3d_metadata_{run_id}.csv"
        
        if merged_file.exists() and csv_metadata_path.exists():
            fasta_ids = {rec.id for rec in SeqIO.parse(str(merged_file), "fasta")}
            
            with open(csv_metadata_path, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                csv_ids = {row["protein_id"] for row in reader}
            
            if fasta_ids == csv_ids:
                print(f"[VALIDATION] ✓ All metadata IDs match FASTA ({len(fasta_ids)} sequences)")
            else:
                print(f"[VALIDATION] ✗ ID mismatch: {len(fasta_ids)} in FASTA, {len(csv_ids)} in metadata")
                validation_passed = False
        
        # Check 2: Sum of per-family sequences equals merged total
        total_family_seqs = 0
        for family in families:
            family_file = DATA_SEQUENCES_DIR / get_family_output_filename(family, run_id)
            if family_file.exists():
                family_count = len(list(SeqIO.parse(str(family_file), "fasta")))
                total_family_seqs += family_count
        
        if merged_file.exists():
            merged_count = len(list(SeqIO.parse(str(merged_file), "fasta")))
            if merged_count == total_family_seqs:
                print(f"[VALIDATION] ✓ Merged file count matches sum of family files ({merged_count} sequences)")
            else:
                print(f"[VALIDATION] ✗ Count mismatch: {merged_count} merged, {total_family_seqs} total from families")
                validation_passed = False
        
        # Check 3: n_sequences_per_family in run metadata matches actual FASTA counts
        metadata_path = DATA_METADATA_DIR / get_metadata_filename(run_id)
        if metadata_path.exists():
            with open(metadata_path, "r", encoding="utf-8") as f:
                meta = json.load(f)
            
            if "n_sequences_per_family" in meta:
                metadata_counts = meta["n_sequences_per_family"]
                counts_match = True
                for family in families:
                    family_file = DATA_SEQUENCES_DIR / get_family_output_filename(family, run_id)
                    if family_file.exists():
                        actual_count = len(list(SeqIO.parse(str(family_file), "fasta")))
                        meta_count = metadata_counts.get(family, 0)
                        if actual_count != meta_count:
                            print(f"[VALIDATION] ✗ {family}: {actual_count} sequences in file, {meta_count} in metadata")
                            counts_match = False
                            validation_passed = False
                
                if counts_match:
                    print(f"[VALIDATION] ✓ Run metadata counts match actual FASTA files")
        
        # Check 4: Family names in metadata match FASTA file associations
        if csv_metadata_path.exists():
            with open(csv_metadata_path, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    family = row["family"]
                    if family not in families:
                        print(f"[VALIDATION] ✗ Unexpected family in metadata: {family}")
                        validation_passed = False
                        break
                else:
                    print(f"[VALIDATION] ✓ All families in metadata are expected")
        
        if validation_passed:
            print("\n[VALIDATION] ✓ All validation checks passed!")
        else:
            print("\n[VALIDATION] ⚠ Some validation checks failed")
            all_steps_successful = False
            
    except Exception as e:
        print(f"[VALIDATION] ✗ Validation failed with error: {e}")
        import traceback
        traceback.print_exc()
        all_steps_successful = False
    
    print()
    
    # ========================================================================
    # NOW: Save Metadata (AFTER all steps and validation)
    # ========================================================================
    print("SAVING METADATA (after all pipeline steps)")
    print("-" * 70)
    
    try:
        metadata_file = save_metadata(
            run_id=run_id,
            families=families,
            uniprot_query=uniprot_query,
            output_dir=DATA_METADATA_DIR,
            additional_info=metadata_info
        )
        print(f"[PIPELINE] ✓ Run metadata saved: {metadata_file.name}")
    except Exception as e:
        print(f"[PIPELINE] ⚠ WARNING: Could not save run metadata: {e}")
    print()
    
    # ========================================================================
    # Summary
    # ========================================================================
    print("=" * 70)
    print("PIPELINE SUMMARY")
    print("=" * 70)
    print(f"Run ID: {run_id}")
    print(f"Output directories:")
    print(f"  - Sequences: {DATA_SEQUENCES_DIR}")
    print(f"  - Metadata: {DATA_METADATA_DIR}")
    print()
    
    # List output files
    print("Sequence files (data/sequences/):")
    
    # Family-specific files
    for family in families:
        family_file = DATA_SEQUENCES_DIR / get_family_output_filename(family, run_id)
        if family_file.exists():
            print(f"  ✓ {family_file.name}")
    
    # Merged file
    merged_file = DATA_SEQUENCES_DIR / get_merged_output_filename(run_id)
    if merged_file.exists():
        print(f"  ✓ {merged_file.name} (ALL FAMILIES COMBINED)")
    
    print()
    print("Metadata files (data/metadata/):")
    
    # Run metadata
    run_meta_file = DATA_METADATA_DIR / get_metadata_filename(run_id)
    if run_meta_file.exists():
        print(f"  ✓ {run_meta_file.name} (run metadata)")
    
    # Sequence metadata
    seq_meta_file = DATA_METADATA_DIR / f"m2_sequence_3d_metadata_{run_id}.csv"
    if seq_meta_file.exists():
        print(f"  ✓ {seq_meta_file.name} (per-sequence 3D info)")
    
    print()
    
    if all_steps_successful:
        print("✓ PIPELINE COMPLETED SUCCESSFULLY")
        return 0
    else:
        print("⚠ PIPELINE COMPLETED WITH WARNINGS/ERRORS")
        print("Check the output above for details")
        return 1


def main():
    """Main entry point for command-line execution."""
    args = parse_arguments()
    
    exit_code = run_module2_pipeline(
        families=args.families,
        ncbi_email=args.email,
        skip_uniprot=args.skip_uniprot,
        skip_cazy=args.skip_cazy,
        skip_merge=args.skip_merge,
        deduplicate_by=args.deduplicate_by
    )
    
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
