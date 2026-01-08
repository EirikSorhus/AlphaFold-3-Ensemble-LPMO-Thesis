"""
Module 2 Refactored – Main entry point for enrichment pipeline.

Orkestrerer:
  1. Iterate CAZy families
  2. For each: fetch CAZy table, enrich, write per-family outputs
  3. Merge all per-family CSVs into single metadata file
  4. Merge all per-family FASTAs into single sequences file
  5. Write run metadata JSON
"""

from __future__ import annotations

import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Optional

# Relative import for package execution
from .family_pipeline import process_families


# ============================================================================
# CONFIGURATION
# ============================================================================

DEFAULT_FAMILIES = [
    "AA1",
    "AA3",
    "AA6",
    "AA9",
    "AA13",
    "AA14",
    "AA15",
    "AA16",
    "AA17",
    "AA18",
    "AA19",
    "AA20",
    "GH61",
]

DEFAULT_OUTPUT_DIR = Path(__file__).parent.parent / "data"
DEFAULT_NCBI_EMAIL = "eirik.melbye@nmbu.no"


# ============================================================================
# LOGGING
# ============================================================================

def setup_logging(
    output_dir: Path,
    run_id: str,
    level: str = "INFO",
) -> None:
    """Configure logging to both console and file."""
    
    log_file = output_dir / "metadata" / f"m2_run_{run_id}.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout),
        ]
    )
    
    logger = logging.getLogger(__name__)
    logger.info(f"Logging to {log_file}")


# ============================================================================
# MERGE OUTPUTS
# ============================================================================

def merge_csv_files(
    csv_files: List[Path],
    output_path: Path,
) -> int:
    """
    Merge multiple CSV files into one.
    
    Args:
        csv_files: List of CSV paths to merge
        output_path: Output merged CSV path
        
    Returns:
        Total rows written
    """
    
    import csv
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Read all CSVs
    all_rows = []
    fieldnames = None
    
    for csv_path in csv_files:
        if not csv_path.exists():
            continue
        
        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            if fieldnames is None:
                fieldnames = reader.fieldnames
            
            all_rows.extend(reader)
    
    # Write merged
    if not fieldnames:
        return 0
    
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)
    
    return len(all_rows)


def merge_fasta_files(
    fasta_files: List[Path],
    output_path: Path,
) -> int:
    """
    Merge multiple FASTA files into one.
    
    Args:
        fasta_files: List of FASTA paths to merge
        output_path: Output merged FASTA path
        
    Returns:
        Total sequences written
    """
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    n_seqs = 0
    
    with open(output_path, "w", encoding="utf-8") as out_f:
        for fasta_path in fasta_files:
            if not fasta_path.exists():
                continue
            
            with open(fasta_path, "r", encoding="utf-8") as in_f:
                for line in in_f:
                    out_f.write(line)
                    if line.startswith(">"):
                        n_seqs += 1
    
    return n_seqs


# ============================================================================
# WRITE METADATA
# ============================================================================

def write_run_metadata(
    output_path: Path,
    run_id: str,
    families: List[str],
    n_sequences: int,
    duration_sec: float,
    errors: List[str],
) -> None:
    """Write run metadata to JSON file."""
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    metadata = {
        "run_id": run_id,
        "timestamp": datetime.now().isoformat(),
        "module": "Module 2 (Refactored)",
        "families": families,
        "n_families": len(families),
        "n_sequences": n_sequences,
        "duration_sec": round(duration_sec, 1),
        "total_errors": len(errors),
        "errors": errors,
    }
    
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)


# ============================================================================
# MAIN ORCHESTRATOR
# ============================================================================

def main(
    families: Optional[List[str]] = None,
    output_dir: Optional[Path] = None,
    ncbi_email: Optional[str] = None,
    run_id: Optional[str] = None,
    log_level: str = "INFO",
) -> int:
    """
    Main entry point for Module 2 refactored pipeline.
    
    Args:
        families: List of CAZy families to process (default: DEFAULT_FAMILIES)
        output_dir: Base output directory (default: DEFAULT_OUTPUT_DIR)
        ncbi_email: Email for NCBI queries (default: DEFAULT_NCBI_EMAIL)
        run_id: Run identifier (default: auto-generated YYYYMMDD_HHMMSS)
        log_level: Logging level
        
    Returns:
        0 if success, 1 if any critical errors
    """
    
    # Defaults
    families = families or DEFAULT_FAMILIES
    output_dir = output_dir or DEFAULT_OUTPUT_DIR
    ncbi_email = ncbi_email or DEFAULT_NCBI_EMAIL
    run_id = run_id or datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Setup logging
    setup_logging(output_dir, run_id, log_level)
    logger = logging.getLogger(__name__)
    
    logger.info("")
    logger.info("=" * 70)
    logger.info("MODULE 2 REFACTORED – LPMO ENRICHMENT PIPELINE")
    logger.info("=" * 70)
    logger.info(f"Run ID: {run_id}")
    logger.info(f"Families: {families}")
    logger.info(f"Output directory: {output_dir}")
    logger.info("")
    
    # ========================================================================
    # PHASE 1: Process all families
    # ========================================================================
    
    stats_list, total_errors, batch_duration = process_families(
        families=families,
        run_id=run_id,
        output_dir=output_dir,
        ncbi_email=ncbi_email,
        log_level=log_level,
    )
    
    # ========================================================================
    # PHASE 2: Merge per-family outputs
    # ========================================================================
    
    logger.info("")
    logger.info("=" * 70)
    logger.info("MERGING PER-FAMILY OUTPUTS")
    logger.info("=" * 70)
    
    # Collect all per-family files
    csv_files = []
    fasta_files = []
    
    for family in families:
        family_dir = output_dir / family
        csv_path = family_dir / f"lpmo_{family}_{run_id}_raw.csv"
        fasta_path = family_dir / f"lpmo_{family}_{run_id}_raw.fasta"
        
        if csv_path.exists():
            csv_files.append(csv_path)
        if fasta_path.exists():
            fasta_files.append(fasta_path)
    
    # Merge CSVs
    logger.info(f"Merging {len(csv_files)} CSV files...")
    merged_csv_path = output_dir / "metadata" / f"m2_sequence_3d_metadata_{run_id}.csv"
    n_csv_rows = merge_csv_files(csv_files, merged_csv_path)
    logger.info(f"  Wrote {n_csv_rows} sequences to {merged_csv_path.name}")
    
    # Merge FASTAs
    logger.info(f"Merging {len(fasta_files)} FASTA files...")
    merged_fasta_path = output_dir / "sequences" / f"lpmo_all_{run_id}_raw.fasta"
    n_sequences = merge_fasta_files(fasta_files, merged_fasta_path)
    logger.info(f"  Wrote {n_sequences} sequences to {merged_fasta_path.name}")
    
    # ========================================================================
    # PHASE 3: Write metadata
    # ========================================================================
    
    logger.info("")
    logger.info("=" * 70)
    logger.info("WRITING RUN METADATA")
    logger.info("=" * 70)
    
    all_errors = []
    for stats in stats_list:
        all_errors.extend(stats.errors)
    
    metadata_path = output_dir / "metadata" / f"m2_run_metadata_{run_id}.json"
    write_run_metadata(
        output_path=metadata_path,
        run_id=run_id,
        families=families,
        n_sequences=n_sequences,
        duration_sec=batch_duration,
        errors=all_errors,
    )
    logger.info(f"Wrote metadata to {metadata_path.name}")
    
    # ========================================================================
    # FINAL SUMMARY
    # ========================================================================
    
    logger.info("")
    logger.info("=" * 70)
    logger.info("FINAL SUMMARY")
    logger.info("=" * 70)
    logger.info(f"Total duration: {batch_duration:.1f}s")
    logger.info(f"Total sequences: {n_sequences}")
    logger.info(f"Total errors: {len(all_errors)}")
    
    if all_errors:
        logger.warning("Errors encountered:")
        for err in all_errors[:10]:
            logger.warning(f"  - {err}")
        if len(all_errors) > 10:
            logger.warning(f"  ... and {len(all_errors) - 10} more")
    
    logger.info("")
    logger.info("Output files:")
    logger.info(f"  CSV: {merged_csv_path}")
    logger.info(f"  FASTA: {merged_fasta_path}")
    logger.info(f"  Metadata: {metadata_path}")
    logger.info("")
    logger.info("✓ Pipeline completed")
    logger.info("=" * 70)
    
    return 1 if total_errors > 0 else 0


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Module 2 Refactored – LPMO enrichment pipeline",
    )
    parser.add_argument(
        "--families",
        nargs="+",
        default=DEFAULT_FAMILIES,
        help=f"CAZy families to process (default: {DEFAULT_FAMILIES})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--ncbi-email",
        default=DEFAULT_NCBI_EMAIL,
        help=f"NCBI email (default: {DEFAULT_NCBI_EMAIL})",
    )
    parser.add_argument(
        "--run-id",
        help="Run ID (default: auto-generated YYYYMMDD_HHMMSS)",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level",
    )
    
    args = parser.parse_args()
    
    exit_code = main(
        families=args.families,
        output_dir=args.output_dir,
        ncbi_email=args.ncbi_email,
        run_id=args.run_id,
        log_level=args.log_level,
    )
    
    sys.exit(exit_code)
