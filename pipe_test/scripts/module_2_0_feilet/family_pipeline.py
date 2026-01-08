"""
Family Pipeline – Per-family orchestration of enrichment.

Orkestrerer komplette per-family workflow:
  1. Fetch CAZy family table
  2. Enrich each sequence (UniProt, NCBI, PDB, taxonomy)
  3. Write CSV + FASTA outputs
"""

from __future__ import annotations

import logging
import time
from typing import List, Dict, Any, Tuple
from pathlib import Path
from dataclasses import dataclass

# Relative imports for package execution
from .cazy_master import get_cazy_family_table
from .enricher import enrich_sequence_metadata
from .metadata_writer import (
    write_enriched_csv,
    write_sequences_fasta,
    validate_rows,
)


# ============================================================================
# LOGGING
# ============================================================================

logger = logging.getLogger(__name__)


# ============================================================================
# STATS DATACLASS
# ============================================================================

@dataclass
class FamilyStats:
    family: str
    n_cazy_entries: int
    n_enriched: int
    n_errors: int
    n_written_csv: int
    n_written_fasta: int
    duration_sec: float
    errors: List[str]
    
    def summary_line(self) -> str:
        """One-line summary for logging."""
        return (
            f"{self.family}: "
            f"{self.n_cazy_entries} CAZy → {self.n_enriched} enriched "
            f"({self.n_errors} errors) → "
            f"{self.n_written_csv} CSV, {self.n_written_fasta} FASTA "
            f"({self.duration_sec:.1f}s)"
        )


# ============================================================================
# MAIN ORCHESTRATOR
# ============================================================================

def process_family(
    family: str,
    run_id: str,
    output_dir: Path,
    ncbi_email: str = "user@example.com",
) -> "FamilyStats":
    """
    Process single CAZy family through full enrichment pipeline.
    
    Args:
        family: CAZy family (e.g., "AA13")
        run_id: Run identifier for output filenames
        output_dir: Base output directory
        ncbi_email: Email for NCBI Entrez queries
        
    Returns:
        FamilyStats with detailed processing results
    """
    t0 = time.time()
    stats = FamilyStats(
        family=family,
        n_cazy_entries=0,
        n_enriched=0,
        n_errors=0,
        n_written_csv=0,
        n_written_fasta=0,
        duration_sec=0.0,
        errors=[],
    )
    
    logger.info(f"Processing family {family}...")
    
    # ========================================================================
    # STEP 1: Fetch CAZy family table
    # ========================================================================
    
    try:
        logger.debug(f"  Fetching CAZy {family} table...")
        cazy_rows = get_cazy_family_table(family)
        stats.n_cazy_entries = len(cazy_rows)
        
        if not cazy_rows:
            logger.warning(f"  No CAZy entries found for {family}")
            stats.duration_sec = time.time() - t0
            return stats
        
        logger.debug(f"  Found {stats.n_cazy_entries} CAZy entries")
        
    except Exception as e:
        err_msg = f"CAZy fetch failed: {e}"
        logger.error(f"  {err_msg}")
        stats.errors.append(err_msg)
        stats.duration_sec = time.time() - t0
        return stats
    
    # ========================================================================
    # STEP 2: Enrich each sequence
    # ========================================================================
    
    enriched_rows = []
    
    for idx, cazy_row in enumerate(cazy_rows, start=1):
        try:
            # Extract identifiers from CAZy row
            uniprot_acc = cazy_row.get("uniprot_acc")
            genbank_acc = cazy_row.get("genbank_acc")
            taxonomy_id = cazy_row.get("taxonomy_id")
            organism_name = cazy_row.get("organism_name")
            pdb_acc = cazy_row.get("pdb_acc")
            
            # Skip if no accession
            if not uniprot_acc and not genbank_acc:
                logger.debug(f"    Skipping row (no accession): {cazy_row}")
                continue
            
            logger.debug(f"  [{idx}/{stats.n_cazy_entries}] Enriching {uniprot_acc or genbank_acc}...")
            
            # Enrich this sequence
            enriched = enrich_sequence_metadata(
                uniprot_acc=uniprot_acc,
                genbank_acc=genbank_acc,
                taxonomy_id=taxonomy_id,
                organism_name=organism_name,
                pdb_acc_from_cazy=pdb_acc,
                ncbi_email=ncbi_email,
            )
            
            # Add CAZy-specific fields
            enriched["family"] = family
            enriched["source_db"] = cazy_row.get("source_db", "CAZy")
            enriched["seq_id"] = uniprot_acc or genbank_acc  # Fallback to GenBank if no UniProt
            
            enriched_rows.append(enriched)
            stats.n_enriched += 1
            
        except Exception as e:
            err_msg = f"Enrichment error (row {idx}): {e}"
            logger.error(f"    {err_msg}")
            stats.errors.append(err_msg)
            stats.n_errors += 1
    
    if not enriched_rows:
        logger.warning(f"  No sequences successfully enriched for {family}")
        stats.duration_sec = time.time() - t0
        return stats
    
    # ========================================================================
    # STEP 3: Validate
    # ========================================================================
    
    validation_warnings = validate_rows(enriched_rows)
    if validation_warnings:
        logger.debug(f"  Validation warnings ({len(validation_warnings)}):")
        for w in validation_warnings[:5]:  # Show first 5
            logger.debug(f"    {w}")
        if len(validation_warnings) > 5:
            logger.debug(f"    ... and {len(validation_warnings) - 5} more")
    
    # ========================================================================
    # STEP 4: Write outputs
    # ========================================================================
    
    family_output_dir = output_dir / family
    
    try:
        # CSV
        csv_path = family_output_dir / f"lpmo_{family}_{run_id}_raw.csv"
        stats.n_written_csv = write_enriched_csv(enriched_rows, csv_path)
        logger.debug(f"  Wrote {stats.n_written_csv} rows to {csv_path.name}")
        
        # FASTA
        fasta_path = family_output_dir / f"lpmo_{family}_{run_id}_raw.fasta"
        stats.n_written_fasta = write_sequences_fasta(enriched_rows, fasta_path)
        logger.debug(f"  Wrote {stats.n_written_fasta} sequences to {fasta_path.name}")
        
    except Exception as e:
        err_msg = f"Output write failed: {e}"
        logger.error(f"  {err_msg}")
        stats.errors.append(err_msg)
    
    # ========================================================================
    # DONE
    # ========================================================================
    
    stats.duration_sec = time.time() - t0
    logger.info(f"  ✓ {stats.summary_line()}")
    
    return stats


# ============================================================================
# BATCH PROCESSING
# ============================================================================

def process_families(
    families: List[str],
    run_id: str,
    output_dir: Path,
    ncbi_email: str = "user@example.com",
    log_level: str = "INFO",
) -> Tuple[List[FamilyStats], int, int]:
    """
    Process multiple families.
    
    Args:
        families: List of CAZy families (e.g., ["AA13", "GH14"])
        run_id: Run identifier
        output_dir: Base output directory
        ncbi_email: Email for NCBI queries
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
        
    Returns:
        Tuple of (stats_list, total_errors, total_duration)
    """
    
    # Setup logging
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    
    logger.info(f"Processing {len(families)} families (run_id={run_id})...")
    logger.info(f"Output directory: {output_dir}")
    
    stats_list = []
    total_errors = 0
    t_batch_start = time.time()
    
    for family in families:
        stats = process_family(
            family=family,
            run_id=run_id,
            output_dir=output_dir,
            ncbi_email=ncbi_email,
        )
        stats_list.append(stats)
        total_errors += stats.n_errors
    
    batch_duration = time.time() - t_batch_start
    
    # Print summary
    logger.info("")
    logger.info("=" * 70)
    logger.info(f"BATCH SUMMARY ({len(families)} families)")
    logger.info("=" * 70)
    
    total_cazy = sum(s.n_cazy_entries for s in stats_list)
    total_enriched = sum(s.n_enriched for s in stats_list)
    total_csv = sum(s.n_written_csv for s in stats_list)
    total_fasta = sum(s.n_written_fasta for s in stats_list)
    
    logger.info(f"Total CAZy entries: {total_cazy}")
    logger.info(f"Total enriched: {total_enriched}")
    logger.info(f"Total written (CSV): {total_csv}")
    logger.info(f"Total written (FASTA): {total_fasta}")
    logger.info(f"Total errors: {total_errors}")
    logger.info(f"Duration: {batch_duration:.1f}s")
    logger.info("=" * 70)
    
    for stats in stats_list:
        status = "✓" if not stats.errors else "⚠"
        logger.info(f"{status} {stats.summary_line()}")
    
    return stats_list, total_errors, batch_duration


if __name__ == "__main__":
    # Example usage
    import tempfile
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir) / "outputs"
        
        stats_list, errors, duration = process_families(
            families=["AA13"],
            run_id="test_20251212_150630",
            output_dir=output_dir,
            ncbi_email="user@example.com",
            log_level="DEBUG",
        )
        
        print(f"\nTest completed with {errors} errors in {duration:.1f}s")
