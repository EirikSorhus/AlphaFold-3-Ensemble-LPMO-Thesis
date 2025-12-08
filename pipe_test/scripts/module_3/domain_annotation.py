# scripts/module_3/domain_annotation.py
"""
Domain annotation pipeline utilities for Module 3.

Responsibilities:
- Run hmmscan against configured HMM databases
- Parse domtblout files and flag weak domain hits
- Select catalytic LPMO domains and CBMs per sequence
- Generate sequence variants and metadata tables
"""

from __future__ import annotations

import csv
import json
import re
import shutil
import subprocess
from collections import defaultdict
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from scripts.module_3.config_m3 import (
    Module3Config,
    M2_FASTA,
    M2_METADATA,
    DOMAINS_PARSED,
    DOMAIN_METADATA,
    RUN_METADATA,
    FULL_LENGTH_FASTA,
    CATALYTIC_FASTA,
    RAW_DOMTBL_DIR,
    ensure_output_directories,
    config_to_dict,
)


DomainHit = Dict[str, object]


# -----------------------------------------------------------------------------
# HMMER helpers
# -----------------------------------------------------------------------------

def hmmer_version(binary: str) -> Optional[str]:
    """Return hmmscan version string if available."""
    if not shutil.which(binary):
        return None
    try:
        result = subprocess.run(
            [binary, "-h"],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        first_line = result.stdout.splitlines()[0] if result.stdout else ""
        return first_line.strip() or None
    except Exception:
        return None


def run_hmmscan(binary: str, hmm_db: Path, fasta: Path, domtblout: Path, cpu: int = 2) -> bool:
    """
    Execute hmmscan and write domtblout.

    Returns:
        bool: True if the command completed successfully, False otherwise.
    """
    if not shutil.which(binary):
        print(f"[HMMER] ERROR: hmmscan binary not found: {binary}")
        return False
    if not hmm_db.exists():
        print(f"[HMMER] ERROR: HMM database missing: {hmm_db}")
        return False
    cmd = [
        binary,
        "--domtblout",
        str(domtblout),
        "--cpu",
        str(cpu),
        str(hmm_db),
        str(fasta),
    ]
    print(f"[HMMER] Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    except Exception as exc:  # pragma: no cover - defensive
        print(f"[HMMER] ERROR: Failed to run hmmscan: {exc}")
        return False

    if result.returncode != 0:
        print(f"[HMMER] ERROR: hmmscan exited with code {result.returncode}")
        if result.stderr:
            print(result.stderr.strip())
        if result.stdout:
            print(result.stdout.strip())
        return False

    print(f"[HMMER] domtblout written to {domtblout}")
    return True


def parse_domtblout(domtblout: Path, source: str) -> List[DomainHit]:
    """Parse a domtblout file into a list of domain hits."""
    hits: List[DomainHit] = []
    if not domtblout.exists():
        print(f"[PARSE] WARNING: domtblout file not found: {domtblout}")
        return hits

    with open(domtblout, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 22:
                continue

            target_name = parts[0]
            target_len = int(parts[2])
            query_name = parts[3]
            query_len = int(parts[5])
            full_evalue = float(parts[6])
            full_score = float(parts[7])
            full_bias = float(parts[8])
            c_evalue = float(parts[11])
            i_evalue = float(parts[12])
            domain_score = float(parts[13])
            domain_bias = float(parts[14])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            ali_from = int(parts[17])
            ali_to = int(parts[18])
            env_from = int(parts[19])
            env_to = int(parts[20])
            acc = float(parts[21])
            description = " ".join(parts[22:]) if len(parts) > 22 else ""

            alignment_length = ali_to - ali_from + 1
            coverage = alignment_length / target_len if target_len else 0.0

            hits.append(
                {
                    "seq_id": query_name.split()[0],
                    "seq_length": query_len,
                    "domain_name": target_name,
                    "domain_family": infer_domain_family(target_name),
                    "source": source,
                    "seq_start": ali_from,
                    "seq_end": ali_to,
                    "hmm_from": hmm_from,
                    "hmm_to": hmm_to,
                    "env_from": env_from,
                    "env_to": env_to,
                    "model_length": target_len,
                    "alignment_length": alignment_length,
                    "coverage": coverage,
                    "evalue": full_evalue,
                    "i_evalue": i_evalue,
                    "score": domain_score,
                    "bias": domain_bias,
                    "acc": acc,
                    "description": description,
                }
            )
    print(f"[PARSE] Parsed {len(hits)} hits from {domtblout.name} ({source})")
    return hits


# -----------------------------------------------------------------------------
# Domain utilities
# -----------------------------------------------------------------------------

def infer_domain_family(domain_name: str) -> Optional[str]:
    """Infer family code from domain name (e.g., AA9, CBM1)."""
    if not domain_name:
        return None
    upper = domain_name.upper()
    match = re.search(r"(AA\d+)", upper)
    if match:
        return match.group(1)
    match = re.search(r"(CBM\d+)", upper)
    if match:
        return match.group(1)
    if upper.startswith("CBM"):
        return "CBM"
    return None


def classify_hit(hit: DomainHit, config: Module3Config) -> Tuple[str, bool]:
    """Return (flag, accepted) for a hit based on thresholds."""
    if float(hit["evalue"]) > config.evalue_cutoff:
        return "high_evalue", False
    if float(hit["coverage"]) < config.coverage_cutoff:
        return "low_coverage", False
    if int(hit["alignment_length"]) < config.min_domain_length:
        return "too_short", False
    return "ok", True


def apply_filters(hits: List[DomainHit], config: Module3Config) -> List[DomainHit]:
    """Attach confidence_flag and accepted to all hits."""
    for hit in hits:
        flag, accepted = classify_hit(hit, config)
        hit["confidence_flag"] = flag
        hit["accepted"] = accepted
    return hits


def overlap_fraction(a: DomainHit, b: DomainHit) -> float:
    """Compute overlap fraction relative to the shorter domain."""
    a_start, a_end = int(a["seq_start"]), int(a["seq_end"])
    b_start, b_end = int(b["seq_start"]), int(b["seq_end"])
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start) + 1)
    if overlap <= 0:
        return 0.0
    a_len = a_end - a_start + 1
    b_len = b_end - b_start + 1
    return overlap / float(min(a_len, b_len))


def prune_overlaps(hits: Sequence[DomainHit], max_fraction: float) -> List[DomainHit]:
    """
    Remove strongly overlapping hits, keeping the best-scoring ones.
    """
    sorted_hits = sorted(hits, key=lambda h: (-float(h["score"]), -float(h["coverage"])))
    kept: List[DomainHit] = []
    for hit in sorted_hits:
        if any(overlap_fraction(hit, kept_hit) > max_fraction for kept_hit in kept):
            continue
        kept.append(hit)
    return kept


def is_lpmo_family(domain_family: Optional[str], config: Module3Config) -> bool:
    if not domain_family:
        return False
    return domain_family.upper() in {fam.upper() for fam in config.lpmo_families}


def select_lpmo_hit(
    hits: Sequence[DomainHit],
    metadata_family: Optional[str],
    config: Module3Config,
) -> Optional[DomainHit]:
    """Pick the best catalytic LPMO domain for a sequence."""
    candidates = [h for h in hits if is_lpmo_family(h.get("domain_family"), config) and h.get("accepted")]
    if config.prefer_family_from_metadata and metadata_family:
        fam_upper = metadata_family.upper()
        preferred = [h for h in candidates if str(h.get("domain_family", "")).upper() == fam_upper]
        if preferred:
            candidates = preferred
    if not candidates:
        return None
    ranked = sorted(
        candidates,
        key=lambda h: (-float(h["coverage"]), -float(h["score"]), float(h["evalue"])),
    )
    return ranked[0]


def summarise_cbms(hits: Sequence[DomainHit]) -> Tuple[bool, str]:
    """Return has_cbm flag and formatted summary string."""
    cbm_hits = [h for h in hits if str(h.get("domain_family", "")).upper().startswith("CBM") and h.get("accepted")]
    cbm_hits = prune_overlaps(cbm_hits, 0.0)  # allow stacked CBMs, but drop exact duplicates
    if not cbm_hits:
        return False, ""
    parts = [
        f"{h.get('domain_family')}({int(h['seq_start'])}-{int(h['seq_end'])})"
        for h in sorted(cbm_hits, key=lambda h: int(h["seq_start"]))
    ]
    return True, ";".join(parts)


# -----------------------------------------------------------------------------
# I/O helpers
# -----------------------------------------------------------------------------

def load_sequences(fasta_path: Path) -> Dict[str, str]:
    """Load sequences into dict keyed by seq_id (first token in header)."""
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    sequences: Dict[str, str] = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_id = record.id.split()[0]
        sequences[seq_id] = str(record.seq)
    print(f"[INPUT] Loaded {len(sequences)} sequences from {fasta_path}")
    return sequences


def load_metadata_table(metadata_path: Path) -> Dict[str, Dict[str, str]]:
    """Load Module 2 metadata into a dict keyed by seq_id."""
    if not metadata_path.exists():
        print(f"[INPUT] WARNING: Metadata file not found: {metadata_path}")
        return {}
    rows: Dict[str, Dict[str, str]] = {}
    with open(metadata_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            seq_id = (row.get("seq_id") or row.get("id") or "").split()[0]
            if seq_id:
                rows[seq_id] = row
    print(f"[INPUT] Loaded metadata for {len(rows)} sequences from {metadata_path}")
    return rows


def write_domain_table(hits: Iterable[DomainHit], output_path: Path) -> None:
    fieldnames = [
        "seq_id",
        "domain_name",
        "domain_family",
        "source",
        "seq_start",
        "seq_end",
        "model_length",
        "alignment_length",
        "coverage",
        "evalue",
        "score",
        "bias",
        "confidence_flag",
        "accepted",
        "description",
    ]
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for hit in hits:
            row = {key: hit.get(key, "") for key in fieldnames}
            writer.writerow(row)
    print(f"[OUTPUT] Domain table written to {output_path}")


def write_metadata_table(rows: Iterable[Dict[str, object]], output_path: Path) -> None:
    fieldnames = [
        "seq_id",
        "family",
        "sequence_length",
        "has_lpmo_domain",
        "lpmo_domain_family",
        "lpmo_domain_name",
        "lpmo_start",
        "lpmo_end",
        "lpmo_source",
        "has_cbm",
        "cbm_domains",
        "use_full_length",
        "use_catalytic",
        "quality_flag",
    ]
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
    print(f"[OUTPUT] Sequence metadata written to {output_path}")


def write_variant_fastas(
    sequences: Dict[str, str],
    seq_metadata: Iterable[Dict[str, object]],
    full_length_path: Path,
    catalytic_path: Path,
) -> Tuple[int, int]:
    """Write full-length and catalytic-domain FASTA files."""
    full_records: List[SeqRecord] = []
    catalytic_records: List[SeqRecord] = []

    for row in seq_metadata:
        seq_id = str(row["seq_id"])
        seq = sequences.get(seq_id)
        if not seq:
            continue

        if row.get("use_full_length"):
            record = SeqRecord(Seq(seq), id=seq_id, description="full_length")
            full_records.append(record)

        if row.get("use_catalytic") and row.get("lpmo_start") and row.get("lpmo_end"):
            try:
                start = int(row["lpmo_start"])
                end = int(row["lpmo_end"])
            except Exception:
                continue
            start_idx = max(0, start - 1)
            end_idx = min(len(seq), end)
            if end_idx > start_idx:
                fragment = seq[start_idx:end_idx]
                desc = f"catalytic lpmo_start={start} lpmo_end={end}"
                catalytic_records.append(SeqRecord(Seq(fragment), id=seq_id, description=desc))

    if full_records:
        SeqIO.write(full_records, full_length_path, "fasta")
        print(f"[OUTPUT] Wrote {len(full_records)} full-length sequences to {full_length_path}")
    else:
        print(f"[OUTPUT] WARNING: No full-length sequences selected for {full_length_path}")

    if catalytic_records:
        SeqIO.write(catalytic_records, catalytic_path, "fasta")
        print(f"[OUTPUT] Wrote {len(catalytic_records)} catalytic domains to {catalytic_path}")
    else:
        print(f"[OUTPUT] WARNING: No catalytic domains selected for {catalytic_path}")

    return len(full_records), len(catalytic_records)


# -----------------------------------------------------------------------------
# Sequence-level annotation logic
# -----------------------------------------------------------------------------

def annotate_sequence(
    seq_id: str,
    seq_len: int,
    family: Optional[str],
    hits: Sequence[DomainHit],
    config: Module3Config,
) -> Dict[str, object]:
    """Annotate a single sequence and determine variant policy."""
    accepted_hits = [h for h in hits if h.get("accepted")]
    accepted_hits = prune_overlaps(accepted_hits, config.max_domain_overlap_fraction)

    lpmo_hit = select_lpmo_hit(accepted_hits, family, config)
    flagged_candidates = [h for h in hits if is_lpmo_family(h.get("domain_family"), config)]
    cbm_flag, cbm_summary = summarise_cbms(accepted_hits)

    has_lpmo_domain = lpmo_hit is not None
    lpmo_length = (int(lpmo_hit["seq_end"]) - int(lpmo_hit["seq_start"]) + 1) if lpmo_hit else 0

    use_full_length = False
    use_catalytic = False
    quality_flag = "ok"

    if not has_lpmo_domain:
        use_full_length = bool(config.fallback_to_full_length_if_no_domain)
        use_catalytic = False
        quality_flag = "no_domain"
        if flagged_candidates:
            quality_flag = str(flagged_candidates[0].get("confidence_flag", "no_domain"))
    else:
        use_full_length = True
        if lpmo_length < config.min_catalytic_domain_length:
            quality_flag = "short_domain"
        else:
            if cbm_flag and config.generate_catalytic_for_cbm_only:
                use_catalytic = True
            elif not cbm_flag and config.generate_catalytic_if_no_cbm:
                use_catalytic = True

    return {
        "seq_id": seq_id,
        "family": family or "",
        "sequence_length": seq_len,
        "has_lpmo_domain": has_lpmo_domain,
        "lpmo_domain_family": lpmo_hit.get("domain_family") if lpmo_hit else "",
        "lpmo_domain_name": lpmo_hit.get("domain_name") if lpmo_hit else "",
        "lpmo_start": int(lpmo_hit["seq_start"]) if lpmo_hit else "",
        "lpmo_end": int(lpmo_hit["seq_end"]) if lpmo_hit else "",
        "lpmo_source": lpmo_hit.get("source") if lpmo_hit else "",
        "has_cbm": cbm_flag,
        "cbm_domains": cbm_summary,
        "use_full_length": use_full_length,
        "use_catalytic": use_catalytic,
        "quality_flag": quality_flag,
    }


# -----------------------------------------------------------------------------
# Run metadata
# -----------------------------------------------------------------------------

def summarise_counts(rows: Sequence[Dict[str, object]]) -> Dict[str, int]:
    """Aggregate useful counts for run metadata."""
    total = len(rows)
    with_lpmo = sum(1 for r in rows if r.get("has_lpmo_domain"))
    with_cbm = sum(1 for r in rows if r.get("has_cbm"))
    use_full = sum(1 for r in rows if r.get("use_full_length"))
    use_cat = sum(1 for r in rows if r.get("use_catalytic"))
    low_quality = sum(1 for r in rows if r.get("quality_flag") != "ok")
    return {
        "num_sequences": total,
        "num_with_lpmo": with_lpmo,
        "num_with_cbm": with_cbm,
        "num_use_full_length": use_full,
        "num_use_catalytic": use_cat,
        "num_low_quality": low_quality,
    }


def write_run_metadata(
    output_path: Path,
    config: Module3Config,
    counts: Dict[str, int],
    hmmer_info: Dict[str, object],
    hmms_used: Dict[str, str],
) -> None:
    """Write run metadata JSON."""
    payload = {
        "run_timestamp": datetime.utcnow().isoformat() + "Z",
        "hmmer_version": hmmer_info.get("version"),
        "hmmscan_binary": hmmer_info.get("binary"),
        "hmm_databases": hmms_used,
        "parameters": config_to_dict(config),
    }
    payload.update(counts)
    with open(output_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    print(f"[OUTPUT] Run metadata written to {output_path}")


# -----------------------------------------------------------------------------
# Pipeline entrypoint
# -----------------------------------------------------------------------------

def run_domain_pipeline(
    config: Module3Config,
    fasta_path: Path = M2_FASTA,
    metadata_path: Path = M2_METADATA,
) -> Tuple[List[DomainHit], List[Dict[str, object]]]:
    """Run the complete domain annotation workflow."""
    ensure_output_directories()

    sequences = load_sequences(fasta_path)
    metadata = load_metadata_table(metadata_path)

    hmms: List[Tuple[str, Path]] = []
    if config.dbcan_hmm:
        hmms.append(("dbcan", Path(config.dbcan_hmm)))
    if config.use_subfamily_hmms and config.dbcan_sub_hmm:
        hmms.append(("dbcan_sub", Path(config.dbcan_sub_hmm)))
    if config.cbm_hmm:
        hmms.append(("cbm", Path(config.cbm_hmm)))

    all_hits: List[DomainHit] = []
    for source, hmm_path in hmms:
        if not hmm_path.exists():
            print(f"[HMMER] WARNING: Skipping {source} (file missing: {hmm_path})")
            continue
        domtblout = RAW_DOMTBL_DIR / f"{source}.domtblout"
        if domtblout.exists() and config.reuse_existing_domtbl:
            print(f"[HMMER] Reusing existing domtblout for {source}: {domtblout}")
        else:
            success = run_hmmscan(config.hmmscan_binary, hmm_path, fasta_path, domtblout, cpu=config.hmmscan_cpu)
            if not success:
                print(f"[HMMER] ERROR: hmmscan failed for {source}, skipping its hits")
                continue
        all_hits.extend(parse_domtblout(domtblout, source))

    if not all_hits:
        print("[PIPELINE] WARNING: No domain hits found; outputs will be minimal.")

    filtered_hits = apply_filters(all_hits, config)
    write_domain_table(filtered_hits, DOMAINS_PARSED)

    hits_by_seq: Dict[str, List[DomainHit]] = defaultdict(list)
    for hit in filtered_hits:
        hits_by_seq[str(hit["seq_id"])].append(hit)

    seq_metadata_rows: List[Dict[str, object]] = []
    for seq_id, sequence in sequences.items():
        seq_len = len(sequence)
        meta_row = metadata.get(seq_id, {})
        family = meta_row.get("family") if meta_row else None
        seq_hits = hits_by_seq.get(seq_id, [])
        seq_metadata_rows.append(annotate_sequence(seq_id, seq_len, family, seq_hits, config))

    write_metadata_table(seq_metadata_rows, DOMAIN_METADATA)
    write_variant_fastas(sequences, seq_metadata_rows, FULL_LENGTH_FASTA, CATALYTIC_FASTA)

    counts = summarise_counts(seq_metadata_rows)
    hmmer_info = {
        "binary": config.hmmscan_binary,
        "version": hmmer_version(config.hmmscan_binary),
    }
    hmms_used = {name: str(path) for name, path in hmms if path and path.exists()}
    write_run_metadata(RUN_METADATA, config, counts, hmmer_info, hmms_used)

    return filtered_hits, seq_metadata_rows


__all__ = [
    "run_domain_pipeline",
    "load_sequences",
    "load_metadata_table",
    "parse_domtblout",
    "apply_filters",
    "annotate_sequence",
]
