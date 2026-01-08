#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime
from typing import Dict, List, Tuple
import requests

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"


def fetch_uniprot_tsv_records(
    ids: List[str],
    fields: List[str],
    batch_size: int,
    user_agent: str,
    timeout_s: int = 60,
) -> Tuple[Dict[str, Dict[str, str]], Dict]:
    """
    Henter UniProt TSV for input-accessions.
    Viktig: setter size=len(chunk) for å unngå default paging.
    Kjører i batches for robusthet.
    """
    session = requests.Session()
    headers = {"User-Agent": user_agent}

    records: Dict[str, Dict[str, str]] = {}
    total_submitted = 0

    for start in range(0, len(ids), batch_size):
        chunk = ids[start : start + batch_size]
        total_submitted += len(chunk)

        query = " OR ".join(f"(accession:{uid})" for uid in chunk)
        params = {
            "format": "tsv",
            "fields": ",".join(fields),
            "query": query,
            "size": str(len(chunk)),  # kritisk: unngå paging-tap
        }

        r = session.get(UNIPROT_SEARCH, params=params, headers=headers, timeout=timeout_s)
        r.raise_for_status()

        lines = [ln for ln in r.text.splitlines() if ln.strip()]
        if len(lines) < 2:
            continue

        header = lines[0].split("\t")
        for ln in lines[1:]:
            cols = ln.split("\t")
            if len(cols) != len(header):
                continue
            row = dict(zip(header, cols))

            # UniProt TSV bruker ofte "Entry" som accession-kolonne
            uid = row.get("Entry") or row.get("accession") or row.get("Entry accession")
            if uid:
                records[uid] = row

    found = len(records)
    missing = total_submitted - found
    stats = {
        "uniprot_ids_submitted": total_submitted,
        "uniprot_ids_found": found,
        "uniprot_ids_missing": missing,
        "timestamp": datetime.now().isoformat(),
        "batch_size": batch_size,
        "fields": fields,
    }
    return records, stats


def main() -> int:
    p = argparse.ArgumentParser(description="Fetch UniProt metadata for a list of UniProt accessions.")
    p.add_argument("--input", "-i", required=True, help="Text file with UniProt IDs (one per line).")
    p.add_argument("--output", "-o", required=True, help="Output TSV/CSV path (writer uses delimiter based on extension).")
    p.add_argument("--stats-out", default=None, help="Write JSON stats for this step.")
    p.add_argument("--batch-size", type=int, default=200)
    p.add_argument("--contact-email", default=None)
    args = p.parse_args()

    ids = [ln.strip() for ln in open(args.input) if ln.strip()]
    if not ids:
        raise SystemExit("No UniProt IDs in input.")

    # Felter: hold dette minimalt og stabilt. Utvid senere når du vet hva som finnes programmatisk.
    fields = [
        "accession",
        "protein_name",
        "gene_names",
        "organism_name",
        "xref_pdb",
        "ec",
        "xref_pfam",
        # ft_signal gir ofte lite/ustabilt i TSV; vi lar det være ute her for nå
        # og heller legger signalp-estimat i en senere module_signalpeptide.
    ]

    ua = "pipe_test-module_2/0.1"
    if args.contact_email:
        ua += f" (contact: {args.contact_email})"

    records, stats = fetch_uniprot_tsv_records(
        ids=ids,
        fields=fields,
        batch_size=args.batch_size,
        user_agent=ua,
    )

    # Bygg output-rader kun for IDs som faktisk ble funnet (dropp NOT_FOUND)
    rows: List[Dict[str, str]] = []
    for uid in ids:
        if uid not in records:
            continue
        entry = records[uid]

        protein_name = entry.get("Protein names") or entry.get("Protein name") or entry.get("protein_name") or ""
        gene_name = entry.get("Gene Names") or entry.get("Gene names") or entry.get("gene_names") or ""
        organism = entry.get("Organism") or entry.get("Organism name") or entry.get("organism_name") or ""
        pdb_field = entry.get("Cross-reference (PDB)") or entry.get("PDB") or entry.get("xref_pdb") or ""
        ec_field = entry.get("EC number") or entry.get("EC") or entry.get("ec") or ""
        pfam_field = entry.get("Cross-reference (Pfam)") or entry.get("Pfam") or entry.get("xref_pfam") or ""

        rows.append(
            {
                "UniProt_ID": uid,
                "Protein_Name": protein_name,
                "Gene_Name": gene_name,
                "Organism": organism,
                "EC_Number": ec_field,
                "PDB_IDs": pdb_field.replace(",", ";").replace("; ", ";").strip(),
                "Pfam_IDs": pfam_field.replace(",", ";").replace("; ", ";").strip(),
            }
        )

    # Dropp kolonner som er helt tomme
    columns_to_drop: List[str] = []
    if rows:
        all_cols = list(rows[0].keys())
        for col in all_cols:
            if all((not r.get(col)) for r in rows):
                columns_to_drop.append(col)
        if columns_to_drop:
            for r in rows:
                for col in columns_to_drop:
                    r.pop(col, None)

    # Velg delimiter ut fra extension
    out_lower = args.output.lower()
    delimiter = "\t"
    if out_lower.endswith(".csv"):
        delimiter = ","

    with open(args.output, "w", newline="") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=list(rows[0].keys()) if rows else ["UniProt_ID"], delimiter=delimiter)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    print(f"Metadata file written to {args.output} ({len(rows)} entries).")
    if columns_to_drop:
        print(f"Dropping columns with no data: {', '.join(columns_to_drop)}")

    if args.stats_out:
        stats["metadata_rows_written"] = len(rows)
        stats["columns_dropped"] = columns_to_drop
        with open(args.stats_out, "w") as jf:
            json.dump(stats, jf, indent=2)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
