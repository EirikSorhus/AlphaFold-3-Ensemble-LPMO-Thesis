#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import sys
import subprocess
from collections import defaultdict
from datetime import datetime
from pathlib import Path
import requests

from cazy_txt_parser import parse_cazy_txt_lines, filter_ncbi_protein_seeds, unique_seed_ids
from uniprot_id_mapping import UniProtIdMapper

CAZY_BASE = "https://www.cazy.org/IMG/cazy_data"


def read_families(args) -> list[str]:
    fams: list[str] = []
    if args.families:
        fams.extend(args.families)
    if args.families_file:
        txt = Path(args.families_file).read_text().splitlines()
        fams.extend([ln.strip() for ln in txt if ln.strip() and not ln.strip().startswith("#")])

    # dedup, preserve order
    seen = set()
    out = []
    for f in fams:
        if f not in seen:
            seen.add(f)
            out.append(f)
    return out


def download_cazy_txt(family: str, cache_dir: Path, user_agent: str, force: bool = False) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    url = f"{CAZY_BASE}/{family}.txt"
    out = cache_dir / f"{family}.txt"

    if out.exists() and not force:
        return out

    r = requests.get(url, headers={"User-Agent": user_agent}, timeout=60)
    r.raise_for_status()
    out.write_text(r.text)
    return out


def write_ids_file(ids: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(ids) + "\n")


def read_tsv_or_csv_as_dict_by_id(path: Path) -> dict[str, dict[str, str]]:
    """
    Leser output fra uniprot_metadata.py (tsv/csv).
    Returnerer dict: UniProt_ID -> row-dict
    """
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with open(path, "r") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        out = {}
        for row in reader:
            uid = row.get("UniProt_ID")
            if uid:
                out[uid] = row
        return out


def write_expanded_metadata(
    base_rows: dict[str, dict[str, str]],
    family_to_uniprot: dict[str, set[str]],
    families: list[str],
    out_path: Path,
) -> int:
    """
    Ekspanderer base metadata (én rad per UniProt_ID) til (UniProt_ID, family).
    Dette gjør at AA13 faktisk vises, fordi family kommer fra CAZy-seed, ikke fra protein_name heuristikk.
    """
    delimiter = "\t" if out_path.suffix.lower() == ".tsv" else ","
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Bygg fieldnames: UniProt_ID + family + resten (bevar rekkefølge)
    any_row = next(iter(base_rows.values())) if base_rows else {"UniProt_ID": ""}
    rest_cols = [c for c in any_row.keys() if c != "UniProt_ID"]
    fieldnames = ["UniProt_ID", "family"] + rest_cols

    rows_written = 0
    with open(out_path, "w", newline="") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        for fam in families:
            for uid in sorted(family_to_uniprot.get(fam, set())):
                if uid not in base_rows:
                    continue
                base = base_rows[uid]
                row = {"UniProt_ID": uid, "family": fam}
                for c in rest_cols:
                    row[c] = base.get(c, "")
                writer.writerow(row)
                rows_written += 1

    return rows_written


def main() -> int:
    p = argparse.ArgumentParser(description="CAZy txt -> NCBI IDs -> UniProt mapping -> UniProt metadata (module_2).")
    p.add_argument("--families", nargs="*", default=None)
    p.add_argument("--families-file", default=None)
    p.add_argument("--project-dir", default=".", help="Pek til pipe_test (default: cwd)")
    p.add_argument("--cache-dir", default=".cache", help="Cache for CAZy downloads")
    p.add_argument("--force-download", action="store_true")
    p.add_argument("--rate-limit", type=float, default=0.25)
    p.add_argument("--contact-email", default=None)
    p.add_argument("--format", choices=["tsv", "csv"], default="tsv")
    p.add_argument("--batch-size-uniprot", type=int, default=200)
    args = p.parse_args()

    families = read_families(args)
    if not families:
        print("ERROR: No families provided.", file=sys.stderr)
        return 2

    project_dir = Path(args.project_dir).resolve()
    data_dir = project_dir / "data"
    run_dir = data_dir / "run"
    meta_dir = data_dir / "metadata"

    cache_dir = project_dir / args.cache_dir

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    ua = "pipe_test-module_2/0.1"
    if args.contact_email:
        ua += f" (contact: {args.contact_email})"

    run_dir.mkdir(parents=True, exist_ok=True)
    meta_dir.mkdir(parents=True, exist_ok=True)

    # ---- 1) CAZy parse per family -> NCBI IDs ----
    family_to_ncbi: dict[str, list[str]] = defaultdict(list)

    for fam in families:
        txt_path = download_cazy_txt(fam, cache_dir=cache_dir, user_agent=ua, force=args.force_download)
        lines = txt_path.read_text().splitlines()

        seeds = parse_cazy_txt_lines(lines, family=fam)
        ncbi_seeds = filter_ncbi_protein_seeds(seeds)
        ncbi_ids = unique_seed_ids(ncbi_seeds)

        print(f"[DEBUG] {fam}: seeds={len(seeds)} ncbi_seeds={len(ncbi_seeds)} unique_ncbi_ids={len(ncbi_ids)}")
        if seeds:
            print("[DEBUG] example seed:", seeds[0].seed_id, seeds[0].seed_source)

        family_to_ncbi[fam].extend(ncbi_ids)

    all_ncbi_ids = []
    for fam in families:
        all_ncbi_ids.extend(family_to_ncbi[fam])

    ncbi_ids_unique = sorted(set(all_ncbi_ids))
    if not ncbi_ids_unique:
        print("ERROR: No NCBI-like protein IDs found in CAZy txt.", file=sys.stderr)
        return 3

    # ---- 2) NCBI -> UniProt mapping (try multiple from_db) ----
    mapper = UniProtIdMapper(user_agent=ua, rate_limit_s=args.rate_limit)

    from_db_candidates = [
        "EMBL-GenBank-DDBJ_CDS",
        "RefSeq_Protein",
    ]

    mapping: dict[str, str] = {}
    mapping_attempts: list[dict] = []

    for from_db in from_db_candidates:
        try:
            m, stats = mapper.map_ids(ncbi_ids_unique, from_db=from_db, to_db="UniProtKB", batch_size=500)
        except Exception as e:
            mapping_attempts.append(
                {"from_db": from_db, "error": str(e), "submitted": len(ncbi_ids_unique), "mapped": 0, "unmapped": len(ncbi_ids_unique)}
            )
            continue

        mapping_attempts.append(
            {"from_db": from_db, "submitted": stats.submitted, "mapped": stats.mapped, "unmapped": stats.unmapped}
        )

        if m:
            mapping = m
            print(f"[INFO] Mapping succeeded with from_db={from_db}: mapped={len(m)}/{len(ncbi_ids_unique)}")
            break

    if not mapping:
        print("ERROR: No UniProt IDs mapped. Try different UniProt 'from' database.", file=sys.stderr)
        # TODO: lagre mapping_attempts i run-json (vi gjør det allerede under run_json)
        return 4

    # family->uniprot
    family_to_uniprot: dict[str, set[str]] = {fam: set() for fam in families}
    for fam in families:
        for ncbi_id in set(family_to_ncbi[fam]):
            uid = mapping.get(ncbi_id)
            if uid:
                family_to_uniprot[fam].add(uid)

    uniprot_ids_unique = sorted(set(mapping.values()))
    if not uniprot_ids_unique:
        print("ERROR: mapping produced no UniProt IDs.", file=sys.stderr)
        return 5

    # ---- 3) Write mapped UniProt IDs file (global) ----
    ids_path = run_dir / f"mapped_uniprot_ids_{timestamp}.txt"
    write_ids_file(uniprot_ids_unique, ids_path)

    # ---- 4) Call uniprot_metadata.py (base metadata, 1 row per UniProt_ID) ----
    base_meta_path = meta_dir / f"metadata_base_{timestamp}.{args.format}"
    uniprot_stats_path = run_dir / f"uniprot_fetch_stats_{timestamp}.json"

    cmd = [
        sys.executable,
        str(project_dir / "scripts" / "module_2" / "uniprot_metadata.py"),
        "--input", str(ids_path),
        "--output", str(base_meta_path),
        "--stats-out", str(uniprot_stats_path),
        "--batch-size", str(args.batch_size_uniprot),
    ]
    if args.contact_email:
        cmd.extend(["--contact-email", args.contact_email])

    print(f"[INFO] Calling: {' '.join(cmd)}")
    proc = subprocess.run(cmd, check=False)
    if proc.returncode != 0:
        print(f"ERROR: uniprot_metadata.py failed with code {proc.returncode}", file=sys.stderr)
        return proc.returncode

    # ---- 5) Expanded metadata (add family column + expand per family) ----
    base_rows = read_tsv_or_csv_as_dict_by_id(base_meta_path)
    expanded_meta_path = meta_dir / f"metadata_expanded_{timestamp}.{args.format}"
    expanded_rows_written = write_expanded_metadata(
        base_rows=base_rows,
        family_to_uniprot=family_to_uniprot,
        families=families,
        out_path=expanded_meta_path,
    )

    # ---- 6) Run JSON ----
    # TODO: senere: lagre per-family “dropped/missing” bedre, og legg til FASTA paths når FASTA er implementert.
    run_json = {
        "run_id": timestamp,  # kan byttes til UUID senere
        "timestamp": timestamp,
        "families": families,
        "cazy": {
            "family_to_ncbi_unique": {fam: len(set(family_to_ncbi[fam])) for fam in families},
            "ncbi_total_unique_global": len(ncbi_ids_unique),
        },
        "uniprot_mapping": {
            "attempts": mapping_attempts,  # TODO: utvid med mer detaljer ved behov
            "mapped_unique_uniprot_global": len(uniprot_ids_unique),
            "family_to_uniprot_unique": {fam: len(family_to_uniprot[fam]) for fam in families},
        },
        "outputs": {
            "mapped_uniprot_ids_file": str(ids_path),
            "metadata_base_file": str(base_meta_path),
            "metadata_expanded_file": str(expanded_meta_path),
            "uniprot_fetch_stats_file": str(uniprot_stats_path),
        },
        "counts": {
            "expanded_metadata_rows_written": expanded_rows_written,
        },
        "fasta": {
            "per_family": "TODO (next step)",
            "all_families": "TODO (next step)",
        },
    }

    run_json_path = run_dir / f"run_metadata_{timestamp}.json"
    with open(run_json_path, "w") as jf:
        json.dump(run_json, jf, indent=2)

    print(f"[OK] Wrote mapped IDs:      {ids_path}")
    print(f"[OK] Wrote base metadata:  {base_meta_path}")
    print(f"[OK] Wrote expanded meta:  {expanded_meta_path}")
    print(f"[OK] Wrote run metadata:   {run_json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
