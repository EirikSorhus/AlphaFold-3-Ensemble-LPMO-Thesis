"""
Signalpeptid module: run SignalP6 on Module 2 outputs, trim signal peptides,
and produce mature sequences + updated metadata.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from textwrap import wrap
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from scripts.module_signalpeptide.config_signalp import (
    SignalPConfig,
    config_to_dict,
    ensure_output_dirs,
    load_config,
)


def kingdom_to_signalp_org(kingdom: str) -> str:
    if kingdom and kingdom.strip().lower() == "eukaryota":
        return "eukarya"
    return "other"


def read_metadata(metadata_path: Path) -> Dict[str, Dict[str, str]]:
    meta: Dict[str, Dict[str, str]] = {}
    if not metadata_path.exists():
        print(f"[signalp] WARNING: Metadata file not found: {metadata_path}")
        return meta

    with open(metadata_path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        headers = [h.lower() for h in reader.fieldnames or []]
        id_fields = ["protein_id", "seq_id", "id"]
        id_field = None
        for cand in id_fields:
            if cand in headers:
                id_field = cand
                break
        if id_field is None:
            raise ValueError("No seq_id/protein_id column found in metadata")

        kingdom_field = "kingdom" if "kingdom" in headers else "organism" if "organism" in headers else None
        species_field = "species" if "species" in headers else None
        length_field = "sequence_length" if "sequence_length" in headers else None

        for row in reader:
            seq_id = row.get(id_field, "").strip()
            if not seq_id:
                continue
            meta[seq_id] = {
                "kingdom": row.get(kingdom_field, "") if kingdom_field else "",
                "species": row.get(species_field, "") if species_field else "",
                "sequence_length": row.get(length_field, "") if length_field else "",
            }
    return meta


def read_fasta(fasta_path: Path) -> List[Tuple[SeqRecord, str]]:
    records: List[Tuple[SeqRecord, str]] = []
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        tokens = rec.description.split()
        seq_id = tokens[0]
        family = tokens[1] if len(tokens) > 1 else ""
        rec.id = seq_id
        rec.description = rec.description
        records.append((rec, family))
    return records


def group_sequences(records: List[Tuple[SeqRecord, str]], meta: Dict[str, Dict[str, str]], treat_missing_as_no_sp: bool):
    groups = {"eukarya": [], "other": []}
    missing = []
    for rec, family in records:
        info = meta.get(rec.id)
        if not info or not info.get("kingdom"):
            missing.append(rec.id)
            continue
        org = kingdom_to_signalp_org(info.get("kingdom", ""))
        groups[org].append((rec, family))
    if missing and not treat_missing_as_no_sp:
        raise ValueError(f"Missing kingdom metadata for {len(missing)} sequences")
    return groups, missing


def write_group_fasta(group_records: List[Tuple[SeqRecord, str]], path: Path) -> None:
    seqs = []
    for rec, _ in group_records:
        seqs.append(rec)
    SeqIO.write(seqs, str(path), "fasta")


def chunk_list(items: List, size: int) -> List[List]:
    return [items[i : i + size] for i in range(0, len(items), size)]


def run_signalp_on_group(group_name: str, records: List[Tuple[SeqRecord, str]], config: SignalPConfig) -> List[Path]:
    output_dir = Path(config.signalp6_output_dir) / group_name
    output_dir.mkdir(parents=True, exist_ok=True)
    result_files: List[Path] = []

    if not records:
        return result_files

    # Remove stale txt outputs to avoid mixing runs
    for old_txt in output_dir.glob("*.txt"):
        try:
            old_txt.unlink()
        except OSError:
            pass

    batches = chunk_list(records, config.signalp6_max_seqs_per_run)
    for idx, batch in enumerate(batches, start=1):
        tmp_fasta = output_dir / f"tmp_{group_name}_{idx}.fasta"
        write_group_fasta(batch, tmp_fasta)

        cmd = [
            config.signalp6_path,
            "--fastafile",
            str(tmp_fasta),
            "--organism",
            group_name,
            "--mode",
            config.signalp6_mode,
            "--format",
            config.signalp6_format,
            "--output_dir",
            str(output_dir),
        ]
        print(f"[signalp] Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

        # Collect newly created txt files in output_dir
        txt_files = sorted(output_dir.glob("*.txt"), key=lambda p: p.stat().st_mtime, reverse=True)
        result_files.extend(txt_files)

    return result_files


def _parse_signalp_txt(txt_path: Path) -> Dict[str, Dict[str, Optional[int]]]:
    results: Dict[str, Dict[str, Optional[int]]] = {}
    with open(txt_path, "r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header: List[str] = []
        for row in reader:
            if not row:
                continue
            if row[0].startswith("#"):
                continue
            if not header:
                header = [col.strip().lower() for col in row]
                continue
            if len(row) < len(header):
                continue
            row_map = {header[i]: row[i].strip() for i in range(len(header))}
            seq_id = row_map.get("id") or row_map.get("sequence id") or row_map.get("sequence") or row_map.get(header[0])
            if not seq_id:
                continue
            prediction = row_map.get("prediction", "").lower()
            cs_val = None
            for key, val in row_map.items():
                key_l = key.lower()
                if "cs" in key_l and "pos" in key_l:
                    token = val.strip().split()[0]
                    if token.isdigit():
                        cs_val = int(token)
                        break
            has_sp = prediction not in {"other", "no_sp", "no signal peptide", "none", "no"}
            results[seq_id] = {
                "has_signal_peptide": has_sp,
                "sp_end": cs_val if has_sp else None,
                "raw_prediction": row_map.get("prediction", ""),
            }
    return results


def parse_signalp_outputs(txt_files: List[Path]) -> Dict[str, Dict[str, Optional[int]]]:
    combined: Dict[str, Dict[str, Optional[int]]] = {}
    for txt in txt_files:
        combined.update(_parse_signalp_txt(txt))
    return combined


def write_parsed_tsv(signalp_calls: Dict[str, Dict[str, Optional[int]]], path: Path) -> None:
    fieldnames = ["seq_id", "has_signal_peptide", "sp_end", "raw_prediction"]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for seq_id, data in sorted(signalp_calls.items()):
            writer.writerow({
                "seq_id": seq_id,
                "has_signal_peptide": data.get("has_signal_peptide", False),
                "sp_end": data.get("sp_end") if data.get("sp_end") is not None else "",
                "raw_prediction": data.get("raw_prediction", ""),
            })


def generate_mature_sequences(
    records: List[Tuple[SeqRecord, str]],
    signalp_calls: Dict[str, Dict[str, Optional[int]]],
    config: SignalPConfig,
    metadata: Dict[str, Dict[str, str]],
) -> Tuple[List[SeqRecord], Dict[str, Dict[str, int]], Dict[str, int]]:
    mature_records: List[SeqRecord] = []
    per_seq_meta: Dict[str, Dict[str, int]] = {}
    counters = {
        "total": 0,
        "with_signal_peptide": 0,
        "without_signal_peptide": 0,
        "missing_prediction": 0,
        "suspicious_length": 0,
    }

    for rec, family in records:
        counters["total"] += 1
        call = signalp_calls.get(rec.id)
        info = metadata.get(rec.id, {})
        seq_len = len(rec.seq)
        declared_len = int(info.get("sequence_length", 0)) if info.get("sequence_length") else None

        has_sp = False
        sp_end = None
        raw_pred = ""
        if call:
            has_sp = bool(call.get("has_signal_peptide", False))
            sp_end = call.get("sp_end")
            raw_pred = call.get("raw_prediction", "")
        else:
            counters["missing_prediction"] += 1

        if has_sp and sp_end:
            if sp_end < 1 or sp_end > seq_len:
                counters["suspicious_length"] += 1
                has_sp = False
                sp_end = None
            elif sp_end < config.min_signal_peptide_length or sp_end > config.max_signal_peptide_length:
                counters["suspicious_length"] += 1
            mature_start = sp_end
            mature_seq = rec.seq[sp_end:]
            signal_peptide_length = sp_end
        else:
            has_sp = False
            signal_peptide_length = 0
            mature_start = 0
            mature_seq = rec.seq

        mature_length = len(mature_seq)
        if declared_len and signal_peptide_length + mature_length != declared_len:
            counters["suspicious_length"] += 1

        counters["with_signal_peptide" if has_sp else "without_signal_peptide"] += 1

        desc_parts = [family] if family else []
        desc_parts.append("mature")
        new_rec = SeqRecord(Seq(str(mature_seq)), id=rec.id, description=" ".join(desc_parts))
        mature_records.append(new_rec)

        per_seq_meta[rec.id] = {
            "has_signal_peptide": int(has_sp),
            "signal_peptide_length": signal_peptide_length,
            "mature_start": mature_start + 1 if mature_start else 1,
            "mature_length": mature_length,
        }

    return mature_records, per_seq_meta, counters


def write_fasta(records: List[SeqRecord], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as handle:
        for rec in records:
            handle.write(f">{rec.id} {rec.description}\n")
            seq_str = str(rec.seq)
            for chunk in wrap(seq_str, 70):
                handle.write(chunk + "\n")


def update_metadata_csv(input_csv: Path, output_csv: Path, per_seq_meta: Dict[str, Dict[str, int]]) -> None:
    rows = []
    if not input_csv.exists():
        print(f"[signalp] WARNING: Input metadata not found: {input_csv}")
        return

    with open(input_csv, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        id_field = None
        for cand in ("protein_id", "seq_id", "id"):
            if cand in fieldnames:
                id_field = cand
                break
        if not id_field:
            raise ValueError("Metadata must contain an id column (protein_id/seq_id/id)")
        extra_fields = [
            "has_signal_peptide",
            "signal_peptide_length",
            "mature_start",
            "mature_length",
        ]
        out_fields = fieldnames + [f for f in extra_fields if f not in fieldnames]

        for row in reader:
            seq_id = row.get(id_field, "")
            if seq_id in per_seq_meta:
                row.update({
                    "has_signal_peptide": per_seq_meta[seq_id]["has_signal_peptide"],
                    "signal_peptide_length": per_seq_meta[seq_id]["signal_peptide_length"],
                    "mature_start": per_seq_meta[seq_id]["mature_start"],
                    "mature_length": per_seq_meta[seq_id]["mature_length"],
                })
            else:
                row.update({
                    "has_signal_peptide": 0,
                    "signal_peptide_length": 0,
                    "mature_start": 1,
                    "mature_length": row.get("sequence_length", ""),
                })
            rows.append(row)

    with open(output_csv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(rows)


def write_run_metadata(config: SignalPConfig, counters: Dict[str, int]) -> None:
    payload = {
        "timestamp": datetime.now().isoformat(),
        "signalp6_path": config.signalp6_path,
        "signalp6_mode": config.signalp6_mode,
        "signalp6_format": config.signalp6_format,
        "input_fasta": str(config.input_fasta),
        "output_fasta": str(config.output_fasta),
        "input_metadata": str(config.input_metadata),
        "output_metadata": str(config.output_metadata),
        "n_sequences_total": counters.get("total", 0),
        "n_with_signal_peptide": counters.get("with_signal_peptide", 0),
        "n_without_signal_peptide": counters.get("without_signal_peptide", 0),
        "n_missing_prediction": counters.get("missing_prediction", 0),
        "n_suspicious_sp_length": counters.get("suspicious_length", 0),
    }
    with open(config.run_metadata_json, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def run_pipeline(args: argparse.Namespace) -> int:
    config = load_config(Path(args.config) if args.config else None)

    # Override paths if provided
    if args.fasta:
        config.input_fasta = Path(args.fasta)
    if args.metadata:
        config.input_metadata = Path(args.metadata)
    if args.output_fasta:
        config.output_fasta = Path(args.output_fasta)
    if args.output_metadata:
        config.output_metadata = Path(args.output_metadata)

    ensure_output_dirs(config)

    print(f"[signalp] Input FASTA: {config.input_fasta}")
    print(f"[signalp] Input metadata: {config.input_metadata}")

    metadata = read_metadata(config.input_metadata)
    records = read_fasta(config.input_fasta)
    groups, missing = group_sequences(records, metadata, config.treat_missing_as_no_sp)

    all_txt_files: List[Path] = []
    for group_name, recs in groups.items():
        txt_files = run_signalp_on_group(group_name, recs, config)
        all_txt_files.extend(txt_files)

    signalp_calls = parse_signalp_outputs(all_txt_files)

    # Add missing predictions as negative if requested
    if config.treat_missing_as_no_sp:
        for rec, _ in records:
            if rec.id not in signalp_calls:
                signalp_calls[rec.id] = {
                    "has_signal_peptide": False,
                    "sp_end": None,
                    "raw_prediction": "",
                }

    write_parsed_tsv(signalp_calls, config.parsed_signalp_tsv)

    mature_records, per_seq_meta, counters = generate_mature_sequences(records, signalp_calls, config, metadata)
    write_fasta(mature_records, config.output_fasta)
    update_metadata_csv(config.input_metadata, config.output_metadata, per_seq_meta)
    write_run_metadata(config, counters)

    print(f"[signalp] Done. Mature FASTA: {config.output_fasta}")
    print(f"[signalp] Updated metadata: {config.output_metadata}")
    print(f"[signalp] Parsed predictions: {config.parsed_signalp_tsv}")
    print(f"[signalp] Run metadata: {config.run_metadata_json}")
    return 0


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="SignalPeptid: SignalP6 + trimming")
    parser.add_argument("--config", help="Optional config file (YAML/JSON)")
    parser.add_argument("--fasta", help="Override input FASTA path")
    parser.add_argument("--metadata", help="Override input metadata path")
    parser.add_argument("--output-fasta", help="Override output mature FASTA path")
    parser.add_argument("--output-metadata", help="Override output metadata path")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    try:
        return run_pipeline(args)
    except subprocess.CalledProcessError as exc:
        print(f"[signalp] ERROR: SignalP6 failed with exit code {exc.returncode}")
        return 1
    except Exception as exc:
        print(f"[signalp] ERROR: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
