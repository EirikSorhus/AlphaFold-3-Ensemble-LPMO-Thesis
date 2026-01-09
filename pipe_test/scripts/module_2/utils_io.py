import csv
import hashlib
import json
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def ensure_parent(path: Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def read_fasta(filepath: Path) -> Dict[str, str]:
    filepath = Path(filepath)
    seqs: Dict[str, str] = {}
    current = None
    with filepath.open("r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:]
                seqs[current] = ""
            elif current:
                seqs[current] += line
    return seqs


def write_json(data, path: Path) -> None:
    ensure_parent(Path(path))
    with Path(path).open("w") as handle:
        json.dump(data, handle, indent=2)


def write_tsv(rows: List[Dict[str, str]], fieldnames: List[str], out_path: Path) -> None:
    ensure_parent(Path(out_path))
    with Path(out_path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_fasta(records: Iterable[Tuple[str, str]], out_path: Path, mode: str = "w") -> None:
    ensure_parent(Path(out_path))
    with Path(out_path).open(mode) as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i+80] + "\n")


def seq_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()
