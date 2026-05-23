#!/bin/bash
#SBATCH --job-name=scan_af3_complex_h_links
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/scan_af3_complex_h_links_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "[ERROR] This script must be submitted with sbatch, not run directly on the login node."
    echo "[INFO] Example: sbatch tests/run_tests_scripts/scan_af3_complex_h_links.sh"
    exit 1
fi

project_root="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse"
python_bin="/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin/python"
out_dir="$project_root/tests/tests_results/complex_h_link_scan_${SLURM_JOB_ID}"

mkdir -p "$out_dir"
cd "$project_root"

"$python_bin" - <<'PY'
from __future__ import annotations

import json
import os
import re
from pathlib import Path

root = Path(".")
job_id = os.environ.get("SLURM_JOB_ID", "manual")
out_dir = root / "tests" / "tests_results" / f"complex_h_link_scan_{job_id}"
summary_path = out_dir / "complex_h_link_scan_summary.json"

files = sorted(root.glob("**/complex_H.pdb"))
atom_re = re.compile(r"^(ATOM  |HETATM)")
conect_re = re.compile(r"^CONECT")
glycan_resnames = {"NAG", "BGC", "GLC"}


def parse_atom(line: str) -> tuple[int, str, str, str, int]:
    serial = int(line[6:11])
    atom = line[12:16].strip()
    resname = line[17:20].strip()
    chain = line[21].strip() or "?"
    resseq = int(line[22:26])
    return serial, atom, resname, chain, resseq


def parse_conect(line: str) -> list[int]:
    serials: list[int] = []
    for start in range(6, len(line), 5):
        field = line[start : start + 5].strip()
        if not field:
            continue
        try:
            serials.append(int(field))
        except ValueError:
            continue
    return serials


def residue_label(key: tuple[str, int, str]) -> str:
    chain, resseq, resname = key
    return f"{resname}{resseq}.{chain}"


results: list[dict[str, object]] = []
for path in files:
    atoms: dict[int, dict[str, object]] = {}
    adjacency: dict[int, set[int]] = {}

    for line in path.read_text(errors="ignore").splitlines():
        if atom_re.match(line):
            try:
                serial, atom, resname, chain, resseq = parse_atom(line)
            except ValueError:
                continue
            atoms[serial] = {
                "atom": atom,
                "resname": resname,
                "chain": chain,
                "resseq": resseq,
            }
            adjacency.setdefault(serial, set())
            continue

        if conect_re.match(line):
            serials = parse_conect(line)
            if len(serials) < 2:
                continue
            src = serials[0]
            adjacency.setdefault(src, set())
            for dst in serials[1:]:
                adjacency.setdefault(dst, set())
                adjacency[src].add(dst)
                adjacency[dst].add(src)

    residue_atoms: dict[tuple[str, int, str], dict[str, int]] = {}
    for serial, atom_info in atoms.items():
        resname = str(atom_info["resname"])
        if resname not in glycan_resnames:
            continue
        key = (str(atom_info["chain"]), int(atom_info["resseq"]), resname)
        residue_atoms.setdefault(key, {})[str(atom_info["atom"])] = serial

    by_chain: dict[str, list[tuple[str, int, str]]] = {}
    for key in residue_atoms:
        by_chain.setdefault(key[0], []).append(key)

    ok_pairs: list[dict[str, str]] = []
    reversed_pairs: list[dict[str, str]] = []
    missing_pairs: list[dict[str, str]] = []
    ambiguous_pairs: list[dict[str, str]] = []

    for chain_id, residue_keys in by_chain.items():
        residue_keys.sort(key=lambda key: key[1])
        for left, right in zip(residue_keys, residue_keys[1:]):
            if right[1] != left[1] + 1:
                continue

            left_o4 = residue_atoms[left].get("O4")
            left_c1 = residue_atoms[left].get("C1")
            right_o4 = residue_atoms[right].get("O4")
            right_c1 = residue_atoms[right].get("C1")
            expected_link = (
                left_o4 is not None
                and right_c1 is not None
                and right_c1 in adjacency.get(left_o4, set())
            )
            reversed_link = (
                left_c1 is not None
                and right_o4 is not None
                and right_o4 in adjacency.get(left_c1, set())
            )

            pair_record = {
                "chain": chain_id,
                "left_residue": residue_label(left),
                "right_residue": residue_label(right),
                "expected_bond": f"{residue_label(left)}:O4 -> {residue_label(right)}:C1",
                "reversed_bond": f"{residue_label(left)}:C1 -> {residue_label(right)}:O4",
            }

            if expected_link and reversed_link:
                ambiguous_pairs.append(pair_record)
            elif expected_link:
                ok_pairs.append(pair_record)
            elif reversed_link:
                reversed_pairs.append(pair_record)
            else:
                missing_pairs.append(pair_record)

    if ok_pairs or reversed_pairs or missing_pairs or ambiguous_pairs:
        results.append(
            {
                "path": str(path),
                "pair_count": len(ok_pairs) + len(reversed_pairs) + len(missing_pairs) + len(ambiguous_pairs),
                "ok_pairs": ok_pairs,
                "reversed_pairs": reversed_pairs,
                "missing_pairs": missing_pairs,
                "ambiguous_pairs": ambiguous_pairs,
            }
        )

summary = {
    "job_id": job_id,
    "root": str(root.resolve()),
    "scanned_file_count": len(files),
    "files_with_glycan_pairs": len(results),
    "files_with_reversed_pairs": sum(1 for item in results if item["reversed_pairs"]),
    "files_with_missing_pairs": sum(1 for item in results if item["missing_pairs"]),
    "files_with_ambiguous_pairs": sum(1 for item in results if item["ambiguous_pairs"]),
    "files": results,
}

out_dir.mkdir(parents=True, exist_ok=True)
summary_path.write_text(json.dumps(summary, indent=2))

print(json.dumps(summary, indent=2))
print(f"[INFO] Summary written to {summary_path}")
PY

echo "[INFO] Results directory: $out_dir"
echo "[INFO] Summary: $out_dir/complex_h_link_scan_summary.json"