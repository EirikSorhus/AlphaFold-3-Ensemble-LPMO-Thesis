#!/bin/bash
#SBATCH --job-name=run_tests_mapping_contracts
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_mapping_contracts_%j.log

set -euo pipefail

# Activate conda environment
export PATH="/cluster/work/projects/nn1003k/eirik/conda/analyse_env/bin:$PATH"

# Paths
project_root="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse"
cd "$project_root"
export PYTHONPATH="$project_root/src:${PYTHONPATH:-}"

test_cif="/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/NAG8/af3/runs/408004/A0A0A1ED04_NAG8/seed-1_sample-0/A0A0A1ED04_NAG8_seed-1_sample-0_model.cif"
out_dir="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/tests_results/"

if [[ ! -f "$test_cif" ]]; then
	echo "[ERROR] test_cif does not exist: $test_cif"
	exit 1
fi

mkdir -p "$out_dir"

job_suffix="${SLURM_JOB_ID:-manual}"
run_dir="$out_dir/mapping_contracts_${job_suffix}"
mkdir -p "$run_dir"

echo "[INFO] Running mapping unit tests"
pytest tests/test_atom_mapping.py -v

echo "[INFO] Running AF3 mapping contract check on: $test_cif"
export TEST_CIF="$test_cif"
export RUN_DIR="$run_dir"
python - <<'PY'
import json
import os
import sys
from pathlib import Path

import gemmi

from lpmo_pipeline.io.normalize_mmcif import NormalizeMMCIFRunner
from lpmo_pipeline.mapping.cross_model_atom_mapping import (
    AtomMapping,
    build_atom_mapping,
    write_atom_map_tsv,
    write_rename_log,
)
from lpmo_pipeline.mapping.rename_atoms import apply_rename, load_atom_map_tsv


def fail(msg: str) -> None:
    print(f"[ERROR] {msg}")
    sys.exit(1)


def find_residue_atom_names(
    structure: gemmi.Structure,
    chain_name: str,
    resnum: int,
) -> list[str]:
    for model in structure:
        for chain in model:
            if chain.name != chain_name:
                continue
            for residue in chain:
                if residue.seqid.num == resnum:
                    return [atom.name for atom in residue]
    return []


def choose_unique_temp_name(old_name: str, existing_names: list[str]) -> str:
    base = old_name.strip() or "X"
    candidates = [
        ("X" + base)[:4],
        ("Q" + base)[:4],
        (base[:3] + "X")[:4],
        "ZZZ1",
        "ZZZ2",
    ]
    for cand in candidates:
        if cand != old_name and cand not in existing_names:
            return cand
    fail(f"Could not find unique temporary atom name for {old_name}")


test_cif = Path(os.environ["TEST_CIF"])
run_dir = Path(os.environ["RUN_DIR"])
normalize_out = run_dir / "normalize_for_mapping"
normalize_out.mkdir(parents=True, exist_ok=True)

# Prepare normalized AF3 CIF so chain IDs are canonical before mapping checks.
ok, normalized_path = NormalizeMMCIFRunner(test_cif, normalize_out).run()
if ok is not True:
    fail("Normalization failed before mapping checks")
if not normalized_path or not Path(normalized_path).exists():
    fail("normalized.cif missing before mapping checks")

normalized_path = Path(normalized_path)
structure = gemmi.read_structure(str(normalized_path))

# AF3-only contract: names are expected to already be canonical.
# Build residue CCD reference keys from observed atom names (identity reference).
ccd_graphs: dict[str, dict[str, list[str]]] = {}
for model in structure:
    for chain in model:
        for residue in chain:
            comp_id = residue.name
            ref = ccd_graphs.setdefault(comp_id, {})
            for atom in residue:
                ref.setdefault(atom.name, [])

result = build_atom_mapping(structure, ccd_graphs)

atom_map_tsv = run_dir / "atom_map.tsv"
rename_log_json = run_dir / "rename_log.json"
write_atom_map_tsv(result, atom_map_tsv)
write_rename_log(result, rename_log_json)

if result.coverage != 1.0:
    fail(f"Expected mapping coverage 1.0, got {result.coverage}")
if len(result.unmapped) != 0:
    fail(f"Expected zero unmapped atoms, got {len(result.unmapped)}")
if len(result.mappings) == 0:
    fail("No mapping entries were generated")

loaded_mappings = load_atom_map_tsv(atom_map_tsv)
if len(loaded_mappings) != len(result.mappings):
    fail("Loaded mapping count does not match generated mapping count")

# Round-trip validation for rename_atoms.py:
# force one reversible rename on a known atom and verify restoration.
target = loaded_mappings[0]
start_names = find_residue_atom_names(structure, target.chain, target.resnum)
if not start_names:
    fail(
        f"Could not find residue for target mapping chain={target.chain} resnum={target.resnum}"
    )

temp_name = choose_unique_temp_name(target.old_atom_name, start_names)

forward_mappings = []
for m in loaded_mappings:
    if (
        m.chain == target.chain
        and m.resnum == target.resnum
        and m.old_atom_name == target.old_atom_name
    ):
        forward_mappings.append(
            AtomMapping(
                chain=m.chain,
                resname=m.resname,
                resnum=m.resnum,
                old_atom_name=m.old_atom_name,
                new_atom_name=temp_name,
                element=m.element,
                match_method="round_trip_forward",
            )
        )
    else:
        forward_mappings.append(m)

structure_fwd = gemmi.read_structure(str(normalized_path))
structure_fwd, n_fwd = apply_rename(structure_fwd, forward_mappings)
if n_fwd < 1:
    fail("Forward rename did not rename any atoms")

names_after_fwd = find_residue_atom_names(structure_fwd, target.chain, target.resnum)
if temp_name not in names_after_fwd:
    fail("Forward rename did not produce expected temporary atom name")

forward_cif = run_dir / "renamed_forward.cif"
structure_fwd.make_mmcif_document().write_file(str(forward_cif))

reverse_mapping = [
    AtomMapping(
        chain=target.chain,
        resname=target.resname,
        resnum=target.resnum,
        old_atom_name=temp_name,
        new_atom_name=target.old_atom_name,
        element=target.element,
        match_method="round_trip_reverse",
    )
]

structure_rev, n_rev = apply_rename(structure_fwd, reverse_mapping)
if n_rev < 1:
    fail("Reverse rename did not rename any atoms")

names_after_rev = find_residue_atom_names(structure_rev, target.chain, target.resnum)
if target.old_atom_name not in names_after_rev:
    fail("Reverse rename did not restore original atom name")
if temp_name in names_after_rev:
    fail("Temporary atom name still present after reverse rename")

roundtrip_cif = run_dir / "renamed_roundtrip.cif"
structure_rev.make_mmcif_document().write_file(str(roundtrip_cif))

summary = {
    "ok": True,
    "af3_only": True,
    "normalized_cif": str(normalized_path),
    "atom_map_tsv": str(atom_map_tsv),
    "rename_log_json": str(rename_log_json),
    "mapping_coverage": result.coverage,
    "mapped_atoms": len(result.mappings),
    "unmapped_atoms": len(result.unmapped),
    "round_trip": {
        "target": {
            "chain": target.chain,
            "resnum": target.resnum,
            "old_atom_name": target.old_atom_name,
            "temp_atom_name": temp_name,
        },
        "forward_renamed_atoms": n_fwd,
        "reverse_renamed_atoms": n_rev,
        "forward_cif": str(forward_cif),
        "roundtrip_cif": str(roundtrip_cif),
    },
}

summary_path = run_dir / "mapping_contract_summary.json"
summary_path.write_text(json.dumps(summary, indent=2))

print("[INFO] Mapping contract summary")
print(json.dumps(summary, indent=2))
PY

echo "[INFO] Mapping contract checks completed successfully"
echo "[INFO] Results directory: $run_dir"
