#!/usr/bin/env python3
"""Fix AF3 inference inputs: replace CIF-based _patch_json.py with BAP (BondAtomPair) version.

Run from the structure_pipeline root directory:
    python fix_bap_inference.py [--dry-run]

What this script does for each oligosaccahride ligand dir in work/:
  1. Reads config/oligo_definitions.yaml to get monomer + bond_atom_pair
  2. Rewrites work/{LIGAND}/af3/input/_patch_json.py with the oligo-mode
     version (multi-monomer ccdCodes + bondedAtomPairs, no userCCDPath)
  3. Removes the 5th CIF-path argument from every _patch_json.py call
     inside work/{LIGAND}/af3/af3_{LIGAND}_inf.sh
  4. Deletes work/{LIGAND}/af3/DONE.ok so the pipeline allows re-submission

MSA data (work/af3_msa/) is NOT touched.

After running, re-submit all inference jobs:
    for sh in work/*/af3/af3_*_inf.sh; do sbatch "$sh"; done
"""

from __future__ import annotations

import re
import sys
import textwrap
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("ERROR: PyYAML not found. Activate the structure_pipeline conda environment.")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).parent
WORK_DIR = SCRIPT_DIR / "work"
OLIGO_DEFS = SCRIPT_DIR / "config" / "oligo_definitions.yaml"

# AF3 chain ID for the main (non-copper) ligand – must match AF3Runner constant
AF3_MAIN_LIGAND_ID = "C"

# Matches the 5th positional argument (CIF path) in _patch_json.py calls, e.g.:
#   /usr/bin/python3 "${INPUT_DIR}/_patch_json.py" ... "CEL4" "/root/af_input/CEL4.cif"
# The pattern intentionally anchors to /root/af_input/ to avoid false positives.
_CIF_ARG_RE = re.compile(r' "/root/af_input/[^"]+\.cif"')


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_oligo_defs(path: Path) -> dict[str, dict]:
    """Return {PREFIX: {monomer, bond_atom_pair}} from the YAML file."""
    with open(path) as fh:
        raw = yaml.safe_load(fh)
    return raw["oligo_definitions"]


def parse_oligo_code(name: str, defs: dict[str, dict]) -> tuple[str, int, dict] | None:
    """Return (prefix, n, spec) if *name* is a known oligo (e.g. 'CEL4'), else None."""
    m = re.fullmatch(r"([A-Za-z]+)(\d+)", name)
    if m is None:
        return None
    prefix = m.group(1).upper()
    n = int(m.group(2))
    spec = defs.get(prefix)
    if spec is None:
        return None
    if n < 2:
        return None
    return prefix, n, spec


def build_bonded_atom_pairs(chain_id: str, spec: dict, n: int) -> list:
    """Generate bondedAtomPairs for n residues in chain_id linked end-to-end."""
    atom_a, atom_b = spec["bond_atom_pair"]
    return [
        [[chain_id, i, atom_a], [chain_id, i + 1, atom_b]]
        for i in range(1, n)          # i = 1 … n-1  (1-based residue indices)
    ]


def write_oligo_patch_script(
    output_path: Path,
    oligo_ccd_codes: list[str],
    oligo_bonded_atom_pairs: list,
) -> None:
    """Write the oligo-mode _patch_json.py helper script to *output_path*.

    This reproduces the logic of AF3Runner._write_patch_script() for the
    oligo case, embedded with the pre-computed codes and BAP pairs so the
    script has no external dependencies at inference time.
    """
    oligo_codes_repr = repr(oligo_ccd_codes)
    oligo_bap_repr = repr(oligo_bonded_atom_pairs)

    script = textwrap.dedent(f"""\
        #!/usr/bin/env python3
        \"\"\"Patch AF3 MSA _data.json for inference (oligo mode) – auto-generated.\"\"\"
        import json, sys

        src, dst, name, ligand = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

        with open(src) as fh:
            data = json.load(fh)

        data["name"] = name

        seqs = [entry for entry in data.get("sequences", []) if "ligand" not in entry]

        # CU is always present
        seqs.append({{"ligand": {{"id": "B", "ccdCodes": ["CU"]}}}})

        # Oligo ligand with multi-monomer ccdCodes
        oligo_codes = {oligo_codes_repr}
        seqs.append({{"ligand": {{"id": "C", "ccdCodes": oligo_codes}}}})

        data["sequences"] = seqs

        # bondedAtomPairs for the oligosaccharide chain
        data["bondedAtomPairs"] = {oligo_bap_repr}

        # Remove CIF-related keys – oligo uses CCD monomers, not a CIF file
        for key in ("userCCD", "userCCDPath"):
            if key in data:
                print(f"Removing {{key}} from AF3 JSON (oligo mode)", file=sys.stderr)
                del data[key]

        with open(dst, "w") as fh:
            json.dump(data, fh, indent=2)
    """)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(script)
    output_path.chmod(0o755)


# ---------------------------------------------------------------------------
# Per-ligand fix
# ---------------------------------------------------------------------------

def fix_ligand(ligand_name: str, defs: dict, dry_run: bool) -> bool:
    """Fix one ligand directory. Returns True if processed, False if skipped."""
    result = parse_oligo_code(ligand_name, defs)
    if result is None:
        print(f"  SKIP  {ligand_name}: not a recognised oligo code – leaving untouched")
        return False

    prefix, n, spec = result
    oligo_codes = [spec["monomer"]] * n
    oligo_bap = build_bonded_atom_pairs(AF3_MAIN_LIGAND_ID, spec, n)

    af3_dir = WORK_DIR / ligand_name / "af3"
    if not af3_dir.exists():
        print(f"  SKIP  {ligand_name}: {af3_dir} not found")
        return False

    input_dir = af3_dir / "input"
    if not input_dir.exists():
        print(f"  SKIP  {ligand_name}: input/ dir not found inside {af3_dir}")
        return False

    patch_script = input_dir / "_patch_json.py"
    inf_sh = af3_dir / f"af3_{ligand_name}_inf.sh"
    done_ok = af3_dir / "DONE.ok"

    if dry_run:
        print(f"  DRY   {ligand_name} → monomer={spec['monomer']}×{n}, "
              f"bond={spec['bond_atom_pair']}")
        print(f"        would rewrite: {patch_script}")
        if inf_sh.exists():
            n_matches = len(_CIF_ARG_RE.findall(inf_sh.read_text()))
            print(f"        would remove {n_matches} CIF arg(s) from: {inf_sh.name}")
        if done_ok.exists():
            print(f"        would delete: {done_ok}")
        return True

    # 1. Rewrite _patch_json.py
    write_oligo_patch_script(patch_script, oligo_codes, oligo_bap)

    # 2. Patch inference script: remove 5th CIF-path argument
    sh_patched_count = 0
    if inf_sh.exists():
        original = inf_sh.read_text()
        sh_patched_count = len(_CIF_ARG_RE.findall(original))
        patched = _CIF_ARG_RE.sub("", original)
        inf_sh.write_text(patched)

    # 3. Delete DONE.ok
    had_done_ok = done_ok.exists()
    if had_done_ok:
        done_ok.unlink()

    print(
        f"  FIXED {ligand_name}: "
        f"monomer={spec['monomer']}×{n}, bond={spec['bond_atom_pair']}, "
        f"{sh_patched_count} CIF arg(s) removed"
        + (" | DONE.ok deleted" if had_done_ok else "")
    )
    return True


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    dry_run = "--dry-run" in sys.argv

    if not OLIGO_DEFS.exists():
        sys.exit(f"ERROR: oligo definitions not found: {OLIGO_DEFS}")
    if not WORK_DIR.exists():
        sys.exit(f"ERROR: work directory not found: {WORK_DIR}")

    defs = load_oligo_defs(OLIGO_DEFS)
    print(f"Loaded {len(defs)} oligo prefix(es): {list(defs.keys())}")
    if dry_run:
        print("DRY-RUN mode – no files will be modified\n")
    else:
        print()

    ligand_dirs = sorted(
        p.name for p in WORK_DIR.iterdir()
        if p.is_dir() and p.name != "af3_msa"
    )

    fixed = 0
    for ligand_name in ligand_dirs:
        if fix_ligand(ligand_name, defs, dry_run):
            fixed += 1

    print()
    print(f"{'Would fix' if dry_run else 'Fixed'}: {fixed}/{len(ligand_dirs)} ligands.")

    if not dry_run and fixed > 0:
        print()
        print("Next steps:")
        print("  1. Inspect a patched file, e.g.:")
        print("       work/CEL4/af3/input/_patch_json.py")
        print("  2. Submit all inference jobs:")
        print('       for sh in work/*/af3/af3_*_inf.sh; do sbatch "$sh"; done')


if __name__ == "__main__":
    main()
