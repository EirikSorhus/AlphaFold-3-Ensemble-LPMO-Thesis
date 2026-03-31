#!/usr/bin/env python3
"""Mini-test: verifiser at _patch_json.py produserer korrekt BAP-JSON for ett protein+ligand.

Kjøres uten AF3 – tester kun JSON-transformasjonen.

Usage:
    python3 test_bap_patch.py [PROTEIN_ID [LIGAND]]

Eksempel:
    python3 test_bap_patch.py P9WEM3 CEL4
"""

import json
import subprocess
import sys
from pathlib import Path

WORK = Path(__file__).parent / "work"

PROTEIN = sys.argv[1] if len(sys.argv) > 1 else "P9WEM3"
LIGAND  = sys.argv[2] if len(sys.argv) > 2 else "CEL4"

src      = WORK / "af3_msa" / PROTEIN / "input" / "produced_msa.json"
patch_py = WORK / LIGAND / "af3" / "input" / "_patch_json.py"
dst      = Path(f"/tmp/test_{PROTEIN}_{LIGAND}.json")


def check(cond: bool, msg: str) -> bool:
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {msg}")
    return cond


print(f"=== BAP-patch test: {PROTEIN} + {LIGAND} ===\n")

# ── Forutsetninger ────────────────────────────────────────────────────────
ok = True
ok &= check(src.exists(),      f"MSA-kildefil finnes:   {src}")
ok &= check(patch_py.exists(), f"_patch_json.py finnes: {patch_py}")
if not ok:
    sys.exit("\nTest avbrutt – mangler kildefiler.")

# ── Vis tilstand i MSA-kildefil ───────────────────────────────────────────
print()
src_data = json.loads(src.read_text())
src_ligs = [s for s in src_data.get("sequences", []) if "ligand" in s]
print(f"MSA-kildefil (før patching):")
print(f"  sequences ligander : {[s['ligand']['ccdCodes'] for s in src_ligs]}")
print(f"  bondedAtomPairs    : {src_data.get('bondedAtomPairs')}")
print(f"  userCCD present    : {'userCCD' in src_data}")
print(f"  userCCDPath present: {'userCCDPath' in src_data}")

# ── Kjør _patch_json.py ───────────────────────────────────────────────────
print(f"\nKjører: {patch_py.name} {src.name} → {dst.name}")
result = subprocess.run(
    ["/usr/bin/python3", str(patch_py), str(src), str(dst), f"{PROTEIN}_{LIGAND}", LIGAND],
    stdout=subprocess.PIPE, stderr=subprocess.PIPE
)
result.stdout = result.stdout.decode()
result.stderr = result.stderr.decode()
if result.returncode != 0:
    print(f"  FEIL: patch-skript krasjet:\n{result.stderr}")
    sys.exit(1)
if result.stderr.strip():
    print(f"  stderr: {result.stderr.strip()}")

# ── Valider output-JSON ───────────────────────────────────────────────────
out = json.loads(dst.read_text())
out_ligs = [s for s in out.get("sequences", []) if "ligand" in s]
bap      = out.get("bondedAtomPairs", [])

print(f"\nOutput-JSON (etter patching):")
print(f"  name               : {out.get('name')}")
print(f"  sequences ligander : {[s['ligand']['ccdCodes'] for s in out_ligs]}")
print(f"  bondedAtomPairs    : {bap}")
print(f"  userCCD present    : {'userCCD' in out}")
print(f"  userCCDPath present: {'userCCDPath' in out}")
print()

# ── Assertions ────────────────────────────────────────────────────────────
passed = 0
total  = 0

def t(cond, msg):
    global passed, total
    total += 1
    if check(cond, msg):
        passed += 1

# Finn forventet monomer fra _patch_json.py-innholdet
patch_src = patch_py.read_text()
import re
m = re.search(r"oligo_codes = (\[.*?\])", patch_src)
expected_codes = eval(m.group(1)) if m else None

main_lig = next((s["ligand"] for s in out_ligs if s["ligand"]["id"] == "C"), None)
cu_lig   = next((s["ligand"] for s in out_ligs if s["ligand"]["id"] == "B"), None)

t(out.get("name") == f"{PROTEIN}_{LIGAND}",    f"name er '{PROTEIN}_{LIGAND}'")
t(cu_lig is not None and cu_lig["ccdCodes"] == ["CU"],
                                                "Kobber (CU) på chain B")
t(main_lig is not None,                        "Ligand på chain C finnes")

if main_lig and expected_codes:
    t(main_lig["ccdCodes"] == expected_codes,
      f"ccdCodes = {expected_codes}  (BAP-monomerer, ikke '{LIGAND}')")
    t("CEL4" not in main_lig["ccdCodes"]
      and "NAG4" not in main_lig["ccdCodes"]
      and "STA4" not in main_lig["ccdCodes"],
      "Ingen enkelt-kode oligo-CCD (ikke CIF-modus)")

n_expected = len(expected_codes) - 1 if expected_codes else 0
t(len(bap) == n_expected,                      f"bondedAtomPairs har {n_expected} bindinger")
t(len(bap) > 0,                                "bondedAtomPairs er ikke tom")
t("userCCD" not in out,                        "userCCD er fjernet")
t("userCCDPath" not in out,                    "userCCDPath er ikke satt")

if bap:
    first = bap[0]
    t(first[0][0] == "C" and first[1][0] == "C",
      f"Første binding kobler chain C til chain C: {first}")
    t(first[0][2] == "C1" and first[1][2] == "O4",
      f"Bindingsatomer er C1–O4: {first[0][2]}–{first[1][2]}")

# ── Resultat ──────────────────────────────────────────────────────────────
print()
print(f"=== {passed}/{total} tester bestått ===")
if passed < total:
    print("FEIL: Noen tester feilet – sjekk output over.")
    sys.exit(1)
else:
    print("OK: _patch_json.py produserer korrekt BAP JSON.")
    print(f"    Output lagret i: {dst}")
