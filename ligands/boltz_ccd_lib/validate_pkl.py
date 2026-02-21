#!/usr/bin/env python3
import sys
import pickle
from rdkit import Chem


def main():
    if len(sys.argv) < 2:
        print("Usage: validate_pkl.py <ligand.pkl>")
        sys.exit(2)
    
    pkl = sys.argv[1]
    try:
        with open(pkl, "rb") as f:
            m = pickle.load(f)
    except FileNotFoundError:
        print(f"ERROR: File not found: {pkl}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to load pickle: {e}")
        sys.exit(1)

    print("Atoms:", m.GetNumAtoms(), "Bonds:", m.GetNumBonds(), "Confs:", m.GetNumConformers())
    print("Mol props:", list(m.GetPropNames()))

    missing_names = [a.GetIdx() for a in m.GetAtoms() if not a.HasProp("name")]
    if missing_names:
        raise SystemExit(f"ERROR: Atoms missing 'name' prop: {missing_names}")

    names = [a.GetProp("name") for a in m.GetAtoms()]
    if len(set(names)) != len(names):
        raise SystemExit("ERROR: Duplicate atom names detected")

    if m.GetNumConformers() == 0:
        raise SystemExit("ERROR: No conformers found")

    # roundtrip SMILES (drops Hs); useful sanity but not required to match CIF exactly
    print("SMILES (no H):", Chem.MolToSmiles(Chem.RemoveHs(m)))

    # print first few atoms
    for a in list(m.GetAtoms())[:10]:
        print(a.GetIdx(), a.GetSymbol(), a.GetProp("name"), "charge", a.GetFormalCharge(), "arom", a.GetIsAromatic())


if __name__ == "__main__":
    main()
