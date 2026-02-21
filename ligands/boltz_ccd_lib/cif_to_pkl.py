#!/usr/bin/env python3
"""
Convert CCD-like chem_comp mmCIF (single component) to Boltz-compatible RDKit Mol pickle.

Reads:
- _chem_comp.id
- loop_ _chem_comp_atom.* with pdbx_model_Cartn_*_ideal
- loop_ _chem_comp_bond.* with value_order + pdbx_aromatic_flag

Writes:
- <CODE>.pkl containing RDKit Mol
  - atom prop "name" = atom_id from CIF
  - one conformer with CIF ideal coordinates
  - mol prop "MOL_NAME" = CODE
"""

import sys
import pickle
from pathlib import Path

import gemmi
from rdkit import Chem
from rdkit.Chem import rdchem


BOND_ORDER_MAP = {
    "SING": rdchem.BondType.SINGLE,
    "DOUB": rdchem.BondType.DOUBLE,
    "TRIP": rdchem.BondType.TRIPLE,
    "AROM": rdchem.BondType.AROMATIC,
}


def _get_single_value(block: gemmi.cif.Block, tag: str) -> str | None:
    # gemmi returns empty string if not present; normalize to None
    try:
        v = block.find_value(tag)
    except Exception:
        return None
    v = (v or "").strip()
    return v if v else None


def cif_to_rdkit_mol(cif_path: Path):
    doc = gemmi.cif.read(str(cif_path))

    # Prefer the first non-empty data block
    blocks = [b for b in doc if b.name and b.name != "global"]
    if not blocks:
        raise ValueError(f"No data blocks found in {cif_path}")
    block = blocks[0]

    comp_id = _get_single_value(block, "_chem_comp.id")
    if not comp_id:
        raise ValueError("Missing _chem_comp.id")
    comp_id = comp_id.strip()

    # --- atoms ---
    atom_col = block.find_loop("_chem_comp_atom.atom_id")
    if atom_col is None:
        raise ValueError("Missing _chem_comp_atom loop")
    
    atom_loop = atom_col.get_loop()

    # Required columns
    req_atom_cols = [
        "_chem_comp_atom.atom_id",
        "_chem_comp_atom.type_symbol",
        "_chem_comp_atom.charge",
        "_chem_comp_atom.pdbx_model_Cartn_x_ideal",
        "_chem_comp_atom.pdbx_model_Cartn_y_ideal",
        "_chem_comp_atom.pdbx_model_Cartn_z_ideal",
    ]
    for c in req_atom_cols:
        if c not in atom_loop.tags:
            raise ValueError(f"Missing atom column: {c}")

    atom_id_col = atom_loop.tags.index("_chem_comp_atom.atom_id")
    elem_col = atom_loop.tags.index("_chem_comp_atom.type_symbol")
    chg_col = atom_loop.tags.index("_chem_comp_atom.charge")
    x_col = atom_loop.tags.index("_chem_comp_atom.pdbx_model_Cartn_x_ideal")
    y_col = atom_loop.tags.index("_chem_comp_atom.pdbx_model_Cartn_y_ideal")
    z_col = atom_loop.tags.index("_chem_comp_atom.pdbx_model_Cartn_z_ideal")

    atoms = []  # list of dicts in CIF order
    atom_index = {}  # atom_id -> idx
    for i in range(atom_loop.length()):
        atom_id = str(atom_loop[i, atom_id_col]).strip()
        elem = str(atom_loop[i, elem_col]).strip()
        charge_str = str(atom_loop[i, chg_col]).strip()
        x_str = str(atom_loop[i, x_col]).strip()
        y_str = str(atom_loop[i, y_col]).strip()
        z_str = str(atom_loop[i, z_col]).strip()
        
        # Handle CIF missing values
        if charge_str in ("?", "."):
            charge = 0
        else:
            charge = int(charge_str)
        
        if any(c in ("?", ".") for c in (x_str, y_str, z_str)):
            raise ValueError(f"Atom {atom_id} has missing ideal coordinates")
        
        x = float(x_str)
        y = float(y_str)
        z = float(z_str)

        if atom_id in atom_index:
            raise ValueError(f"Duplicate atom_id in CIF: {atom_id}")

        atom_index[atom_id] = i
        atoms.append({"atom_id": atom_id, "elem": elem, "charge": charge, "x": x, "y": y, "z": z})

    if not atoms:
        raise ValueError("No atoms parsed from _chem_comp_atom loop")

    # --- bonds ---
    bond_col = block.find_loop("_chem_comp_bond.atom_id_1")
    if bond_col is None:
        raise ValueError("Missing _chem_comp_bond loop")
    
    bond_loop = bond_col.get_loop()

    req_bond_cols = [
        "_chem_comp_bond.atom_id_1",
        "_chem_comp_bond.atom_id_2",
        "_chem_comp_bond.value_order",
        "_chem_comp_bond.pdbx_aromatic_flag",
    ]
    for c in req_bond_cols:
        if c not in bond_loop.tags:
            raise ValueError(f"Missing bond column: {c}")

    a1_col = bond_loop.tags.index("_chem_comp_bond.atom_id_1")
    a2_col = bond_loop.tags.index("_chem_comp_bond.atom_id_2")
    ord_col = bond_loop.tags.index("_chem_comp_bond.value_order")
    arom_col = bond_loop.tags.index("_chem_comp_bond.pdbx_aromatic_flag")

    bonds = []
    for i in range(bond_loop.length()):
        a1 = str(bond_loop[i, a1_col]).strip()
        a2 = str(bond_loop[i, a2_col]).strip()
        order = str(bond_loop[i, ord_col]).strip().upper()
        arom = str(bond_loop[i, arom_col]).strip().upper()

        if a1 not in atom_index or a2 not in atom_index:
            raise ValueError(f"Bond references unknown atom_id: {a1}-{a2}")

        if order not in BOND_ORDER_MAP:
            raise ValueError(f"Unknown bond value_order: {order}")

        bonds.append((a1, a2, order, arom))

    if not bonds:
        raise ValueError("No bonds parsed from _chem_comp_bond loop")

    # --- build RDKit mol ---
    rw = Chem.RWMol()
    for a in atoms:
        rd_atom = Chem.Atom(a["elem"])
        rd_atom.SetFormalCharge(a["charge"])
        rw.AddAtom(rd_atom)

    # add bonds
    for a1, a2, order, arom in bonds:
        i1 = atom_index[a1]
        i2 = atom_index[a2]
        bt = BOND_ORDER_MAP[order]
        rw.AddBond(i1, i2, bt)
        b = rw.GetBondBetweenAtoms(i1, i2)
        if bt == rdchem.BondType.AROMATIC or arom == "Y":
            b.SetIsAromatic(True)
            rw.GetAtomWithIdx(i1).SetIsAromatic(True)
            rw.GetAtomWithIdx(i2).SetIsAromatic(True)

    mol = rw.GetMol()

    # atom names from CIF atom_id
    for a in atoms:
        mol.GetAtomWithIdx(atom_index[a["atom_id"]]).SetProp("name", a["atom_id"])

    # conformer with CIF coords (ideal)
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    for a in atoms:
        idx = atom_index[a["atom_id"]]
        conf.SetAtomPosition(idx, (a["x"], a["y"], a["z"]))
    mol.AddConformer(conf, assignId=True)

    # minimal identity
    mol.SetProp("MOL_NAME", comp_id)

    # conservative sanitize
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        raise ValueError(f"RDKit sanitization failed: {e}")

    # ensure props get pickled
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

    return comp_id, mol


def main():
    if len(sys.argv) < 2:
        print("Usage: cif_to_pkl.py <ligand.cif> [out_dir]")
        sys.exit(2)

    cif_path = Path(sys.argv[1]).resolve()
    if not cif_path.exists():
        print(f"ERROR: Input file not found: {cif_path}")
        sys.exit(1)
    
    out_dir = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else Path(".").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        comp_id, mol = cif_to_rdkit_mol(cif_path)
    except Exception as e:
        print(f"ERROR: Failed to convert {cif_path.name}: {e}")
        sys.exit(1)
    
    out_path = out_dir / f"{comp_id}.pkl"

    try:
        with open(out_path, "wb") as f:
            pickle.dump(mol, f)
    except Exception as e:
        print(f"ERROR: Failed to write pickle file: {e}")
        sys.exit(1)

    print(f"Wrote {out_path}")
    print(f"  atoms={mol.GetNumAtoms()} bonds={mol.GetNumBonds()} conformers={mol.GetNumConformers()}")
    print(f"  MOL_NAME={mol.GetProp('MOL_NAME') if mol.HasProp('MOL_NAME') else '(missing)'}")


if __name__ == "__main__":
    main()
