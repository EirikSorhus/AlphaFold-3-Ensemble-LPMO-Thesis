"""Convert MOL (V2000/V3000) files into CIF components for AlphaFold 3 custom CCDs."""

from __future__ import annotations

import argparse
import re
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

ATOMIC_WEIGHTS = {
    "H": 1.00784,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.999,
    "P": 30.973762,
    "S": 32.065,
    "F": 18.998403,
    "CL": 35.453,
    "BR": 79.904,
    "I": 126.90447,
}

BOND_ORDER_MAP = {
    "1": "SING",
    "S": "SING",
    "2": "DOUB",
    "D": "DOUB",
    "3": "TRIP",
    "T": "TRIP",
    "4": "AROM",
    "A": "AROM",
    "AR": "AROM",
}

CHARGE_CODE_MAP = {
    0: 0,
    1: 3,
    2: 2,
    3: 1,
    4: 0,
    5: -1,
    6: -2,
    7: -3,
}

ATOM_ID_MAX_LEN = 4


def sanitize_code(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]", "", value.upper())
    return cleaned or "COMP"


def sanitize_atom_id(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]", "", value.upper())
    return cleaned or "ATOM"


def parse_mol(text: str) -> Tuple[List[dict], List[dict]]:
    if "V3000" in text.upper():
        return parse_v3000_mol(text)
    return parse_v2000_mol(text)


def parse_v3000_mol(text: str) -> Tuple[List[dict], List[dict]]:
    atoms: List[dict] = []
    bonds: List[dict] = []
    section: str | None = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line.startswith("M  V30"):
            continue
        payload = line[6:].strip()
        parts = payload.split()
        if not parts:
            continue
        keyword = parts[0]
        if keyword == "BEGIN" and len(parts) > 1:
            section = parts[1]
            continue
        if keyword == "END":
            section = None
            continue
        if keyword == "COUNTS":
            continue
        if section == "ATOM":
            atoms.append(_parse_v3000_atom(parts))
        elif section == "BOND":
            bonds.append(_parse_v3000_bond(parts))
    if not atoms:
        raise ValueError("No atoms parsed from MOL content; ensure V3000 format")
    return atoms, bonds


def _parse_v3000_atom(parts: List[str]) -> dict:
    if len(parts) < 6:
        raise ValueError(f"Malformed atom line: {' '.join(parts)}")
    idx = int(parts[0])
    element = parts[1]
    x, y, z = map(float, parts[2:5])
    props = _parse_key_values(parts[5:])
    charge = int(props.get("CHG", "0"))
    label = props.get("ATOMNAME") or props.get("LAB") or ""
    return {"idx": idx, "element": element, "x": x, "y": y, "z": z, "charge": charge, "label": label}


def _parse_v3000_bond(parts: List[str]) -> dict:
    if len(parts) < 4:
        raise ValueError(f"Malformed bond line: {' '.join(parts)}")
    order = parts[1].upper()
    a1 = int(parts[2])
    a2 = int(parts[3])
    props = _parse_key_values(parts[4:])
    return {"order": order, "a1": a1, "a2": a2, "props": props}


def _parse_key_values(tokens: List[str]) -> Dict[str, str]:
    props: Dict[str, str] = {}
    for token in tokens:
        if token in {"0", "-", "."}:
            continue
        if "=" in token:
            key, value = token.split("=", 1)
            props[key.upper()] = value
    return props


def parse_v2000_mol(text: str) -> Tuple[List[dict], List[dict]]:
    lines = text.splitlines()
    counts_idx = None
    for idx, line in enumerate(lines):
        if "V2000" in line:
            counts_idx = idx
            counts_line = line
            break
    if counts_idx is None:
        raise ValueError("Unable to locate V2000 counts line in MOL file")
    counts_prefix = counts_line.split("V2000", 1)[0]
    count_tokens = counts_prefix.split()
    if len(count_tokens) < 2:
        raise ValueError("Malformed V2000 counts line")
    
    # Handle cases where atom and bond counts are concatenated (e.g., "108112" instead of "108 112")
    first_token = count_tokens[0]
    if len(first_token) > 3 and first_token.isdigit():
        # Try to split concatenated counts (assume 3 digits each for atom/bond counts)
        mid = len(first_token) // 2
        num_atoms = int(first_token[:mid])
        num_bonds = int(first_token[mid:])
    else:
        num_atoms = int(count_tokens[0])
        num_bonds = int(count_tokens[1])
    atom_start = counts_idx + 1
    atom_end = atom_start + num_atoms
    bond_end = atom_end + num_bonds
    atom_lines = lines[atom_start:atom_end]
    bond_lines = lines[atom_end:bond_end]
    if len(atom_lines) != num_atoms or len(bond_lines) != num_bonds:
        raise ValueError("Unexpected end of file while parsing V2000 MOL")
    atoms = [_parse_v2000_atom(i + 1, atom_lines[i]) for i in range(num_atoms)]
    bonds = [_parse_v2000_bond(line) for line in bond_lines]
    _apply_m_chg_records(atoms, lines[bond_end:])
    return atoms, bonds


def _parse_v2000_atom(idx: int, line: str) -> dict:
    parts = line.split()
    if len(parts) < 4:
        raise ValueError(f"Malformed V2000 atom line: {line}")
    x, y, z = map(float, parts[:3])
    element = parts[3]
    charge_code = _safe_int(line[36:39])
    if charge_code is None and len(parts) > 6:
        charge_code = _safe_int(parts[6])
    charge = CHARGE_CODE_MAP.get(charge_code, 0) if charge_code is not None else 0
    label = line[60:66].strip()
    return {"idx": idx, "element": element, "x": x, "y": y, "z": z, "charge": charge, "label": label}


def _parse_v2000_bond(line: str) -> dict:
    parts = line.split()
    
    # Handle concatenated atom indices (e.g., "91105" instead of "91 105")
    # Normal bond line has 7 tokens: a1 a2 order flag flag flag flag
    # If we have 6 tokens and first token is long digits, it's likely concatenated
    if len(parts) == 6 and len(parts[0]) > 3 and parts[0].isdigit():
        # Likely concatenated: split in half
        first_token = parts[0]
        mid = len(first_token) // 2
        a1 = int(first_token[:mid])
        a2 = int(first_token[mid:])
        order = parts[1]
    elif len(parts) < 3:
        raise ValueError(f"Malformed V2000 bond line: {line}")
    else:
        a1 = int(parts[0])
        a2 = int(parts[1])
        order = parts[2]
    return {"order": order, "a1": a1, "a2": a2, "props": {}}


def _apply_m_chg_records(atoms: List[dict], trailer_lines: List[str]) -> None:
    atoms_by_idx = {atom["idx"]: atom for atom in atoms}
    for raw in trailer_lines:
        stripped = raw.strip()
        if not stripped.startswith("M  CHG"):
            continue
        tokens = stripped.split()
        if len(tokens) < 4:
            continue
        count = _safe_int(tokens[2]) or 0
        values = tokens[3:]
        for i in range(count):
            pos = 2 * i
            if pos + 1 >= len(values):
                break
            atom_idx = _safe_int(values[pos])
            charge = _safe_int(values[pos + 1])
            if atom_idx is None or charge is None:
                continue
            if atom_idx in atoms_by_idx:
                atoms_by_idx[atom_idx]["charge"] = charge


def _safe_int(value: str | None) -> int | None:
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def build_formula(atoms: List[dict]) -> Tuple[str, float]:
    counts = Counter(atom["element"].upper() for atom in atoms)
    ordered: List[Tuple[str, int]] = []
    for symbol in ("C", "H"):
        if symbol in counts:
            ordered.append((symbol, counts.pop(symbol)))
    ordered.extend(sorted(counts.items()))
    parts = [f"{sym}{cnt if cnt != 1 else ''}" for sym, cnt in ordered]
    formula = " ".join(parts) if parts else "?"
    weight = 0.0
    for sym, cnt in Counter(atom["element"].upper() for atom in atoms).items():
        weight += ATOMIC_WEIGHTS.get(sym, 0.0) * cnt
    return formula, weight if weight else float("nan")


def quote_if_needed(value: str) -> str:
    if value in {"?", "."}:
        return value
    if re.search(r"\s", value):
        return f"'{value}'"
    return value


def mol_to_cif_text(path: Path) -> Tuple[str, str]:
    name = path.stem
    comp_id = sanitize_code(name)
    atoms, bonds = parse_mol(path.read_text())

    # assign atom ids
    element_counts: Counter[str] = Counter()
    atom_ids: Dict[int, str] = {}
    atom_rows: List[str] = []
    total_charge = 0
    for ordinal, atom in enumerate(atoms, start=1):
        element = atom["element"].upper()
        total_charge += atom["charge"]
        atom_id = _assign_atom_id(atom, element_counts, atom_ids.values())
        atom_ids[atom["idx"]] = atom_id
        atom_rows.append(
            f"{comp_id} {atom_id} {element} {atom['charge']} N "
            f"{atom['x']:.4f} {atom['y']:.4f} {atom['z']:.4f} {ordinal}"
        )

    bond_rows: List[str] = []
    for ordinal, bond in enumerate(bonds, start=1):
        atom_id_1 = atom_ids[bond["a1"]]
        atom_id_2 = atom_ids[bond["a2"]]
        order = BOND_ORDER_MAP.get(bond["order"], "SING")
        arom_flag = "Y" if order == "AROM" else "N"
        bond_rows.append(f"{comp_id} {atom_id_1} {atom_id_2} {order} {arom_flag} {ordinal}")

    formula, formula_weight = build_formula(atoms)
    if formula_weight != formula_weight:  # NaN check
        formula_weight_value = "?"
    else:
        formula_weight_value = f"{formula_weight:.4f}"

    lines: List[str] = []
    lines.append(f"data_{comp_id}")
    lines.append(f"_chem_comp.id {comp_id}")
    lines.append(f"_chem_comp.name {quote_if_needed(name or '?')}")
    lines.append("_chem_comp.type 'non-polymer'")
    lines.append(f"_chem_comp.formula {quote_if_needed(formula or '?')}")
    lines.append("_chem_comp.mon_nstd_parent_comp_id ?")
    lines.append("_chem_comp.pdbx_synonyms ?")
    lines.append(f"_chem_comp.formula_weight {formula_weight_value}")
    lines.append(f"_chem_comp.pdbx_formal_charge {total_charge}")
    lines.append("_chem_comp.pdbx_release_status 'PRD'")
    lines.append("#")

    lines.append("loop_")
    lines.append("_chem_comp_atom.comp_id")
    lines.append("_chem_comp_atom.atom_id")
    lines.append("_chem_comp_atom.type_symbol")
    lines.append("_chem_comp_atom.charge")
    lines.append("_chem_comp_atom.pdbx_leaving_atom_flag")
    lines.append("_chem_comp_atom.pdbx_model_Cartn_x_ideal")
    lines.append("_chem_comp_atom.pdbx_model_Cartn_y_ideal")
    lines.append("_chem_comp_atom.pdbx_model_Cartn_z_ideal")
    lines.append("_chem_comp_atom.pdbx_ordinal")
    lines.extend(atom_rows)
    lines.append("#")

    if bond_rows:
        lines.append("loop_")
        lines.append("_chem_comp_bond.comp_id")
        lines.append("_chem_comp_bond.atom_id_1")
        lines.append("_chem_comp_bond.atom_id_2")
        lines.append("_chem_comp_bond.value_order")
        lines.append("_chem_comp_bond.pdbx_aromatic_flag")
        lines.append("_chem_comp_bond.pdbx_ordinal")
        lines.extend(bond_rows)
        lines.append("#")

    return comp_id, "\n".join(lines) + "\n"


def _assign_atom_id(atom: dict, element_counts: Counter[str], used_ids) -> str:
    element = atom["element"].upper()
    element_counts[element] += 1
    ordinal = element_counts[element]
    label = atom.get("label", "").strip()
    candidate = _normalize_atom_id(label)
    if not candidate:
        candidate = _generate_default_atom_id(element, ordinal)
    candidate = _ensure_length(candidate, element, ordinal)
    existing = set(used_ids)
    if candidate not in existing:
        return candidate
    suffix = 1
    while True:
        alternate = _generate_default_atom_id(element, ordinal + suffix)
        if alternate not in existing:
            return alternate
        suffix += 1


def _normalize_atom_id(raw: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]", "", raw.upper())
    if not cleaned:
        return ""
    if not re.search(r"[A-Z]", cleaned):
        return ""
    return cleaned


def _generate_default_atom_id(element: str, ordinal: int) -> str:
    digits = str(ordinal)
    if len(digits) >= ATOM_ID_MAX_LEN:
        return digits[-ATOM_ID_MAX_LEN:]
    prefix_len = ATOM_ID_MAX_LEN - len(digits)
    prefix = element[:prefix_len] or element[0]
    return f"{prefix}{digits}"


def _ensure_length(atom_id: str, element: str, ordinal: int) -> str:
    trimmed = atom_id[:ATOM_ID_MAX_LEN]
    if trimmed:
        return trimmed
    return _generate_default_atom_id(element, ordinal)


def convert_directory(input_dir: Path, output_dir: Path) -> None:
    if not input_dir.is_dir():
        raise SystemExit(f"Input directory {input_dir} does not exist")
    output_dir.mkdir(parents=True, exist_ok=True)
    mol_files = sorted(input_dir.glob("*.mol"))
    if not mol_files:
        raise SystemExit(f"No .mol files found in {input_dir}")
    for mol_file in mol_files:
        comp_id, cif_text = mol_to_cif_text(mol_file)
        output_path = output_dir / f"{comp_id}.cif"
        output_path.write_text(cif_text)
        print(f"Wrote {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, help="Directory with MOL files")
    parser.add_argument("--output-dir", type=Path, default=Path("cif"), help="Where CIF files will be written")
    args = parser.parse_args()
    convert_directory(args.input_dir, args.output_dir)


if __name__ == "__main__":
    main()
