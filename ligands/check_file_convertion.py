"""
Generic conversion checker for folders like:
  base/
    pdb/
    mol/
    cif/

Matches triplets by oligomer length n, regardless of family (amylose/cellulose/nag_chitin/...).

Handles broken MOL V2000 spacing by fixed-width parsing:
- Counts line: first 3 chars = natoms, next 3 = nbonds
- Bond line: first 3 = a1, next 3 = a2, next 3 = order

Usage:
  python3 check_file_convertion.py /path/to/base --tol 1e-5
"""

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

Vec3 = Tuple[float, float, float]


@dataclass
class MolLike:
    name: str
    elements: List[str]
    coords: List[Vec3]
    bonds: Set[Tuple[int, int, int]]  # (i,j,order), i<j, 0-based


# -------------------------
# Helpers
# -------------------------

def dist(a: Vec3, b: Vec3) -> float:
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

def normalize_element(e: str) -> str:
    e = (e or "").strip()
    if not e:
        return e
    return e[0].upper() + (e[1:].lower() if len(e) > 1 else "")

def trailing_int(s: str) -> Optional[int]:
    m = re.search(r"(\d+)\s*$", s)
    return int(m.group(1)) if m else None


# -------------------------
# Robust MOL V2000 parsing (fixed width)
# -------------------------

def _fw_int(s: str) -> Optional[int]:
    s = s.strip()
    if not s:
        return None
    try:
        return int(s)
    except ValueError:
        return None

def parse_mol_v2000(path: Path) -> MolLike:
    lines = path.read_text(errors="replace").splitlines()
    if len(lines) < 4:
        raise ValueError(f"{path}: too few lines for MOL")

    counts = lines[3].rstrip("\n")

    # V2000 counts line is fixed width:
    # columns 1-3: natoms, 4-6: nbonds (1-indexed in spec; 0-index slicing below)
    # This survives missing spaces like "171178".
    natoms = _fw_int(counts[0:3])
    nbonds = _fw_int(counts[3:6])

    # Fallback if file is even more mangled
    if natoms is None or nbonds is None:
        nums = re.findall(r"\d+", counts)
        if len(nums) >= 2:
            natoms, nbonds = int(nums[0]), int(nums[1])
        else:
            raise ValueError(f"{path}: cannot parse counts line: {counts!r}")

    atom_lines = lines[4:4+natoms]
    bond_lines = lines[4+natoms:4+natoms+nbonds]

    if len(atom_lines) != natoms:
        raise ValueError(f"{path}: expected {natoms} atom lines, got {len(atom_lines)}")
    if len(bond_lines) != nbonds:
        raise ValueError(f"{path}: expected {nbonds} bond lines, got {len(bond_lines)}")

    elements: List[str] = []
    coords: List[Vec3] = []

    for idx, ln in enumerate(atom_lines, 1):
        # Atom lines are usually whitespace-separated; fixed-width exists but split is OK here.
        parts = ln.split()
        if len(parts) < 4:
            raise ValueError(f"{path}: bad atom line {idx}: {ln!r}")
        x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
        elem = normalize_element(parts[3])
        elements.append(elem)
        coords.append((x, y, z))

    bonds: Set[Tuple[int, int, int]] = set()
    for idx, ln in enumerate(bond_lines, 1):
        # Bond lines are fixed width (3 chars each). This survives missing spaces: "91104  1 ..."
        a1 = _fw_int(ln[0:3])
        a2 = _fw_int(ln[3:6])
        order = _fw_int(ln[6:9])

        # Fallback to split if needed
        if a1 is None or a2 is None or order is None:
            parts = ln.split()
            if len(parts) >= 3:
                a1, a2, order = int(parts[0]), int(parts[1]), int(parts[2])
            else:
                raise ValueError(f"{path}: bad bond line {idx}: {ln!r}")

        i, j = a1-1, a2-1
        if i == j:
            continue
        if i > j:
            i, j = j, i
        bonds.add((i, j, int(order)))

    return MolLike(name=path.stem, elements=elements, coords=coords, bonds=bonds)


# -------------------------
# PDB parsing
# -------------------------

def parse_pdb(path: Path) -> MolLike:
    elements: List[str] = []
    coords: List[Vec3] = []
    serials: List[Optional[int]] = []
    lines = path.read_text(errors="replace").splitlines()

    for ln in lines:
        if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
            continue
        if len(ln) < 54:
            continue
        x = float(ln[30:38])
        y = float(ln[38:46])
        z = float(ln[46:54])
        elem = ln[76:78].strip() if len(ln) >= 78 else ""
        if not elem:
            aname = ln[12:16].strip()
            elem = re.sub(r"[^A-Za-z]", "", aname)[:2]
        elem = normalize_element(elem)
        elements.append(elem)
        coords.append((x, y, z))
        try:
            serials.append(int(ln[6:11]))
        except Exception:
            serials.append(None)

    serial_to_idx: Dict[int, int] = {}
    for i, s in enumerate(serials):
        if s is not None:
            serial_to_idx[s] = i

    bonds: Set[Tuple[int, int, int]] = set()
    for ln in lines:
        if not ln.startswith("CONECT"):
            continue
        nums = re.findall(r"\d+", ln)
        if len(nums) < 2:
            continue
        src = int(nums[0])
        for tgt_s in nums[1:]:
            tgt = int(tgt_s)
            if src in serial_to_idx and tgt in serial_to_idx:
                i, j = serial_to_idx[src], serial_to_idx[tgt]
                if i == j:
                    continue
                if i > j:
                    i, j = j, i
                bonds.add((i, j, 1))  # order unknown in PDB
    return MolLike(name=path.stem, elements=elements, coords=coords, bonds=bonds)


# -------------------------
# CIF chem_comp parsing (minimal)
# -------------------------

def parse_cif_chemcomp(path: Path, expected_comp_id: Optional[str] = None) -> Tuple[MolLike, str]:
    """
    Returns (MolLike, comp_id_detected).
    If expected_comp_id is None, detect it and use that.
    """
    text = path.read_text(errors="replace").splitlines()

    def parse_loop(start_idx: int):
        tags: List[str] = []
        rows: List[List[str]] = []
        i = start_idx
        assert text[i].strip() == "loop_"
        i += 1
        while i < len(text) and text[i].strip().startswith("_"):
            tags.append(text[i].strip())
            i += 1
        while i < len(text):
            ln = text[i].strip()
            if not ln:
                i += 1
                continue
            if ln.startswith("#") or ln == "loop_" or ln.startswith("_"):
                break
            rows.append(ln.split())
            i += 1
        return tags, rows, i

    atom_tags = None
    atom_rows: List[List[str]] = []
    bond_tags = None
    bond_rows: List[List[str]] = []

    i = 0
    while i < len(text):
        if text[i].strip() == "loop_":
            tags, rows, nxt = parse_loop(i)
            tagset = set(tags)
            if "_chem_comp_atom.atom_id" in tagset and "_chem_comp_atom.pdbx_ordinal" in tagset:
                atom_tags, atom_rows = tags, rows
            if "_chem_comp_bond.atom_id_1" in tagset and "_chem_comp_bond.atom_id_2" in tagset:
                bond_tags, bond_rows = tags, rows
            i = nxt
        else:
            i += 1

    if atom_tags is None:
        raise ValueError(f"{path}: did not find _chem_comp_atom loop")
    if bond_tags is None:
        raise ValueError(f"{path}: did not find _chem_comp_bond loop")

    at = {t: k for k, t in enumerate(atom_tags)}
    bt = {t: k for k, t in enumerate(bond_tags)}

    comp_col_a = at.get("_chem_comp_atom.comp_id")
    comp_col_b = bt.get("_chem_comp_bond.comp_id")
    if comp_col_a is None or comp_col_b is None:
        raise ValueError(f"{path}: missing comp_id in loops")

    # Detect comp_id if not supplied
    comp_ids = set()
    for r in atom_rows:
        if len(r) >= len(atom_tags):
            comp_ids.add(r[comp_col_a])
    if not comp_ids:
        raise ValueError(f"{path}: no comp_id values found in atom loop")

    comp_id = expected_comp_id if expected_comp_id else sorted(comp_ids)[0]

    atoms = []
    for r in atom_rows:
        if len(r) < len(atom_tags):
            continue
        if r[comp_col_a] != comp_id:
            continue
        atom_id = r[at["_chem_comp_atom.atom_id"]]
        elem = normalize_element(r[at["_chem_comp_atom.type_symbol"]])
        x = float(r[at["_chem_comp_atom.pdbx_model_Cartn_x_ideal"]])
        y = float(r[at["_chem_comp_atom.pdbx_model_Cartn_y_ideal"]])
        z = float(r[at["_chem_comp_atom.pdbx_model_Cartn_z_ideal"]])
        ordn = int(r[at["_chem_comp_atom.pdbx_ordinal"]])
        atoms.append((ordn, atom_id, elem, (x, y, z)))

    if not atoms:
        raise ValueError(f"{path}: no atoms found for comp_id={comp_id}")

    atoms.sort(key=lambda t: t[0])
    elements = [a[2] for a in atoms]
    coords = [a[3] for a in atoms]
    atom_ids = [a[1] for a in atoms]
    atom_id_to_idx = {aid: i for i, aid in enumerate(atom_ids)}

    bonds: Set[Tuple[int, int, int]] = set()
    for r in bond_rows:
        if len(r) < len(bond_tags):
            continue
        if r[comp_col_b] != comp_id:
            continue
        a1 = r[bt["_chem_comp_bond.atom_id_1"]]
        a2 = r[bt["_chem_comp_bond.atom_id_2"]]
        vo = r[bt["_chem_comp_bond.value_order"]].upper()
        order = {"SING": 1, "DOUB": 2, "TRIP": 3, "AROM": 1}.get(vo, 1)
        if a1 not in atom_id_to_idx or a2 not in atom_id_to_idx:
            continue
        i, j = atom_id_to_idx[a1], atom_id_to_idx[a2]
        if i == j:
            continue
        if i > j:
            i, j = j, i
        bonds.add((i, j, int(order)))

    return MolLike(name=path.stem, elements=elements, coords=coords, bonds=bonds), comp_id


# -------------------------
# Comparisons
# -------------------------

def compare_by_order(a: MolLike, b: MolLike, tol: float) -> Dict[str, object]:
    out = {
        "same_natoms": len(a.elements) == len(b.elements),
        "natoms_a": len(a.elements),
        "natoms_b": len(b.elements),
        "max_coord_diff": None,
        "coord_mismatches": [],
        "element_mismatches": [],
        "bond_diff_count": None,
        "bond_only_in_a": 0,
        "bond_only_in_b": 0,
    }
    if len(a.elements) != len(b.elements):
        return out

    maxd = 0.0
    for i in range(len(a.elements)):
        if a.elements[i] != b.elements[i]:
            out["element_mismatches"].append((i, a.elements[i], b.elements[i]))
        d = dist(a.coords[i], b.coords[i])
        maxd = max(maxd, d)
        if d > tol:
            out["coord_mismatches"].append((i, a.coords[i], b.coords[i], d))
    out["max_coord_diff"] = maxd

    only_a = a.bonds - b.bonds
    only_b = b.bonds - a.bonds
    out["bond_only_in_a"] = len(only_a)
    out["bond_only_in_b"] = len(only_b)
    out["bond_diff_count"] = len(only_a) + len(only_b)
    return out


def greedy_map_by_coords(ref: MolLike, mov: MolLike, tol: float):
    """
    Map mov -> ref by nearest neighbor within tol, element-respecting.
    Requires equal natoms.
    """
    if len(ref.elements) != len(mov.elements):
        return None, float("inf"), [f"natoms differ (ref={len(ref.elements)} mov={len(mov.elements)})"]

    ref_by_elem: Dict[str, List[int]] = {}
    for i, e in enumerate(ref.elements):
        ref_by_elem.setdefault(e, []).append(i)

    pairs: List[Tuple[float, int, int]] = []
    for i, e in enumerate(mov.elements):
        for j in ref_by_elem.get(e, []):
            d = dist(mov.coords[i], ref.coords[j])
            if d <= tol:
                pairs.append((d, i, j))

    pairs.sort(key=lambda t: t[0])
    mov_to_ref = [-1] * len(mov.elements)
    used_ref: Set[int] = set()

    for d, i, j in pairs:
        if mov_to_ref[i] != -1:
            continue
        if j in used_ref:
            continue
        mov_to_ref[i] = j
        used_ref.add(j)

    warnings: List[str] = []
    if any(x == -1 for x in mov_to_ref):
        unmatched = [i for i, x in enumerate(mov_to_ref) if x == -1]
        for i in unmatched[:10]:
            e = mov.elements[i]
            cands = ref_by_elem.get(e, [])
            if not cands:
                warnings.append(f"unmatched mov atom {i} element {e}: no ref atoms with that element")
            else:
                dmin = min(dist(mov.coords[i], ref.coords[j]) for j in cands)
                warnings.append(f"unmatched mov atom {i} element {e}: closest same-element distance {dmin:.6g} > tol")
        return None, float("inf"), warnings

    maxd = 0.0
    for i, j in enumerate(mov_to_ref):
        maxd = max(maxd, dist(mov.coords[i], ref.coords[j]))
    return mov_to_ref, maxd, warnings


def chirality_signatures(m: MolLike) -> Dict[int, int]:
    """
    Simple chirality signature for centers with degree 4 (explicit Hs required).
    Signature = sign( a · (b × d) ) for 3 neighbors chosen deterministically.
    """
    adj = {i: [] for i in range(len(m.elements))}
    for i, j, _o in m.bonds:
        adj[i].append(j)
        adj[j].append(i)

    sig: Dict[int, int] = {}
    for c, neigh in adj.items():
        if len(neigh) != 4:
            continue
        ns = sorted(neigh)
        a, b, d = ns[0], ns[1], ns[2]
        rc = m.coords[c]
        ra, rb, rd = m.coords[a], m.coords[b], m.coords[d]

        ax, ay, az = ra[0]-rc[0], ra[1]-rc[1], ra[2]-rc[2]
        bx, by, bz = rb[0]-rc[0], rb[1]-rc[1], rb[2]-rc[2]
        dx, dy, dz = rd[0]-rc[0], rd[1]-rc[1], rd[2]-rc[2]

        cx = by*dz - bz*dy
        cy = bz*dx - bx*dz
        cz = bx*dy - by*dx
        tp = ax*cx + ay*cy + az*cz
        if abs(tp) < 1e-12:
            continue
        sig[c] = 1 if tp > 0 else -1
    return sig


def compare_chirality(a: MolLike, b: MolLike, mapping_b_to_a=None) -> Dict[str, object]:
    sa = chirality_signatures(a)
    sb_raw = chirality_signatures(b)

    if mapping_b_to_a is None:
        common = sorted(set(sa) & set(sb_raw))
        mism = [(c, sa[c], sb_raw[c]) for c in common if sa[c] != sb_raw[c]]
        return {"centers_a": len(sa), "centers_b": len(sb_raw), "common_centers": len(common), "mismatches": mism}

    sb = {}
    for cb, sign in sb_raw.items():
        ca = mapping_b_to_a[cb]
        sb[ca] = sign
    common = sorted(set(sa) & set(sb))
    mism = [(c, sa[c], sb[c]) for c in common if sa[c] != sb[c]]
    return {"centers_a": len(sa), "centers_b": len(sb_raw), "common_centers": len(common), "mismatches": mism}


# -------------------------
# Matching triplets by n
# -------------------------

def infer_n_from_mol_filename(path: Path) -> Optional[int]:
    # Prefer trailing _<n> or -<n> or ...<n>
    return trailing_int(re.sub(r"\.(mol|sdf)$", "", path.name, flags=re.I))

def infer_n_from_cif_file(path: Path) -> Optional[int]:
    # Try to find _chem_comp.id or data_ line and extract trailing digits
    txt = path.read_text(errors="replace")
    m = re.search(r"^\s*_chem_comp\.id\s+([A-Za-z0-9_\-]+)\s*$", txt, flags=re.M)
    if m:
        return trailing_int(m.group(1))
    m = re.search(r"^\s*data_([A-Za-z0-9_\-]+)\s*$", txt, flags=re.M)
    if m:
        return trailing_int(m.group(1))
    # fallback to filename
    return trailing_int(path.stem)

def infer_monomer_token_from_pdb_filenames(pdb_files: List[Path]) -> Optional[str]:
    """
    Guess a monomer token in filenames like:
      DGlcpNAcb1-4DGlcpNAcb1-4... etc
    We pick the most frequent token matching D....p... across filenames.
    """
    token_counts: Dict[str, int] = {}
    pat = re.compile(r"(D[A-Za-z0-9]+p[A-Za-z0-9]*)")
    for p in pdb_files:
        toks = pat.findall(p.name)
        for t in toks:
            token_counts[t] = token_counts.get(t, 0) + 1
    if not token_counts:
        return None
    # choose most frequent
    return sorted(token_counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]

def infer_n_from_pdb_filename(path: Path, monomer_token: Optional[str]) -> Optional[int]:
    if monomer_token and monomer_token in path.name:
        return path.name.count(monomer_token)
    # fallback: count all D...p... tokens and use the max count token
    pat = re.compile(r"(D[A-Za-z0-9]+p[A-Za-z0-9]*)")
    toks = pat.findall(path.name)
    if not toks:
        return None
    # count each token in name, take max (best guess)
    best = 0
    for t in set(toks):
        best = max(best, path.name.count(t))
    return best if best > 0 else None


def index_files_by_n(base: Path, subdir: str, exts: Tuple[str, ...], infer_n_func) -> Dict[int, Path]:
    d = base / subdir
    out: Dict[int, Path] = {}
    for ext in exts:
        for p in d.glob(f"*{ext}"):
            n = infer_n_func(p)
            if n is None:
                continue
            # If multiple, keep shortest filename (usually the “main” one)
            if n not in out or len(p.name) < len(out[n].name):
                out[n] = p
    return out


# -------------------------
# Main
# -------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base_dir", type=str)
    ap.add_argument("--pdb_subdir", default="pdb")
    ap.add_argument("--mol_subdir", default="mol")
    ap.add_argument("--cif_subdir", default="cif")
    ap.add_argument("--tol", type=float, default=1e-5)
    ap.add_argument("--map_tol_pdb", type=float, default=0.05, help="tolerance for PDB->MOL mapping")
    args = ap.parse_args()

    base = Path(args.base_dir)
    pdb_dir = base / args.pdb_subdir
    mol_dir = base / args.mol_subdir
    cif_dir = base / args.cif_subdir

    if not pdb_dir.exists() or not mol_dir.exists() or not cif_dir.exists():
        raise SystemExit(f"Expected subdirs not found under {base} (have: {pdb_dir}, {mol_dir}, {cif_dir})")

    pdb_files = list(pdb_dir.glob("*.pdb")) + list(pdb_dir.glob("*.ent"))
    monomer_token = infer_monomer_token_from_pdb_filenames(pdb_files)

    print(f"Base: {base}")
    print(f"tol (MOL<->CIF): {args.tol:g} Å")
    print(f"map_tol_pdb (PDB->MOL): {args.map_tol_pdb:g} Å")
    print(f"Detected PDB monomer token: {monomer_token!r}")
    print("")

    # Index all files by n
    mol_by_n = index_files_by_n(base, args.mol_subdir, (".mol", ".sdf"), infer_n_from_mol_filename)
    cif_by_n = index_files_by_n(base, args.cif_subdir, (".cif", ".mmcif"), infer_n_from_cif_file)

    pdb_by_n: Dict[int, Path] = {}
    for p in pdb_files:
        n = infer_n_from_pdb_filename(p, monomer_token)
        if n is None:
            continue
        if n not in pdb_by_n or len(p.name) < len(pdb_by_n[n].name):
            pdb_by_n[n] = p

    all_n = sorted(set(mol_by_n) | set(cif_by_n) | set(pdb_by_n))
    if not all_n:
        raise SystemExit("No files matched any n. Check filenames and folders.")

    any_fail = False

    for n in all_n:
        print("=" * 80)
        print(f"n={n}")
        pdb_path = pdb_by_n.get(n)
        mol_path = mol_by_n.get(n)
        cif_path = cif_by_n.get(n)

        print(f"PDB: {pdb_path.name if pdb_path else 'NOT FOUND'}")
        print(f"MOL: {mol_path.name if mol_path else 'NOT FOUND'}")
        print(f"CIF: {cif_path.name if cif_path else 'NOT FOUND'}")

        if not (pdb_path and mol_path and cif_path):
            print("-> SKIP (missing file)")
            any_fail = True
            continue

        try:
            mol = parse_mol_v2000(mol_path)
            cif, cif_comp_id = parse_cif_chemcomp(cif_path, expected_comp_id=None)
            pdb = parse_pdb(pdb_path)
        except Exception as e:
            print(f"-> ERROR parsing: {e}")
            any_fail = True
            continue

        # 1) MOL vs CIF (order-based)
        print("\n[1] MOL vs CIF (by atom order / CIF pdbx_ordinal)")
        r = compare_by_order(mol, cif, tol=args.tol)
        print(f"  natoms: MOL={r['natoms_a']} CIF={r['natoms_b']} same={r['same_natoms']}")
        if r["same_natoms"]:
            print(f"  max |Δr|: {r['max_coord_diff']:.6g} Å")
            print(f"  element mismatches: {len(r['element_mismatches'])}")
            print(f"  coord mismatches (>tol): {len(r['coord_mismatches'])}")
            print(f"  bond diffs: {r['bond_diff_count']} (only_in_MOL={r['bond_only_in_a']}, only_in_CIF={r['bond_only_in_b']})")
            if r["element_mismatches"]:
                print("    element mismatch examples:", r["element_mismatches"][:5])
            if r["coord_mismatches"]:
                print("    coord mismatch examples:", r["coord_mismatches"][:3])
        else:
            any_fail = True

        ch_mc = compare_chirality(mol, cif, mapping_b_to_a=None)
        print("\n  chirality (MOL vs CIF):")
        print(f"    centers: MOL={ch_mc['centers_a']} CIF={ch_mc['centers_b']} common={ch_mc['common_centers']}")
        print(f"    mismatches: {len(ch_mc['mismatches'])}")
        if ch_mc["mismatches"]:
            print("    mismatch examples:", ch_mc["mismatches"][:5])
            any_fail = True

        # 2) PDB vs MOL mapping (coord-based)
        print("\n[2] PDB vs MOL (coord-based mapping)")
        mapping, maxd_pm, warns = greedy_map_by_coords(ref=mol, mov=pdb, tol=args.map_tol_pdb)
        if mapping is None:
            print("  mapping FAILED")
            for w in warns[:10]:
                print("   -", w)
            any_fail = True
        else:
            print(f"  mapping OK (PDB->MOL). max |Δr|: {maxd_pm:.6g} Å")
            ch_pm = compare_chirality(mol, pdb, mapping_b_to_a=mapping)
            print("  chirality (PDB vs MOL):")
            print(f"    centers: MOL={ch_pm['centers_a']} PDB={ch_pm['centers_b']} common={ch_pm['common_centers']}")
            print(f"    mismatches: {len(ch_pm['mismatches'])}")
            if ch_pm["mismatches"]:
                print("    mismatch examples:", ch_pm["mismatches"][:5])
                any_fail = True

        print("")

    print("=" * 80)
    if any_fail:
        print("DONE: Some checks FAILED or files were missing. See details above.")
        raise SystemExit(2)
    else:
        print("DONE: All checks PASSED.")
        raise SystemExit(0)


if __name__ == "__main__":
    main()
