"""
LPMO Pipeline: Step 7b - Convert normalized mmCIF to PDB.
Responsibility: Prepare a PoseBusters-compatible PDB artifact from normalized.cif.

Backend selection (adapted from PoseBench; see ATTRIBUTION.md):
    1. PDBFixer/OpenMM (preferred) — handles AF3 mmCIF correctly; preserves
       chain identity and HETATM records.  Biopython does NOT handle AF3 mmCIF
       correctly.  No hydrogens are added here — that is step 7's job.
    2. gemmi (fallback) — used if pdbfixer/openmm are not importable.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Optional, Tuple

from lpmo_pipeline.io.gemmi_compat import gemmi
from lpmo_pipeline.utils.logging import FailureLog, StructuredLogger


@dataclass
class CIFToPDBReport:
    """Summary for one mmCIF to PDB conversion."""

    backend: str
    input_cif: str
    output_pdb: str
    atom_count: int
    model_count: int
    chain_count: int
    has_conect_records: bool
    has_link_records: bool
    stripped_metal_atom_count: int = 0
    stripped_metal_elements: list[str] = field(default_factory=list)
    backend_fallback_reason: str = ""


@dataclass(frozen=True)
class PoseBustersPDBSanitization:
    """Summary of records removed from a PoseBusters-specific PDB export."""

    stripped_metal_atom_count: int = 0
    stripped_metal_elements: tuple[str, ...] = ()


_POSEBUSTERS_INCOMPATIBLE_METALS: frozenset[str] = frozenset(
    {
        "LI", "NA", "K", "RB", "CS",
        "MG", "CA", "SR", "BA",
        "AL", "MN", "FE", "CO", "NI", "CU", "ZN", "CD", "HG",
    }
)

_GLYCAN_LINK_RESNAMES: frozenset[str] = frozenset(["NAG", "BGC", "GLC"])

_GLYCAN_CORE_BONDS: tuple[tuple[str, str], ...] = (
    ("C1", "C2"),
    ("C1", "O5"),
    ("C2", "C3"),
    ("C3", "C4"),
    ("C3", "O3"),
    ("C4", "C5"),
    ("C4", "O4"),
    ("C5", "C6"),
    ("C5", "O5"),
    ("C6", "O6"),
)

_GLYCAN_OPTIONAL_BONDS: tuple[tuple[str, str], ...] = (
    ("C2", "N2"),
    ("N2", "C7"),
    ("C7", "C8"),
    ("C7", "O7"),
    ("C2", "O2"),
)


class CIFToPDBRunner:
    """Convert normalized mmCIF to a PoseBusters-ready PDB artifact."""

    def __init__(
        self,
        input_cif: Path,
        output_dir: Path,
        *,
        strip_metals_for_posebusters: bool = True,
    ):
        self.input_cif = Path(input_cif)
        self.output_dir = Path(output_dir)
        self.strip_metals_for_posebusters = strip_metals_for_posebusters
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.logger = StructuredLogger("cif_to_pdb", self.output_dir)
        self.failures = FailureLog(self.output_dir / "cif_to_pdb_failures.json")

    def run(self) -> Tuple[bool, Optional[Path]]:
        """Convert ``input_cif`` into ``for_posebusters.pdb``.

        Tries PDBFixer/OpenMM first (preferred for AF3 mmCIF), then falls back
        to gemmi if PDBFixer is not importable.

        Returns:
            Tuple of success flag and output PDB path.
        """
        start = datetime.utcnow()
        self.logger.log_step_start("cif_to_pdb", {"input": str(self.input_cif)})

        try:
            report, output_path = self._convert_auto()
            if self.strip_metals_for_posebusters:
                sanitization = _write_posebusters_ready_pdb(output_path, output_path)
                summary = _summarize_pdb_file(output_path)
                report.atom_count = summary.atom_count
                report.chain_count = summary.chain_count
                report.has_conect_records = summary.has_conect_records
                report.has_link_records = summary.has_link_records
                report.stripped_metal_atom_count = sanitization.stripped_metal_atom_count
                report.stripped_metal_elements = list(sanitization.stripped_metal_elements)
        except Exception as exc:
            self.logger.log_failure("cif_to_pdb_exception", {"error": str(exc)})
            self.failures.record("cif_to_pdb", None, "exception", {"error": str(exc)})
            self.failures.write()
            elapsed = (datetime.utcnow() - start).total_seconds()
            self.logger.log_step_end("cif_to_pdb", "failure", {}, elapsed)
            self.logger.close()
            return False, None

        report_path = self.output_dir / "cif_to_pdb_report.json"
        report_path.write_text(json.dumps(asdict(report), indent=2))
        self.logger.log_artifact("for_posebusters_pdb", output_path, "PDB export for PoseBusters")
        self.logger.log_artifact("cif_to_pdb_report", report_path, "mmCIF to PDB conversion report")

        elapsed = (datetime.utcnow() - start).total_seconds()
        self.logger.log_step_end(
            "cif_to_pdb",
            "success",
            {
                "backend": report.backend,
                "atom_count": report.atom_count,
                "chain_count": report.chain_count,
                "has_conect_records": report.has_conect_records,
                "has_link_records": report.has_link_records,
            },
            elapsed,
        )
        self.failures.write()
        self.logger.close()
        return True, output_path

    def _convert_auto(self) -> tuple[CIFToPDBReport, Path]:
        """Try PDBFixer first; fall back to gemmi."""
        try:
            return self._convert_with_pdbfixer()
        except ImportError as exc:
            return self._convert_with_gemmi(fallback_reason=f"pdbfixer_not_importable: {exc}")
        except Exception as exc:
            return self._convert_with_gemmi(fallback_reason=f"pdbfixer_failed: {exc}")

    def _convert_with_pdbfixer(self) -> tuple[CIFToPDBReport, Path]:
        """Convert using PDBFixer/OpenMM (preferred for AF3 mmCIF).

        No hydrogens are added — that is step 7's responsibility.
        Chain IDs are preserved via keepIds=True.

        Adapted from PoseBench (MIT); see ATTRIBUTION.md.
        """
        from pdbfixer import PDBFixer  # type: ignore[import]
        from openmm.app import PDBFile  # type: ignore[import]

        fixer = PDBFixer(filename=str(self.input_cif))
        _add_name_based_glycosidic_bonds(fixer.topology, fixer.positions)
        output_path = self.output_dir / "for_posebusters.pdb"

        with open(output_path, "w") as out_f:
            PDBFile.writeFile(fixer.topology, fixer.positions, out_f, keepIds=True)

        _rewrite_conect_from_topology(output_path, fixer.topology)

        pdb_text = output_path.read_text()
        atom_lines = [
            ln for ln in pdb_text.splitlines() if ln.startswith(("ATOM  ", "HETATM"))
        ]
        chain_ids = {ln[21] for ln in atom_lines if len(ln) > 21}

        report = CIFToPDBReport(
            backend="pdbfixer",
            input_cif=str(self.input_cif),
            output_pdb=str(output_path),
            atom_count=len(atom_lines),
            model_count=1,
            chain_count=len(chain_ids),
            has_conect_records="CONECT" in pdb_text,
            has_link_records="LINK" in pdb_text,
            backend_fallback_reason="",
        )
        return report, output_path

    def _convert_with_gemmi(self, fallback_reason: str = "") -> tuple[CIFToPDBReport, Path]:
        if gemmi is None:
            raise ImportError("Gemmi is required for cif_to_pdb conversion")

        structure = gemmi.read_structure(str(self.input_cif))
        output_path = self.output_dir / "for_posebusters.pdb"

        model = structure[0] if len(structure) > 0 else None
        atom_count = model.count_atom_sites() if model is not None else 0

        if atom_count == 0:
            return self._convert_from_atom_site_block(output_path, fallback_reason=fallback_reason)

        options = gemmi.PdbWriteOptions()
        options.minimal_file = False
        options.atom_records = True
        options.seqres_records = False
        options.ssbond_records = False
        options.link_records = True
        options.conect_records = True
        options.ter_records = True
        options.end_record = True
        options.preserve_serial = True

        structure.write_pdb(str(output_path), options)

        pdb_text = output_path.read_text()
        report = CIFToPDBReport(
            backend="gemmi",
            input_cif=str(self.input_cif),
            output_pdb=str(output_path),
            atom_count=atom_count,
            model_count=len(structure),
            chain_count=len(model) if model is not None else 0,
            has_conect_records="CONECT" in pdb_text,
            has_link_records="LINK" in pdb_text,
            backend_fallback_reason=fallback_reason,
        )
        return report, output_path

    def _convert_from_atom_site_block(self, output_path: Path, fallback_reason: str = "") -> tuple[CIFToPDBReport, Path]:
        """Fallback conversion for minimal mmCIF files with only an _atom_site loop."""
        if gemmi is None:
            raise ImportError("Gemmi is required for cif_to_pdb conversion")

        doc = gemmi.cif.read(str(self.input_cif))
        block = doc.sole_block()
        table = block.find(
            "_atom_site.",
            [
                "group_PDB",
                "id",
                "type_symbol",
                "label_atom_id",
                "label_comp_id",
                "label_asym_id",
                "label_seq_id",
                "Cartn_x",
                "Cartn_y",
                "Cartn_z",
                "B_iso_or_equiv",
            ],
        )
        if len(table) == 0:
            raise ValueError("No atom records found in _atom_site loop")

        lines = [
            f"HEADER    {'MMCIF TO PDB CONVERSION':<40}{block.name[:12]:>12}",
            "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1",
        ]

        current_chain = None
        current_resseq = None
        atom_count = 0
        chain_ids: set[str] = set()
        for row in table:
            record_name = row[0] if row[0] in {"ATOM", "HETATM"} else "HETATM"
            serial = int(row[1])
            element = row[2].strip() or row[3].strip()[:1]
            atom_name = row[3].strip()
            resname = row[4].strip()
            chain_id = (row[5].strip() or "A")[:1]
            resseq = int(float(row[6]))
            x = float(row[7])
            y = float(row[8])
            z = float(row[9])
            b_iso = float(row[10])

            if current_chain is not None and (chain_id != current_chain or resseq != current_resseq):
                lines.append(f"TER   {serial:>5}      {resname:>3} {current_chain:1}{current_resseq:>4}")

            lines.append(
                self._format_pdb_atom_line(
                    record_name=record_name,
                    serial=serial,
                    atom_name=atom_name,
                    resname=resname,
                    chain_id=chain_id,
                    resseq=resseq,
                    x=x,
                    y=y,
                    z=z,
                    b_iso=b_iso,
                    element=element,
                )
            )
            current_chain = chain_id
            current_resseq = resseq
            chain_ids.add(chain_id)
            atom_count += 1

        lines.append("END")
        output_path.write_text("\n".join(lines) + "\n")

        report = CIFToPDBReport(
            backend="gemmi_atom_site_fallback",
            input_cif=str(self.input_cif),
            output_pdb=str(output_path),
            atom_count=atom_count,
            model_count=1,
            chain_count=len(chain_ids),
            has_conect_records=False,
            has_link_records=False,
            backend_fallback_reason=fallback_reason,
        )
        return report, output_path

    @staticmethod
    def _format_pdb_atom_line(
        *,
        record_name: str,
        serial: int,
        atom_name: str,
        resname: str,
        chain_id: str,
        resseq: int,
        x: float,
        y: float,
        z: float,
        b_iso: float,
        element: str,
    ) -> str:
        atom_field = atom_name[:4].rjust(4)
        element_field = element[:2].rjust(2)
        return (
            f"{record_name:<6}{serial:>5} {atom_field} {'':1}{resname:>3} {chain_id:1}"
            f"{resseq:>4}{'':1}   {x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{b_iso:>6.2f}"
            f"          {element_field}"
        )


@dataclass(frozen=True)
class _PDBFileSummary:
    atom_count: int
    chain_count: int
    has_conect_records: bool
    has_link_records: bool


def _summarize_pdb_file(pdb_path: Path) -> _PDBFileSummary:
    pdb_text = pdb_path.read_text()
    atom_lines = [
        line for line in pdb_text.splitlines() if line.startswith(("ATOM  ", "HETATM"))
    ]
    chain_ids = {line[21] for line in atom_lines if len(line) > 21 and line[21].strip()}
    return _PDBFileSummary(
        atom_count=len(atom_lines),
        chain_count=len(chain_ids),
        has_conect_records="CONECT" in pdb_text,
        has_link_records="LINK" in pdb_text,
    )


def _parse_pdb_xyz(line: str) -> tuple[float, float, float] | None:
    try:
        return (
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54]),
        )
    except ValueError:
        return None


def _distance_squared(
    point_a: tuple[float, float, float],
    point_b: tuple[float, float, float],
) -> float:
    dx = point_a[0] - point_b[0]
    dy = point_a[1] - point_b[1]
    dz = point_a[2] - point_b[2]
    return dx * dx + dy * dy + dz * dz


def _candidate_distance_by_serials(
    serial_a: int | None,
    serial_b: int | None,
    atom_xyz: dict[int, tuple[float, float, float]],
) -> float | None:
    if serial_a is None or serial_b is None:
        return None
    point_a = atom_xyz.get(serial_a)
    point_b = atom_xyz.get(serial_b)
    if point_a is None or point_b is None:
        return None
    return _distance_squared(point_a, point_b)


def _choose_glycosidic_serial_pair(
    left_atoms: dict[str, int],
    right_atoms: dict[str, int],
    atom_xyz: dict[int, tuple[float, float, float]],
) -> tuple[int | None, int | None]:
    forward_pair = (left_atoms.get("O4"), right_atoms.get("C1"))
    reverse_pair = (left_atoms.get("C1"), right_atoms.get("O4"))

    forward_distance = _candidate_distance_by_serials(*forward_pair, atom_xyz)
    reverse_distance = _candidate_distance_by_serials(*reverse_pair, atom_xyz)

    if forward_distance is not None and reverse_distance is not None:
        return forward_pair if forward_distance <= reverse_distance else reverse_pair
    if forward_distance is not None:
        return forward_pair
    if reverse_distance is not None:
        return reverse_pair
    if forward_pair[0] is not None and forward_pair[1] is not None:
        return forward_pair
    if reverse_pair[0] is not None and reverse_pair[1] is not None:
        return reverse_pair
    return None, None


def _position_to_xyz(position: Any) -> tuple[float, float, float] | None:
    try:
        from openmm import unit  # type: ignore[import]

        values = position.value_in_unit(unit.angstrom)
        return (float(values[0]), float(values[1]), float(values[2]))
    except Exception:
        pass

    if hasattr(position, "x") and hasattr(position, "y") and hasattr(position, "z"):
        try:
            return (float(position.x), float(position.y), float(position.z))
        except Exception:
            return None

    try:
        return (float(position[0]), float(position[1]), float(position[2]))
    except Exception:
        return None


def _build_topology_atom_xyz(topology, positions) -> dict[int, tuple[float, float, float]]:
    if positions is None:
        return {}

    atom_xyz: dict[int, tuple[float, float, float]] = {}
    for atom, position in zip(topology.atoms(), positions):
        xyz = _position_to_xyz(position)
        if xyz is None:
            continue
        atom_xyz[atom.index] = xyz
    return atom_xyz


def _choose_glycosidic_topology_atoms(left_residue, right_residue, atom_xyz: dict[int, tuple[float, float, float]]):
    left_atoms = {atom.name: atom for atom in left_residue.atoms() if atom.name in {"C1", "O4"}}
    right_atoms = {atom.name: atom for atom in right_residue.atoms() if atom.name in {"C1", "O4"}}

    forward_pair = (left_atoms.get("O4"), right_atoms.get("C1"))
    reverse_pair = (left_atoms.get("C1"), right_atoms.get("O4"))

    forward_distance = _candidate_distance_by_serials(
        None if forward_pair[0] is None else forward_pair[0].index,
        None if forward_pair[1] is None else forward_pair[1].index,
        atom_xyz,
    )
    reverse_distance = _candidate_distance_by_serials(
        None if reverse_pair[0] is None else reverse_pair[0].index,
        None if reverse_pair[1] is None else reverse_pair[1].index,
        atom_xyz,
    )

    if forward_distance is not None and reverse_distance is not None:
        return forward_pair if forward_distance <= reverse_distance else reverse_pair
    if forward_distance is not None:
        return forward_pair
    if reverse_distance is not None:
        return reverse_pair
    if forward_pair[0] is not None and forward_pair[1] is not None:
        return forward_pair
    if reverse_pair[0] is not None and reverse_pair[1] is not None:
        return reverse_pair
    return None, None


def _rewrite_conect_from_topology(pdb_path: Path, topology) -> None:
    """Rewrite PDB CONECT records from full topology bonds.

    PDBFixer/OpenMM can emit sparse CONECT lines for glycan-heavy AF3 exports.
    For downstream ligand-aware tools like PoseBusters, replace the emitted
    connectivity with the full bond graph and patch in known glycan bonds by
    residue/atom naming where the topology is still incomplete.
    """
    lines = pdb_path.read_text().splitlines(keepends=True)

    atom_serials: list[int] = []
    kept_lines: list[str] = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            serial_txt = line[6:11].strip()
            try:
                atom_serials.append(int(serial_txt))
            except ValueError:
                atom_serials.append(len(atom_serials) + 1)
            kept_lines.append(line)
            continue

        if line.startswith("CONECT"):
            continue
        kept_lines.append(line)

    n_top_atoms = sum(1 for _ in topology.atoms())
    if n_top_atoms != len(atom_serials):
        return

    adjacency: dict[int, set[int]] = {serial: set() for serial in atom_serials}
    for atom_a, atom_b in topology.bonds():
        serial_a = atom_serials[atom_a.index]
        serial_b = atom_serials[atom_b.index]
        adjacency[serial_a].add(serial_b)
        adjacency[serial_b].add(serial_a)

    _augment_glycan_bonds_from_names(kept_lines, adjacency)

    conect_lines: list[str] = []
    for src in sorted(adjacency):
        neighbors = sorted(adjacency[src])
        if not neighbors:
            continue
        for start in range(0, len(neighbors), 4):
            chunk = neighbors[start : start + 4]
            conect_lines.append(
                f"CONECT{src:5d}" + "".join(f"{serial:5d}" for serial in chunk) + "\n"
            )

    if kept_lines and kept_lines[-1].startswith("END"):
        end_line = kept_lines.pop()
        kept_lines.extend(conect_lines)
        kept_lines.append(end_line)
    else:
        kept_lines.extend(conect_lines)

    pdb_path.write_text("".join(kept_lines))


def _augment_glycan_bonds_from_names(
    pdb_lines: list[str],
    adjacency: dict[int, set[int]],
) -> None:
    """Add missing glycan bonds using explicit atom-name templates."""
    atoms: list[dict[str, object]] = []
    atom_xyz: dict[int, tuple[float, float, float]] = {}
    for line in pdb_lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        serial_txt = line[6:11].strip()
        if not serial_txt:
            continue
        try:
            serial = int(serial_txt)
        except ValueError:
            continue

        xyz = _parse_pdb_xyz(line)
        if xyz is not None:
            atom_xyz[serial] = xyz

        resname = line[17:20].strip()
        if resname not in _GLYCAN_LINK_RESNAMES:
            continue

        atom_name = line[12:16].strip()
        chain_id = line[21].strip() or "?"
        resseq_txt = line[22:26].strip()
        try:
            resseq = int(resseq_txt)
        except ValueError:
            continue

        atoms.append(
            {
                "serial": serial,
                "atom_name": atom_name,
                "resname": resname,
                "chain": chain_id,
                "resseq": resseq,
            }
        )

    if not atoms:
        return

    by_residue: dict[tuple[str, int, str], dict[str, int]] = {}
    for atom in atoms:
        key = (str(atom["chain"]), int(atom["resseq"]), str(atom["resname"]))
        by_residue.setdefault(key, {})[str(atom["atom_name"])] = int(atom["serial"])

    for atom_map in by_residue.values():
        serial_c1 = atom_map.get("C1")
        serial_o4 = atom_map.get("O4")
        if serial_c1 is not None and serial_o4 is not None:
            _remove_adjacency_edge(adjacency, serial_c1, serial_o4)

    for atom_map in by_residue.values():
        for atom_a_name, atom_b_name in (_GLYCAN_CORE_BONDS + _GLYCAN_OPTIONAL_BONDS):
            serial_a = atom_map.get(atom_a_name)
            serial_b = atom_map.get(atom_b_name)
            if serial_a is None or serial_b is None:
                continue
            _add_adjacency_edge(adjacency, serial_a, serial_b)

    by_chain: dict[str, list[tuple[str, int, str]]] = {}
    for key in by_residue.keys():
        chain_id, resseq, _resname = key
        by_chain.setdefault(chain_id, []).append(key)

    for residue_keys in by_chain.values():
        residue_keys.sort(key=lambda key: key[1])
        for index in range(len(residue_keys) - 1):
            left = residue_keys[index]
            right = residue_keys[index + 1]
            if right[1] != left[1] + 1:
                continue

            bond_a, bond_b = _choose_glycosidic_serial_pair(
                by_residue[left],
                by_residue[right],
                atom_xyz,
            )
            if bond_a is None or bond_b is None:
                continue
            _add_adjacency_edge(adjacency, bond_a, bond_b)


def _add_adjacency_edge(adjacency: dict[int, set[int]], serial_a: int, serial_b: int) -> None:
    if serial_a == serial_b:
        return
    if serial_a not in adjacency or serial_b not in adjacency:
        return
    adjacency[serial_a].add(serial_b)
    adjacency[serial_b].add(serial_a)


def _remove_adjacency_edge(adjacency: dict[int, set[int]], serial_a: int, serial_b: int) -> None:
    if serial_a not in adjacency or serial_b not in adjacency:
        return
    adjacency[serial_a].discard(serial_b)
    adjacency[serial_b].discard(serial_a)


def _residue_seq_id(residue) -> int | None:
    try:
        return int(str(residue.id).strip())
    except (TypeError, ValueError):
        return None


def _add_name_based_glycosidic_bonds(topology, positions=None) -> int:
    """Add glycosidic bonds for adjacent NAG/BGC/GLC residues.

    Residue numbering is not guaranteed to follow the physical left-to-right
    order in crystal-derived subsets, so choose the shorter cross-residue
    C1/O4 pair instead of assuming O4(i)-C1(i+1) from sequence ids alone.
    """
    existing: set[tuple[int, int]] = set()
    for atom_a, atom_b in topology.bonds():
        idx_a, idx_b = atom_a.index, atom_b.index
        if idx_a > idx_b:
            idx_a, idx_b = idx_b, idx_a
        existing.add((idx_a, idx_b))

    atom_xyz = _build_topology_atom_xyz(topology, positions)

    by_chain: dict[str, list] = {}
    for chain in topology.chains():
        residues = []
        for residue in chain.residues():
            if residue.name not in _GLYCAN_LINK_RESNAMES:
                continue
            seq_id = _residue_seq_id(residue)
            if seq_id is None:
                continue
            residues.append(residue)
        residues.sort(key=lambda residue: _residue_seq_id(residue) or 10**9)
        by_chain[chain.id] = residues

    added = 0
    for residues in by_chain.values():
        for index in range(len(residues) - 1):
            left = residues[index]
            right = residues[index + 1]
            left_seq = _residue_seq_id(left)
            right_seq = _residue_seq_id(right)
            if left_seq is None or right_seq is None or right_seq != left_seq + 1:
                continue

            bond_a, bond_b = _choose_glycosidic_topology_atoms(left, right, atom_xyz)
            if bond_a is None or bond_b is None:
                continue

            idx_a, idx_b = bond_a.index, bond_b.index
            if idx_a > idx_b:
                idx_a, idx_b = idx_b, idx_a
            if (idx_a, idx_b) in existing:
                continue

            topology.addBond(bond_a, bond_b)
            existing.add((idx_a, idx_b))
            added += 1

    return added


def _write_posebusters_ready_pdb(input_pdb: Path, output_pdb: Path) -> PoseBustersPDBSanitization:
    """Remove metal records that cause PoseBusters to crash, preserving other bonds."""
    lines = input_pdb.read_text().splitlines()

    stripped_atom_serials: set[int] = set()
    stripped_residues: set[tuple[str, str, str, str]] = set()
    stripped_elements: set[str] = set()
    kept_lines: list[str] = []

    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and _is_posebusters_incompatible_metal(line):
            serial = _parse_atom_serial(line)
            if serial is not None:
                stripped_atom_serials.add(serial)
            residue_key = _parse_residue_key(line)
            if residue_key is not None:
                stripped_residues.add(residue_key)
            element = _parse_element_symbol(line)
            if element:
                stripped_elements.add(element)
            continue
        kept_lines.append(line)

    filtered_lines: list[str] = []
    for line in kept_lines:
        if line.startswith("TER"):
            residue_key = _parse_residue_key(line)
            if residue_key is not None and residue_key in stripped_residues:
                continue
        elif line.startswith("LINK"):
            first_residue, second_residue = _parse_link_residue_keys(line)
            if first_residue in stripped_residues or second_residue in stripped_residues:
                continue
        elif line.startswith("CONECT"):
            serials = _parse_conect_serials(line)
            kept_serials = [serial for serial in serials if serial not in stripped_atom_serials]
            if len(kept_serials) < 2:
                continue
            line = "CONECT" + "".join(f"{serial:>5}" for serial in kept_serials)
        filtered_lines.append(line)

    output_pdb.write_text("\n".join(filtered_lines) + "\n")

    return PoseBustersPDBSanitization(
        stripped_metal_atom_count=len(stripped_atom_serials),
        stripped_metal_elements=tuple(sorted(stripped_elements)),
    )


def _is_posebusters_incompatible_metal(line: str) -> bool:
    return _parse_element_symbol(line) in _POSEBUSTERS_INCOMPATIBLE_METALS


def _parse_atom_serial(line: str) -> int | None:
    try:
        return int(line[6:11].strip())
    except ValueError:
        return None


def _parse_element_symbol(line: str) -> str:
    element = line[76:78].strip().upper()
    if element:
        return element

    atom_name = "".join(ch for ch in line[12:16].strip() if ch.isalpha())
    if not atom_name:
        return ""
    if len(atom_name) >= 2:
        first_two = atom_name[:2].upper()
        if first_two in _POSEBUSTERS_INCOMPATIBLE_METALS:
            return first_two
    return atom_name[:1].upper()


def _parse_residue_key(line: str) -> tuple[str, str, str, str] | None:
    if line.startswith("TER"):
        tokens = line.split()
        if len(tokens) >= 5:
            return (tokens[2], tokens[3], tokens[4], "")
        return None
    if len(line) < 27:
        return None
    return (
        line[17:20].strip(),
        line[21].strip(),
        line[22:26].strip(),
        line[26].strip(),
    )


def _parse_link_residue_keys(
    line: str,
) -> tuple[tuple[str, str, str, str] | None, tuple[str, str, str, str] | None]:
    if len(line) < 57:
        return None, None
    first = (
        line[17:20].strip(),
        line[21].strip(),
        line[22:26].strip(),
        line[26].strip(),
    )
    second = (
        line[47:50].strip(),
        line[51].strip(),
        line[52:56].strip(),
        line[56].strip(),
    )
    return first, second


def _parse_conect_serials(line: str) -> list[int]:
    serials: list[int] = []
    for token in line[6:].split():
        try:
            serials.append(int(token))
        except ValueError:
            continue
    return serials


def convert_cif_to_pdb(
    input_cif: Path,
    output_dir: Path,
    *,
    strip_metals_for_posebusters: bool = True,
) -> Tuple[bool, Optional[Path]]:
    """Convenience wrapper for step 7b conversion."""
    return CIFToPDBRunner(
        input_cif=input_cif,
        output_dir=output_dir,
        strip_metals_for_posebusters=strip_metals_for_posebusters,
    ).run()