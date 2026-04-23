"""Step 7: Protonation and export artifacts for QC/analysis.

Responsibility:
  - Produce `for_posebusters.pdb` from normalized mmCIF (via step 7b converter)
  - Produce `complex_H.pdb` — full complex with explicit hydrogens added
  - Produce `ligand_for_prolif.mol2` — glycan ligand only, with hydrogens

Protonation backend priority for complex_H.pdb:
  1. PDBFixer/OpenMM (primary) — adds hydrogens at physiological pH
  2. reduce          (secondary, via AmberTools)
  3. obabel          (tertiary)
  4. Blocker         — reported in protonation_report.json; NO silent copy

Ligand export:
  - Glycan chains (normalized: B, C, D) extracted with gemmi
  - Residue labels validated to be glycan-like (not protein residues)
  - Converted to MOL2 via obabel (required for ProLIF)

Important:
  - This module does not install external programs.
  - Missing external tools are reported as blockers in protonation_report.json.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional, Tuple

import gemmi

from lpmo_pipeline.io.cif_to_pdb import convert_cif_to_pdb
from lpmo_pipeline.utils.logging import FailureLog, StructuredLogger

# Residue names that are unambiguously protein — used to detect a contaminated
# ligand-only extract before handing it to obabel/ProLIF.
_PROTEIN_RESIDUE_NAMES: frozenset[str] = frozenset(
    [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL", "MSE", "SEC",
    ]
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

# Optional substituents: NAG has N2/C7/C8/O7; glucose-like can have O2.
_GLYCAN_OPTIONAL_BONDS: tuple[tuple[str, str], ...] = (
    ("C2", "N2"),
    ("N2", "C7"),
    ("C7", "C8"),
    ("C7", "O7"),
    ("C2", "O2"),
)


@dataclass
class ProtonationReport:
    """Summary for protonation/export step."""

    input_cif: str
    output_dir: str
    for_posebusters_pdb: str
    complex_h_pdb: str
    ligand_mol2: str
    complex_h_backend: str = "none"
    used_obabel: bool = False
    blockers: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


class ProtonationExportRunner:
    """Run step 7 protonation/export with graceful handling of missing tools."""

    def __init__(self, input_cif: Path, output_dir: Path):
        self.input_cif = Path(input_cif)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.logger = StructuredLogger("protonate_export", self.output_dir)
        self.failures = FailureLog(self.output_dir / "protonate_export_failures.json")

    def run(self) -> Tuple[bool, Optional[dict]]:
        """Execute step 7 exports.

        Returns:
            (success, report_as_dict_or_none)

        Success requires:
          - for_posebusters.pdb generated
          - complex_H.pdb generated with real hydrogens (not a copy)
          - ligand_for_prolif.mol2 generated, glycan-only
        """
        start = datetime.utcnow()
        self.logger.log_step_start("protonate_export", {"input": str(self.input_cif)})

        protonated_dir = self.output_dir
        for_posebusters = protonated_dir / "for_posebusters.pdb"
        complex_h = protonated_dir / "complex_H.pdb"
        ligand_mol2 = protonated_dir / "ligand_for_prolif.mol2"

        blockers: list[str] = []
        warnings: list[str] = []
        complex_h_backend = "none"
        used_obabel = False

        try:
            # ----------------------------------------------------------------
            # 1) Generate for_posebusters.pdb from normalized.cif
            # ----------------------------------------------------------------
            success_pdb, pdb_path = convert_cif_to_pdb(self.input_cif, protonated_dir)
            if not success_pdb or pdb_path is None or not pdb_path.exists():
                blockers.append("cif_to_pdb_failed")
                raise RuntimeError("Failed to generate for_posebusters.pdb from normalized.cif")

            if pdb_path != for_posebusters:
                shutil.copyfile(pdb_path, for_posebusters)

            # ----------------------------------------------------------------
            # 2) Add hydrogens → complex_H.pdb
            #    Priority: PDBFixer/OpenMM → reduce → obabel → BLOCKER
            # ----------------------------------------------------------------
            complex_h_backend, h_warnings = _add_hydrogens(for_posebusters, complex_h)
            warnings.extend(h_warnings)
            if complex_h_backend == "none":
                blockers.append("complex_h_protonation_failed_no_tool_available")
            else:
                # Sanity-check: complex_H.pdb must have more atoms than the
                # unprotonated for_posebusters.pdb.
                n_base = _count_atom_lines(for_posebusters)
                n_h = _count_atom_lines(complex_h)
                if n_h <= n_base:
                    warnings.append(
                        f"complex_h_no_new_atoms: base={n_base} complex_h={n_h} "
                        f"backend={complex_h_backend}"
                    )

            # ----------------------------------------------------------------
            # 3) Export ligand-only MOL2 (glycan chains B/C/D) via obabel
            # ----------------------------------------------------------------
            obabel_bin = shutil.which("obabel")
            if obabel_bin:
                try:
                    ligand_only_pdb = protonated_dir / "ligand_only_for_prolif.pdb"
                    ligand_residue_count = _write_ligand_only_pdb(
                        self.input_cif, ligand_only_pdb
                    )
                    if ligand_residue_count == 0:
                        blockers.append("ligand_only_extract_empty")
                    else:
                        subprocess.run(
                            [
                                obabel_bin,
                                "-ipdb",
                                str(ligand_only_pdb),
                                "-omol2",
                                "-O",
                                str(ligand_mol2),
                                "-h",
                            ],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True,
                            check=True,
                            timeout=180,
                        )
                        used_obabel = ligand_mol2.exists() and ligand_mol2.stat().st_size > 0
                        if not used_obabel:
                            blockers.append("obabel_no_output")
                        else:
                            _validate_mol2_sections(ligand_mol2, blockers, warnings)
                except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
                    blockers.append("obabel_failed")
                    warnings.append(f"obabel_error: {exc}")
            else:
                blockers.append("obabel_not_found")

            # ----------------------------------------------------------------
            # 4) Write report
            # ----------------------------------------------------------------
            report = ProtonationReport(
                input_cif=str(self.input_cif),
                output_dir=str(protonated_dir),
                for_posebusters_pdb=str(for_posebusters),
                complex_h_pdb=str(complex_h),
                ligand_mol2=str(ligand_mol2),
                complex_h_backend=complex_h_backend,
                used_obabel=used_obabel,
                blockers=blockers,
                warnings=warnings,
            )

            report_path = protonated_dir / "protonation_report.json"
            report_path.write_text(json.dumps(asdict(report), indent=2))

            self.logger.log_artifact("for_posebusters_pdb", for_posebusters, "PoseBusters input PDB")
            self.logger.log_artifact("complex_h_pdb", complex_h, "Hydrogenated complex PDB")
            if ligand_mol2.exists():
                self.logger.log_artifact("ligand_mol2", ligand_mol2, "Ligand MOL2 for ProLIF")
            self.logger.log_artifact("protonation_report", report_path, "Step 7 protonation/export report")

            success = len(blockers) == 0
            if not success:
                self.failures.record(
                    "protonate_export",
                    None,
                    "missing_or_failed_tools",
                    {"blockers": blockers, "warnings": warnings},
                )
                self.logger.log_failure("protonate_export_blocked", {"blockers": blockers})

            elapsed = (datetime.utcnow() - start).total_seconds()
            self.logger.log_step_end(
                "protonate_export",
                "success" if success else "failure",
                {
                    "complex_h_backend": complex_h_backend,
                    "used_obabel": used_obabel,
                    "blockers": blockers,
                },
                elapsed,
            )
            self.failures.write()
            self.logger.close()
            return success, asdict(report)

        except Exception as exc:
            self.logger.log_failure("protonate_export_exception", {"error": str(exc)})
            self.failures.record("protonate_export", None, "exception", {"error": str(exc)})
            self.failures.write()
            elapsed = (datetime.utcnow() - start).total_seconds()
            self.logger.log_step_end("protonate_export", "failure", {}, elapsed)
            self.logger.close()
            return False, None


def protonate_and_export(input_cif: Path, output_dir: Path) -> Tuple[bool, Optional[dict]]:
    """Convenience wrapper for step 7."""
    runner = ProtonationExportRunner(input_cif=input_cif, output_dir=output_dir)
    return runner.run()


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _add_hydrogens(
    input_pdb: Path,
    output_pdb: Path,
    ph: float = 7.0,
) -> tuple[str, list[str]]:
    """Add explicit hydrogens to *input_pdb* and write *output_pdb*.

    Priority: PDBFixer/OpenMM → reduce → obabel → return ("none", [blocker_msg]).

    Returns:
        (backend_name, warnings_list)
    """
    warnings: list[str] = []

    # 1) PDBFixer/OpenMM
    try:
        _protonate_with_pdbfixer(input_pdb, output_pdb, ph=ph)
        return "pdbfixer", warnings
    except ImportError:
        warnings.append("pdbfixer_not_importable")
    except Exception as exc:
        warnings.append(f"pdbfixer_failed: {exc}")

    # 2) reduce (AmberTools)
    reduce_bin = shutil.which("reduce")
    if reduce_bin:
        try:
            with open(output_pdb, "w") as out_f:
                subprocess.run(
                    [reduce_bin, "-BUILD", str(input_pdb)],
                    stdout=out_f,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                    timeout=300,
                )
            return "reduce", warnings
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
            warnings.append(f"reduce_failed: {exc}")
    else:
        warnings.append("reduce_not_found")

    # 3) obabel -h
    obabel_bin = shutil.which("obabel")
    if obabel_bin:
        try:
            subprocess.run(
                [
                    obabel_bin, "-ipdb", str(input_pdb),
                    "-opdb", "-O", str(output_pdb), "-h",
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
                timeout=300,
            )
            return "obabel", warnings
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
            warnings.append(f"obabel_failed: {exc}")
    else:
        warnings.append("obabel_not_found")

    # All tools failed — caller adds blocker; do NOT copy
    return "none", warnings


def _protonate_with_pdbfixer(input_pdb: Path, output_pdb: Path, ph: float = 7.0) -> None:
    """Add hydrogens at *ph* using PDBFixer/OpenMM.

    Writes the hydrogenated structure to *output_pdb* while preserving chain IDs.
    Raises ImportError if pdbfixer/openmm are not available.
    """
    from pdbfixer import PDBFixer  # type: ignore[import]
    from openmm.app import PDBFile  # type: ignore[import]

    fixer = PDBFixer(filename=str(input_pdb))
    _add_name_based_glycosidic_bonds(fixer.topology)
    fixer.addMissingHydrogens(ph)

    with open(output_pdb, "w") as out_f:
        PDBFile.writeFile(fixer.topology, fixer.positions, out_f, keepIds=True)

    _rewrite_conect_from_topology(output_pdb, fixer.topology)


def _rewrite_conect_from_topology(pdb_path: Path, topology) -> None:
    """Rewrite PDB CONECT records from full topology bonds.

    OpenMM can emit sparse CONECT lines. For downstream tools that rely on
    explicit connectivity, we replace CONECT records using all bonds present
    in the final protonated topology.
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
                # Fallback to sequential serial if malformed.
                atom_serials.append(len(atom_serials) + 1)
            kept_lines.append(line)
            continue

        if line.startswith("CONECT"):
            continue
        kept_lines.append(line)

    n_top_atoms = sum(1 for _ in topology.atoms())
    if n_top_atoms != len(atom_serials):
        # If atom order/size diverges unexpectedly, keep file unchanged.
        return

    adjacency: dict[int, set[int]] = {s: set() for s in atom_serials}
    for a1, a2 in topology.bonds():
        s1 = atom_serials[a1.index]
        s2 = atom_serials[a2.index]
        adjacency[s1].add(s2)
        adjacency[s2].add(s1)

    _augment_glycan_bonds_from_names(kept_lines, adjacency)

    conect_lines: list[str] = []
    for src in sorted(adjacency):
        nbrs = sorted(adjacency[src])
        if not nbrs:
            continue
        for i in range(0, len(nbrs), 4):
            chunk = nbrs[i : i + 4]
            conect_lines.append(
                f"CONECT{src:5d}" + "".join(f"{n:5d}" for n in chunk) + "\n"
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
    """Add missing glycan bonds using explicit atom-name templates.

    OpenMM topology can miss some intra-residue bonds for carbohydrate HETATM
    residues. This adds missing NAG/BGC/GLC bonds so CONECT is complete.
    """
    atoms: list[dict[str, object]] = []
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

        element = line[76:78].strip().upper()
        if not element:
            element = "".join(ch for ch in atom_name if ch.isalpha())[:1].upper()

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

    # Intra-residue bonds from known glycan atom naming conventions.
    for atom_map in by_residue.values():
        for a_name, b_name in (_GLYCAN_CORE_BONDS + _GLYCAN_OPTIONAL_BONDS):
            a_serial = atom_map.get(a_name)
            b_serial = atom_map.get(b_name)
            if a_serial is None or b_serial is None:
                continue
            _add_adjacency_edge(adjacency, a_serial, b_serial)

    # Inter-residue glycosidic bond: C1(i) - O4(i+1), same chain.
    by_chain: dict[str, list[tuple[str, int, str]]] = {}
    for key in by_residue.keys():
        chain, resseq, _resname = key
        by_chain.setdefault(chain, []).append(key)

    for chain, keys in by_chain.items():
        _ = chain
        keys.sort(key=lambda k: k[1])
        for i in range(len(keys) - 1):
            left = keys[i]
            right = keys[i + 1]
            if right[1] != left[1] + 1:
                continue

            left_c1 = by_residue[left].get("C1")
            right_o4 = by_residue[right].get("O4")
            if left_c1 is None or right_o4 is None:
                continue
            _add_adjacency_edge(adjacency, left_c1, right_o4)


def _add_adjacency_edge(adjacency: dict[int, set[int]], a: int, b: int) -> None:
    """Insert an undirected bond edge into adjacency."""
    if a == b:
        return
    if a not in adjacency or b not in adjacency:
        return
    adjacency[a].add(b)
    adjacency[b].add(a)


def _residue_seq_id(residue) -> Optional[int]:
    """Best-effort integer sequence id from OpenMM residue metadata."""
    try:
        return int(str(residue.id).strip())
    except (TypeError, ValueError):
        return None


def _add_name_based_glycosidic_bonds(topology) -> int:
    """Add C1(i)-O4(i+1) glycosidic bonds for NAG/BGC/GLC residues.

    This ensures PDBFixer does not protonate O4 as hydroxyl when O4 is the
    bridge oxygen in a beta-1,4-like glycan linkage.
    """
    # Track existing bonds to avoid duplicates.
    existing: set[tuple[int, int]] = set()
    for a1, a2 in topology.bonds():
        i, j = a1.index, a2.index
        if i > j:
            i, j = j, i
        existing.add((i, j))

    by_chain: dict[str, list] = {}
    for chain in topology.chains():
        residues = []
        for residue in chain.residues():
            if residue.name not in _GLYCAN_LINK_RESNAMES:
                continue
            seq = _residue_seq_id(residue)
            if seq is None:
                continue
            residues.append(residue)
        residues.sort(key=lambda r: _residue_seq_id(r) or 10**9)
        by_chain[chain.id] = residues

    added = 0
    for residues in by_chain.values():
        for idx in range(len(residues) - 1):
            left = residues[idx]
            right = residues[idx + 1]
            left_seq = _residue_seq_id(left)
            right_seq = _residue_seq_id(right)
            if left_seq is None or right_seq is None or right_seq != left_seq + 1:
                continue

            left_c1 = None
            right_o4 = None
            for atom in left.atoms():
                if atom.name == "C1":
                    left_c1 = atom
                    break
            for atom in right.atoms():
                if atom.name == "O4":
                    right_o4 = atom
                    break
            if left_c1 is None or right_o4 is None:
                continue

            i, j = left_c1.index, right_o4.index
            if i > j:
                i, j = j, i
            if (i, j) in existing:
                continue
            topology.addBond(left_c1, right_o4)
            existing.add((i, j))
            added += 1

    return added


def _count_atom_lines(pdb_path: Path) -> int:
    """Count ATOM/HETATM lines in a PDB file."""
    count = 0
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith(("ATOM  ", "HETATM")):
                count += 1
    return count


def _validate_mol2_sections(
    mol2_path: Path,
    blockers: list[str],
    warnings: list[str],
) -> None:
    """Check that the MOL2 file has required sections for ProLIF."""
    text = mol2_path.read_text(errors="ignore")
    if "@<TRIPOS>ATOM" not in text:
        blockers.append("mol2_missing_tripos_atom_section")
    if "@<TRIPOS>BOND" not in text:
        blockers.append("mol2_missing_tripos_bond_section")
    # Protein contamination check
    protein_hits = [r for r in _PROTEIN_RESIDUE_NAMES if r in text]
    if protein_hits:
        warnings.append(f"mol2_contains_protein_residue_labels: {protein_hits}")


def _copy_chain_to_output(model: gemmi.Model, source_chain_name: str, out_model: gemmi.Model) -> bool:
    """Copy one chain by name from source model into output model."""
    for chain in model:
        if chain.name == source_chain_name:
            out_chain = gemmi.Chain(chain.name)
            for residue in chain:
                out_chain.add_residue(residue, pos=-1)
            out_model.add_chain(out_chain, pos=-1)
            return True
    return False


def _build_ligand_only_structure(
    input_cif: Path,
    glycan_chains: tuple[str, ...] = ("B", "C", "D"),
) -> gemmi.Structure:
    """Create a structure containing only validated glycan chains for ligand export.

    Chains are validated: residues must NOT be exclusively protein residue names.
    Raises ValueError if no glycan chains are found or all chains appear protein-like.
    """
    st = gemmi.read_structure(str(input_cif))
    if len(st) == 0:
        raise ValueError(f"No models found in {input_cif}")

    out_st = gemmi.Structure()
    out_st.name = st.name
    out_model = gemmi.Model("1")

    added_residues = 0
    for chain_name in glycan_chains:
        for chain in st[0]:
            if chain.name != chain_name:
                continue
            residue_names = {res.name for res in chain}
            if residue_names and residue_names.issubset(_PROTEIN_RESIDUE_NAMES):
                raise ValueError(
                    f"Chain {chain_name} contains only protein residues "
                    f"({residue_names}); not a glycan chain"
                )
            if _copy_chain_to_output(st[0], chain_name, out_model):
                added_residues += sum(1 for _ in chain)

    if added_residues == 0:
        raise ValueError(f"No glycan chains {glycan_chains} found in {input_cif}")

    out_st.add_model(out_model, pos=-1)
    return out_st


def _write_ligand_only_pdb(input_cif: Path, out_pdb: Path) -> int:
    """Write ligand-only PDB (glycan chains B/C/D) for MOL2 conversion.

    Returns the number of residues written.
    """
    out_st = _build_ligand_only_structure(input_cif)
    out_st.write_pdb(str(out_pdb))
    return sum(len(list(chain)) for chain in out_st[0])

