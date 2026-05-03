from __future__ import annotations

import copy
import shlex
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
from typing import Any

try:
    import gemmi as _native_gemmi  # type: ignore[import-not-found]
except ImportError:
    _native_gemmi = None


def _has_native_gemmi(module: Any) -> bool:
    return bool(module) and all(
        hasattr(module, attr)
        for attr in ("read_structure", "cif", "Structure", "Model", "Chain")
    )


_REAL_GEMMI = _native_gemmi if _has_native_gemmi(_native_gemmi) else None


def gemmi_is_functional() -> bool:
    return _REAL_GEMMI is not None


def gemmi_version() -> str:
    if _REAL_GEMMI is None:
        return "not_functional"
    return getattr(_REAL_GEMMI, "__version__", "unknown")


def remap_cif_block_values(block: Any, remap: dict[str, str], tags: list[str]) -> None:
    if hasattr(block, "remap_values"):
        for tag in tags:
            block.remap_values(tag, remap)
        return

    for tag in tags:
        pair = block.find_pair(tag)
        if pair:
            value = pair[1]
            new_value = remap.get(str(value).strip(), value)
            if new_value != value:
                block.set_pair(tag, str(new_value))

        column = block.find_loop(tag)
        if not column:
            continue

        loop = column.get_loop()
        tags_in_loop = list(loop.tags)
        if tag not in tags_in_loop:
            continue

        column_index = tags_in_loop.index(tag)
        width = loop.width()
        columns: list[list[str]] = []
        for current_column_index in range(width):
            current_column: list[str] = []
            for row_index in range(loop.length()):
                offset = row_index * width + current_column_index
                current_column.append(str(loop.values[offset]))
            columns.append(current_column)

        columns[column_index] = [
            str(remap.get(str(value).strip(), value))
            for value in columns[column_index]
        ]
        loop.set_all_values(columns)


if _REAL_GEMMI is not None:
    gemmi = _REAL_GEMMI
else:

    @dataclass
    class _CompatSeqId:
        num: int


    @dataclass
    class _CompatElement:
        name: str


    @dataclass
    class _CompatAtom:
        name: str
        element: _CompatElement
        b_iso: float = 0.0
        serial: int = 0
        x: float = 0.0
        y: float = 0.0
        z: float = 0.0
        record_name: str = "HETATM"


    @dataclass
    class Residue:
        name: str
        seqid: _CompatSeqId
        subchain: str = ""
        atoms: list[_CompatAtom] = field(default_factory=list)

        def __iter__(self):
            return iter(self.atoms)

        def __len__(self) -> int:
            return len(self.atoms)

        def add_atom(self, atom: _CompatAtom) -> None:
            self.atoms.append(atom)


    @dataclass
    class Chain:
        name: str
        residues: list[Residue] = field(default_factory=list)

        def __iter__(self):
            return iter(self.residues)

        def __len__(self) -> int:
            return len(self.residues)

        def add_residue(self, residue: Residue, pos: int = -1) -> None:
            copied = copy.deepcopy(residue)
            if pos < 0 or pos >= len(self.residues):
                self.residues.append(copied)
            else:
                self.residues.insert(pos, copied)


    @dataclass
    class Model:
        name: str
        chains: list[Chain] = field(default_factory=list)

        def __iter__(self):
            return iter(self.chains)

        def __len__(self) -> int:
            return len(self.chains)

        def add_chain(self, chain: Chain, pos: int = -1) -> None:
            copied = copy.deepcopy(chain)
            if pos < 0 or pos >= len(self.chains):
                self.chains.append(copied)
            else:
                self.chains.insert(pos, copied)

        def count_atom_sites(self) -> int:
            return sum(len(residue.atoms) for chain in self.chains for residue in chain.residues)


    @dataclass
    class _CompatLoop:
        tags: list[str]
        rows: list[list[str]] = field(default_factory=list)


    class Block:
        def __init__(self, name: str) -> None:
            self.name = name
            self._scalars: dict[str, str] = {}
            self._loops: list[_CompatLoop] = []
            self._records: list[tuple[str, str] | _CompatLoop] = []
            self._tag_to_loop: dict[str, _CompatLoop] = {}

        def add_scalar(self, tag: str, value: str) -> None:
            self._scalars[tag] = value
            self._records.append((tag, value))

        def add_loop(self, loop: _CompatLoop) -> None:
            self._loops.append(loop)
            self._records.append(loop)
            for tag in loop.tags:
                self._tag_to_loop[tag] = loop

        def find_values(self, tag: str) -> list[str]:
            if tag in self._scalars:
                return [self._scalars[tag]]
            loop = self._tag_to_loop.get(tag)
            if loop is None:
                return []
            index = loop.tags.index(tag)
            return [row[index] for row in loop.rows]

        def find(self, tags_or_prefix: list[str] | str, columns: list[str] | None = None) -> list[list[str]]:
            if columns is None:
                tags = list(tags_or_prefix)
            else:
                prefix = str(tags_or_prefix)
                tags = [f"{prefix}{column}" for column in columns]

            for loop in self._loops:
                if all(tag in loop.tags for tag in tags):
                    indexes = [loop.tags.index(tag) for tag in tags]
                    return [[row[index] for index in indexes] for row in loop.rows]
            return []

        def remap_values(self, tag: str, remap: dict[str, str]) -> None:
            if tag in self._scalars:
                old_value = self._scalars[tag]
                new_value = remap.get(old_value.strip(), old_value)
                if new_value == old_value:
                    return
                self._scalars[tag] = new_value
                for index, record in enumerate(self._records):
                    if isinstance(record, tuple) and record[0] == tag:
                        self._records[index] = (tag, new_value)
                        break
                return

            loop = self._tag_to_loop.get(tag)
            if loop is None:
                return
            index = loop.tags.index(tag)
            for row in loop.rows:
                row[index] = remap.get(row[index].strip(), row[index])

        def to_mmcif_text(self) -> str:
            lines = [f"data_{self.name or 'compat'}", "#"]
            for record in self._records:
                if isinstance(record, tuple):
                    tag, value = record
                    lines.append(f"{tag} {_format_mmcif_value(value)}")
                    lines.append("#")
                    continue

                lines.append("loop_")
                lines.extend(record.tags)
                for row in record.rows:
                    lines.append(" ".join(_format_mmcif_value(value) for value in row))
                lines.append("#")
            return "\n".join(lines) + "\n"


    class _CompatDocument:
        def __init__(self, block: Block) -> None:
            self._block = block

        def sole_block(self) -> Block:
            return self._block

        def write_file(self, path: str) -> None:
            Path(path).write_text(self._block.to_mmcif_text())


    @dataclass
    class PdbWriteOptions:
        minimal_file: bool = False
        atom_records: bool = True
        seqres_records: bool = False
        ssbond_records: bool = False
        link_records: bool = True
        conect_records: bool = True
        ter_records: bool = True
        end_record: bool = True
        preserve_serial: bool = True


    @dataclass
    class Structure:
        name: str = ""
        models: list[Model] = field(default_factory=list)
        block: Block | None = None

        def __iter__(self):
            return iter(self.models)

        def __len__(self) -> int:
            return len(self.models)

        def __getitem__(self, index: int) -> Model:
            return self.models[index]

        def add_model(self, model: Model, pos: int = -1) -> None:
            copied = copy.deepcopy(model)
            if pos < 0 or pos >= len(self.models):
                self.models.append(copied)
            else:
                self.models.insert(pos, copied)

        def make_mmcif_document(self) -> _CompatDocument:
            block = self.block if self.block is not None else _build_block_from_structure(self)
            return _CompatDocument(block)

        def write_pdb(self, path: str, options: PdbWriteOptions | None = None) -> None:
            lines: list[str] = []
            serial = 1
            for model in self.models:
                current_chain = None
                current_resseq = None
                current_resname = None
                for chain in model:
                    for residue in chain:
                        if current_chain is not None and (
                            chain.name != current_chain or residue.seqid.num != current_resseq
                        ):
                            lines.append(
                                f"TER   {serial:>5}      {current_resname:>3} {current_chain:1}{current_resseq:>4}"
                            )
                            serial += 1

                        for atom in residue:
                            atom_serial = atom.serial or serial
                            lines.append(
                                _format_pdb_atom_line(
                                    record_name=atom.record_name,
                                    serial=atom_serial,
                                    atom_name=atom.name,
                                    resname=residue.name,
                                    chain_id=chain.name,
                                    resseq=residue.seqid.num,
                                    x=atom.x,
                                    y=atom.y,
                                    z=atom.z,
                                    b_iso=atom.b_iso,
                                    element=atom.element.name,
                                )
                            )
                            serial = atom_serial + 1

                        current_chain = chain.name
                        current_resseq = residue.seqid.num
                        current_resname = residue.name

            lines.append("END")
            Path(path).write_text("\n".join(lines) + "\n")


    class _CompatCifModule:
        Block = Block

        @staticmethod
        def read(path: str) -> _CompatDocument:
            return _CompatDocument(_parse_mmcif_block(Path(path).read_text(), Path(path).stem))


    def read_structure(path: str) -> Structure:
        input_path = Path(path)
        if input_path.suffix.lower() in {".cif", ".mmcif"}:
            block = _parse_mmcif_block(input_path.read_text(), input_path.stem)
            return _build_structure_from_block(block)
        return _build_structure_from_pdb(input_path.read_text(), input_path.stem)


    gemmi = SimpleNamespace(
        read_structure=read_structure,
        cif=_CompatCifModule(),
        Structure=Structure,
        Model=Model,
        Chain=Chain,
        PdbWriteOptions=PdbWriteOptions,
        __version__="compat-fallback",
    )


    def _parse_mmcif_block(text: str, default_name: str) -> Block:
        lines = text.splitlines()
        block_name = default_name
        block = Block(default_name)
        index = 0
        while index < len(lines):
            stripped = lines[index].strip()
            if not stripped or stripped == "#":
                index += 1
                continue
            if stripped.startswith("data_"):
                block_name = stripped[5:].strip() or default_name
                block = Block(block_name)
                index += 1
                continue
            if stripped == "loop_":
                index += 1
                tags: list[str] = []
                while index < len(lines) and lines[index].strip().startswith("_"):
                    tags.append(lines[index].strip())
                    index += 1
                loop = _CompatLoop(tags=tags)
                while index < len(lines):
                    stripped = lines[index].strip()
                    if not stripped:
                        index += 1
                        continue
                    if stripped == "#":
                        index += 1
                        break
                    if stripped == "loop_" or stripped.startswith("_") or stripped.startswith("data_"):
                        break
                    fields = shlex.split(lines[index], posix=True)
                    if len(fields) >= len(tags):
                        loop.rows.append(fields[: len(tags)])
                    index += 1
                block.add_loop(loop)
                continue
            if stripped.startswith("_"):
                fields = shlex.split(stripped, posix=True)
                tag = fields[0]
                if len(fields) >= 2:
                    value = fields[1]
                    index += 1
                else:
                    index += 1
                    value = lines[index].strip() if index < len(lines) else "?"
                    index += 1
                block.add_scalar(tag, value)
                continue
            index += 1

        return block


    def _build_structure_from_block(block: Block) -> Structure:
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
        if not table:
            alt_table = block.find(
                [
                    "_atom_site.group_PDB",
                    "_atom_site.id",
                    "_atom_site.type_symbol",
                    "_atom_site.label_atom_id",
                    "_atom_site.label_comp_id",
                    "_atom_site.auth_asym_id",
                    "_atom_site.auth_seq_id",
                    "_atom_site.Cartn_x",
                    "_atom_site.Cartn_y",
                    "_atom_site.Cartn_z",
                    "_atom_site.B_iso_or_equiv",
                ]
            )
            table = alt_table

        structure = Structure(name=block.name, block=block)
        model = Model("1")
        structure.add_model(model)

        chain_map: dict[str, Chain] = {}
        residue_map: dict[tuple[str, str, int], Residue] = {}

        # Prefer auth identifiers when present.
        auth_chain = block.find_values("_atom_site.auth_asym_id")
        auth_seq = block.find_values("_atom_site.auth_seq_id")
        model_nums = block.find_values("_atom_site.pdbx_PDB_model_num")

        for row_index, row in enumerate(table):
            record_name, serial_text, element_text, atom_name, resname, label_chain, label_seq, x, y, z, b_iso = row
            chain_name = auth_chain[row_index].strip() if row_index < len(auth_chain) and auth_chain[row_index].strip() not in {"", ".", "?"} else label_chain.strip()
            seq_text = auth_seq[row_index].strip() if row_index < len(auth_seq) and auth_seq[row_index].strip() not in {"", ".", "?"} else label_seq.strip()
            model_name = model_nums[row_index].strip() if row_index < len(model_nums) else "1"

            while not structure.models or structure.models[-1].name != model_name:
                structure.add_model(Model(model_name))
            model = structure.models[-1]

            chain = chain_map.get(f"{model_name}:{chain_name}")
            if chain is None:
                chain = Chain(chain_name)
                model.add_chain(chain)
                chain = model.chains[-1]
                chain_map[f"{model_name}:{chain_name}"] = chain

            resseq = _safe_int(seq_text)
            residue_key = (f"{model_name}:{chain_name}", resname.strip(), resseq)
            residue = residue_map.get(residue_key)
            if residue is None:
                residue = Residue(name=resname.strip(), seqid=_CompatSeqId(resseq), subchain=chain_name)
                chain.add_residue(residue)
                residue = chain.residues[-1]
                residue_map[residue_key] = residue

            residue.add_atom(
                _CompatAtom(
                    name=atom_name.strip(),
                    element=_CompatElement(element_text.strip() or atom_name.strip()[:1]),
                    b_iso=_safe_float(b_iso),
                    serial=_safe_int(serial_text),
                    x=_safe_float(x),
                    y=_safe_float(y),
                    z=_safe_float(z),
                    record_name="ATOM" if record_name.strip() == "ATOM" else "HETATM",
                )
            )

        return structure


    def _build_structure_from_pdb(text: str, name: str) -> Structure:
        structure = Structure(name=name)
        model = Model("1")
        structure.add_model(model)

        chain_map: dict[str, Chain] = {}
        residue_map: dict[tuple[str, str, int], Residue] = {}

        for line in text.splitlines():
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            record_name = line[:6].strip() or "HETATM"
            serial = _safe_int(line[6:11].strip())
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            chain_name = line[21].strip() or "A"
            resseq = _safe_int(line[22:26].strip())
            x = _safe_float(line[30:38].strip())
            y = _safe_float(line[38:46].strip())
            z = _safe_float(line[46:54].strip())
            b_iso = _safe_float(line[60:66].strip())
            element = line[76:78].strip() or atom_name[:1]

            chain = chain_map.get(chain_name)
            if chain is None:
                chain = Chain(chain_name)
                model.add_chain(chain)
                chain = model.chains[-1]
                chain_map[chain_name] = chain

            residue_key = (chain_name, resname, resseq)
            residue = residue_map.get(residue_key)
            if residue is None:
                residue = Residue(name=resname, seqid=_CompatSeqId(resseq), subchain=chain_name)
                chain.add_residue(residue)
                residue = chain.residues[-1]
                residue_map[residue_key] = residue

            residue.add_atom(
                _CompatAtom(
                    name=atom_name,
                    element=_CompatElement(element),
                    b_iso=b_iso,
                    serial=serial,
                    x=x,
                    y=y,
                    z=z,
                    record_name=record_name,
                )
            )

        return structure


    def _build_block_from_structure(structure: Structure) -> Block:
        block = Block(structure.name or "compat")
        atom_loop = _CompatLoop(
            tags=[
                "_atom_site.group_PDB",
                "_atom_site.id",
                "_atom_site.type_symbol",
                "_atom_site.label_atom_id",
                "_atom_site.label_alt_id",
                "_atom_site.label_comp_id",
                "_atom_site.label_asym_id",
                "_atom_site.label_entity_id",
                "_atom_site.label_seq_id",
                "_atom_site.pdbx_PDB_ins_code",
                "_atom_site.Cartn_x",
                "_atom_site.Cartn_y",
                "_atom_site.Cartn_z",
                "_atom_site.occupancy",
                "_atom_site.B_iso_or_equiv",
                "_atom_site.auth_seq_id",
                "_atom_site.auth_asym_id",
                "_atom_site.pdbx_PDB_model_num",
            ]
        )
        serial = 1
        for model_index, model in enumerate(structure.models, start=1):
            for entity_index, chain in enumerate(model, start=1):
                for residue in chain:
                    for atom in residue:
                        atom_loop.rows.append(
                            [
                                "ATOM" if atom.record_name == "ATOM" else "HETATM",
                                str(atom.serial or serial),
                                atom.element.name or atom.name[:1],
                                atom.name,
                                ".",
                                residue.name,
                                chain.name,
                                str(entity_index),
                                str(residue.seqid.num),
                                "?",
                                f"{atom.x:.3f}",
                                f"{atom.y:.3f}",
                                f"{atom.z:.3f}",
                                "1.00",
                                f"{atom.b_iso:.2f}",
                                str(residue.seqid.num),
                                chain.name,
                                str(model_index),
                            ]
                        )
                        serial += 1
        block.add_loop(atom_loop)
        return block


    def _format_mmcif_value(value: str) -> str:
        if value in {"", None}:
            return "?"
        text = str(value)
        if any(ch.isspace() for ch in text) or any(ch in text for ch in ("'", '"', "#")):
            escaped = text.replace("'", "''")
            return f"'{escaped}'"
        return text


    def _safe_int(value: str) -> int:
        try:
            return int(float(value))
        except (TypeError, ValueError):
            return 0


    def _safe_float(value: str) -> float:
        try:
            return float(value)
        except (TypeError, ValueError):
            return 0.0


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
        return (
            f"{record_name:<6}{serial:>5} {atom_field} {resname:>3} {chain_id[:1]:1}"
            f"{resseq:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{b_iso:>6.2f}          {element[:2].rjust(2)}"
        )