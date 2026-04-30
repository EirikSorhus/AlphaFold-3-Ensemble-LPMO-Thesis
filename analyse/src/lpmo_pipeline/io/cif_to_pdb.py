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
from typing import Optional, Tuple

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
    backend_fallback_reason: str = ""


class CIFToPDBRunner:
    """Convert normalized mmCIF to a PoseBusters-ready PDB artifact."""

    def __init__(self, input_cif: Path, output_dir: Path):
        self.input_cif = Path(input_cif)
        self.output_dir = Path(output_dir)
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
        output_path = self.output_dir / "for_posebusters.pdb"

        with open(output_path, "w") as out_f:
            PDBFile.writeFile(fixer.topology, fixer.positions, out_f, keepIds=True)

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


def convert_cif_to_pdb(input_cif: Path, output_dir: Path) -> Tuple[bool, Optional[Path]]:
    """Convenience wrapper for step 7b conversion."""
    return CIFToPDBRunner(input_cif=input_cif, output_dir=output_dir).run()