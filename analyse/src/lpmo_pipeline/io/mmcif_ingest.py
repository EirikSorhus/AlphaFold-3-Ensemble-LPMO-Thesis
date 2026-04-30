"""
LPMO Pipeline: Step 1 — Ingest & QC
Responsibility: Parse raw mmCIF from prediction models, validate schema

Step 1 contract:
  Input: raw prediction output (PDB or mmCIF)
  Output: ingest_qc.log, pass/fail flag
  Tools: Gemmi, custom validation
  Checks: mmcif_parse_ok, required_categories_present, atom_count > 0,
          entity_types_present, chain_summary
"""

from pathlib import Path
from typing import Tuple, Optional, Dict, List
import json
from datetime import datetime

from lpmo_pipeline.io.gemmi_compat import gemmi
from lpmo_pipeline.utils.logging import StructuredLogger, FailureLog
from lpmo_pipeline.utils.data_models import QCFlag, QCStatus, FailureReason


# Categories that must be present in all AF3 output CIFs
REQUIRED_MMCIF_CATEGORIES = [
    "_atom_site",
    "_chem_comp",
    "_entity",
]


class IngestQCRunner:
    """Step 1: Ingest raw prediction, parse, validate schema."""

    def __init__(self, input_path: Path, output_dir: Path):
        """
        Args:
            input_path: raw PDB or mmCIF file
            output_dir: directory for logs + outputs
        """
        self.input_path = Path(input_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.logger = StructuredLogger("ingest_qc", self.output_dir)
        self.failures = FailureLog(self.output_dir / "ingest_failures.json")

        self.qc_flags: list[QCFlag] = []
        # Populated after successful parse
        self.structure: Optional[gemmi.Structure] = None
        self.cif_block: Optional[gemmi.cif.Block] = None

    def run(self) -> Tuple[bool, Optional[str]]:
        """
        Execute ingestion and QC.

        Returns:
            (pass_flag, input_mmcif_path_or_none)
        """
        start = datetime.utcnow()

        # Check: file exists and readable (before stat)
        if not self.input_path.exists():
            self.logger.log_step_start(
                "ingest_qc",
                {"input_file": str(self.input_path), "size_bytes": 0},
            )
            self._check_fail("file_exists", f"Input file not found: {self.input_path}")
            self.logger.log_step_end("ingest_qc", "failure", {}, 0.0)
            self.failures.write()
            return False, None

        self.logger.log_step_start(
            "ingest_qc",
            {"input_file": str(self.input_path), "size_bytes": self.input_path.stat().st_size},
        )

        # Check: parse as mmCIF
        try:
            self._parse_with_gemmi()
        except Exception as e:
            self._check_fail("mmcif_parse_ok", f"Parse error: {e}")
            self.logger.log_step_end("ingest_qc", "failure", {}, 0.0)
            self.failures.write()
            return False, None

        # Check: required categories present in CIF block
        self._check_categories()

        # Check: atom count > 0
        self._check_atom_count()

        # Check: entity types make sense for LPMO complex
        self._check_entity_types()

        # Collect all checks
        hard_fail = any(f.status == QCStatus.HARD_FAIL for f in self.qc_flags)

        # Build chain summary
        chain_summary = self._build_chain_summary()

        # Write QC report
        qc_report_path = self.output_dir / "ingest_qc_report.json"
        self._write_qc_report(qc_report_path, chain_summary)
        self.logger.log_artifact("ingest_qc_report", qc_report_path, "QC checks for ingest step")

        # Write log file (pass/fail)
        log_dict = {
            "timestamp": datetime.utcnow().isoformat(),
            "step": "ingest_qc",
            "status": "failure" if hard_fail else "pass",
            "checks": [
                {
                    "name": f.check_name,
                    "status": f.status.value,
                    "message": f.message,
                }
                for f in self.qc_flags
            ],
            "chain_summary": chain_summary,
        }
        log_path = self.output_dir / "ingest_qc.log"
        with open(log_path, "w") as f:
            json.dump(log_dict, f, indent=2)

        elapsed = (datetime.utcnow() - start).total_seconds()
        self.logger.log_step_end(
            "ingest_qc",
            "failure" if hard_fail else "success",
            {"qc_flags": len(self.qc_flags), "hard_fail": hard_fail},
            elapsed,
        )

        self.failures.write()
        self.logger.close()

        return not hard_fail, str(self.input_path)

    def _parse_with_gemmi(self) -> None:
        """Parse input file using Gemmi. Stores structure and CIF block."""
        if gemmi is None:
            raise ImportError("Gemmi not installed")

        # Read as structure (handles both mmCIF and PDB)
        self.structure = gemmi.read_structure(str(self.input_path))
        self._check_pass("mmcif_parse_ok", "Parsed successfully with gemmi")

        # Also read the raw CIF document for category checks
        if self.input_path.suffix in (".cif", ".mmcif"):
            doc = gemmi.cif.read(str(self.input_path))
            self.cif_block = doc.sole_block()

    def _check_categories(self) -> None:
        """Verify required mmCIF categories are present in the CIF block."""
        if self.cif_block is None:
            self._check_pass(
                "required_categories",
                "Non-CIF input (PDB) — category check skipped",
            )
            return

        missing: List[str] = []
        for cat in REQUIRED_MMCIF_CATEGORIES:
            # gemmi Block.find() returns a Table; check if it has rows
            tag = f"{cat}.id"
            col = self.cif_block.find_values(tag)
            if len(col) == 0:
                # Try a common first column instead
                alt_tags = {
                    "_atom_site": "_atom_site.id",
                    "_chem_comp": "_chem_comp.id",
                    "_entity": "_entity.id",
                }
                alt = alt_tags.get(cat, tag)
                col = self.cif_block.find_values(alt)
                if len(col) == 0:
                    missing.append(cat)

        if missing:
            self._check_fail(
                "required_categories",
                f"Missing mmCIF categories: {', '.join(missing)}",
            )
        else:
            self._check_pass(
                "required_categories",
                f"All required categories present: {', '.join(REQUIRED_MMCIF_CATEGORIES)}",
            )

    def _check_atom_count(self) -> None:
        """Ensure we have atoms in the structure."""
        if self.structure is None:
            self._check_fail("atom_count_positive", "No structure parsed")
            return

        try:
            model = self.structure[0]
            n_atoms = model.count_atom_sites()
            if n_atoms > 0:
                self._check_pass("atom_count_positive", f"Atom count: {n_atoms}")
            else:
                self._check_fail("atom_count_positive", "No atoms in structure")
        except Exception as e:
            self._check_fail("atom_count_positive", str(e))

    def _check_entity_types(self) -> None:
        """Check that expected entity types are present (polymer + non-polymer)."""
        if self.cif_block is None:
            return

        entity_types = list(self.cif_block.find_values("_entity.type"))
        if not entity_types:
            self._check_fail("entity_types", "No _entity.type found")
            return

        has_polymer = any(t.strip() == "polymer" for t in entity_types)
        has_nonpolymer = any(t.strip() == "non-polymer" for t in entity_types)

        if has_polymer and has_nonpolymer:
            self._check_pass(
                "entity_types",
                f"Entity types: {[t.strip() for t in entity_types]}",
            )
        elif has_polymer:
            self._check_pass(
                "entity_types",
                f"Entity types: {[t.strip() for t in entity_types]} (no non-polymer — Cu may be missing)",
            )
        else:
            self._check_fail("entity_types", f"Missing polymer entity; found: {entity_types}")

    def _build_chain_summary(self) -> Dict:
        """Build a summary of chain composition."""
        if self.structure is None:
            return {}

        model = self.structure[0]
        summary = {}
        for chain in model:
            residues = list(chain)
            comp_ids = sorted({r.name for r in residues})
            n_atoms = sum(len(list(r)) for r in residues)
            summary[chain.name] = {
                "n_residues": len(residues),
                "n_atoms": n_atoms,
                "comp_ids": comp_ids,
            }
        return summary
    
    def _check_pass(self, check_name: str, message: str):
        """Record a passed check."""
        flag = QCFlag(
            check_name=check_name,
            status=QCStatus.PASS,
            message=message,
        )
        self.qc_flags.append(flag)
        self.logger.log_check(check_name, "pass", message)
    
    def _check_fail(self, check_name: str, message: str):
        """Record a hard-fail check."""
        flag = QCFlag(
            check_name=check_name,
            status=QCStatus.HARD_FAIL,
            message=message,
        )
        self.qc_flags.append(flag)
        self.logger.log_check(check_name, "hard_fail", message)
        self.failures.record("ingest_qc", None, f"check_{check_name}", {"message": message})
    
    def _write_qc_report(self, path: Path, chain_summary: Optional[Dict] = None):
        """Write structured QC report."""
        report = {
            "timestamp": datetime.utcnow().isoformat(),
            "step": "ingest_qc",
            "input_file": str(self.input_path),
            "checks": [
                {
                    "name": f.check_name,
                    "status": f.status.value,
                    "message": f.message,
                    "data": f.data,
                }
                for f in self.qc_flags
            ],
            "chain_summary": chain_summary or {},
        }
        with open(path, "w") as f:
            json.dump(report, f, indent=2)
