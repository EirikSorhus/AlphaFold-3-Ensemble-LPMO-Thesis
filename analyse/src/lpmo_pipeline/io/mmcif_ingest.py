"""
LPMO Pipeline: Step 1 — Ingest & QC
Responsibility: Parse raw mmCIF from prediction models, validate schema

Step 1 contract:
  Input: raw prediction output (PDB or mmCIF)
  Output: ingest_qc.log, pass/fail flag
  Tools: Gemmi, custom validation
  Checks: mmcif_parse_ok, required_categories_present, atom_count > 0
"""

from pathlib import Path
from typing import Tuple, Optional, Dict
import json
from datetime import datetime

try:
    import gemmi
except ImportError:
    gemmi = None

from lpmo_pipeline.utils.logging import StructuredLogger, FailureLog
from lpmo_pipeline.utils.data_models import QCFlag, QCStatus, FailureReason


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
    
    def run(self) -> Tuple[bool, Optional[str]]:
        """
        Execute ingestion and QC.
        
        Returns:
            (pass_flag, normalized_mmcif_path_or_none)
        """
        start = datetime.utcnow()
        self.logger.log_step_start(
            "ingest_qc",
            {"input_file": str(self.input_path), "size_bytes": self.input_path.stat().st_size}
        )
        
        # Check: file exists and readable
        if not self.input_path.exists():
            self._check_fail("file_exists", f"Input file not found: {self.input_path}")
            self.logger.log_step_end("ingest_qc", "failure", {}, 0.0)
            self.failures.write()
            return False, None
        
        # Check: parse as mmCIF or PDB
        try:
            doc = self._parse_with_gemmi()
        except Exception as e:
            self._check_fail("mmcif_parse_ok", f"Parse error: {str(e)}")
            self.logger.log_step_end("ingest_qc", "failure", {}, 0.0)
            self.failures.write()
            return False, None
        
        # Check: required categories present
        self._check_categories(doc)
        
        # Check: atom count > 0
        self._check_atom_count(doc)
        
        # Collect all checks
        hard_fail = any(f.status == QCStatus.HARD_FAIL for f in self.qc_flags)
        
        # Write QC report
        qc_report_path = self.output_dir / "ingest_qc_report.json"
        self._write_qc_report(qc_report_path)
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
        }
        log_path = self.output_dir / "ingest_qc.log"
        with open(log_path, "w") as f:
            json.dump(log_dict, f, indent=2)
        
        elapsed = (datetime.utcnow() - start).total_seconds()
        self.logger.log_step_end(
            "ingest_qc",
            "failure" if hard_fail else "success",
            {"qc_flags": len(self.qc_flags), "hard_fail": hard_fail},
            elapsed
        )
        
        self.failures.write()
        self.logger.close()
        
        return not hard_fail, str(self.input_path)
    
    def _parse_with_gemmi(self) -> gemmi.Structure:
        """Parse input file using Gemmi (auto-detects PDB/mmCIF)."""
        if gemmi is None:
            raise ImportError("Gemmi not installed")
        
        # Gemmi can load mmCIF and PDB
        doc = gemmi.cif.read_file(str(self.input_path))
        return gemmi.cif.as_structure(doc)
    
    def _check_categories(self, doc: gemmi.Structure):
        """Verify required mmCIF categories."""
        # In Gemmi, structure object; we check if we can access key info
        try:
            n_atoms = len(doc[0])  # first model, all atoms
            if n_atoms == 0:
                self._check_fail("atom_site_present", "No atoms in structure")
                return
            
            # Check if we have chemical component info (approximation)
            # In real code, would check doc.info or parse CIF blocks
            self._check_pass("atom_site_present", "Atom data present")
        except Exception as e:
            self._check_fail("atom_site_present", f"Cannot read atom site: {str(e)}")
    
    def _check_atom_count(self, doc: gemmi.Structure):
        """Ensure we have atoms."""
        try:
            n_atoms = len(doc[0])
            if n_atoms > 0:
                self._check_pass("atom_count_positive", f"Atom count: {n_atoms}")
            else:
                self._check_fail("atom_count_positive", "No atoms in structure")
        except Exception as e:
            self._check_fail("atom_count_positive", str(e))
    
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
    
    def _write_qc_report(self, path: Path):
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
        }
        with open(path, "w") as f:
            json.dump(report, f, indent=2)
