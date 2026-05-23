"""Non-protonating analysis exports for ProLIF and convergence inputs."""
from __future__ import annotations

import json
import shutil
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional, Tuple

from lpmo_pipeline.io.cif_to_pdb import convert_cif_to_pdb
from lpmo_pipeline.io.protonate_export import _write_ligand_only_pdb
from lpmo_pipeline.utils.logging import FailureLog, StructuredLogger


@dataclass
class AnalysisExportReport:
    """Summary for non-protonated downstream analysis exports."""

    input_cif: str
    output_dir: str
    complex_for_prolif_pdb: str
    ligand_pdb: str
    ligand_only_residue_count: int
    ligand_only_skipped_residue_count: int
    ligand_only_skipped_residue_names: list[str]
    blockers: list[str]
    warnings: list[str]


class AnalysisExportRunner:
    """Write ProLIF-ready PDB artifacts without protonating the structure."""

    def __init__(self, input_cif: Path, output_dir: Path):
        self.input_cif = Path(input_cif)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logger = StructuredLogger("analysis_export", self.output_dir)
        self.failures = FailureLog(self.output_dir / "analysis_export_failures.json")

    def run(self) -> Tuple[bool, Optional[dict]]:
        """Write complex and ligand PDB artifacts for implicit-H ProLIF."""
        start = datetime.utcnow()
        self.logger.log_step_start("analysis_export", {"input": str(self.input_cif)})

        complex_pdb = self.output_dir / "complex_for_prolif.pdb"
        ligand_pdb = self.output_dir / "ligand_only_for_prolif.pdb"
        blockers: list[str] = []
        warnings: list[str] = []
        ligand_only_residue_count = 0
        ligand_only_skipped_residue_count = 0
        ligand_only_skipped_residue_names: list[str] = []

        try:
            with tempfile.TemporaryDirectory(prefix="complex_for_prolif_", dir=self.output_dir) as tmp_dir:
                ok_pdb, tmp_pdb = convert_cif_to_pdb(
                    self.input_cif,
                    Path(tmp_dir),
                    strip_metals_for_posebusters=False,
                )
                if not ok_pdb or tmp_pdb is None or not tmp_pdb.exists():
                    blockers.append("complex_pdb_export_failed")
                else:
                    shutil.copyfile(tmp_pdb, complex_pdb)

            try:
                ligand_only_residue_count = _write_ligand_only_pdb(
                    self.input_cif,
                    ligand_pdb,
                    warnings=warnings,
                )
                skipped_summary = getattr(_write_ligand_only_pdb, "last_skipped_residue_summary", {})
                ligand_only_skipped_residue_count = int(skipped_summary.get("count", 0))
                ligand_only_skipped_residue_names = list(skipped_summary.get("residue_names", []))
                if ligand_only_residue_count == 0:
                    blockers.append("ligand_only_extract_empty")
            except Exception as exc:
                blockers.append("ligand_only_extract_failed")
                warnings.append(f"ligand_only_extract_error: {exc}")

            report = AnalysisExportReport(
                input_cif=str(self.input_cif),
                output_dir=str(self.output_dir),
                complex_for_prolif_pdb=str(complex_pdb),
                ligand_pdb=str(ligand_pdb),
                ligand_only_residue_count=ligand_only_residue_count,
                ligand_only_skipped_residue_count=ligand_only_skipped_residue_count,
                ligand_only_skipped_residue_names=ligand_only_skipped_residue_names,
                blockers=blockers,
                warnings=warnings,
            )
            report_path = self.output_dir / "analysis_export_report.json"
            report_path.write_text(json.dumps(asdict(report), indent=2))

            if complex_pdb.exists():
                self.logger.log_artifact("complex_for_prolif_pdb", complex_pdb, "Non-protonated complex PDB for ProLIF")
            if ligand_pdb.exists():
                self.logger.log_artifact("ligand_pdb", ligand_pdb, "Non-protonated ligand PDB for ProLIF")
            self.logger.log_artifact("analysis_export_report", report_path, "Non-protonated analysis export report")

            success = len(blockers) == 0
            if not success:
                self.failures.record(
                    "analysis_export",
                    None,
                    "missing_or_failed_artifacts",
                    {"blockers": blockers, "warnings": warnings},
                )
                self.logger.log_failure("analysis_export_blocked", {"blockers": blockers})

            elapsed = (datetime.utcnow() - start).total_seconds()
            self.logger.log_step_end(
                "analysis_export",
                "success" if success else "failure",
                {"blockers": blockers},
                elapsed,
            )
            self.failures.write()
            self.logger.close()
            return success, asdict(report)
        except Exception as exc:
            self.logger.log_failure("analysis_export_exception", {"error": str(exc)})
            self.failures.record("analysis_export", None, "exception", {"error": str(exc)})
            self.failures.write()
            elapsed = (datetime.utcnow() - start).total_seconds()
            self.logger.log_step_end("analysis_export", "failure", {}, elapsed)
            self.logger.close()
            return False, None


def export_analysis_artifacts(input_cif: Path, output_dir: Path) -> Tuple[bool, Optional[dict]]:
    """Convenience wrapper for non-protonated downstream analysis exports."""
    return AnalysisExportRunner(input_cif=input_cif, output_dir=output_dir).run()
