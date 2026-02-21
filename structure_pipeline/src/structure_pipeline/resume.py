"""Resume logic and status tracking.

Directory convention: work/{ligand_ccd}/{model}/
  - status.jsonl  (append-only event log)
  - DONE.ok       (completion marker)
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Any

from .cases import Case, CaseStatus


@dataclass
class StatusEntry:
    """A single status entry in the status log."""

    timestamp: str
    case_id: str
    protein_id: str
    ligand_id: str
    ligand_ccd_code: str
    model: str
    event: str  # started, completed, failed
    message: str = ""
    metadata: dict[str, Any] | None = None

    def to_json_line(self) -> str:
        """Convert to JSON line."""
        return json.dumps(asdict(self))

    @classmethod
    def from_json_line(cls, line: str) -> "StatusEntry":
        """Parse from JSON line."""
        data = json.loads(line)
        return cls(**data)


class StatusTracker:
    """Track job status with append-only log files.

    Path layout::

        work_dir/{ligand_ccd}/{model}/status.jsonl
        work_dir/{ligand_ccd}/{model}/DONE.ok
    """

    def __init__(self, work_dir: Path):
        self.work_dir = work_dir
        self.work_dir.mkdir(parents=True, exist_ok=True)

    def _get_status_file(self, ligand_ccd: str, model: str) -> Path:
        """Get path to status file for a ligand/model combination."""
        return self.work_dir / ligand_ccd / model / "status.jsonl"

    def _get_done_marker(self, ligand_ccd: str, model: str) -> Path:
        """Get path to DONE.ok marker."""
        return self.work_dir / ligand_ccd / model / "DONE.ok"

    def log_event(
        self,
        case: Case,
        event: str,
        message: str = "",
        metadata: dict[str, Any] | None = None,
    ) -> None:
        """Log a status event for a case."""
        status_file = self._get_status_file(case.ligand_ccd_code, case.model)
        status_file.parent.mkdir(parents=True, exist_ok=True)

        entry = StatusEntry(
            timestamp=datetime.now().isoformat(),
            case_id=case.case_id,
            protein_id=case.protein_id,
            ligand_id=case.ligand_id,
            ligand_ccd_code=case.ligand_ccd_code,
            model=case.model,
            event=event,
            message=message,
            metadata=metadata,
        )

        with open(status_file, "a") as f:
            f.write(entry.to_json_line() + "\n")

    def get_completed_cases(self, ligand_ccd: str, model: str) -> set[str]:
        """Get set of completed case IDs for a ligand/model."""
        completed = set()
        status_file = self._get_status_file(ligand_ccd, model)

        if not status_file.exists():
            return completed

        with open(status_file) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    entry = StatusEntry.from_json_line(line)
                    if entry.event == "completed":
                        completed.add(entry.case_id)
                    elif entry.event == "failed":
                        completed.discard(entry.case_id)
                except Exception:
                    continue

        return completed

    def is_ligand_model_done(self, ligand_ccd: str, model: str) -> bool:
        """Check if all cases for a ligand/model are done."""
        return self._get_done_marker(ligand_ccd, model).exists()

    def mark_ligand_model_done(self, ligand_ccd: str, model: str) -> None:
        """Mark a ligand/model group as fully completed."""
        done_marker = self._get_done_marker(ligand_ccd, model)
        done_marker.parent.mkdir(parents=True, exist_ok=True)
        done_marker.touch()

    def get_pending_cases(
        self,
        cases: list[Case],
        resume: bool = True,
    ) -> list[Case]:
        """Filter cases to only those that need to be run."""
        if not resume:
            return [c for c in cases if c.status != CaseStatus.SKIPPED.value]

        pending = []
        for case in cases:
            if case.status == CaseStatus.SKIPPED.value:
                continue

            # Check if this ligand/model group is fully done
            if self.is_ligand_model_done(case.ligand_ccd_code, case.model):
                continue

            # Check if this specific case is completed
            completed = self.get_completed_cases(case.ligand_ccd_code, case.model)
            if case.case_id in completed:
                continue

            pending.append(case)

        return pending

    def generate_summary(self) -> dict[str, Any]:
        """Generate a summary of all job statuses."""
        summary: dict[str, Any] = {
            "total_ligands": 0,
            "by_model": {},
            "by_status": {},
        }

        ligands = set()

        for ligand_dir in self.work_dir.iterdir():
            if not ligand_dir.is_dir():
                continue

            ligands.add(ligand_dir.name)

            for model_dir in ligand_dir.iterdir():
                if not model_dir.is_dir():
                    continue

                model = model_dir.name
                if model not in summary["by_model"]:
                    summary["by_model"][model] = {
                        "completed": 0,
                        "failed": 0,
                        "pending": 0,
                    }

                if (model_dir / "DONE.ok").exists():
                    summary["by_model"][model]["completed"] += 1
                elif (model_dir / "status.jsonl").exists():
                    summary["by_model"][model]["pending"] += 1

        summary["total_ligands"] = len(ligands)
        return summary
