"""
LPMO Pipeline Logging Utilities
Responsibility: Structured logging for all pipeline steps, artifact writing

Logging policy:
  - Every step writes a JSON log with status, timing, errors, outputs
  - All failures go to dedicated failure logs (not lost in general logs)
  - Artifact outputs explicitly registered
"""

import logging
import json
from pathlib import Path
from typing import Dict, Optional, Any, List
from datetime import datetime
import sys


class StructuredLogger:
    """Wrapper for structured JSON logging."""
    
    def __init__(self, name: str, log_dir: Path):
        """
        Args:
            name: logger identifier (e.g., "normalize_mmcif", "placer_runner")
            log_dir: output directory for logs
        """
        self.name = name
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        
        # Set up dual logging: console + JSON file
        self.logger = logging.getLogger(name)
        self.logger.setLevel(logging.DEBUG)
        
        # Console handler (INFO+)
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_formatter = logging.Formatter(
            '%(asctime)s [%(name)s] %(levelname)s: %(message)s'
        )
        console_handler.setFormatter(console_formatter)
        self.logger.addHandler(console_handler)
        
        # File handler (JSON)
        self.log_file = self.log_dir / f"{name}.jsonl"
        self.file_handler = open(self.log_file, "w")
    
    def log_step_start(self, step_name: str, input_summary: Dict):
        """Log the beginning of a pipeline step."""
        record = {
            "timestamp": datetime.utcnow().isoformat(),
            "event": "step_start",
            "step": step_name,
            "input": input_summary,
        }
        self._write_json_record(record)
        self.logger.info(f"[{step_name}] Starting...")
    
    def log_step_end(self, step_name: str, status: str, output_summary: Dict, elapsed_sec: float):
        """Log the completion of a pipeline step."""
        record = {
            "timestamp": datetime.utcnow().isoformat(),
            "event": "step_end",
            "step": step_name,
            "status": status,  # "success", "failure", "skip"
            "output": output_summary,
            "elapsed_sec": elapsed_sec,
        }
        self._write_json_record(record)
        level = logging.ERROR if status == "failure" else logging.INFO
        self.logger.log(level, f"[{step_name}] {status.upper()} (elapsed: {elapsed_sec:.2f}s)")
    
    def log_check(self, check_name: str, status: str, message: str, data: Optional[Dict] = None):
        """Log a single QC/validation check."""
        record = {
            "timestamp": datetime.utcnow().isoformat(),
            "event": "check",
            "check": check_name,
            "status": status,  # "pass", "soft_flag", "hard_fail"
            "message": message,
            "data": data,
        }
        self._write_json_record(record)
        level = logging.ERROR if status == "hard_fail" else (logging.WARNING if status == "soft_flag" else logging.INFO)
        self.logger.log(level, f"Check '{check_name}': {status} — {message}")
    
    def log_failure(self, reason: str, details: Dict):
        """Log a hard failure (dedicate to failure log)."""
        record = {
            "timestamp": datetime.utcnow().isoformat(),
            "event": "failure",
            "reason": reason,
            "details": details,
        }
        self._write_json_record(record)
        self.logger.error(f"FAILURE: {reason}\n  Details: {json.dumps(details, indent=2)}")
    
    def log_artifact(self, artifact_name: str, path: Path, description: str):
        """Register an output artifact."""
        record = {
            "timestamp": datetime.utcnow().isoformat(),
            "event": "artifact",
            "artifact": artifact_name,
            "path": str(path),
            "size_bytes": path.stat().st_size if path.exists() else 0,
            "description": description,
        }
        self._write_json_record(record)
        self.logger.info(f"Artifact '{artifact_name}' → {path}")
    
    def _write_json_record(self, record: Dict):
        """Write a single JSON record (one per line)."""
        self.file_handler.write(json.dumps(record) + "\n")
        self.file_handler.flush()
    
    def close(self):
        """Clean up file handles."""
        if self.file_handler:
            self.file_handler.close()


class FailureLog:
    """Dedicated log for all hard failures (for post-run analysis)."""
    
    def __init__(self, output_path: Path):
        """
        Args:
            output_path: path to failures JSON file
        """
        self.output_path = Path(output_path)
        self.failures: List[Dict] = []
    
    def record(self, step_name: str, pose_id: Optional[str], reason: str, details: Dict):
        """Record one failure."""
        failure_entry = {
            "timestamp": datetime.utcnow().isoformat(),
            "step": step_name,
            "pose_id": pose_id,
            "reason": reason,
            "details": details,
        }
        self.failures.append(failure_entry)
    
    def write(self):
        """Flush all failures to JSON."""
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.output_path, "w") as f:
            json.dump(self.failures, f, indent=2, default=str)


class ArtifactRegistry:
    """Track all output artifacts for reproducibility."""
    
    def __init__(self):
        self.artifacts: Dict[str, Dict] = {}
    
    def register(self, artifact_id: str, path: Path, description: str, artifact_type: str):
        """
        Args:
            artifact_id: unique identifier
            path: file path
            description: human-readable description
            artifact_type: "normalized.cif", "atom_map.tsv", "qc_report.json", etc.
        """
        self.artifacts[artifact_id] = {
            "path": str(path),
            "description": description,
            "type": artifact_type,
            "exists": path.exists(),
            "size_bytes": path.stat().st_size if path.exists() else None,
            "timestamp": datetime.utcnow().isoformat(),
        }
    
    def to_json(self) -> str:
        """Serialize registry to JSON."""
        return json.dumps(self.artifacts, indent=2, default=str)
