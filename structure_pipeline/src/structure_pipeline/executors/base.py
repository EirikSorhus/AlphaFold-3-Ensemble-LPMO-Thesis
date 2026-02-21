"""Base executor interface."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any


class JobStatus(str, Enum):
    """Status of a submitted job."""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"
    UNKNOWN = "unknown"


@dataclass
class SubmittedJob:
    """A submitted job."""

    job_id: str  # SLURM job ID or local process ID
    job_name: str
    work_dir: Path
    script_path: Path
    submitted_at: datetime = field(default_factory=datetime.now)
    status: JobStatus = JobStatus.PENDING
    exit_code: int | None = None
    error_message: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def done_marker(self) -> Path:
        """Path to DONE.ok marker file."""
        return self.work_dir / "DONE.ok"

    @property
    def is_done(self) -> bool:
        """Check if job completed successfully."""
        return self.done_marker.exists()


class ExecutorInterface(ABC):
    """Abstract base class for job executors."""

    @abstractmethod
    def submit(
        self,
        script_content: str,
        work_dir: Path,
        job_name: str,
        dependency: str | None = None,
    ) -> SubmittedJob:
        """Submit a job for execution.

        Args:
            script_content: Full script content to execute
            work_dir: Working directory for the job
            job_name: Name for the job
            dependency: Optional dependency (e.g., afterok:12345)

        Returns:
            SubmittedJob object
        """
        pass

    @abstractmethod
    def get_status(self, job: SubmittedJob) -> JobStatus:
        """Get current status of a submitted job.

        Args:
            job: Previously submitted job

        Returns:
            Current job status
        """
        pass

    @abstractmethod
    def cancel(self, job: SubmittedJob) -> bool:
        """Cancel a running job.

        Args:
            job: Job to cancel

        Returns:
            True if cancellation was successful
        """
        pass

    def wait(
        self,
        job: SubmittedJob,
        timeout: float | None = None,
        poll_interval: float = 5.0,
    ) -> JobStatus:
        """Wait for a job to complete.

        Args:
            job: Job to wait for
            timeout: Maximum time to wait (None for infinite)
            poll_interval: Seconds between status checks

        Returns:
            Final job status
        """
        import time

        start = time.time()
        while True:
            status = self.get_status(job)
            if status in (
                JobStatus.COMPLETED,
                JobStatus.FAILED,
                JobStatus.CANCELLED,
            ):
                return status

            if timeout and (time.time() - start) > timeout:
                return JobStatus.UNKNOWN

            time.sleep(poll_interval)
