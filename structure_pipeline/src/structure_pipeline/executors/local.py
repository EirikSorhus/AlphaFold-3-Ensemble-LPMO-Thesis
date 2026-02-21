"""Local executor for running jobs as subprocesses."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

from .base import ExecutorInterface, JobStatus, SubmittedJob


class LocalExecutor(ExecutorInterface):
    """Execute jobs locally as subprocesses.

    Useful for testing and debugging without SLURM.
    Note: Does not support dependencies between jobs.
    """

    def __init__(self, dry_run: bool = False):
        """Initialize local executor.

        Args:
            dry_run: If True, don't actually run jobs
        """
        self.dry_run = dry_run
        self._processes: dict[str, subprocess.Popen] = {}
        self._job_counter = 0

    def submit(
        self,
        script_content: str,
        work_dir: Path,
        job_name: str,
        dependency: str | None = None,
    ) -> SubmittedJob:
        """Submit a job for local execution.

        Args:
            script_content: Full script content to execute
            work_dir: Working directory for the job
            job_name: Name for the job
            dependency: Ignored for local executor

        Returns:
            SubmittedJob object
        """
        work_dir.mkdir(parents=True, exist_ok=True)

        # Write script to file
        script_path = work_dir / f"{job_name}.sh"
        with open(script_path, "w") as f:
            # Remove SBATCH directives for local execution
            lines = script_content.split("\n")
            filtered = [
                line for line in lines if not line.strip().startswith("#SBATCH")
            ]
            f.write("\n".join(filtered))
        script_path.chmod(0o755)

        # Generate job ID
        self._job_counter += 1
        job_id = f"local_{self._job_counter}"

        job = SubmittedJob(
            job_id=job_id,
            job_name=job_name,
            work_dir=work_dir,
            script_path=script_path,
        )

        if self.dry_run:
            job.status = JobStatus.PENDING
            job.metadata["dry_run"] = True
            return job

        # Execute script
        stdout_path = work_dir / f"{job_name}_stdout.log"
        stderr_path = work_dir / f"{job_name}_stderr.log"

        with open(stdout_path, "w") as stdout_f, open(stderr_path, "w") as stderr_f:
            process = subprocess.Popen(
                ["bash", str(script_path)],
                cwd=work_dir,
                stdout=stdout_f,
                stderr=stderr_f,
                env={**os.environ, "SLURM_JOB_ID": job_id},
            )

        self._processes[job_id] = process
        job.status = JobStatus.RUNNING

        return job

    def get_status(self, job: SubmittedJob) -> JobStatus:
        """Get current status of a local job."""
        if job.metadata.get("dry_run"):
            return JobStatus.PENDING

        process = self._processes.get(job.job_id)
        if process is None:
            # Check if DONE.ok exists (job completed before we tracked it)
            if job.is_done:
                return JobStatus.COMPLETED
            return JobStatus.UNKNOWN

        poll_result = process.poll()
        if poll_result is None:
            return JobStatus.RUNNING

        job.exit_code = poll_result
        if poll_result == 0:
            return JobStatus.COMPLETED
        return JobStatus.FAILED

    def cancel(self, job: SubmittedJob) -> bool:
        """Cancel a local job."""
        process = self._processes.get(job.job_id)
        if process is None:
            return False

        try:
            process.terminate()
            process.wait(timeout=5)
            return True
        except subprocess.TimeoutExpired:
            process.kill()
            return True
        except Exception:
            return False
