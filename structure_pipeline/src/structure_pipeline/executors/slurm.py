"""SLURM executor for HPC job submission."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

from .base import ExecutorInterface, JobStatus, SubmittedJob


class SlurmExecutor(ExecutorInterface):
    """Execute jobs via SLURM sbatch.

    Supports job dependencies for multi-stage workflows.
    """

    def __init__(self, dry_run: bool = False):
        """Initialize SLURM executor.

        Args:
            dry_run: If True, don't actually submit jobs
        """
        self.dry_run = dry_run

    def submit(
        self,
        script_content: str,
        work_dir: Path,
        job_name: str,
        dependency: str | None = None,
    ) -> SubmittedJob:
        """Submit a job via sbatch.

        Args:
            script_content: Full SLURM script content
            work_dir: Working directory for the job
            job_name: Name for the job
            dependency: Optional dependency (e.g., afterok:12345)

        Returns:
            SubmittedJob object with SLURM job ID
        """
        work_dir.mkdir(parents=True, exist_ok=True)

        # Write script to file
        script_path = work_dir / f"{job_name}.sh"
        with open(script_path, "w") as f:
            f.write(script_content)
        script_path.chmod(0o755)

        job = SubmittedJob(
            job_id="",  # Will be filled after submission
            job_name=job_name,
            work_dir=work_dir,
            script_path=script_path,
        )

        if self.dry_run:
            job.job_id = f"dry_run_{job_name}"
            job.status = JobStatus.PENDING
            job.metadata["dry_run"] = True
            return job

        # Build sbatch command
        cmd = ["sbatch"]
        if dependency:
            cmd.extend(["--dependency", dependency])
        cmd.append(str(script_path))

        # Submit job
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                cwd=work_dir,
            )

            # Parse job ID from output: "Submitted batch job 12345"
            match = re.search(r"Submitted batch job (\d+)", result.stdout)
            if match:
                job.job_id = match.group(1)
                job.status = JobStatus.PENDING
            else:
                job.status = JobStatus.FAILED
                job.error_message = f"Could not parse job ID from: {result.stdout}"

        except subprocess.CalledProcessError as e:
            job.status = JobStatus.FAILED
            job.error_message = f"sbatch failed: {e.stderr}"

        return job

    def get_status(self, job: SubmittedJob) -> JobStatus:
        """Get current status of a SLURM job.

        Uses sacct for completed jobs, squeue for running jobs.
        """
        if job.metadata.get("dry_run"):
            return JobStatus.PENDING

        # First check DONE.ok marker (faster than sacct)
        if job.is_done:
            return JobStatus.COMPLETED

        # Try squeue first (for running/pending jobs)
        try:
            result = subprocess.run(
                [
                    "squeue",
                    "-j", job.job_id,
                    "-h",  # No header
                    "-o", "%T",  # State only
                ],
                capture_output=True,
                text=True,
            )

            if result.returncode == 0 and result.stdout.strip():
                state = result.stdout.strip().upper()
                state_map = {
                    "PENDING": JobStatus.PENDING,
                    "RUNNING": JobStatus.RUNNING,
                    "COMPLETING": JobStatus.RUNNING,
                    "COMPLETED": JobStatus.COMPLETED,
                    "FAILED": JobStatus.FAILED,
                    "CANCELLED": JobStatus.CANCELLED,
                    "TIMEOUT": JobStatus.FAILED,
                    "NODE_FAIL": JobStatus.FAILED,
                }
                return state_map.get(state, JobStatus.UNKNOWN)

        except Exception:
            pass

        # Job not in queue, try sacct for historical info
        try:
            result = subprocess.run(
                [
                    "sacct",
                    "-j", job.job_id,
                    "-n",  # No header
                    "-o", "State",
                    "-X",  # No sub-jobs
                ],
                capture_output=True,
                text=True,
            )

            if result.returncode == 0 and result.stdout.strip():
                state = result.stdout.strip().split()[0].upper()
                if "COMPLETED" in state:
                    return JobStatus.COMPLETED
                elif "FAILED" in state or "TIMEOUT" in state:
                    return JobStatus.FAILED
                elif "CANCELLED" in state:
                    return JobStatus.CANCELLED

        except Exception:
            pass

        return JobStatus.UNKNOWN

    def cancel(self, job: SubmittedJob) -> bool:
        """Cancel a SLURM job."""
        if job.metadata.get("dry_run"):
            return True

        try:
            subprocess.run(
                ["scancel", job.job_id],
                capture_output=True,
                check=True,
            )
            return True
        except Exception:
            return False

    def submit_with_dependency_chain(
        self,
        scripts: list[tuple[str, Path, str]],
    ) -> list[SubmittedJob]:
        """Submit multiple jobs with sequential dependencies.

        Args:
            scripts: List of (script_content, work_dir, job_name) tuples

        Returns:
            List of SubmittedJob objects
        """
        jobs = []
        prev_job_id = None

        for script_content, work_dir, job_name in scripts:
            dependency = f"afterok:{prev_job_id}" if prev_job_id else None
            job = self.submit(script_content, work_dir, job_name, dependency)
            jobs.append(job)
            prev_job_id = job.job_id

        return jobs
