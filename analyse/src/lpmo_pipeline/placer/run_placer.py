# src/lpmo_pipeline/placer/run_placer.py
"""
Responsibility: Run PLACER refinement to generate an ensemble of ligand poses.
Input:  normalized.cif (protein + ligand + Cu)
Output: placer_ensemble/ directory with 50–200 ranked poses

HARD RULE: PLACER must run BEFORE any final scoring / reporting.
GATE:     ≥1 pose returned; if 0 → PlacerZeroPosesError (hard fail)
"""
from __future__ import annotations

import json
import logging
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from lpmo_pipeline.utils.exceptions import PlacerZeroPosesError

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Result model
# ---------------------------------------------------------------------------
@dataclass
class PlacerPose:
    """Single PLACER-generated pose."""

    pose_id: str
    pdb_path: Path
    placer_score: float = 0.0
    clash_count: int = 0
    cu_geom_sanity: bool = True
    rank: int = 0


@dataclass
class PlacerResult:
    """Full PLACER run result."""

    run_id: str = ""
    n_modes: int = 0
    poses: list[PlacerPose] = field(default_factory=list)
    best_pose: PlacerPose | None = None
    output_dir: Path = Path(".")


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------
def run_placer(
    input_cif: Path,
    output_dir: Path,
    n_modes: int = 200,
    run_id: str = "",
    placer_executable: str = "placer",
    extra_args: dict[str, Any] | None = None,
) -> PlacerResult:
    """Run PLACER on a normalized complex.

    Args:
        input_cif: Path to normalized.cif (protein + ligand + Cu).
        output_dir: Directory to write ensemble poses into.
        n_modes: Number of PLACER modes to sample (50–200).
        run_id: Identifier for logging context.
        placer_executable: Path or name of PLACER binary.
        extra_args: Additional CLI args as key:value dict.

    Returns:
        PlacerResult with all poses and references.

    Raises:
        PlacerZeroPosesError: If PLACER returns 0 poses.
    """
    logger.info("PLACER: running %d modes for %s → %s", n_modes, input_cif, output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Step 1: Build PLACER command ---
    cmd = [
        placer_executable,
        "--input", str(input_cif),
        "--output", str(output_dir),
        "--n_modes", str(n_modes),
    ]
    if extra_args:
        for k, v in extra_args.items():
            cmd.extend([f"--{k}", str(v)])

    # --- Step 2: Execute ---
    logger.info("PLACER cmd: %s", " ".join(cmd))
    # PSEUDOCODE:
    # proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    # if proc.returncode != 0:
    #     logger.error("PLACER failed: %s", proc.stderr)
    #     raise PlacerZeroPosesError(f"PLACER exited with code {proc.returncode}",
    #                                 run_id=run_id, step="placer_refine")

    # --- Step 3: Collect output PDBs ---
    pose_pdbs = sorted(output_dir.glob("pose_*.pdb"))
    if len(pose_pdbs) == 0:
        raise PlacerZeroPosesError(
            "PLACER returned 0 poses",
            run_id=run_id,
            step="placer_refine",
        )

    # --- Step 4: Parse PLACER scores ---
    scores_file = output_dir / "placer_scores.txt"
    scores = _parse_placer_scores(scores_file) if scores_file.exists() else {}

    # --- Step 5: Build PlacerPose list ---
    poses: list[PlacerPose] = []
    for rank, pdb_path in enumerate(pose_pdbs, 1):
        pose_id = pdb_path.stem
        score = scores.get(pose_id, 0.0)
        poses.append(PlacerPose(
            pose_id=pose_id,
            pdb_path=pdb_path,
            placer_score=score,
            rank=rank,
        ))

    # Sort by score (ascending = better)
    poses.sort(key=lambda p: p.placer_score)
    for rank, p in enumerate(poses, 1):
        p.rank = rank

    result = PlacerResult(
        run_id=run_id,
        n_modes=n_modes,
        poses=poses,
        best_pose=poses[0] if poses else None,
        output_dir=output_dir,
    )

    logger.info("PLACER: %d poses generated, best score=%.3f", len(poses),
                result.best_pose.placer_score if result.best_pose else float("nan"))

    # Write structured scores JSON
    _write_placer_scores_json(result, output_dir / "placer_scores.json")

    return result


def _parse_placer_scores(scores_file: Path) -> dict[str, float]:
    """Parse PLACER score output file.

    PSEUDOCODE: Format depends on actual PLACER output.
    Expected: tab-separated pose_id\\tscore per line.
    """
    scores: dict[str, float] = {}
    with open(scores_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                scores[parts[0]] = float(parts[1])
    return scores


def _write_placer_scores_json(result: PlacerResult, output_path: Path) -> None:
    """Write structured PLACER scores to JSON."""
    data = {
        "run_id": result.run_id,
        "n_modes": result.n_modes,
        "n_poses": len(result.poses),
        "poses": {
            p.pose_id: {
                "rank": p.rank,
                "placer_score": p.placer_score,
                "clash_count": p.clash_count,
                "cu_geom_sanity": p.cu_geom_sanity,
                "pdb_path": str(p.pdb_path),
            }
            for p in result.poses
        },
    }
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    logger.info("Wrote PLACER scores JSON to %s", output_path)
