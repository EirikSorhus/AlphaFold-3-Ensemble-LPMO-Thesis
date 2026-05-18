"""
LPMO Pipeline: Work-Root Discovery
Responsibility: Traverse a structure_pipeline work/ directory and build a
WorkRootManifest describing every AF3/RF3 prediction artifact found.

Directory convention (discovered, not assumed):
    work/
      {TARGET}/                        # e.g. CEL6, STA8, NAG4
        {model}/                       # af3 | rf3  (boltz is excluded)
          runs/
            {run_id}/                  # numeric SLURM job ID
              {UNIPROT}_{TARGET}/      # e.g. B6EQJ6_CEL6
                {UNIPROT}_{TARGET}_model.cif   (top-level consolidated CIF)
                *.json                          (confidence JSON, if present)
                seed-{N}_sample-{M}/            (per-sample directories)
                  *.cif
                  *.json

The module discovers everything dynamically — no hardcoded target/uniprot lists.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import List, Optional, Tuple

from lpmo_pipeline.utils.data_models import (
    ALLOWED_MODELS,
    ModelEntry,
    RunEntry,
    SampleEntry,
    SampleStatus,
    TargetEntry,
    UniprotTargetEntry,
    WorkRootManifest,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Regex patterns
# ---------------------------------------------------------------------------
# Target names: uppercase letter prefix + numeric suffix (e.g. CEL6, STA8, NAG4)
TARGET_RE = re.compile(r"^[A-Z]{2,10}\d{1,4}$")

# Seed-sample directory: seed-{int}_sample-{int}
SEED_SAMPLE_RE = re.compile(r"^seed-(\d+)_sample-(\d+)$")

# Numeric run ID (SLURM job ID)
RUN_ID_RE = re.compile(r"^\d+$")

# Uniprot-target directory: {UNIPROT}_{TARGET}
# UniProt accession: 6 or 10 alphanumeric characters
UNIPROT_TARGET_RE = re.compile(r"^([A-Za-z0-9]{6,10})_([A-Z]{2,10}\d{1,4})$")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def discover_work_root(
    work_root: Path,
    *,
    af3_only: bool = False,
    latest_only: bool = False,
    include_targets: Tuple[str, ...] = (),
) -> WorkRootManifest:
    """Scan a work/ directory and return a fully-populated manifest.

    Args:
        work_root: Path to the ``work/`` directory
            (e.g. ``/cluster/.../structure_pipeline/work``).
        af3_only: If True, only discover AF3 predictions (skip RF3).
        latest_only: If True, only discover the run pointed to by the
            ``latest`` symlink in each model directory. If the symlink is
            absent or unusable, the highest numeric run ID is used.
        include_targets: Optional target-name allowlist. If provided,
            only matching target directories are traversed.

    Returns:
        WorkRootManifest with all discovered targets, models, runs,
        uniprot-target directories, seed/sample slots, and file paths.
    """
    work_root = Path(work_root)
    manifest = WorkRootManifest(work_root=work_root)

    if not work_root.is_dir():
        manifest.errors.append(f"work_root is not a directory: {work_root}")
        return manifest

    models_filter = frozenset({"af3"}) if af3_only else ALLOWED_MODELS
    targets_filter = frozenset(include_targets)

    for child in sorted(work_root.iterdir()):
        if not child.is_dir():
            continue
        if not TARGET_RE.match(child.name):
            logger.debug("Skipping non-target directory: %s", child.name)
            continue
        if targets_filter and child.name not in targets_filter:
            continue
        target_entry = _discover_target(
            child, manifest.errors,
            models_filter=models_filter, latest_only=latest_only,
        )
        if target_entry.models:  # only include targets that have models
            manifest.targets.append(target_entry)

    return manifest


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------
def _discover_target(
    target_dir: Path,
    errors: List[str],
    *,
    models_filter: frozenset = ALLOWED_MODELS,
    latest_only: bool = False,
) -> TargetEntry:
    """Discover model directories under a target."""
    entry = TargetEntry(target=target_dir.name, directory=target_dir)

    for model_name in sorted(models_filter):
        model_dir = target_dir / model_name
        if not model_dir.is_dir():
            continue
        model_entry = _discover_model(
            model_dir, target_dir.name, errors, latest_only=latest_only,
        )
        entry.models.append(model_entry)

    return entry


def _discover_model(
    model_dir: Path,
    target: str,
    errors: List[str],
    *,
    latest_only: bool = False,
) -> ModelEntry:
    """Discover runs/ under a model directory.

    When *latest_only* is True the ``latest`` symlink in *model_dir* is
    resolved and only that single run is scanned. If the symlink does
    not exist or cannot be mapped onto the local runs/ tree, the highest
    numeric run directory is scanned instead.
    """
    entry = ModelEntry(model=model_dir.name, directory=model_dir)

    runs_dir = model_dir / "runs"
    if not runs_dir.is_dir():
        logger.debug("No runs/ directory in %s", model_dir)
        return entry

    # Resolve the ``latest`` symlink if requested.
    if latest_only:
        latest_link = model_dir / "latest"
        if latest_link.is_symlink() or latest_link.is_dir():
            resolved = latest_link.resolve()
            candidate_run_dir = _resolve_latest_run_dir(model_dir, resolved)
            if candidate_run_dir is not None:
                run_entry = _discover_run(candidate_run_dir, target, errors)
                entry.runs.append(run_entry)
                return entry
            else:
                logger.warning(
                    "latest link in %s resolves to %s which is not a "
                    "valid run directory — falling back to full scan",
                    model_dir, resolved,
                )
        else:
            logger.debug(
                "No 'latest' symlink in %s — scanning all runs", model_dir,
            )

    if latest_only:
        numeric_runs = [
            run_child
            for run_child in sorted(runs_dir.iterdir())
            if run_child.is_dir() and RUN_ID_RE.match(run_child.name)
        ]
        if numeric_runs:
            latest_run = max(numeric_runs, key=lambda path: int(path.name))
            entry.runs.append(_discover_run(latest_run, target, errors))
        return entry

    for run_child in sorted(runs_dir.iterdir()):
        if not run_child.is_dir():
            continue
        if not RUN_ID_RE.match(run_child.name):
            logger.debug("Skipping non-numeric run dir: %s", run_child.name)
            continue
        run_entry = _discover_run(run_child, target, errors)
        entry.runs.append(run_entry)

    return entry


def _resolve_latest_run_dir(model_dir: Path, resolved: Path) -> Optional[Path]:
    """Map a latest symlink target onto the local runs/ directory when possible."""
    if RUN_ID_RE.match(resolved.name):
        local_candidate = model_dir / "runs" / resolved.name
        if local_candidate.is_dir():
            return local_candidate
    return None


def _discover_run(
    run_dir: Path, target: str, errors: List[str]
) -> RunEntry:
    """Discover uniprot-target directories under a run."""
    entry = RunEntry(run_id=run_dir.name, directory=run_dir)

    for child in sorted(run_dir.iterdir()):
        if not child.is_dir():
            continue
        m = UNIPROT_TARGET_RE.match(child.name)
        if not m:
            # Skip non-matching dirs (e.g. "cache")
            logger.debug("Skipping non-uniprot-target dir: %s", child.name)
            continue
        uniprot_id, parsed_target = m.group(1), m.group(2)
        if parsed_target != target:
            errors.append(
                f"Target mismatch: directory {child} has target "
                f"{parsed_target!r} but parent target is {target!r}"
            )
        ut_entry = _discover_uniprot_target(child, uniprot_id, parsed_target, errors)
        entry.uniprot_targets.append(ut_entry)

    return entry


def _discover_uniprot_target(
    ut_dir: Path,
    uniprot_id: str,
    target: str,
    errors: List[str],
) -> UniprotTargetEntry:
    """Discover CIF, JSON, and seed/sample directories inside a uniprot-target dir."""
    entry = UniprotTargetEntry(
        uniprot_id=uniprot_id,
        target=target,
        directory=ut_dir,
    )

    # Scan top-level files for model CIF and confidence JSONs
    try:
        children = sorted(ut_dir.iterdir())
    except PermissionError:
        errors.append(f"Permission denied reading {ut_dir}")
        return entry

    for child in children:
        if child.is_file():
            if child.suffix == ".cif":
                # Take the first/only CIF at top level as the model CIF
                if entry.model_cif_path is None:
                    entry.model_cif_path = child
                else:
                    # Multiple CIFs at top level — unusual but record both
                    logger.warning(
                        "Multiple top-level CIFs in %s: %s (keeping first: %s)",
                        ut_dir, child.name, entry.model_cif_path.name,
                    )
            elif child.suffix == ".json":
                entry.confidence_json_paths.append(child)
        elif child.is_dir():
            m = SEED_SAMPLE_RE.match(child.name)
            if m:
                seed, sample = int(m.group(1)), int(m.group(2))
                sample_entry = _discover_sample(child, seed, sample)
                entry.samples.append(sample_entry)

    # Sort samples for deterministic ordering
    entry.samples.sort(key=lambda s: (s.seed, s.sample))

    return entry


def _discover_sample(sample_dir: Path, seed: int, sample: int) -> SampleEntry:
    """Classify a seed-sample directory based on its file contents."""
    entry = SampleEntry(seed=seed, sample=sample, directory=sample_dir)

    try:
        children = list(sample_dir.iterdir())
    except PermissionError:
        entry.status = SampleStatus.EMPTY_DIRECTORY
        return entry

    if not children:
        entry.status = SampleStatus.EMPTY_DIRECTORY
        return entry

    for child in sorted(children):
        if not child.is_file():
            continue
        if child.suffix == ".cif":
            entry.cif_paths.append(child)
        elif child.suffix == ".json":
            entry.confidence_json_paths.append(child)

    # Determine status
    has_cif = len(entry.cif_paths) > 0
    has_json = len(entry.confidence_json_paths) > 0

    if has_cif and has_json:
        entry.status = SampleStatus.COMPLETE
    elif has_cif and not has_json:
        entry.status = SampleStatus.MISSING_CONFIDENCE_JSON
    elif not has_cif and has_json:
        entry.status = SampleStatus.MISSING_CIF
    else:
        entry.status = SampleStatus.MISSING_ALL

    return entry


# ---------------------------------------------------------------------------
# Convenience: parse a single seed-sample dirname
# ---------------------------------------------------------------------------
def parse_seed_sample(name: str) -> Optional[Tuple[int, int]]:
    """Parse ``seed-{N}_sample-{M}`` → ``(N, M)`` or ``None``."""
    m = SEED_SAMPLE_RE.match(name)
    if m:
        return int(m.group(1)), int(m.group(2))
    return None
