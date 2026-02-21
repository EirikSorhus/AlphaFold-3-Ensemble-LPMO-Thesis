"""Bind-mount helpers for Apptainer/Singularity containers.

Computes minimal common-parent bind-mount prefixes from all referenced
file and directory paths, avoiding redundant overlapping mounts.
"""

from __future__ import annotations

import os
from pathlib import Path, PurePosixPath


def _common_prefix(a: PurePosixPath, b: PurePosixPath) -> PurePosixPath:
    """Return the longest common ancestor of two absolute paths."""
    parts_a = a.parts
    parts_b = b.parts
    common = []
    for pa, pb in zip(parts_a, parts_b):
        if pa == pb:
            common.append(pa)
        else:
            break
    if not common:
        return PurePosixPath("/")
    return PurePosixPath(*common)


def compute_bind_mounts(
    paths: list[str | Path],
    *,
    min_depth: int = 3,
) -> list[str]:
    """Compute minimal bind-mount prefixes for a list of paths.

    Groups paths that share common parent directories and returns the
    fewest ``--bind src:src`` entries needed to make all paths accessible.

    Args:
        paths: Absolute file / directory paths that the job needs to access.
        min_depth: Never collapse above this directory depth (default 3,
            which prevents mounting ``/`` or ``/cluster``).

    Returns:
        Sorted list of ``"src:src"`` bind-mount strings, one per unique
        prefix, ready for ``apptainer exec --bind ...``.
    """
    if not paths:
        return []

    # Resolve and deduplicate paths
    resolved: list[PurePosixPath] = []
    for p in paths:
        pp = PurePosixPath(os.path.realpath(str(p)))
        if not pp.is_absolute():
            continue
        resolved.append(pp)

    if not resolved:
        return []

    # Get the parent directory for each path (we bind dirs, not files)
    dirs: set[PurePosixPath] = set()
    for pp in resolved:
        # If the path ends with a suffix (looks like a file), use parent
        if pp.suffix:
            dirs.add(pp.parent)
        else:
            dirs.add(pp)

    # Sort dirs to process in order
    sorted_dirs = sorted(dirs, key=lambda d: d.parts)

    # Greedily merge dirs that share common prefixes at or above min_depth
    prefixes: list[PurePosixPath] = []
    for d in sorted_dirs:
        merged = False
        for i, existing in enumerate(prefixes):
            cp = _common_prefix(existing, d)
            if len(cp.parts) >= min_depth:
                # Merge: replace with the common prefix
                prefixes[i] = cp
                merged = True
                break
        if not merged:
            prefixes.append(d)

    # Deduplicate: remove prefixes that are children of other prefixes
    final: list[PurePosixPath] = []
    prefixes_sorted = sorted(prefixes, key=lambda p: len(p.parts))
    for p in prefixes_sorted:
        if not any(p != f and _is_subpath(p, f) for f in final):
            final.append(p)

    return sorted(f"{p}:{p}" for p in final)


def _is_subpath(child: PurePosixPath, parent: PurePosixPath) -> bool:
    """Check if child is under parent."""
    try:
        child.relative_to(parent)
        return True
    except ValueError:
        return False


def format_bind_flags(mounts: list[str]) -> str:
    """Format bind-mounts as apptainer --bind flags.

    Args:
        mounts: List of ``"src:dst"`` strings from :func:`compute_bind_mounts`.

    Returns:
        Multi-line string of ``--bind "src:dst"`` flags joined with `` \\\\\\n``.
    """
    if not mounts:
        return ""
    lines = [f'  --bind "{m}"' for m in mounts]
    return " \\\n".join(lines)


def log_mount_plan(mounts: list[str], paths: list[str | Path]) -> str:
    """Generate a human-readable mount plan for dry-run logging.

    Args:
        mounts: List of bind-mount strings.
        paths: Original paths that were analysed.

    Returns:
        Multi-line description string.
    """
    lines = ["Mount plan:"]
    lines.append(f"  Input paths ({len(paths)}):")
    for p in sorted(str(pp) for pp in paths):
        lines.append(f"    - {p}")
    lines.append(f"  Bind mounts ({len(mounts)}):")
    for m in mounts:
        lines.append(f"    --bind {m}")
    return "\n".join(lines)
