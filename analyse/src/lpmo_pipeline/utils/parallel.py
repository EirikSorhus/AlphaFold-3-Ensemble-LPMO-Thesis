# src/lpmo_pipeline/utils/parallel.py
"""
Responsibility: Thin wrapper around concurrent.futures for embarrassingly-parallel jobs.
Input:  A callable + list of arguments
Output: List of results (preserving order)
"""
from __future__ import annotations

import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, Callable, TypeVar

logger = logging.getLogger(__name__)

T = TypeVar("T")
R = TypeVar("R")


def parallel_map(
    func: Callable[[T], R],
    items: list[T],
    max_workers: int = 4,
    description: str = "parallel job",
) -> list[R | None]:
    """Run *func* over *items* in parallel using ProcessPoolExecutor.

    Returns results in the same order as *items*.  On per-item failure
    the slot is filled with ``None`` and the error is logged.

    Args:
        func: Picklable callable accepting one positional arg.
        items: Ordered list of inputs.
        max_workers: Number of parallel workers.
        description: For logging context.

    Returns:
        List of results (same length as *items*).
    """
    results: list[R | None] = [None] * len(items)
    if not items:
        return results

    logger.info("%s: submitting %d items to %d workers", description, len(items), max_workers)

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit with index so we can restore order
        future_to_idx = {
            executor.submit(func, item): idx for idx, item in enumerate(items)
        }

        for future in as_completed(future_to_idx):
            idx = future_to_idx[future]
            try:
                results[idx] = future.result()
            except Exception:
                logger.exception(
                    "%s: item %d failed (input=%r)", description, idx, items[idx]
                )
                results[idx] = None

    n_ok = sum(1 for r in results if r is not None)
    logger.info("%s: %d/%d succeeded", description, n_ok, len(items))
    return results
