"""Thin joblib helpers for per-group parallel work.

Workers must not mutate a shared AnnData; pass subsets / arrays and collect
results on the main process. ``n_jobs=1`` stays fully sequential.
"""

from __future__ import annotations

import os
from collections.abc import Callable, Sequence
from typing import TypeVar

T = TypeVar("T")
R = TypeVar("R")


def resolve_n_jobs(n_jobs: int = 1) -> int:
    """Normalize n_jobs: 1 = sequential, -1 = all CPUs, 0 treated as 1."""
    if n_jobs is None or n_jobs == 0:
        return 1
    n = int(n_jobs)
    if n < 0:
        return max(1, os.cpu_count() or 1)
    return n


def parallel_map(
    fn: Callable[[T], R],
    items: Sequence[T],
    *,
    n_jobs: int = 1,
    backend: str = "loky",
) -> list[R]:
    """Map ``fn`` over ``items``; preserves input order.

    With ``n_jobs==1`` (or a single item) runs in-process sequentially so
    behavior matches pre-parallel code and pickling is avoided.
    """
    seq = list(items)
    if not seq:
        return []
    n = resolve_n_jobs(n_jobs)
    if n == 1 or len(seq) == 1:
        return [fn(x) for x in seq]
    from joblib import Parallel, delayed

    return Parallel(n_jobs=n, backend=backend)(delayed(fn)(x) for x in seq)


__all__ = ["parallel_map", "resolve_n_jobs"]
