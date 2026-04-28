"""Bounded process-based parallel execution."""

from __future__ import annotations

import os
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, TypeVar

DEFAULT_AUTO_WORKERS: int = 4

T = TypeVar("T")


def resolve_worker_count(workers: int | None) -> int:
    if workers is not None and workers > 0:
        return workers
    cpu = os.cpu_count() or 1
    return min(cpu, DEFAULT_AUTO_WORKERS)


def build_chunks(items: list[Any], n_chunks: int) -> list[list[Any]]:
    if n_chunks <= 1 or not items:
        return [items]
    size = max(1, len(items) // n_chunks)
    chunks = [items[i : i + size] for i in range(0, len(items), size)]
    return chunks


def iter_chunk_results(
    fn: Callable[..., T],
    chunks: list[list[Any]],
    workers: int,
    *,
    extra_kwargs: dict[str, Any] | None = None,
) -> Iterator[T]:
    kwargs = extra_kwargs or {}
    if workers <= 1 or len(chunks) <= 1:
        for chunk in chunks:
            yield fn(chunk, **kwargs)
        return
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(fn, chunk, **kwargs): chunk for chunk in chunks}
        for future in as_completed(futures):
            yield future.result()
