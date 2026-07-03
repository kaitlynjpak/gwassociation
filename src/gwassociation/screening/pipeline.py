"""All-pairs overlap computation over a directory of sky maps.

Generates every unique event pair from the cached sky maps, computes the
overlap integral for each (optionally in parallel), and persists the results to
JSON so the expensive step need not be repeated.
"""

from __future__ import annotations

import json
import math
import os
from itertools import combinations

from .download import SKYMAP_EXT, list_cached_skymaps
from .overlap import overlap_from_files

# Pair names join the two event IDs with a comma, matching the published JSON.
PAIR_SEP = ","


def _event_name(filename: str) -> str:
    """Strip sky-map extensions from a filename to recover the event ID."""
    for ext in (".fits.gz", ".fits"):
        if filename.endswith(ext):
            return filename[: -len(ext)]
    return filename


def generate_pairs(skymap_dir: str):
    """Return ``(pair_paths, pair_names)`` for every unique pair in a directory.

    ``pair_paths`` is a list of ``[path1, path2]`` and ``pair_names`` the
    matching ``"event1,event2"`` strings. Pairs are unordered and exclude
    self-pairs.
    """
    files = list_cached_skymaps(skymap_dir)
    pair_paths = []
    pair_names = []
    for map_i, map_j in combinations(files, 2):
        pair_paths.append([
            os.path.join(skymap_dir, map_i),
            os.path.join(skymap_dir, map_j),
        ])
        pair_names.append(f"{_event_name(map_i)}{PAIR_SEP}{_event_name(map_j)}")
    return pair_paths, pair_names


def compute_overlaps(skymap_dir: str, workers: int = 1, progress=None) -> dict:
    """Compute overlap integrals for all unique pairs in ``skymap_dir``.

    Parameters
    ----------
    skymap_dir : str
        Directory of cached MOC sky maps.
    workers : int
        Number of worker processes. ``1`` runs serially (no pool overhead).
    progress : callable, optional
        Wraps an iterable for a progress bar (e.g. ``tqdm``); identity if None.

    Returns
    -------
    dict
        Mapping ``"event1,event2" -> overlap`` (``nan`` for failed pairs).
    """
    if progress is None:
        def progress(iterable, **_kwargs):
            return iterable

    pair_paths, pair_names = generate_pairs(skymap_dir)
    if not pair_paths:
        return {}

    if workers and workers > 1:
        values = _compute_parallel(pair_paths, workers)
    else:
        values = [
            overlap_from_files(paths)
            for paths in progress(pair_paths, desc="Computing overlaps")
        ]

    return dict(zip(pair_names, values))


def _compute_parallel(pair_paths, workers):
    """Compute overlaps with a process pool, preserving input order."""
    from multiprocessing import Pool

    with Pool(processes=workers) as pool:
        return pool.map(overlap_from_files, pair_paths)


def save_overlaps(overlaps: dict, path: str) -> None:
    """Write the overlaps mapping to JSON, encoding ``nan`` as null."""
    safe = {
        name: (None if value is None or math.isnan(value) else float(value))
        for name, value in overlaps.items()
    }
    with open(path, "w") as handle:
        json.dump(safe, handle)


def load_overlaps(path: str) -> dict:
    """Load an overlaps JSON, decoding null back to ``nan``."""
    with open(path) as handle:
        raw = json.load(handle)
    return {
        name: (float("nan") if value is None else float(value))
        for name, value in raw.items()
    }


def expected_pair_count(skymap_dir: str) -> int:
    """Return n*(n-1)/2 for the number of cached sky maps in a directory."""
    n = len([
        f for f in os.listdir(skymap_dir) if f.endswith(SKYMAP_EXT)
    ]) if os.path.isdir(skymap_dir) else 0
    return n * (n - 1) // 2
