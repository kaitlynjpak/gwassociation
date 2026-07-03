"""Empirical statistics and ranking for the overlap distribution.

Given the dictionary of pairwise overlaps, this module builds the empirical null
distribution, assigns each pair a p-value (the fraction of pairs with overlap at
least as large), applies multiple-comparison corrections (Bonferroni and
Benjamini--Hochberg FDR), and ranks candidate pairs -- optionally excluding
same-day pairs, which are typically duplicate pipeline detections of one trigger
rather than independent events.
"""

from __future__ import annotations

import datetime as dt
import math
from dataclasses import dataclass

import numpy as np


def valid_overlap_values(overlaps: dict) -> np.ndarray:
    """Return the strictly-positive, non-NaN overlap values as an array."""
    values = np.array([
        v for v in overlaps.values()
        if v is not None and not math.isnan(v) and v > 0
    ], dtype=float)
    return values


def describe(overlaps: dict) -> dict:
    """Return summary statistics of the overlap distribution.

    Keys: ``n``, ``mean``, ``median``, ``std``, ``min``, ``max`` and a
    ``percentiles`` sub-dict for the 90/95/99/99.9th percentiles.
    """
    values = valid_overlap_values(overlaps)
    if values.size == 0:
        return {"n": 0}
    return {
        "n": int(values.size),
        "mean": float(values.mean()),
        "median": float(np.median(values)),
        "std": float(values.std()),
        "min": float(values.min()),
        "max": float(values.max()),
        "percentiles": {
            str(p): float(np.percentile(values, p)) for p in (90, 95, 99, 99.9)
        },
    }


def empirical_pvalues(overlaps: dict) -> dict:
    """Return per-pair empirical p-values.

    The p-value of a pair is the fraction of all valid pairs whose overlap is
    greater than or equal to that pair's overlap. Invalid pairs map to ``nan``.
    """
    values = valid_overlap_values(overlaps)
    n_total = values.size
    if n_total == 0:
        return {name: float("nan") for name in overlaps}

    ordered = np.sort(values)
    names = list(overlaps.keys())
    raw = np.array(
        [overlaps[name] if overlaps[name] is not None else np.nan for name in names],
        dtype=float,
    )
    valid = np.isfinite(raw) & (raw > 0)

    # p-value = fraction of pairs with overlap >= this one. Computed for all
    # pairs at once via the insertion point on the sorted value array.
    pvals = np.full(raw.shape, np.nan)
    idx = np.searchsorted(ordered, raw[valid], side="left")
    pvals[valid] = (n_total - idx) / n_total
    return dict(zip(names, (float(p) for p in pvals)))


@dataclass
class Corrections:
    """Multiple-comparison correction summary."""

    n_tests: int
    alpha: float
    bonferroni_threshold: float
    n_bonferroni: int
    fdr_threshold: float
    n_fdr: int


def multiple_comparison_corrections(pvalues: dict, alpha: float = 0.05) -> Corrections:
    """Apply Bonferroni and Benjamini--Hochberg FDR corrections.

    Parameters
    ----------
    pvalues : dict
        Mapping of pair name to empirical p-value (NaNs ignored).
    alpha : float
        Family-wise / false-discovery target.
    """
    pval_array = np.array([
        p for p in pvalues.values() if p is not None and not math.isnan(p)
    ], dtype=float)
    n_tests = pval_array.size
    if n_tests == 0:
        return Corrections(0, alpha, 0.0, 0, 0.0, 0)

    bonferroni_threshold = alpha / n_tests
    n_bonferroni = int(np.sum(pval_array < bonferroni_threshold))

    sorted_pvals = np.sort(pval_array)
    m = sorted_pvals.size
    bh_thresholds = alpha * np.arange(1, m + 1) / m
    below = sorted_pvals <= bh_thresholds
    if np.any(below):
        k_max = int(np.max(np.where(below)[0]))
        n_fdr = k_max + 1
        fdr_threshold = float(bh_thresholds[k_max])
    else:
        n_fdr = 0
        fdr_threshold = 0.0

    return Corrections(
        n_tests=n_tests,
        alpha=alpha,
        bonferroni_threshold=bonferroni_threshold,
        n_bonferroni=n_bonferroni,
        fdr_threshold=fdr_threshold,
        n_fdr=n_fdr,
    )


def _event_date(event_id: str) -> dt.datetime:
    """Parse the YYMMDD date encoded in a GraceDB superevent ID (e.g. S250611f)."""
    return dt.datetime.strptime(event_id[1:7], "%y%m%d")


def days_apart(pair_name: str) -> int:
    """Return the absolute number of days between the two events in a pair."""
    a, b = pair_name.split(",")
    return abs((_event_date(a) - _event_date(b)).days)


def rank_pairs(
    overlaps: dict,
    pvalues: dict | None = None,
    top_n: int | None = None,
    min_days_apart: int = 0,
) -> list[dict]:
    """Rank pairs by overlap (descending).

    Parameters
    ----------
    overlaps : dict
        Mapping of pair name to overlap value.
    pvalues : dict, optional
        Per-pair empirical p-values to attach to each row.
    top_n : int, optional
        Truncate to this many pairs (after filtering). None keeps all.
    min_days_apart : int
        Drop pairs whose events are fewer than this many days apart. Use ``1``
        to exclude same-day duplicate detections of a single trigger.

    Returns
    -------
    list of dict
        Each row has ``rank``, ``event1``, ``event2``, ``overlap``,
        ``days_apart`` and (if available) ``pvalue``.
    """
    rows = []
    for name, ov in overlaps.items():
        if ov is None or math.isnan(ov) or ov <= 0:
            continue
        try:
            gap = days_apart(name)
        except (ValueError, IndexError):
            gap = None
        if min_days_apart and (gap is None or gap < min_days_apart):
            continue
        event1, event2 = name.split(",")
        row = {
            "event1": event1,
            "event2": event2,
            "overlap": float(ov),
            "days_apart": gap,
        }
        if pvalues is not None:
            row["pvalue"] = pvalues.get(name, float("nan"))
        rows.append(row)

    rows.sort(key=lambda r: r["overlap"], reverse=True)
    if top_n is not None:
        rows = rows[:top_n]
    for i, row in enumerate(rows, 1):
        row["rank"] = i
    return rows


def write_top_pairs_csv(rows: list[dict], path: str) -> None:
    """Write ranked pairs to CSV (stdlib only; no pandas dependency)."""
    import csv

    fieldnames = ["rank", "event1", "event2", "overlap", "days_apart", "pvalue"]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
