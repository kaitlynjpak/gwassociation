"""Query GraceDB for superevents to screen.

Wraps the ``ligo-gracedb`` client to return a list of superevent IDs matching a
time window and false-alarm-rate cut, optionally filtered to likely binary
black holes (BBH) by their ``p_astro`` classification.

For lensing screening the BBH filter is deliberately loose: events whose
classification files cannot be fetched are *kept*, because excluding them would
needlessly shrink the sample (they may be real BBHs simply missing a
``p_astro`` file).
"""

from __future__ import annotations

import ast
import importlib.util
from dataclasses import dataclass

# Defaults mirror the published analysis.
DEFAULT_TIME_WINDOW = "32 months ago .. now"
DEFAULT_FAR_THRESHOLD = 2.3e-5  # events / second
DEFAULT_BBH_THRESHOLD = 0.1     # loose p_astro(BBH) cut for screening


def _require_gracedb():
    if importlib.util.find_spec("ligo.gracedb") is None:  # pragma: no cover
        raise ImportError(
            "Querying GraceDB requires 'ligo-gracedb'. Install the screening "
            'extra with:\n    pip install -e ".[screen]"'
        )
    from ligo.gracedb.rest import GraceDb

    return GraceDb


@dataclass
class QueryResult:
    """Outcome of a GraceDB superevent query."""

    event_ids: list[str]
    n_before_classification: int
    n_kept_unclassified: int


def query_superevents(
    time_window: str = DEFAULT_TIME_WINDOW,
    far_threshold: float = DEFAULT_FAR_THRESHOLD,
    bbh_threshold: float | None = DEFAULT_BBH_THRESHOLD,
    client=None,
    progress=None,
) -> QueryResult:
    """Return superevent IDs matching the cuts.

    Parameters
    ----------
    time_window : str
        GraceDB ``created:`` expression, e.g. ``'32 months ago .. now'``.
    far_threshold : float
        Keep events with false alarm rate below this value (events / second).
    bbh_threshold : float or None
        If given, keep events with ``p_astro(BBH) >`` this value (events lacking
        a classification file are kept). If None, skip BBH filtering entirely.
    client : optional
        A pre-built ``GraceDb`` client (mainly for testing). One is created if
        not supplied.
    progress : callable, optional
        Wraps an iterable for a progress bar (e.g. ``tqdm``); identity if None.

    Returns
    -------
    QueryResult
    """
    if progress is None:
        def progress(iterable, **_kwargs):
            return iterable

    if client is None:
        client = _require_gracedb()()

    query_string = f"created: {time_window} far < {far_threshold}"
    event_ids = [se["superevent_id"] for se in client.superevents(query_string)]
    n_before = len(event_ids)

    if bbh_threshold is None:
        return QueryResult(event_ids, n_before, 0)

    kept: list[str] = []
    n_kept_unclassified = 0
    for se in progress(event_ids, desc="Filtering BBH events"):
        classification = _fetch_classification(client, se)
        if classification is None:
            kept.append(se)
            n_kept_unclassified += 1
            continue
        if classification.get("BBH", 0) > bbh_threshold:
            kept.append(se)

    return QueryResult(kept, n_before, n_kept_unclassified)


def _fetch_classification(client, superevent_id):
    """Return the parsed p_astro dict for an event, or None if unavailable."""
    for pipeline in ("gstlal.p_astro.json", "pycbc.p_astro.json"):
        try:
            response = client.files(superevent_id, pipeline)
            return ast.literal_eval(response.read().decode("UTF-8"))
        except Exception:
            continue
    return None
