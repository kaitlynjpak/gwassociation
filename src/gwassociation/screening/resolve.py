"""Resolve a GW event identifier to a local sky-map file.

An identifier can be any of:

1. an existing sky-map file path (used as-is);
2. an event ID found in one of the local search directories;
3. a well-known legacy event name (e.g. ``GW170817``) fetched from its public
   archive URL;
4. a GraceDB superevent ID (e.g. ``S250727dc``) fetched via the GraceDB API.

Downloads are cached under ``cache_dir`` so repeat calls are fast.
"""

from __future__ import annotations

import os

from .download import (
    _require_requests,
    download_skymaps,
    find_skymap_file,
)

DEFAULT_CACHE = os.path.expanduser("~/.cache/gwassociation/skymaps")

# Legacy events that predate GraceDB superevents, with public sky-map URLs.
KNOWN_EVENT_URLS = {
    "GW170817": "https://dcc.ligo.org/public/0157/P1800381/007/GW170817_skymap.fits.gz",
}


def resolve_skymap(event, search_dirs=(), cache_dir=DEFAULT_CACHE, session=None, log=print):
    """Return a local sky-map path for ``event``, downloading + caching if needed.

    Parameters
    ----------
    event : str
        A file path, local event ID, legacy event name, or GraceDB superevent ID.
    search_dirs : iterable of str
        Local directories to look in before downloading.
    cache_dir : str
        Where downloads are written (and also searched).
    session : optional
        A ``requests.Session``-like object (mainly for testing).
    log : callable
        Progress sink.

    Raises
    ------
    FileNotFoundError
        If the event cannot be resolved by any route.
    """
    # 1. Direct file path.
    if os.path.isfile(event):
        return event

    # 2. Local search directories, then the cache.
    for directory in list(search_dirs) + [cache_dir]:
        if directory and os.path.isdir(directory):
            try:
                return find_skymap_file(event, directory)
            except FileNotFoundError:
                pass

    os.makedirs(cache_dir, exist_ok=True)

    # 3. Known legacy event (fetched from a public URL).
    if event in KNOWN_EVENT_URLS:
        dest = os.path.join(cache_dir, f"{event}.fits.gz")
        log(f"downloading {event} from public archive ...")
        _download_url(KNOWN_EVENT_URLS[event], dest, session)
        return dest

    # 4. GraceDB superevent.
    log(f"downloading {event} from GraceDB ...")
    result = download_skymaps([event], cache_dir, session=session)
    if event in result.downloaded:
        return find_skymap_file(event, cache_dir)

    raise FileNotFoundError(
        f"Could not resolve '{event}': not a local file, not found in the search "
        f"directories, not a known legacy event, and no GraceDB sky map was "
        f"available. Pass a file path or a valid GraceDB superevent ID."
    )


def _download_url(url, dest, session=None):
    """Download ``url`` to ``dest`` (raises FileNotFoundError on failure)."""
    if session is None:
        session = _require_requests().Session()
    response = session.get(url)
    if getattr(response, "status_code", None) != 200:
        raise FileNotFoundError(
            f"download failed for {url}: HTTP {getattr(response, 'status_code', '?')}"
        )
    with open(dest, "wb") as handle:
        handle.write(response.content)
    return dest
