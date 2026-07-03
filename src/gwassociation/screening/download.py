"""Download the best-available sky map for each superevent from GraceDB.

Sky maps are fetched over the public GraceDB HTTP API and cached on disk. The
preferred product is the multi-order ``bayestar.multiorder.fits`` (used by the
overlap integral); a flattened ``bayestar.fits.gz`` is accepted as a fallback.
Already-cached files are skipped so the step is resumable.
"""

from __future__ import annotations

import importlib.util
import os
from dataclasses import dataclass, field

GRACEDB_API = "https://gracedb.ligo.org/api"

# Preferred sky-map products, in priority order.
SKYMAP_PRODUCTS = ("bayestar.multiorder.fits", "bayestar.fits.gz")

SKYMAP_EXT = (".fits", ".fits.gz")


def _require_requests():
    if importlib.util.find_spec("requests") is None:  # pragma: no cover
        raise ImportError(
            "Downloading sky maps requires 'requests'. Install the screening "
            'extra with:\n    pip install -e ".[screen]"'
        )
    import requests

    return requests


@dataclass
class DownloadResult:
    """Outcome of a sky-map download batch."""

    downloaded: list[str] = field(default_factory=list)
    failed: list[str] = field(default_factory=list)
    skymap_dir: str = "."


def find_skymap_file(event_id: str, skymap_dir: str) -> str:
    """Return the cached sky-map path for an event, trying known extensions.

    Raises
    ------
    FileNotFoundError
        If no sky map for ``event_id`` exists in ``skymap_dir``.
    """
    for ext in SKYMAP_EXT:
        path = os.path.join(skymap_dir, f"{event_id}{ext}")
        if os.path.exists(path):
            return path
    raise FileNotFoundError(f"No skymap found for {event_id} in {skymap_dir}")


def list_cached_skymaps(skymap_dir: str) -> list[str]:
    """Return the sky-map filenames cached in ``skymap_dir`` (sorted)."""
    if not os.path.isdir(skymap_dir):
        return []
    return sorted(f for f in os.listdir(skymap_dir) if f.endswith(SKYMAP_EXT))


def download_skymaps(
    event_ids,
    skymap_dir: str,
    products=SKYMAP_PRODUCTS,
    session=None,
    progress=None,
) -> DownloadResult:
    """Download the best sky map for each event into ``skymap_dir``.

    Parameters
    ----------
    event_ids : iterable of str
        Superevent IDs to fetch.
    skymap_dir : str
        Destination directory (created if needed).
    products : sequence of str
        Candidate sky-map filenames to try per event, in priority order.
    session : optional
        A ``requests.Session``-like object (mainly for testing). One is created
        if not supplied.
    progress : callable, optional
        Wraps an iterable for a progress bar (e.g. ``tqdm``); identity if None.

    Returns
    -------
    DownloadResult
    """
    if progress is None:
        def progress(iterable, **_kwargs):
            return iterable

    if session is None:
        session = _require_requests().Session()

    os.makedirs(skymap_dir, exist_ok=True)
    result = DownloadResult(skymap_dir=skymap_dir)

    for se_id in progress(list(event_ids), desc="Downloading skymaps"):
        # Cached already (under any known extension)?
        try:
            find_skymap_file(se_id, skymap_dir)
            result.downloaded.append(se_id)
            continue
        except FileNotFoundError:
            pass

        output_path = os.path.join(skymap_dir, f"{se_id}.fits.gz")
        if _download_one(session, se_id, products, output_path):
            result.downloaded.append(se_id)
        else:
            result.failed.append(se_id)

    return result


def _download_one(session, se_id, products, output_path) -> bool:
    """Fetch the first available product for one event; return success."""
    for fname in products:
        url = f"{GRACEDB_API}/superevents/{se_id}/files/{fname}"
        try:
            response = session.get(url)
        except Exception:
            continue
        if getattr(response, "status_code", None) == 200:
            with open(output_path, "wb") as handle:
                handle.write(response.content)
            return True
    return False
