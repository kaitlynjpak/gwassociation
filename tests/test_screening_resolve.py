"""Tests for event-ID -> sky-map resolution (no real network)."""
import pytest

from gwassociation.screening import resolve


class FakeResponse:
    def __init__(self, content=b"FITS", status_code=200):
        self.content = content
        self.status_code = status_code


class FakeSession:
    """Records URLs and serves canned content for known ones."""

    def __init__(self, available_urls=()):
        self.available = set(available_urls)
        self.requested = []

    def get(self, url):
        self.requested.append(url)
        if url in self.available:
            return FakeResponse(b"FITSDATA", 200)
        # GraceDB-style file endpoints used by download_skymaps.
        return FakeResponse(b"", 404)


def test_resolve_direct_file_path(tmp_path):
    f = tmp_path / "S1a.fits.gz"
    f.write_bytes(b"x")
    assert resolve.resolve_skymap(str(f)) == str(f)


def test_resolve_from_search_dir(tmp_path):
    (tmp_path / "S250727dc.fits.gz").write_bytes(b"x")
    got = resolve.resolve_skymap("S250727dc", search_dirs=[str(tmp_path)],
                                 cache_dir=str(tmp_path / "cache"))
    assert got.endswith("S250727dc.fits.gz")


def test_resolve_known_legacy_event_downloads(tmp_path):
    url = resolve.KNOWN_EVENT_URLS["GW170817"]
    session = FakeSession(available_urls=[url])
    got = resolve.resolve_skymap("GW170817", cache_dir=str(tmp_path), session=session)
    assert got.endswith("GW170817.fits.gz")
    assert (tmp_path / "GW170817.fits.gz").read_bytes() == b"FITSDATA"
    assert url in session.requested


def test_resolve_uses_cache_before_download(tmp_path):
    (tmp_path / "GW170817.fits.gz").write_bytes(b"CACHED")
    session = FakeSession()  # nothing available
    got = resolve.resolve_skymap("GW170817", cache_dir=str(tmp_path), session=session)
    assert got.endswith("GW170817.fits.gz")
    assert session.requested == []  # served from cache, no network


def test_resolve_gracedb_superevent(tmp_path):
    # download_skymaps tries the GraceDB API URL for the multiorder product.
    from gwassociation.screening.download import GRACEDB_API, SKYMAP_PRODUCTS
    url = f"{GRACEDB_API}/superevents/S250727dc/files/{SKYMAP_PRODUCTS[0]}"
    session = FakeSession(available_urls=[url])
    got = resolve.resolve_skymap("S250727dc", cache_dir=str(tmp_path), session=session)
    assert got.endswith("S250727dc.fits.gz")


def test_resolve_unknown_raises(tmp_path):
    session = FakeSession()  # nothing resolves
    with pytest.raises(FileNotFoundError):
        resolve.resolve_skymap("NotAnEvent", cache_dir=str(tmp_path), session=session)
