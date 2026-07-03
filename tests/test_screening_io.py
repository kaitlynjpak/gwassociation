"""Tests for the GraceDB query and sky-map download with fakes (no network)."""
from gwassociation.screening import download as dl
from gwassociation.screening import gracedb


class FakeResponse:
    def __init__(self, payload, status_code=200):
        self._payload = payload
        self.status_code = status_code
        self.content = payload if isinstance(payload, bytes) else payload.encode()

    def read(self):
        return self.content


class FakeGraceClient:
    """Minimal stand-in for ligo.gracedb.rest.GraceDb."""

    def __init__(self, events, classifications):
        self._events = events
        self._classifications = classifications

    def superevents(self, query_string):  # noqa: D401 - mimic client API
        assert "far <" in query_string
        return ({"superevent_id": e} for e in self._events)

    def files(self, superevent_id, pipeline):
        key = (superevent_id, pipeline)
        if key not in self._classifications:
            raise FileNotFoundError(pipeline)
        return FakeResponse(str(self._classifications[key]))


def test_query_keeps_unclassified_and_filters_bbh():
    client = FakeGraceClient(
        events=["S1a", "S2b", "S3c"],
        classifications={
            ("S1a", "gstlal.p_astro.json"): {"BBH": 0.9, "Terrestrial": 0.1},
            ("S2b", "gstlal.p_astro.json"): {"BBH": 0.01, "BNS": 0.99},
            # S3c has no classification file -> kept as a possible BBH.
        },
    )
    res = gracedb.query_superevents(bbh_threshold=0.1, client=client)
    assert set(res.event_ids) == {"S1a", "S3c"}
    assert res.n_kept_unclassified == 1
    assert res.n_before_classification == 3


def test_query_can_disable_bbh_filter():
    client = FakeGraceClient(events=["S1a", "S2b"], classifications={})
    res = gracedb.query_superevents(bbh_threshold=None, client=client)
    assert res.event_ids == ["S1a", "S2b"]


class FakeSession:
    def __init__(self, available):
        # available: dict event_id -> filename that returns 200
        self._available = available
        self.requested = []

    def get(self, url):
        self.requested.append(url)
        parts = url.split("/")
        se_id, fname = parts[-3], parts[-1]
        if self._available.get(se_id) == fname:
            return FakeResponse(b"FITSDATA", status_code=200)
        return FakeResponse(b"", status_code=404)


def test_download_prefers_multiorder_then_falls_back(tmp_path):
    session = FakeSession({
        "S1a": "bayestar.multiorder.fits",
        "S2b": "bayestar.fits.gz",
    })
    res = dl.download_skymaps(["S1a", "S2b"], str(tmp_path), session=session)
    assert set(res.downloaded) == {"S1a", "S2b"}
    assert res.failed == []
    assert (tmp_path / "S1a.fits.gz").read_bytes() == b"FITSDATA"


def test_download_skips_cached(tmp_path):
    (tmp_path / "S1a.fits.gz").write_bytes(b"CACHED")
    session = FakeSession({})  # no successful downloads available
    res = dl.download_skymaps(["S1a"], str(tmp_path), session=session)
    assert res.downloaded == ["S1a"]
    assert session.requested == []  # never hit the network


def test_download_records_failures(tmp_path):
    session = FakeSession({})
    res = dl.download_skymaps(["S9z"], str(tmp_path), session=session)
    assert res.failed == ["S9z"]


def test_find_skymap_file_and_list(tmp_path):
    (tmp_path / "S1a.fits.gz").write_bytes(b"x")
    (tmp_path / "S2b.fits").write_bytes(b"y")
    assert dl.find_skymap_file("S1a", str(tmp_path)).endswith("S1a.fits.gz")
    assert dl.find_skymap_file("S2b", str(tmp_path)).endswith("S2b.fits")
    assert dl.list_cached_skymaps(str(tmp_path)) == ["S1a.fits.gz", "S2b.fits"]
