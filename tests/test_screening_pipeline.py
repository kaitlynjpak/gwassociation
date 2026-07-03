"""Tests for pair generation and overlaps persistence (no network/FITS)."""
import math

from gwassociation.screening import pipeline


def _make_skymaps(tmp_path, names):
    for name in names:
        (tmp_path / f"{name}.fits.gz").write_bytes(b"stub")
    return str(tmp_path)


def test_generate_pairs_unique(tmp_path):
    skymap_dir = _make_skymaps(tmp_path, ["S1a", "S2b", "S3c"])
    paths, pair_names = pipeline.generate_pairs(skymap_dir)
    assert len(paths) == 3  # 3 choose 2
    assert set(pair_names) == {"S1a,S2b", "S1a,S3c", "S2b,S3c"}
    # No self-pairs, no reversed duplicates.
    assert all(a != b for a, b in (p.split(",") for p in pair_names))


def test_expected_pair_count(tmp_path):
    skymap_dir = _make_skymaps(tmp_path, ["S1a", "S2b", "S3c", "S4d"])
    assert pipeline.expected_pair_count(skymap_dir) == 6


def test_save_load_overlaps_roundtrip(tmp_path):
    overlaps = {"A,B": 12.5, "A,C": float("nan"), "B,C": 0.0}
    path = tmp_path / "ov.json"
    pipeline.save_overlaps(overlaps, str(path))
    loaded = pipeline.load_overlaps(str(path))
    assert loaded["A,B"] == 12.5
    assert loaded["B,C"] == 0.0
    assert math.isnan(loaded["A,C"])  # nan survives as null -> nan


def test_event_name_strips_extensions():
    assert pipeline._event_name("S250611f.fits.gz") == "S250611f"
    assert pipeline._event_name("S250611f.fits") == "S250611f"
