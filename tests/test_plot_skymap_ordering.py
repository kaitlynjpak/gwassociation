"""Regression tests for ``plot_skymap``: ordering consistency and backends.

The GW probability map is stored in NESTED HEALPix ordering. The contour
renderer must be told that ordering, otherwise the localization is scrambled and
no longer coincides with the transient marker (which is projected straight from
RA/Dec through the same WCS). These tests lock in the ordering plumbing and the
``ligo.skymap`` -> healpy -> basic fallback chain without requiring the heavy
plotting dependencies to be installed.
"""
from types import SimpleNamespace

import numpy as np
import matplotlib

matplotlib.use("Agg")

from gwassociation.plotting import skymap as skymap_mod
from gwassociation.plotting.skymap import (
    _greedy_credible_levels,
    _skymap_is_nested,
    plot_skymap,
)


def test_ligo_backend_receives_nested_ordering(monkeypatch, tmp_path):
    seen = {}
    monkeypatch.setattr(
        skymap_mod, "_plot_ligo_style",
        lambda fig, prob, nest, transient, levels: seen.__setitem__("nest", nest),
    )
    gw = SimpleNamespace(skymap={"prob": np.ones(12), "nside": 1, "nest": True})
    plot_skymap(gw, SimpleNamespace(ra=197.45, dec=-23.38), out_file=str(tmp_path / "s.png"))
    assert seen["nest"] is True


def test_ligo_backend_receives_ring_ordering(monkeypatch, tmp_path):
    seen = {}
    monkeypatch.setattr(
        skymap_mod, "_plot_ligo_style",
        lambda fig, prob, nest, transient, levels: seen.__setitem__("nest", nest),
    )
    gw = SimpleNamespace(skymap={"prob": np.ones(12), "nside": 1, "nest": False})
    plot_skymap(gw, SimpleNamespace(ra=10.0, dec=20.0), out_file=str(tmp_path / "s.png"))
    assert seen["nest"] is False


def test_falls_back_to_healpy_when_ligo_missing(monkeypatch, tmp_path):
    seen = {}

    def raise_import(*a, **k):
        raise ImportError("no ligo.skymap")

    monkeypatch.setattr(skymap_mod, "_plot_ligo_style", raise_import)
    monkeypatch.setattr(
        skymap_mod, "_plot_healpy_style",
        lambda fig, prob, nest, transient: seen.__setitem__("nest", nest),
    )
    gw = SimpleNamespace(skymap={"prob": np.ones(12), "nside": 1, "nest": True})
    plot_skymap(gw, SimpleNamespace(ra=1.0, dec=2.0), out_file=str(tmp_path / "s.png"))
    assert seen["nest"] is True


def test_basic_fallback_when_no_backend(monkeypatch, tmp_path):
    def raise_import(*a, **k):
        raise ImportError("no ligo.skymap")

    def raise_runtime(*a, **k):
        raise RuntimeError("no healpy")

    monkeypatch.setattr(skymap_mod, "_plot_ligo_style", raise_import)
    monkeypatch.setattr(skymap_mod, "_plot_healpy_style", raise_runtime)
    gw = SimpleNamespace(skymap={"prob": np.ones(12), "nside": 1, "nest": True})
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plot_skymap(gw, SimpleNamespace(ra=1.0, dec=2.0), out_file=str(tmp_path / "s.png"))
    assert (tmp_path / "s.png").exists()


def test_greedy_credible_levels_orders_by_probability():
    cl = _greedy_credible_levels(np.array([0.7, 0.2, 0.1]))
    # Highest-probability pixel is enclosed by the smallest credible level.
    assert np.argmin(cl) == 0
    np.testing.assert_allclose(np.sort(cl), [0.7, 0.9, 1.0])


def test_greedy_credible_levels_handles_empty_map():
    cl = _greedy_credible_levels(np.zeros(12))
    np.testing.assert_array_equal(cl, np.ones(12))


def test_skymap_is_nested_prefers_dict_flag():
    gw = SimpleNamespace(skymap={"prob": np.ones(12), "nest": False}, nest=True)
    assert _skymap_is_nested(gw) is False


def test_skymap_is_nested_falls_back_to_event_attribute():
    gw = SimpleNamespace(skymap=None, nest=False)
    assert _skymap_is_nested(gw) is False


def test_skymap_is_nested_defaults_to_nested():
    gw = SimpleNamespace()
    assert _skymap_is_nested(gw) is True
