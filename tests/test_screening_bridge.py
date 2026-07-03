"""Tests for the GW--GW point-source bridge (event 2 as a pseudo transient)."""
import numpy as np
import pytest

pytest.importorskip("ligo.skymap")
pytest.importorskip("astropy_healpix")

from gwassociation.screening import bridge  # noqa: E402


def test_skymap_peak_returns_valid_coords(gw_fits_path):
    ra, dec, dist = bridge.skymap_peak(gw_fits_path)
    assert 0.0 <= ra <= 360.0
    assert -90.0 <= dec <= 90.0
    # The bundled fixture is a 3D map, so a positive distance is expected.
    assert dist is None or dist > 0


def test_event_to_transient_has_point_fields(gw_fits_path):
    t = bridge.event_to_transient(gw_fits_path, name="ev2")
    assert t["name"] == "ev2"
    assert "ra" in t and "dec" in t
    # 3D fixture -> redshift is derived by default and activates the distance term.
    assert "z" in t and t["z"] > 0
    # Only fields the Transient container accepts are present.
    assert "distance_mpc" not in t


def test_event_to_transient_distance_opt_out(gw_fits_path):
    t = bridge.event_to_transient(gw_fits_path, include_distance=False)
    assert "z" not in t


def test_event_to_transient_optional_times(gw_fits_path):
    t = bridge.event_to_transient(gw_fits_path, gw_time=100.0, event_time=110.0)
    assert t["gw_time"] == 100.0 and t["time"] == 110.0
    t2 = bridge.event_to_transient(gw_fits_path)
    assert "time" not in t2 and "gw_time" not in t2


def test_pair_point_odds_self_pair(gw_fits_path):
    """Event 2 at its own peak lands in a high-probability region of event 1."""
    r = bridge.pair_point_odds(gw_fits_path, gw_fits_path, name1="A", name2="B")
    assert r["event1"] == "A" and r["event2"] == "B"
    for key in ("I_omega", "I_dl", "log_posterior_odds", "confidence"):
        assert key in r
    # Peak sits in an above-isotropic-density region -> spatial Bayes factor > 1.
    assert r["I_omega"] > 1.0
    assert np.isfinite(r["I_omega"])
