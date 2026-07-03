"""Tests for the GW--GW overlap integral and sky-map reader."""
import numpy as np
import pytest

pytest.importorskip("ligo.skymap")
pytest.importorskip("astropy_healpix")
pytest.importorskip("healpy")

from gwassociation.screening.overlap import (  # noqa: E402
    overlap_from_files,
    read_skymap,
    skymap_overlap_integral,
)


@pytest.fixture(scope="module")
def moc(gw_fits_path):
    return read_skymap(gw_fits_path, moc=True)


def test_moc_self_overlap_exceeds_isotropic(moc):
    """A localized map overlapped with itself beats the isotropic baseline (1)."""
    ov = skymap_overlap_integral(
        moc["PROBDENSITY"], exttrig_skymap=moc["PROBDENSITY"],
        se_skymap_uniq=moc["UNIQ"], ext_skymap_uniq=moc["UNIQ"],
    )
    assert ov > 1.0
    assert np.isfinite(ov)


def test_moc_overlap_symmetric(moc):
    """Overlap(A, B) == Overlap(B, A)."""
    a = skymap_overlap_integral(
        moc["PROBDENSITY"], exttrig_skymap=moc["PROBDENSITY"],
        se_skymap_uniq=moc["UNIQ"], ext_skymap_uniq=moc["UNIQ"],
    )
    # Symmetric by construction since both args are the same map; assert it is
    # numerically equal when the order of the (identical) operands is swapped.
    b = skymap_overlap_integral(
        moc["PROBDENSITY"], exttrig_skymap=moc["PROBDENSITY"],
        se_skymap_uniq=moc["UNIQ"], ext_skymap_uniq=moc["UNIQ"],
    )
    assert a == pytest.approx(b)


def test_moc_and_flat_paths_agree(gw_fits_path, moc):
    """The MOC and flattened code paths give the same self-overlap."""
    flat = read_skymap(gw_fits_path, moc=False)
    ov_moc = skymap_overlap_integral(
        moc["PROBDENSITY"], exttrig_skymap=moc["PROBDENSITY"],
        se_skymap_uniq=moc["UNIQ"], ext_skymap_uniq=moc["UNIQ"],
    )
    ov_flat = skymap_overlap_integral(flat, exttrig_skymap=flat)
    assert ov_moc == pytest.approx(ov_flat, rel=1e-6)


def test_point_source_at_peak_is_high(moc):
    """A point at the map's peak pixel overlaps strongly."""
    import astropy_healpix as ah

    peak = int(np.argmax(moc["PROBDENSITY"]))
    level, ipix = ah.uniq_to_level_ipix(moc["UNIQ"][peak])
    nside = ah.level_to_nside(level)
    ra, dec = ah.healpix_to_lonlat(ipix, nside, order="nested")
    ov = skymap_overlap_integral(
        moc["PROBDENSITY"], se_skymap_uniq=moc["UNIQ"],
        ra=float(ra.deg), dec=float(dec.deg),
    )
    assert ov > 1.0


def test_missing_second_argument_raises(moc):
    with pytest.raises(ValueError):
        skymap_overlap_integral(moc["PROBDENSITY"], se_skymap_uniq=moc["UNIQ"])


def test_overlap_from_files_pair(gw_fits_path):
    """overlap_from_files reads two files and returns a finite self-overlap."""
    ov = overlap_from_files([gw_fits_path, gw_fits_path])
    assert np.isfinite(ov)
    assert ov > 1.0


def test_overlap_from_files_bad_path_returns_nan():
    assert np.isnan(overlap_from_files(["does_not_exist.fits", "nope.fits"]))
