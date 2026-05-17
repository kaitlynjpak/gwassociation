"""
Lightweight healpy wrapper with fallbacks using astropy-healpix so that the
package can run even when healpy is not installed in the current interpreter.
"""
from __future__ import annotations

import importlib
import importlib.util
import math
from pathlib import Path
from typing import Tuple

import numpy as np
import astropy.units as u
from astropy.io import fits

HEALPY_AVAILABLE = importlib.util.find_spec("healpy") is not None
_hp = importlib.import_module("healpy") if HEALPY_AVAILABLE else None
ASTROPY_HEALPIX_AVAILABLE = importlib.util.find_spec("astropy_healpix") is not None
HEALPix = (
    importlib.import_module("astropy_healpix").HEALPix
    if ASTROPY_HEALPIX_AVAILABLE and not HEALPY_AVAILABLE
    else None
)


def _ensure_astropy_healpix():
    """Ensure the astropy-healpix fallback is importable.

    The package never installs dependencies at runtime; callers should install
    either ``healpy`` or ``astropy-healpix`` if HEALPix coordinate transforms
    are needed and healpy is unavailable.
    """
    global HEALPix
    if HEALPix is not None:
        return
    if importlib.util.find_spec("astropy_healpix") is None:  # pragma: no cover
        raise ModuleNotFoundError(
            "Install either 'healpy' or 'astropy-healpix>=0.7' to use HEALPix "
            "coordinate transforms without healpy."
        )
    HEALPix = importlib.import_module("astropy_healpix").HEALPix  # type: ignore


def _healpix_obj(nside: int, nest: bool = True) -> "HEALPix":
    if HEALPY_AVAILABLE:
        raise RuntimeError("HEALPix fallback requested even though healpy is available.")
    _ensure_astropy_healpix()
    order = "nested" if nest else "ring"
    return HEALPix(nside=nside, order=order, frame="icrs")  # type: ignore


def npix2nside(npix: int) -> int:
    if HEALPY_AVAILABLE:
        return int(_hp.npix2nside(npix))
    return int(math.sqrt(npix / 12.0))


def nside2npix(nside: int) -> int:
    if HEALPY_AVAILABLE:
        return int(_hp.nside2npix(nside))
    return 12 * int(nside) * int(nside)


def nside2pixarea(nside: int) -> float:
    if HEALPY_AVAILABLE:
        return float(_hp.nside2pixarea(nside))
    return 4.0 * math.pi / nside2npix(nside)


def ang2pix(nside: int, theta, phi, nest: bool = True):
    if HEALPY_AVAILABLE:
        return _hp.ang2pix(nside, theta, phi, nest=nest)
    hp = _healpix_obj(nside, nest)
    lon = np.asarray(phi) * u.rad
    lat = np.asarray(np.pi / 2 - np.asarray(theta)) * u.rad
    return hp.lonlat_to_healpix(lon, lat)


def pix2ang(nside: int, ipix, nest: bool = True):
    if HEALPY_AVAILABLE:
        return _hp.pix2ang(nside, ipix, nest=nest)
    hp = _healpix_obj(nside, nest)
    lon, lat = hp.healpix_to_lonlat(np.asarray(ipix))
    theta = (np.pi / 2) - lat.to_value(u.rad)
    phi = lon.to_value(u.rad)
    return theta, phi


def ang2vec(theta, phi):
    if HEALPY_AVAILABLE:
        return _hp.ang2vec(theta, phi)
    theta = np.asarray(theta, dtype=float)
    phi = np.asarray(phi, dtype=float)
    sin_theta = np.sin(theta)
    x = sin_theta * np.cos(phi)
    y = sin_theta * np.sin(phi)
    z = np.cos(theta)
    return np.array([x, y, z])


def pix2vec(nside: int, ipix, nest: bool = True):
    if HEALPY_AVAILABLE:
        return _hp.pix2vec(nside, ipix, nest=nest)
    theta, phi = pix2ang(nside, ipix, nest=nest)
    return ang2vec(theta, phi)


def read_map(filename: str, h: bool = False, verbose: bool = False):
    if HEALPY_AVAILABLE:
        return _hp.read_map(filename, h=h, verbose=verbose)
    data, header = _read_map_fallback(filename)
    return (data, header) if h else data


def write_map(filename: str, m, overwrite: bool = True):
    if HEALPY_AVAILABLE:
        return _hp.write_map(filename, m, overwrite=overwrite)
    data = np.asarray(m, dtype=float)
    primary = fits.PrimaryHDU(data)
    hdul = fits.HDUList([primary])
    hdul.writeto(filename, overwrite=overwrite)
    return filename


def mollview(*args, **kwargs):
    if not HEALPY_AVAILABLE:
        raise RuntimeError("healpy is required for sky-map plotting. Install healpy to enable mollview.")
    return _hp.mollview(*args, **kwargs)


def projscatter(*args, **kwargs):
    if not HEALPY_AVAILABLE:
        raise RuntimeError("healpy is required for projscatter plotting. Install healpy to enable this feature.")
    return _hp.projscatter(*args, **kwargs)


def graticule(*args, **kwargs):
    if not HEALPY_AVAILABLE:
        raise RuntimeError("healpy is required for graticule plotting. Install healpy to enable this feature.")
    return _hp.graticule(*args, **kwargs)


def _read_map_fallback(path: str) -> Tuple[np.ndarray, dict]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Skymap file not found: {path}")

    with fits.open(path) as hdul:
        header = dict(hdul[0].header)
        if len(hdul) > 1 and isinstance(hdul[1], fits.BinTableHDU):
            table = hdul[1].data
            colname = table.columns.names[0]
            data = np.asarray(table[colname]).ravel()
            header.update(dict(hdul[1].header))
        else:
            data = np.asarray(hdul[0].data, dtype=float).ravel()
    return data, header

