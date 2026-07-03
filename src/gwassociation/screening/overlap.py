"""Sky-map reading and the GW--GW overlap integral.

The overlap integral quantifies the spatial overlap between two probability
distributions on the celestial sphere. It was introduced for GW--EM coincidence
in Singer et al. (2014), https://doi.org/10.3847/1538-4357/aabfd2, and is used
here for GW--GW lensing screening.

These helpers operate on multi-order (MOC/UNIQ) sky maps as served by GraceDB
(``bayestar.multiorder.fits``), which is what the published analysis used. They
require the optional ``[screen]`` dependencies (``ligo.skymap``,
``astropy-healpix``, ``healpy``); importing this module without them raises a
clear error pointing at the extra.
"""

from __future__ import annotations

import importlib.util

import numpy as np

_MISSING = [
    name
    for name in ("ligo.skymap", "astropy_healpix", "healpy")
    if importlib.util.find_spec(name) is None
]
if _MISSING:  # pragma: no cover - exercised only without the extra installed
    raise ImportError(
        "gwassociation.screening requires the optional screening dependencies "
        f"({', '.join(_MISSING)}). Install them with:\n"
        '    pip install -e ".[screen]"'
    )

import astropy_healpix as ah  # noqa: E402
import astropy.units as u  # noqa: E402
import healpy as hp  # noqa: E402
import ligo.skymap.io  # noqa: E402
from astropy.coordinates import SkyCoord  # noqa: E402


def read_skymap(filepath, moc=True):
    """Read a FITS sky map.

    Parameters
    ----------
    filepath : str
        Path to the FITS sky-map file (``.fits`` or ``.fits.gz``).
    moc : bool, optional
        If True (default), return the multi-order table with ``PROBDENSITY`` and
        ``UNIQ`` columns. If False, return a flat HEALPix probability array
        (the metadata returned by ``ligo.skymap`` is stripped for convenience).

    Returns
    -------
    table or ndarray
        MOC table when ``moc`` is True, otherwise a flat probability array.
    """
    result = ligo.skymap.io.read_sky_map(filepath, moc=moc)
    if not moc and isinstance(result, tuple):
        # ligo.skymap returns (array, metadata) for flattened maps.
        return result[0]
    return result


def skymap_overlap_integral(
    se_skymap,
    exttrig_skymap=(),
    se_skymap_uniq=(),
    ext_skymap_uniq=(),
    ra=None,
    dec=None,
    se_nested=True,
    ext_nested=True,
):
    """Compute the sky-map overlap integral between two localizations.

    Supports multi-order (UNIQ) sky maps, flattened (nested/ring) sky maps, and
    point sources specified by ``ra``/``dec``. Higher values indicate greater
    spatial overlap; the value is dimensionless.

    Parameters
    ----------
    se_skymap : array
        Probability-density array of the first sky map.
    exttrig_skymap : array, optional
        Probability-density array of the second sky map.
    se_skymap_uniq, ext_skymap_uniq : array, optional
        UNIQ pixel indices, supplied when the corresponding map is multi-order.
    ra, dec : float, optional
        Point-source coordinates in degrees (alternative to a second sky map).
    se_nested, ext_nested : bool, optional
        Pixel ordering of flattened maps (nested when True, ring when False).

    Returns
    -------
    float
        The overlap integral.
    """
    se_order = "nested" if se_nested or any(se_skymap_uniq) else "ring"
    ext_order = "nested" if ext_nested or any(ext_skymap_uniq) else "ring"

    se_skymap = np.abs(se_skymap)
    exttrig_skymap = np.abs(exttrig_skymap)

    if not any(exttrig_skymap) and not (ra is not None and dec is not None):
        raise ValueError("Please provide either external sky map or ra/dec coordinates")

    # --- Case 1: multi-order first sky map (UNIQ ordering) ---
    if any(se_skymap_uniq):
        level, ipix = ah.uniq_to_level_ipix(se_skymap_uniq)
        nsides = ah.level_to_nside(level)
        areas = ah.nside_to_pixel_area(nsides)
        ra_gw, dec_gw = ah.healpix_to_lonlat(ipix, nsides, order="nested")

        sky_prior = 1 / (4 * np.pi * u.sr)
        se_norm = np.sum(se_skymap * areas)

        # Case 1a: both maps multi-order.
        if any(ext_skymap_uniq):
            level_ext, ipix_ext = ah.uniq_to_level_ipix(ext_skymap_uniq)
            nsides_ext = ah.level_to_nside(level_ext)
            ra_ext, dec_ext = ah.healpix_to_lonlat(ipix_ext, nsides_ext, order=ext_order)

            c = SkyCoord(ra=ra_gw, dec=dec_gw)
            catalog = SkyCoord(ra=ra_ext, dec=dec_ext)
            ext_ind, _d2d, _d3d = c.match_to_catalog_sky(catalog)

            ext_norm = np.sum(exttrig_skymap * ah.nside_to_pixel_area(nsides_ext))

            return np.sum(
                se_skymap * areas * exttrig_skymap[ext_ind]
                / sky_prior / se_norm / ext_norm
            ).to(1).value

        # Case 1b: multi-order first map, point source.
        if ra is not None and dec is not None:
            c = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
            catalog = SkyCoord(ra=ra_gw, dec=dec_gw)
            ind, _d2d, _d3d = c.match_to_catalog_sky(catalog)
            return (se_skymap[ind] / u.sr / sky_prior / se_norm).to(1).value

        # Case 1c: multi-order first map, flattened second map.
        if any(exttrig_skymap):
            ext_nside = ah.npix_to_nside(len(exttrig_skymap))
            ext_ind = ah.lonlat_to_healpix(ra_gw, dec_gw, ext_nside, order=ext_order)
            ext_norm = np.sum(exttrig_skymap)

            return np.sum(
                se_skymap * areas * exttrig_skymap[ext_ind]
                / ah.nside_to_pixel_area(ext_nside)
                / sky_prior / se_norm / ext_norm
            ).to(1).value

    # --- Case 2: flattened first sky map ---
    else:
        # Case 2a: flattened map, point source.
        if ra is not None and dec is not None:
            se_nside = ah.npix_to_nside(len(se_skymap))
            ind = ah.lonlat_to_healpix(ra * u.deg, dec * u.deg, se_nside, order=se_order)
            se_norm = sum(se_skymap)
            return se_skymap[ind] * len(se_skymap) / se_norm

        # Case 2b: both maps flattened.
        if any(exttrig_skymap):
            if se_nested != ext_nested:
                raise ValueError("Sky maps must both use nested or ring ordering")

            nside_s = hp.npix2nside(len(se_skymap))
            nside_e = hp.npix2nside(len(exttrig_skymap))

            if nside_s > nside_e:
                exttrig_skymap = hp.ud_grade(
                    exttrig_skymap, nside_out=nside_s,
                    order_in=("NESTED" if ext_nested else "RING"),
                )
            else:
                se_skymap = hp.ud_grade(
                    se_skymap, nside_out=nside_e,
                    order_in=("NESTED" if se_nested else "RING"),
                )

            se_norm = se_skymap.sum()
            exttrig_norm = exttrig_skymap.sum()

            if se_norm > 0 and exttrig_norm > 0:
                return (
                    np.dot(se_skymap, exttrig_skymap)
                    / se_norm / exttrig_norm * len(se_skymap)
                )
            raise ValueError("At least one sky map has zero total probability")

    raise ValueError("Please provide valid GW and external sky map information")


def overlap_from_files(skymap_files):
    """Compute the overlap integral for a ``(path1, path2)`` pair of MOC maps.

    Returns ``nan`` if either map cannot be read or the integral fails, so this
    is safe to map across many pairs in parallel.
    """
    try:
        s1 = read_skymap(skymap_files[0], moc=True)
        s2 = read_skymap(skymap_files[1], moc=True)
        return skymap_overlap_integral(
            s1["PROBDENSITY"],
            exttrig_skymap=s2["PROBDENSITY"],
            se_skymap_uniq=s1["UNIQ"],
            ext_skymap_uniq=s2["UNIQ"],
        )
    except Exception:
        return float("nan")
