"""Spatial overlap calculations for point and map localizations."""
from __future__ import annotations

import numpy as np

from ..utils import healpix as hp_utils


def _prob_array(skymap_or_event) -> np.ndarray:
    if isinstance(skymap_or_event, dict):
        return np.asarray(skymap_or_event.get("prob", skymap_or_event.get("data")), dtype=float)
    if hasattr(skymap_or_event, "prob") and skymap_or_event.prob is not None:
        return np.asarray(skymap_or_event.prob, dtype=float)
    if hasattr(skymap_or_event, "skymap"):
        if isinstance(skymap_or_event.skymap, dict):
            return np.asarray(skymap_or_event.skymap.get("prob", skymap_or_event.skymap.get("data")), dtype=float)
        if skymap_or_event.skymap is not None:
            return np.asarray(skymap_or_event.skymap, dtype=float)
    return np.asarray(skymap_or_event, dtype=float)


def _nside(skymap_or_event, prob: np.ndarray) -> int:
    if isinstance(skymap_or_event, dict) and skymap_or_event.get("nside") is not None:
        return int(skymap_or_event["nside"])
    if hasattr(skymap_or_event, "nside") and skymap_or_event.nside is not None:
        return int(skymap_or_event.nside)
    return hp_utils.npix2nside(prob.size)


def _as_skymap_dict(skymap_or_dict):
    prob = _prob_array(skymap_or_dict)
    return {"prob": prob, "data": prob, "nside": _nside(skymap_or_dict, prob)}


def _normalize_probabilities(prob_array):
    arr = np.asarray(prob_array, dtype=float)
    total = np.sum(arr)
    if total <= 0:
        raise ValueError("Skymap probabilities must sum to a positive value.")
    return arr / total


def IOmega_point(ra_rad, dec_rad, gw_prob, nside, nest=True):
    """Point-source sky Bayes factor relative to an isotropic sky prior."""
    prob = _normalize_probabilities(gw_prob)
    theta = np.pi / 2 - dec_rad
    phi = ra_rad
    ipix = hp_utils.ang2pix(nside, theta, phi, nest=nest)
    pix_area = hp_utils.nside2pixarea(nside)
    p_density = prob[ipix] / pix_area
    sky_prior = 1.0 / (4.0 * np.pi)
    return float(p_density / sky_prior)


def IOmega_maps(gw_prob, ext_prob, nside):
    """Map-vs-map sky Bayes factor relative to isotropic sky priors."""
    gw = _normalize_probabilities(gw_prob)
    ext = _normalize_probabilities(ext_prob)
    return float(gw.size * np.sum(gw * ext))


def _point_overlap(gw_map, ra, dec, gw_nested=True):
    return IOmega_point(np.radians(ra), np.radians(dec), gw_map["prob"], gw_map["nside"], nest=gw_nested)


def _map_overlap(gw_map, ext_map):
    if gw_map["nside"] != ext_map["nside"]:
        raise ValueError("GW and external skymaps must share the same NSIDE for overlap computation.")
    return IOmega_maps(gw_map["prob"], ext_map["prob"], gw_map["nside"])


def skymap_overlap_integral(gw_skymap, ext_skymap=None, ra=None, dec=None, gw_nested=True, ext_nested=True):
    """Compute spatial overlap against a point source or another sky map."""
    gw_map = _as_skymap_dict(gw_skymap)
    if ext_skymap is not None:
        return _map_overlap(gw_map, _as_skymap_dict(ext_skymap))
    if ra is not None and dec is not None:
        return _point_overlap(gw_map, ra, dec, gw_nested=gw_nested)
    raise ValueError("Provide either (ra, dec) for point sources or ext_skymap for map overlap.")


class SpatialOverlap:
    """Calculate spatial overlap integral I_Ω."""

    @staticmethod
    def compute(gw_event, em_transient, search_radius: float = 0.0) -> float:
        if getattr(gw_event, "skymap", None) is None and hasattr(gw_event, "load_skymap"):
            gw_event.load_skymap()
        return skymap_overlap_integral(
            gw_event,
            ra=em_transient.ra,
            dec=em_transient.dec,
            gw_nested=getattr(gw_event, "nest", True),
        )

    @staticmethod
    def compute_map_overlap(primary_event, secondary_event) -> float:
        if getattr(primary_event, "skymap", None) is None and hasattr(primary_event, "load_skymap"):
            primary_event.load_skymap()
        if getattr(secondary_event, "skymap", None) is None and hasattr(secondary_event, "load_skymap"):
            secondary_event.load_skymap()
        return skymap_overlap_integral(primary_event, ext_skymap=secondary_event)
