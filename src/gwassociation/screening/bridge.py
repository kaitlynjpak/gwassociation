"""Treat a second GW event as a point source for the GW--EM odds path.

The point-source association path (a GW sky map vs. a point-like counterpart) is
the validated, physically-scaled odds calculation. This module lets it run on
*any pair of GW events* by turning the second event's sky map into a pseudo
point-source "transient":

* its **peak sky position** (RA, Dec), and
* for 3D maps, a **redshift** derived from the conditional distance at that peak,
  which activates the distance-overlap term.

This is the GW--GW "lensing bridge" described in the analysis notebooks: it
compares event 2's most-likely position/distance against event 1's full 3D
localization, reusing the point-source machinery rather than the (separate,
un-normalized) skymap-vs-skymap path.
"""

from __future__ import annotations

import numpy as np

from .overlap import read_skymap


def skymap_peak(skymap_path):
    """Return ``(ra_deg, dec_deg, dist_mpc)`` at the peak of a MOC sky map.

    ``dist_mpc`` is the conditional distance mean (``DISTMU``) at the peak pixel
    for 3D maps, or ``None`` for 2D maps / non-finite entries.
    """
    import astropy_healpix as ah

    m = read_skymap(skymap_path, moc=True)
    peak = int(np.argmax(m["PROBDENSITY"]))
    level, ipix = ah.uniq_to_level_ipix(m["UNIQ"][peak])
    nside = ah.level_to_nside(level)
    ra, dec = ah.healpix_to_lonlat(ipix, nside, order="nested")

    dist = None
    if "DISTMU" in m.colnames:
        d = float(m["DISTMU"][peak])
        if np.isfinite(d) and d > 0:
            dist = d
    return float(ra.deg), float(dec.deg), dist


def distance_to_redshift(dist_mpc, cosmo=None):
    """Convert a luminosity distance (Mpc) to redshift via an Astropy cosmology."""
    import astropy.units as u
    from astropy.cosmology import Planck15, z_at_value

    cosmo = cosmo or Planck15
    return float(z_at_value(cosmo.luminosity_distance, dist_mpc * u.Mpc))


def event_to_transient(skymap_path, name=None, gw_time=None, event_time=None,
                       include_distance=True, cosmo=None):
    """Build a point-source ``transient_info`` dict from a GW event's sky map.

    ``include_distance`` (default True) adds a redshift derived from the peak
    distance so the distance-overlap term is active, matching the full GW170817
    breakdown. Set it False to score on sky position (and timing, if given)
    alone. ``time``/``gw_time`` are included only when supplied, so the temporal
    term is otherwise left inactive (I_t = 1) rather than forced with a
    meaningless delay.
    """
    ra, dec, dist = skymap_peak(skymap_path)
    transient = {"name": name or "event2", "ra": ra, "dec": dec}
    if include_distance and dist is not None:
        transient["z"] = distance_to_redshift(dist, cosmo)
    if event_time is not None:
        transient["time"] = event_time
    if gw_time is not None:
        transient["gw_time"] = gw_time
    return transient


def pair_point_odds(
    skymap1_path,
    skymap2_path,
    name1="event1",
    name2="event2",
    gw_time=None,
    event2_time=None,
    em_model="kilonova",
    include_distance=True,
    cosmo=None,
    **compute_kwargs,
):
    """Point-source association odds for a pair of GW events.

    Event 1's 3D sky map is the "GW event"; event 2 is reduced to a point source
    (peak position + distance-derived redshift) and scored with the validated
    point-source path. Returns the ``compute_odds`` result dict, augmented with
    the pseudo-transient under ``_transient`` and the event names.

    Provide ``gw_time`` (event 1) and ``event2_time`` to activate the temporal
    term; omit them to score on sky + distance overlap alone (recommended for
    long-delay lensing candidates, where EM light-curve timing does not apply).
    """
    from .. import Association

    transient = event_to_transient(
        skymap2_path, name=name2, gw_time=gw_time, event_time=event2_time,
        include_distance=include_distance, cosmo=cosmo,
    )
    assoc = Association(skymap1_path, transient)
    result = assoc.compute_odds(em_model=em_model, **compute_kwargs)
    result["event1"] = name1
    result["event2"] = name2
    result["_transient"] = transient
    return result
