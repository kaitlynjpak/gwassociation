"""Sky-map plotting helper for GW localizations and point-like EM candidates.

The preferred rendering uses ``ligo.skymap``'s ``astro hours mollweide``
projection with credible-region contours (the publication style used by the
overlap-screening plots). The EM transient is drawn from an ICRS ``SkyCoord``,
so it passes through the same WCS transform as the contours and always lands on
the correct sky position. If ``ligo.skymap`` is unavailable the code falls back
to a healpy ``mollview`` render, and finally to a bare RA/Dec scatter.
"""

import warnings

import numpy as np
import matplotlib.pyplot as plt

from ..utils import healpix as hp_utils

# Contour colors from the innermost (smallest) credible level outward, matching
# the blue palette used elsewhere in the package.
_CREDIBLE_COLORS = {
    1: ["#2563eb"],
    2: ["#1d4ed8", "#60a5fa"],
    3: ["#1e3a8a", "#2563eb", "#60a5fa"],
}


def _skymap_is_nested(gw):
    """Return the HEALPix ordering flag for ``gw``'s probability map.

    GW sky maps in this package default to NESTED ordering (see
    :func:`gwassociation.io.skymap.load_gw_skymap`). The ordering may be exposed
    either on the sky-map dict (``gw.skymap['nest']``) or as an attribute on the
    event (``gw.nest``); fall back to NESTED when neither is present.
    """
    if isinstance(getattr(gw, "skymap", None), dict) and "nest" in gw.skymap:
        return bool(gw.skymap["nest"])
    if hasattr(gw, "nest") and gw.nest is not None:
        return bool(gw.nest)
    return True


def _extract_prob(gw):
    """Return the normalized HEALPix probability array for ``gw`` (or ``None``)."""
    skymap = None
    if getattr(gw, "skymap", None) is not None:
        skymap = gw.skymap
    elif getattr(gw, "data", None) is not None:
        return np.asarray(gw.data, dtype=float)
    elif hasattr(gw, "load_skymap"):
        gw.load_skymap()
        skymap = getattr(gw, "skymap", None)

    if skymap is None:
        return None
    if isinstance(skymap, dict):
        data = skymap.get("prob", skymap.get("data"))
        return np.asarray(data, dtype=float) if data is not None else None
    return np.asarray(skymap, dtype=float)


def _greedy_credible_levels(prob):
    """Map each pixel to the smallest credible level (0-1) that encloses it.

    Contours drawn at ``level`` then bound the smallest region containing that
    fraction of the total probability (e.g. the 90% credible region).
    """
    prob = np.asarray(prob, dtype=float).ravel()
    prob = np.where(np.isfinite(prob), prob, 0.0)
    total = prob.sum()
    if total <= 0:
        return np.ones_like(prob)
    prob = prob / total
    order = np.argsort(-prob)
    credible = np.empty_like(prob)
    credible[order] = np.cumsum(prob[order])
    return credible


def _contour_colors(n):
    if n in _CREDIBLE_COLORS:
        return _CREDIBLE_COLORS[n]
    cmap = plt.get_cmap("Blues")
    return [cmap(x) for x in np.linspace(0.9, 0.4, n)]


def _plot_ligo_style(fig, prob, nest, transient, levels):
    """Render the ligo.skymap astro-Mollweide contour plot onto ``fig``."""
    import astropy.units as u
    import matplotlib.lines as mlines
    import matplotlib.patches as mpatches
    import ligo.skymap.plot  # noqa: F401 -- registers the astro-mollweide projection
    from astropy.coordinates import SkyCoord

    credible = _greedy_credible_levels(prob)
    lvls = sorted(levels)
    colors = _contour_colors(len(lvls))
    linewidths = np.linspace(2.4, 1.4, len(lvls)).tolist()

    ax = fig.add_subplot(111, projection="astro hours mollweide")
    ax.grid()
    # ``nested=nest`` keeps the map's ordering aligned with the marker's WCS
    # transform; a mismatch here displaces the localization from the transient.
    ax.contour_hpx((credible, "ICRS"), levels=lvls, nested=nest,
                   colors=colors, linewidths=linewidths)

    handles = [
        mpatches.Patch(color=c, label=f"{int(round(level * 100))}% credible region")
        for c, level in zip(colors, lvls)
    ]

    if transient is not None and hasattr(transient, "ra") and hasattr(transient, "dec"):
        coord = SkyCoord(float(transient.ra) * u.deg, float(transient.dec) * u.deg)
        ax.plot_coord(coord, marker="o", color="gold", markersize=10,
                      markeredgecolor="black", markeredgewidth=1.0, linestyle="none")
        handles.append(
            mlines.Line2D([], [], marker="o", color="gold", markeredgecolor="black",
                          markersize=9, linestyle="none", label="EM transient")
        )

    ax.legend(handles=handles, loc="upper right", fontsize=11)
    ax.set_title("GW Localization with EM Candidate")


def _plot_healpy_style(fig, prob, nest, transient):
    """Fallback render using healpy ``mollview`` + ``projscatter``."""
    hp_utils.mollview(prob, title="GW Skymap with EM Candidate",
                      unit="Probability", fig=fig.number, nest=nest, cmap="YlOrRd")
    if transient is not None and hasattr(transient, "ra") and hasattr(transient, "dec"):
        hp_utils.projscatter(transient.ra, transient.dec, lonlat=True, marker="o",
                             s=120, color="blue", edgecolor="white", linewidth=2,
                             label="EM Transient")
    hp_utils.graticule(dpar=30, dmer=30, alpha=0.3)


def _plot_basic(fig, transient):
    """Last-resort RA/Dec scatter used when no HEALPix backend is available."""
    ax = fig.add_subplot(111)
    ax.set_title("Skymap (basic view)")
    if transient is not None and hasattr(transient, "ra") and hasattr(transient, "dec"):
        ax.plot(transient.ra, transient.dec, "o", markersize=12, color="blue",
                label="EM Transient")
        ax.set_xlabel("RA [deg]")
        ax.set_ylabel("Dec [deg]")
        ax.legend()
        ax.grid(True, alpha=0.3)


def plot_skymap(gw, transient, out_file="skymap.png", levels=(0.5, 0.9)):
    """Plot a GW localization with an EM candidate marker.

    Uses the ligo.skymap astro-Mollweide contour style when available, falling
    back to healpy ``mollview`` and then a bare scatter. ``levels`` sets the
    credible-region contours (defaults to the 50% and 90% regions).
    """
    fig = plt.figure(figsize=(12, 6))
    prob = _extract_prob(gw)
    nest = _skymap_is_nested(gw)

    try:
        if prob is None:
            raise ValueError("Sky-map probability data unavailable.")
        try:
            _plot_ligo_style(fig, prob, nest, transient, levels)
        except ImportError:
            # ligo.skymap (or a transitive dependency) is not installed.
            _plot_healpy_style(fig, prob, nest, transient)
    except Exception as e:
        warnings.warn(
            f"Detailed skymap plot unavailable ({e}); using basic fallback plot. "
            "Install the 'screen' extra (ligo.skymap + healpy) for full rendering.",
            RuntimeWarning,
        )
        fig.clf()
        _plot_basic(fig, transient)

    fig.savefig(out_file, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig
