"""Publication-style plots for the overlap screening pipeline.

Provides the overlap-distribution histogram, the empirical survival function,
and the joint-probability sky-map plot for a high-overlap pair. The sky-map plot
requires the optional ``[screen]`` dependencies (``ligo.skymap``, ``healpy``);
the histogram and survival function need only matplotlib + numpy.
"""

from __future__ import annotations

import numpy as np

from .statistics import valid_overlap_values


def plot_overlap_distribution(overlaps: dict, save_path: str, high_threshold: float = 100.0):
    """Two-panel histogram: full log distribution and the high-overlap tail.

    Mirrors Figure 1 of the poster. Returns the matplotlib Figure.
    """
    import matplotlib.pyplot as plt

    values = valid_overlap_values(overlaps)
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Left: full distribution on a log10 axis (drop float-noise zeros).
    overlap_for_plot = values[values > 1e-6]
    log10_overlaps = np.log10(overlap_for_plot)
    axes[0].hist(log10_overlaps, bins=60, alpha=0.7, color="steelblue",
                 edgecolor="black", linewidth=0.5)
    axes[0].set_xlim(-6, 3)
    axes[0].set_xlabel(r"$\log_{10}$(Overlap Integral)", fontsize=14)
    axes[0].set_ylabel("Number of Pairs", fontsize=14)
    axes[0].set_title("Full Distribution of Overlap Integrals", fontsize=16)
    axes[0].grid(True, alpha=0.3)

    # Right: high-overlap tail on a linear axis.
    high_overlaps = values[values > high_threshold]
    axes[1].hist(high_overlaps, bins=30, alpha=0.7, color="crimson",
                 edgecolor="black", linewidth=0.5)
    axes[1].set_xlabel("Overlap Integral", fontsize=14)
    axes[1].set_ylabel("Number of Pairs", fontsize=14)
    axes[1].set_title(
        f"High-Overlap Pairs (> {high_threshold:g})\n{len(high_overlaps):,} pairs",
        fontsize=16,
    )
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig


def plot_survival_function(overlaps: dict, save_path: str):
    """Empirical survival function (complementary CDF) on a log-log axis.

    Shows, for each overlap value, the fraction of pairs with overlap at least
    that large -- the empirical p-value curve. Returns the matplotlib Figure.
    """
    import matplotlib.pyplot as plt

    values = valid_overlap_values(overlaps)
    ordered = np.sort(values)
    n = ordered.size
    # Survival fraction at each observed value: P(overlap >= x).
    survival = 1.0 - np.arange(n) / n

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.loglog(ordered, survival, color="steelblue", lw=2)
    for alpha, label in [(0.05, "p = 0.05"), (0.01, "p = 0.01"), (0.001, "p = 0.001")]:
        ax.axhline(alpha, ls="--", color="gray", alpha=0.7)
        ax.text(ordered.min(), alpha * 1.1, label, fontsize=10, color="gray")
    ax.set_xlabel("Overlap Integral", fontsize=14)
    ax.set_ylabel("Survival Fraction  P(overlap $\\geq$ x)", fontsize=14)
    ax.set_title("Empirical Survival Function", fontsize=16)
    ax.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig


def plot_joint_skymap(
    skymap1_path: str,
    skymap2_path: str,
    save_path: str,
    event1_id: str = "Event 1",
    event2_id: str = "Event 2",
    overlap_value: float | None = None,
    nside_out: int = 512,
    levels=(0.05, 0.5, 0.9),
):
    """Plot two events and their joint ``P1*P2`` sky map on one Mollweide panel.

    Mirrors Figure 6 of the poster. Returns a dict of 90% credible-region areas
    (deg^2) for event 1, event 2, and the joint -- a small joint area relative to
    the individuals indicates the events are consistent with a shared position.
    """
    import healpy as hp
    import ligo.skymap.plot  # noqa: F401 -- registers the 'astro ... mollweide' projection
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    from ligo.skymap.moc import rasterize

    from .overlap import read_skymap

    s1_moc = read_skymap(skymap1_path, moc=True)
    s2_moc = read_skymap(skymap2_path, moc=True)

    order_out = int(np.log2(nside_out))
    r1 = rasterize(s1_moc, order=order_out)
    r2 = rasterize(s2_moc, order=order_out)
    p1 = np.asarray(r1["PROBDENSITY"], dtype=float)
    p2 = np.asarray(r2["PROBDENSITY"], dtype=float)

    pixel_area = hp.nside2pixarea(nside_out)
    if p1.sum() > 0:
        p1 = p1 / (p1.sum() * pixel_area)
    if p2.sum() > 0:
        p2 = p2 / (p2.sum() * pixel_area)

    joint = p1 * p2
    total = joint.sum() * pixel_area
    if total > 0:
        joint = joint / total

    def credible_from_density(prob_density):
        pmass = prob_density * pixel_area
        if pmass.sum() <= 0:
            return np.ones_like(pmass)
        pmass = pmass / pmass.sum()
        order = np.argsort(-pmass)
        cumulative = np.empty_like(pmass)
        cumulative[order] = np.cumsum(pmass[order])
        return cumulative

    cl1 = credible_from_density(p1)
    cl2 = credible_from_density(p2)
    cljoint = credible_from_density(joint)

    fig = plt.figure(figsize=(16, 10))
    ax = plt.axes(projection="astro hours mollweide")
    ax.grid()

    lvls = sorted(levels)
    blue_shades = ["#1e3a8a", "#2563eb", "#60a5fa"]
    red_shades = ["#7f1d1d", "#dc2626", "#fb7185"]
    gold_shades = ["#78350f", "#d97706", "#fbbf24"]
    line_widths = [2.6, 2.0, 1.5]

    ax.contour_hpx((cl1, "ICRS"), levels=lvls, nested=True,
                   colors=blue_shades, linewidths=line_widths)
    ax.contour_hpx((cl2, "ICRS"), levels=lvls, nested=True,
                   colors=red_shades, linewidths=line_widths)
    ax.contour_hpx((cljoint, "ICRS"), levels=lvls, nested=True,
                   colors=gold_shades, linewidths=line_widths)

    def _patches(shades, name):
        return [
            mpatches.Patch(color=c, label=f"{name} {int(level * 100)}%")
            for c, level in zip(shades, lvls)
        ]

    handles = (
        _patches(blue_shades, event1_id)
        + _patches(red_shades, event2_id)
        + _patches(gold_shades, "Joint (P1·P2)")
    )
    ax.legend(handles=handles, loc="upper right", fontsize=16, ncol=1)

    title = f"{event1_id} & {event2_id}  +  Joint Probability Skymap"
    if overlap_value is not None:
        fmt = ".2e" if abs(overlap_value) < 0.01 else ".2f"
        title += rf"     $\mathcal{{O}}={overlap_value:{fmt}}$"
    plt.title(title, fontsize=20, pad=20)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    sr_to_deg2 = (180.0 / np.pi) ** 2

    def cr_area_deg2(credible_map, cr):
        return float((credible_map <= cr).sum() * pixel_area * sr_to_deg2)

    return {
        "area1_90_deg2": cr_area_deg2(cl1, 0.9),
        "area2_90_deg2": cr_area_deg2(cl2, 0.9),
        "areaj_90_deg2": cr_area_deg2(cljoint, 0.9),
    }
