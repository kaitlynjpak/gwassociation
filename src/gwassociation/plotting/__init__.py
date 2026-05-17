"""
Plotting modules for GW-EM association
"""

def __getattr__(name):
    if name == "plot_skymap":
        from .skymap import plot_skymap
        return plot_skymap
    elif name == "plot_distance_posteriors":
        from .distributions import plot_distance_posteriors
        return plot_distance_posteriors
    elif name == "plot_temporal_distribution":
        from .distributions import plot_temporal_distribution
        return plot_temporal_distribution
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    "plot_skymap",
    "plot_distance_posteriors",
    "plot_temporal_distribution"
]