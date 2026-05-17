"""
Analysis modules for GW-EM association
"""

# Import main functions when needed
def __getattr__(name):
    if name == "compute_posterior_odds":
        from .odds import compute_posterior_odds
        return compute_posterior_odds
    elif name == "compute_coincident_odds":
        from .odds import compute_coincident_odds
        return compute_coincident_odds
    elif name == "SpatialOverlap":
        from .spatial import SpatialOverlap
        return SpatialOverlap
    elif name == "DistanceOverlap":
        from .los import DistanceOverlap
        return DistanceOverlap
    elif name == "TemporalOverlap":
        from .temporal import TemporalOverlap
        return TemporalOverlap
    elif name == "RadialOverlap":
        from .radial import RadialOverlap
        return RadialOverlap
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    "compute_posterior_odds",
    "compute_coincident_odds",
    "SpatialOverlap",
    "DistanceOverlap",
    "TemporalOverlap",
    "RadialOverlap",
]