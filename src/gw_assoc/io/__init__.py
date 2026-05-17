"""
Input/Output modules for GW-EM association
"""

def __getattr__(name):
    if name == "load_gw_skymap":
        from .skymap import load_gw_skymap
        return load_gw_skymap
    elif name == "GWEvent":
        from .skymap import GWEvent
        return GWEvent
    elif name == "Transient":
        from .transient import Transient
        return Transient
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = ["load_gw_skymap", "GWEvent", "Transient"]