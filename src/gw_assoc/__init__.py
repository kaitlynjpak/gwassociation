"""
GW-EM Association Framework

Main package initialization
"""

__version__ = '0.2.0'

# Lazy loading of main class
def __getattr__(name):
    if name == "Association":
        from .association import Association
        return Association
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = ["Association", "__version__"]