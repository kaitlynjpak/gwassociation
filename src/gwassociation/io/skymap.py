"""Sky-map loading and event containers."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any
import importlib
import importlib.util

import numpy as np
from astropy.io import fits

from ..utils import healpix as hp_utils

LIGO_SKYMAP_AVAILABLE = (
    importlib.util.find_spec("ligo") is not None
    and importlib.util.find_spec("ligo.skymap") is not None
    and importlib.util.find_spec("ligo.skymap.io.fits") is not None
)
ligo_fits = importlib.import_module("ligo.skymap.io.fits") if LIGO_SKYMAP_AVAILABLE else None


def _normalize_prob(probability: np.ndarray) -> np.ndarray:
    """Return a finite HEALPix probability array normalized to sum to one."""
    prob = np.asarray(probability, dtype=float).ravel()
    prob = np.where(np.isfinite(prob), prob, 0.0)
    total = prob.sum()
    if total <= 0:
        raise ValueError("Sky-map probability column has no positive finite values.")
    return prob / total


def _empty_distance(prob: np.ndarray) -> dict[str, np.ndarray]:
    """Distance placeholders for 2D maps so downstream code can branch cleanly."""
    shape = prob.shape
    return {
        "distmu": np.full(shape, np.nan),
        "distsigma": np.full(shape, np.nan),
        "distnorm": np.zeros(shape),
    }


def load_gw_skymap(path: str) -> dict[str, Any]:
    """Load a 2D or 3D GW HEALPix sky map from FITS.

    The returned dictionary intentionally exposes both ``prob`` and legacy
    ``data`` keys, each pointing to the normalized sky probability map. 3D
    distance columns are available at top level and inside ``distances``.
    """
    if LIGO_SKYMAP_AVAILABLE:
        try:
            table = ligo_fits.read_sky_map(path, distances=True, nest=None)
            if isinstance(table, tuple):
                # ligo.skymap has varied return conventions across versions.
                prob = _normalize_prob(table[0])
                distmu = np.asarray(table[1], dtype=float) if len(table) > 1 else None
                distsigma = np.asarray(table[2], dtype=float) if len(table) > 2 else None
                distnorm = np.asarray(table[3], dtype=float) if len(table) > 3 else None
                metadata = table[4] if len(table) > 4 and isinstance(table[4], dict) else {}
            else:
                names = [name.lower() for name in table.colnames]
                prob = _normalize_prob(table[table.colnames[names.index("prob")]])
                distmu = np.asarray(table[table.colnames[names.index("distmu")]], dtype=float) if "distmu" in names else None
                distsigma = np.asarray(table[table.colnames[names.index("distsigma")]], dtype=float) if "distsigma" in names else None
                distnorm = np.asarray(table[table.colnames[names.index("distnorm")]], dtype=float) if "distnorm" in names else None
                metadata = dict(getattr(table, "meta", {}))
            nest = bool(metadata.get("nest", metadata.get("ORDERING", "NESTED") == "NESTED"))
            nside = hp_utils.npix2nside(prob.size)
            has_3d = distmu is not None and distsigma is not None and distnorm is not None
            distances = (
                {"distmu": distmu, "distsigma": distsigma, "distnorm": distnorm}
                if has_3d
                else _empty_distance(prob)
            )
            return {
                "file": path,
                "kind": "gw_skymap_3d" if has_3d else "gw_skymap_2d",
                "prob": prob,
                "data": prob,
                "nside": nside,
                "nest": nest,
                "is_3d": has_3d,
                "header": metadata,
                "distances": distances if has_3d else None,
                **distances,
            }
        except Exception:
            # Fall through to the local FITS parser; it is more predictable for tests.
            pass

    with fits.open(path) as hdul:
        header: dict[str, Any] = {}
        for hdu in hdul:
            header.update(dict(hdu.header))
            if getattr(hdu, "data", None) is None:
                continue
            data = hdu.data
            if hasattr(data, "columns"):
                names = list(data.columns.names)
                lower = [name.lower() for name in names]
                prob_col = names[lower.index("prob")] if "prob" in lower else names[0]
                prob = _normalize_prob(data[prob_col])
                distmu = np.asarray(data[names[lower.index("distmu")]], dtype=float) if "distmu" in lower else None
                distsigma = np.asarray(data[names[lower.index("distsigma")]], dtype=float) if "distsigma" in lower else None
                distnorm = np.asarray(data[names[lower.index("distnorm")]], dtype=float) if "distnorm" in lower else None
                break
            else:
                prob = _normalize_prob(np.asarray(data, dtype=float))
                distmu = distsigma = distnorm = None
                break
        else:
            raise ValueError(f"No probability map found in {path!r}.")

    nside = hp_utils.npix2nside(prob.size)
    nest = str(header.get("ORDERING", "NESTED")).upper() == "NESTED"
    has_3d = distmu is not None and distsigma is not None and distnorm is not None
    distances = (
        {"distmu": distmu, "distsigma": distsigma, "distnorm": distnorm}
        if has_3d
        else _empty_distance(prob)
    )
    return {
        "file": path,
        "kind": "gw_skymap_3d" if has_3d else "gw_skymap_2d",
        "prob": prob,
        "data": prob,
        "nside": nside,
        "nest": nest,
        "is_3d": has_3d,
        "header": header,
        "distances": distances if has_3d else None,
        **distances,
    }


@dataclass
class GWEvent:
    """Container for a GW event and its loaded sky-map metadata."""

    skymap_path: str
    event_time: float
    event_name: str | None = None
    skymap: dict[str, Any] | None = None
    distances: dict[str, np.ndarray] | None = None
    nside: int | None = None
    is_3d: bool = False
    nest: bool = True
    prob: np.ndarray | None = field(default=None, init=False)

    def load_skymap(self) -> dict[str, Any]:
        """Load the sky map and cache commonly used aliases on the event."""
        skymap_dict = load_gw_skymap(self.skymap_path)
        self.skymap = skymap_dict
        self.prob = skymap_dict["prob"]
        self.distances = skymap_dict.get("distances")
        self.nside = skymap_dict.get("nside")
        self.is_3d = bool(skymap_dict.get("is_3d", False))
        self.nest = bool(skymap_dict.get("nest", True))
        return skymap_dict
