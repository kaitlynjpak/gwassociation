# tests/test_distance_norm.py
import numpy as np
from scipy.integrate import quad
from gw_assoc.analysis.los import line_of_sight_pdf

def test_los_integrates_to_one(gwevent, tol=5e-2):
    gw = gwevent  # fixture: load one map (or construct a small dummy dict)
    sm = gw.skymap
    ra, dec = np.deg2rad(120.0), np.deg2rad(-30.0)  # pick any point in support
    gw_map = dict(prob=sm["prob"], distmu=sm["distmu"], distsigma=sm["distsigma"],
                  distnorm=sm["distnorm"], nside=sm["nside"], nest=sm.get("nest", True))
    f = lambda r: line_of_sight_pdf(ra, dec, r, gw_map)
    val = quad(f, 1.0, 1.0e4, limit=200)[0]  # Mpc bounds match your prior
    assert abs(val - 1.0) < tol
