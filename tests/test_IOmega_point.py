# tests/test_IOmega_point.py
import numpy as np
from gw_assoc.analysis.spatial import IOmega_point

def test_IOmega_uniform_expectation_is_one(gwevent, N=2000, tol=0.1):
    sm = gwevent.skymap
    nside, prob = sm["nside"], sm["prob"]
    # Random uniform sky samples
    u = np.random.rand(N)
    v = np.random.rand(N)
    ra  = 2*np.pi*u
    dec = np.arcsin(2*v - 1)

    vals = [IOmega_point(ra[i], dec[i], prob, nside, nest=sm.get("nest", True)) for i in range(N)]
    mean = np.mean(vals)
    assert abs(mean - 1.0) < tol  # law of large numbers

def test_IOmega_high_pixel_is_large(gwevent):
    from gw_assoc.utils import healpix as hp_utils

    sm = gwevent.skymap
    ip = np.argmax(sm["prob"])
    # invert pix->ang to pick that pixel’s center
    theta, phi = hp_utils.pix2ang(sm["nside"], ip, nest=sm.get("nest", True))
    ra, dec = phi, (np.pi/2 - theta)
    val = IOmega_point(ra, dec, sm["prob"], sm["nside"], nest=sm.get("nest", True))
    assert val > 1.0
