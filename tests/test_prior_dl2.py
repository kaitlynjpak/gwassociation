# tests/test_prior_dl2.py
import numpy as np
from scipy.integrate import quad
from gw_assoc.stats import prior_dl2, DL_MIN_MPC, DL_MAX_MPC

def test_prior_normalizes_to_one():
    val = quad(lambda r: prior_dl2(r), DL_MIN_MPC, DL_MAX_MPC, limit=200)[0]
    assert abs(val - 1.0) < 1e-6

def test_prior_shape_increasing():
    xs = np.array([10, 100, 1000])
    p = prior_dl2(xs)
    assert p[0] < p[1] < p[2]  # ∝ r^2
