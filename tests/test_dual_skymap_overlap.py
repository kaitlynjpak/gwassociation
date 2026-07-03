import numpy as np

from gwassociation.analysis.spatial import SpatialOverlap
from gwassociation.analysis.radial import RadialOverlap


class DummyEvent:
    def __init__(self, prob, mu=None, sigma=None):
        self.skymap = np.array(prob, dtype=float)
        self.nside = 1
        self.nest = True
        self.is_3d = mu is not None and sigma is not None
        if self.is_3d:
            self.distances = {
                'distmu': np.array(mu, dtype=float),
                'distsigma': np.array(sigma, dtype=float),
                'distnorm': np.ones_like(mu, dtype=float)
            }
        else:
            self.distances = None

    def load_skymap(self):
        return {'data': self.skymap, 'nside': self.nside, 'is_3d': self.is_3d}


def test_spatial_overlap_map_vs_map_uniform():
    prob_primary = np.ones(12) / 12
    prob_secondary = np.ones(12) / 12

    event_a = DummyEvent(prob_primary)
    event_b = DummyEvent(prob_secondary)

    spatial = SpatialOverlap()
    overlap = spatial.compute_map_overlap(event_a, event_b)

    assert np.isclose(overlap, 1.0)


def _radial(mu_a, mu_b, sig_a=10.0, sig_b=15.0, npix=12):
    a = DummyEvent(np.ones(npix) / npix, np.full(npix, mu_a), np.full(npix, sig_a))
    b = DummyEvent(np.ones(npix) / npix, np.full(npix, mu_b), np.full(npix, sig_b))
    return RadialOverlap().compute_map_overlap(a, b)


def test_radial_overlap_positive_and_finite():
    overlap = _radial(100.0, 110.0)
    assert np.isfinite(overlap)
    assert overlap > 0.0


def test_radial_overlap_decreases_with_distance_separation():
    # Consistent distances overlap more than mismatched ones.
    close = _radial(100.0, 100.0)
    near = _radial(100.0, 130.0)
    far = _radial(100.0, 400.0)
    assert close > near > far


def test_radial_overlap_masks_pathological_pixels():
    """Low-probability pixels with garbage distance params must not blow up.

    Real 3D sky maps store distmu ~ -1e5 and distnorm ~ 1e80 in near-zero
    pixels; these must be masked so the result stays finite and bounded.
    """
    good_mu, good_sig = 100.0, 10.0
    prob = np.ones(12) / 12
    clean = DummyEvent(prob, np.full(12, good_mu), np.full(12, good_sig))

    mu = np.full(12, good_mu)
    sig = np.full(12, good_sig)
    mu[0] = -1.0e5        # pathological distmu (negative)
    sig[0] = 1.0e6        # pathological distsigma
    poisoned = DummyEvent(prob, mu, sig)
    poisoned.distances["distnorm"] = np.full(12, 1.0)
    poisoned.distances["distnorm"][0] = 1.0e80  # pathological distnorm

    overlap = RadialOverlap().compute_map_overlap(clean, poisoned)
    assert np.isfinite(overlap)
    assert 0.0 < overlap < 1.0e6  # bounded, not the old ~1e47 blow-up
