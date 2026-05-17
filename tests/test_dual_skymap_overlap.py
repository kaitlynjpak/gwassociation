import numpy as np

from gw_assoc.analysis.spatial import SpatialOverlap
from gw_assoc.analysis.radial import RadialOverlap
from gw_assoc.stats import prior_dl2


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


def test_radial_overlap_gaussian_products():
    prob_primary = np.ones(12) / 12
    prob_secondary = np.ones(12) / 12
    mu_primary = np.full(12, 100.0)
    mu_secondary = np.full(12, 110.0)
    sigma_primary = np.full(12, 10.0)
    sigma_secondary = np.full(12, 15.0)

    event_a = DummyEvent(prob_primary, mu_primary, sigma_primary)
    event_b = DummyEvent(prob_secondary, mu_secondary, sigma_secondary)

    radial = RadialOverlap()
    overlap = radial.compute_map_overlap(event_a, event_b)

    var = sigma_primary[0] ** 2 + sigma_secondary[0] ** 2
    gaussian_overlap = (1 / np.sqrt(2 * np.pi * var)) * np.exp(-0.5 * (10 ** 2) / var)
    prior = prior_dl2(105.0)
    expected = gaussian_overlap / prior

    assert np.isclose(overlap, expected)


