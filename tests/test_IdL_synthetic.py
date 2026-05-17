# tests/test_IdL_synthetic.py
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from gw_assoc.analysis.los import line_of_sight_pdf
from gw_assoc.stats import prior_dl2, IdL_at_host

def test_IdL_peaks_near_mu():
    mu, sigma = 1000.0, 100.0
    gw_map = {
        "prob":      np.array([1.0]),      # not used in LOS ratio
        "distmu":    np.array([mu]),
        "distsigma": np.array([sigma]),
        "distnorm":  np.array([1.0]),      # choose 1; normalization is in gaussian + r^2
        "nside":     1, "nest": True
    }
    # Fake ang2pix→0 by forcing the function to read index 0: pass any ra/dec, nside=1 makes all → 0
    f = lambda ra,dec,dL: line_of_sight_pdf(ra, dec, dL, gw_map)
    cosmo = FlatLambdaCDM(H0=73, Om0=0.315)

    ra, dec = 0.0, 0.0
    # evaluate IdL around mu
    def IdL_at(dL):
        p_los = f(ra, dec, dL)
        return p_los / prior_dl2(dL)

    r_grid = np.array([mu - 3*sigma, mu - 2*sigma, mu, mu + 2*sigma, mu + 3*sigma])
    vals = [IdL_at(r) for r in r_grid]
    # peak near mu
    assert vals[2] == max(vals)
