# tests/test_end_to_end_parity.py
import numpy as np
from gw_assoc.analysis.odds import compute_posterior_odds

def test_parity_with_notebook(gwevent, host, notebook):
    """
    Fixtures:
      - gwevent: loaded GWEvent
      - host: object with ra, dec, z
      - notebook: dict with keys 'IOmega', 'IdL', 'overlap'
    """
    res = compute_posterior_odds(gwevent, host, prior_odds=1.0, H0=73.0, Om0=0.315)
    assert np.isclose(res["IOmega"],  notebook["IOmega"],  rtol=0.1)
    assert np.isclose(res["IdL"],     notebook["IdL"],     rtol=0.1)
    assert np.isclose(res["posterior_odds"], notebook["overlap"], rtol=0.1)
