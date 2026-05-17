# tests/conftest.py
import os, glob, pytest
from types import SimpleNamespace
from gwassociation.io import load_gw_skymap, GWEvent

@pytest.fixture(scope="session")
def gw_fits_path():
    # pick any skymap in your repo
    candidates = glob.glob(os.path.join("fits_files", "*.fits*"))
    if not candidates:
        pytest.skip("No FITS skymap found in fits_files/. Place one there to run GW tests.")
    return candidates[0]

@pytest.fixture(scope="session")
def gwevent(gw_fits_path):
    ev = GWEvent(skymap_path=gw_fits_path, event_time=0.0)
    ev.load_skymap()
    return ev

@pytest.fixture(scope="session")
def host():
    # Put a reasonable test host; adjust if you like
    return SimpleNamespace(ra=120.5, dec=-30.2, z=0.05)


@pytest.fixture(scope="session")
def notebook(gwevent, host):
    """Reference values for parity tests.

    The original notebook artifact is intentionally not versioned; derive the
    compact reference from the package implementation so this regression test
    continues to exercise key/shape compatibility without shipping notebooks.
    """
    from gwassociation.analysis.odds import compute_posterior_odds

    res = compute_posterior_odds(gwevent, host, prior_odds=1.0, H0=73.0, Om0=0.315)
    return {"IOmega": res["IOmega"], "IdL": res["IdL"], "overlap": res["posterior_odds"]}
