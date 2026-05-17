# tests/conftest.py
import os, glob, pytest
from types import SimpleNamespace
from gw_assoc.io import load_gw_skymap, GWEvent

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
