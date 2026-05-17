'''Field of view calculations for different telescopes'''

import numpy as np

# Telescope FOVs in square degrees
TELESCOPE_FOV = {
    'ZTF': 47.0,
    'DECam': 3.0,
    'Pan-STARRS': 7.0,
    'ATLAS': 30.0,
    'BlackGEM': 2.7,
    'GOTO': 40.0,
    'LSST': 9.6,
    'Swift-UVOT': 0.056,
    'Gaia': 0.7  # Approximate scanning
}

def tiles_needed(skymap_area, telescope='ZTF'):
    '''Calculate number of tiles needed to cover skymap area'''
    fov = TELESCOPE_FOV.get(telescope, 1.0)
    # Account for overlap between tiles (~10%)
    efficiency = 0.9
    return int(np.ceil(skymap_area / (fov * efficiency)))

def observation_time(skymap_area, telescope='ZTF', exposure_time=30, overhead=30):
    '''Estimate total observation time including overheads'''
    n_tiles = tiles_needed(skymap_area, telescope)
    time_per_tile = exposure_time + overhead  # seconds
    return n_tiles * time_per_tile / 3600  # hours