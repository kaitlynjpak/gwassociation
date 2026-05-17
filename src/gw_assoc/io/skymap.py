import numpy as np
from typing import Dict, Optional, Tuple
try:
    import ligo.skymap.io
    import ligo.skymap.distance
    LIGO_SKYMAP_AVAILABLE = True
except ImportError:
    LIGO_SKYMAP_AVAILABLE = False

from ..utils import healpix as hp_utils

def load_gw_skymap(path: str) -> Dict:
    '''
    Load GW skymap from FITS file
    
    Returns dict with skymap data and metadata
    '''
    if LIGO_SKYMAP_AVAILABLE:
        try:
            skymap_data = ligo.skymap.io.read_sky_map(path)
            prob_map = skymap_data[0]
            nside = hp_utils.npix2nside(len(prob_map))
            
            if len(skymap_data) == 4:  # 2D skymap
                return {
                    'file': path,
                    'kind': 'gw_skymap_2d',
                    'data': prob_map,
                    'header': skymap_data[1],
                    'nside': nside,
                    'is_3d': False,
                    'nest': True
                }
            else:  # 3D skymap
                distances = ligo.skymap.distance.parameters_to_marginal_moments(
                    skymap_data[0], skymap_data[1]
                )
                return {
                    'file': path,
                    'kind': 'gw_skymap_3d',
                    'data': prob_map,
                    'distances': distances,
                    'nside': nside,
                    'is_3d': True,
                    'nest': True
                }
        except Exception as e:
            print(f"Error loading with ligo.skymap: {e}")
    
    # Fallback to simple healpy loading
    try:
        skymap, header = hp_utils.read_map(path, h=True, verbose=False)
        return {
            'file': path,
            'kind': 'gw_skymap_2d',
            'data': skymap,
            'header': dict(header),
            'nside': hp_utils.npix2nside(len(skymap)),
            'is_3d': False,
            'nest': False
        }
    except Exception as e:
        print(f"Error loading skymap: {e}")
        # Return minimal placeholder
        return {'file': path, 'kind': 'gw_skymap', 'data': None}

class GWEvent:
    '''Container for GW event data'''
    def __init__(self, skymap_path: str, event_time: float, event_name: str = None):
        self.skymap_path = skymap_path
        self.event_time = event_time  # GPS time
        self.event_name = event_name
        self.skymap = None
        self.distances = None
        self.nside = None
        self.is_3d = False
        self.nest = True
        
    def load_skymap(self):
        '''Load skymap data'''
        skymap_dict = load_gw_skymap(self.skymap_path)
        self.skymap = skymap_dict.get('data')
        self.distances = skymap_dict.get('distances')
        self.nside = skymap_dict.get('nside')
        self.is_3d = skymap_dict.get('is_3d', False)
        self.nest = skymap_dict.get('nest', True)
        return skymap_dict