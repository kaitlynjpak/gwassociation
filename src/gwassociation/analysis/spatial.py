import numpy as np

from ..utils import healpix as hp_utils


def _as_skymap_dict(skymap_or_dict):
    if isinstance(skymap_or_dict, dict):
        return skymap_or_dict
    return {'data': np.asarray(skymap_or_dict), 'nside': hp_utils.npix2nside(len(skymap_or_dict))}


def _normalize_probabilities(prob_array):
    arr = np.asarray(prob_array, dtype=float)
    total = np.sum(arr)
    if total <= 0:
        raise ValueError("Skymap probabilities must sum to a positive value.")
    return arr / total


def _point_overlap(gw_map, ra, dec, gw_nested=True):
    theta = np.radians(90 - dec)
    phi = np.radians(ra)
    nside = gw_map['nside']
    data = np.asarray(gw_map['data'], dtype=float)

    ipix = hp_utils.ang2pix(nside, theta, phi, nest=gw_nested)
    norm = data.sum()
    if norm <= 0:
        return 0.0
    probability_mass = data[ipix] / norm
    npix = data.size
    return probability_mass * npix


def _map_overlap(gw_map, ext_map):
    if gw_map['nside'] != ext_map['nside']:
        raise ValueError("GW and external skymaps must share the same NSIDE for overlap computation.")

    gw_prob = _normalize_probabilities(gw_map['data'])
    ext_prob = _normalize_probabilities(ext_map['data'])

    overlap_mass = np.sum(gw_prob * ext_prob)
    npix = gw_prob.size
    return npix * overlap_mass


def skymap_overlap_integral(gw_skymap, ext_skymap=None,
                            ra=None, dec=None,
                            gw_nested=True, ext_nested=True):
    '''
    Compute the spatial overlap integral between a GW skymap and either a point source
    or another skymap.
    '''
    gw_map = _as_skymap_dict(gw_skymap)

    if ext_skymap is not None:
        ext_map = _as_skymap_dict(ext_skymap)
        return _map_overlap(gw_map, ext_map)

    if ra is not None and dec is not None:
        return _point_overlap(gw_map, ra, dec, gw_nested=gw_nested)

    raise ValueError("Provide either (ra, dec) for point sources or ext_skymap for map overlap.")

class SpatialOverlap:
    '''Calculate spatial overlap integral I_Ω between GW skymap and EM position'''
    
    @staticmethod
    def compute(gw_event, em_transient, search_radius: float = 0.0) -> float:
        '''
        Calculate spatial overlap integral I_Ω
        
        Parameters:
        -----------
        gw_event: GWEvent object with loaded skymap
        em_transient: Transient object with position
        search_radius: Error radius in degrees (default 0 for point source)
        
        Returns:
        --------
        I_omega: Spatial overlap probability
        '''
        if gw_event.skymap is None:
            gw_event.load_skymap()
        
        # Use the existing skymap_overlap_integral function
        return skymap_overlap_integral(
            {'data': gw_event.skymap, 'nside': gw_event.nside},
            ra=em_transient.ra,
            dec=em_transient.dec,
            gw_nested=True
        )

    @staticmethod
    def compute_map_overlap(primary_event, secondary_event) -> float:
        """
        Calculate spatial overlap for two full skymaps (skymap-vs-skymap).
        """
        if primary_event.skymap is None:
            primary_event.load_skymap()
        if secondary_event.skymap is None:
            secondary_event.load_skymap()

        primary_map = {'data': primary_event.skymap, 'nside': primary_event.nside}
        secondary_map = {'data': secondary_event.skymap, 'nside': secondary_event.nside}

        return skymap_overlap_integral(
            primary_map,
            ext_skymap=secondary_map,
            gw_nested=getattr(primary_event, "nest", True),
            ext_nested=getattr(secondary_event, "nest", True)
        )