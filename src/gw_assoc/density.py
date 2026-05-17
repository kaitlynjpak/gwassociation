'''Galaxy density and rate calculations'''
import numpy as np
from astropy.cosmology import FlatLambdaCDM

from .utils import healpix as hp_utils

cosmo = FlatLambdaCDM(H0=73, Om0=0.315)  # keep this if you want a default; or pass a cosmo around

def galaxy_density(z, galaxy_type='all'):
    '''
    Number density of galaxies at redshift z
    Returns: galaxies per Mpc^3
    '''
    # Simplified Schechter function parameters
    if galaxy_type == 'all':
        phi_star = 1e-2  # Mpc^-3
        alpha = -1.0
    elif galaxy_type == 'star_forming':
        phi_star = 5e-3
        alpha = -1.2
    else:  # elliptical
        phi_star = 3e-3
        alpha = -0.8
    
    # Evolution with redshift (simplified)
    evolution = (1 + z)**(-1)
    
    return phi_star * evolution

def expected_transients(volume, rate=1e-4):
    '''
    Expected number of transients in a volume
    
    Parameters:
    volume: Comoving volume in Mpc^3
    rate: Transient rate per Mpc^3 per year
    '''
    observation_time = 1.0  # year
    return volume * rate * observation_time

def false_alarm_probability(n_candidates, search_area, search_time):
    '''
    Probability of chance coincidence
    
    Parameters:
    n_candidates: Number of candidates found
    search_area: Sky area searched (square degrees)
    search_time: Time window (days)
    '''
    # All-sky transient rate (rough estimate)
    all_sky_rate = 100  # per day for magnitude < 20
    
    # Scale by search area (all sky = 41253 sq deg)
    area_fraction = search_area / 41253
    
    # Expected false alarms
    expected_false = all_sky_rate * area_fraction * search_time
    
    # Poisson probability
    from scipy.stats import poisson
    return 1 - poisson.cdf(n_candidates - 1, expected_false)

def line_of_sight_pdf(ra_rad, dec_rad, dL_mpc, gw_map):
    """
    LOS distance density p_LOS(dL | Ω) evaluated at (ra, dec) for given luminosity distance.

    Parameters
    ----------
    ra_rad, dec_rad : float
        Sky position in radians.
    dL_mpc : float or ndarray
        Luminosity distance(s) in Mpc.
    gw_map : dict-like
        Must provide keys:
          - 'prob'      : HEALPix probability map (sum ~ 1 over sphere)
          - 'distmu'    : per-pixel distance mean (Mpc)
          - 'distsigma' : per-pixel distance std (Mpc)
          - 'distnorm'  : per-pixel normalization factor from the 3D skymap
          - 'nside'     : HEALPix NSIDE
          - 'nest'      : bool, whether map is NESTED

    Returns
    -------
    p_LOS : float or ndarray
        Distance density in 1/Mpc at the given (ra, dec).
        Formula matches the notebook: dL^2 * Normal(dL; mu, sigma) * distnorm(Ω).
        (Optionally you may multiply by the sky prob if your code expects p(Ω, dL).)
    """
    nside = int(gw_map['nside'])
    nest  = bool(gw_map.get('nest', True))

    theta = np.pi/2 - dec_rad
    phi   = ra_rad
    ipix  = hp_utils.ang2pix(nside, theta, phi, nest=nest)

    mu      = float(gw_map['distmu'][ipix])
    sigma   = float(gw_map['distsigma'][ipix])
    distnor = float(gw_map['distnorm'][ipix])

    dL = np.asarray(dL_mpc, dtype=float)
    # 1D Gaussian in distance
    norm  = 1.0 / (np.sqrt(2.0*np.pi) * sigma)
    gauss = norm * np.exp(-0.5 * ((dL - mu) / sigma)**2)

    # === The required r^2 factor lives RIGHT HERE ===
    p_los = (dL**2) * gauss * distnor

    # If you want p(Ω, dL) instead of p(dL | Ω), uncomment next line to include sky prob:
    # p_los = p_los * float(gw_map['prob'][ipix])

    return p_los if p_los.ndim else float(p_los)
