"""Line-of-sight distance calculations for GW--EM association tests.

This module contains helper functions for converting between redshift and
luminosity distance under a flat :class:`astropy.cosmology.FlatLambdaCDM`
cosmology, plus :class:`DistanceOverlap`, the public calculator used by the
high-level association API.  Distances are expressed in megaparsecs unless noted
otherwise.
"""

import numpy as np
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM, Planck15
from astropy import constants
import astropy.units as u

from ..density import line_of_sight_pdf

# [Your existing code from los.py]
zinterp = np.linspace(0, 0.5, 5000)
dz = zinterp[1] - zinterp[0]
speed_of_light = constants.c.to('km/s').value
H0Planck = Planck15.H0.value
Om0Planck = Planck15.Om0

def dL_at_z_H0(z, h0, Om0):
    """Return luminosity distance in Mpc for redshift ``z`` and cosmology parameters."""
    cosmo = FlatLambdaCDM(H0=h0, Om0=Om0)
    dLs = cosmo.luminosity_distance(z).to(u.Mpc).value
    return dLs

def z_at_dL_H0(dL, h0, Om0):
    """Interpolate the redshift corresponding to luminosity distance ``dL`` in Mpc."""
    cosmo = FlatLambdaCDM(H0=h0, Om0=Om0)
    dLs = cosmo.luminosity_distance(zinterp).to(u.Mpc).value  
    z_at_dL = interp1d(dLs, zinterp)
    return z_at_dL(dL)

def E(z, Om):
    """Return the dimensionless Hubble expansion factor for matter density ``Om``."""
    return np.sqrt(Om*(1+z)**3 + (1.0-Om))

def dz_by_dL_H0(z, dL, h0, Om0):
    """Return the Jacobian ``dz / dD_L`` for a flat cosmology."""
    return 1/(dL/(1+z) + speed_of_light*(1+z)/(h0*E(z,Om0)))

def dL_by_z_H0(z, dL, h0, Om0):
    """Return the Jacobian ``dD_L / dz`` for a flat cosmology."""
    return dL/(1+z) + speed_of_light*(1+z)/(h0*E(z,Om0))

def dvdz(z, H0, Om0):
    """Return all-sky differential comoving volume in Gpc^3 per unit redshift."""
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    dvdz = 4*np.pi*cosmo.differential_comoving_volume(z).to(u.Gpc**3/u.sr).value
    return dvdz

class DistanceOverlap:
    '''Calculate distance overlap integral I_D_L between GW and EM distances'''
    
    def __init__(self, cosmo=Planck15):
        """Create a distance-overlap calculator using an Astropy cosmology."""
        self.cosmo = cosmo
        self.setup_interpolation()
    
    def setup_interpolation(self):
        '''Setup redshift-distance interpolation'''
        self.z_grid = np.linspace(0, 0.5, 5000)
        self.dL_grid = self.cosmo.luminosity_distance(self.z_grid).to(u.Mpc).value
        self.z_at_dL = interp1d(self.dL_grid, self.z_grid, 
                                bounds_error=False, fill_value=np.nan)
        self.dL_at_z = interp1d(self.z_grid, self.dL_grid,
                                bounds_error=False, fill_value=np.nan)
    
    def compute(self, gw_event, em_transient, H0_uncertainty: float = 7.0) -> float:
        '''
        Calculate distance overlap integral I_D_L
        '''
        if em_transient.z is None:
            return 1.0  # No distance constraint
        
        # Get EM luminosity distance with uncertainties. Accept both the
        # Transient dataclass and light-weight objects exposing only ``z``.
        if hasattr(em_transient, "get_luminosity_distance"):
            dL_em = em_transient.get_luminosity_distance(self.cosmo)
        else:
            dL_em = self.cosmo.luminosity_distance(em_transient.z).to(u.Mpc).value
        
        # Add peculiar velocity uncertainty (~300 km/s typical)
        v_pec = 300  # km/s
        z_pec = v_pec / speed_of_light
        z_total = em_transient.z
        z_err = getattr(em_transient, "z_err", None) or 0.01
        
        # Approximate uncertainty in distance
        dL_em_err = dL_em * np.sqrt((z_err/z_total)**2 + 
                                    (z_pec/z_total)**2 + 
                                    (H0_uncertainty/self.cosmo.H0.value)**2)
        
        if gw_event.is_3d and gw_event.distances is not None:
            # Use GW distance posterior
            dL_gw_mean = gw_event.distances.get('distmean', 100)
            dL_gw_std = gw_event.distances.get('diststd', 50)
            
            # Calculate overlap integral (simplified Gaussian approximation)
            var_sum = dL_gw_std**2 + dL_em_err**2
            I_DL = np.exp(-0.5 * (dL_em - dL_gw_mean)**2 / var_sum) / np.sqrt(2*np.pi*var_sum)
            
            # Normalize by prior (uniform in comoving volume)
            I_DL *= dL_em**2
        else:
            # No GW distance information
            I_DL = 1.0
            
        return I_DL