"""Container objects and convenience methods for electromagnetic transients.

The :class:`Transient` dataclass stores the minimum candidate information used by
association calculations: sky position, optional redshift, optional timing, and
optional classification metadata.
"""

from dataclasses import dataclass
from typing import Optional
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck15
import astropy.units as u

@dataclass
class Transient:
    '''Enhanced container for a transient event with full functionality'''
    ra: float  # degrees
    dec: float  # degrees
    z: Optional[float] = None  # redshift
    z_err: Optional[float] = None  # redshift error
    time: Optional[float] = None  # detection time (GPS or MJD)
    time_err: Optional[float] = None  # time uncertainty
    magnitude: Optional[float] = None
    filter_band: Optional[str] = None
    name: Optional[str] = None
    
    def __repr__(self):
        """Return a concise debug representation with sky position and redshift."""
        return f"<Transient {self.name if self.name else ''} ra={self.ra:.3f}, dec={self.dec:.3f}, z={self.z}>"
    
    def get_skycoord(self):
        '''Return SkyCoord object'''
        return SkyCoord(ra=self.ra*u.deg, dec=self.dec*u.deg)
    
    def get_luminosity_distance(self, cosmo=Planck15):
        '''Convert redshift to luminosity distance'''
        if self.z is not None:
            return cosmo.luminosity_distance(self.z).to(u.Mpc).value
        return None