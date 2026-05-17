# src/gw_assoc/stats.py
from __future__ import annotations
import numpy as np
from typing import Literal

from .density import line_of_sight_pdf
from .utils import healpix as hp_utils

# --- Distance-prior (match notebook): p0(dL) ∝ dL^2 on [1, 1e4] Mpc, normalized ---
DL_MIN_MPC = 1.0
DL_MAX_MPC = 1.0e4

def prior_dl2(dL, dmin=DL_MIN_MPC, dmax=DL_MAX_MPC):
    """Normalized p0(dL) ∝ dL^2 on [dmin, dmax] (Mpc)."""
    dL = np.asarray(dL, dtype=float)
    p = np.zeros_like(dL, dtype=float)
    mask = (dL >= dmin) & (dL <= dmax)
    Z = (dmax**3 - dmin**3) / 3.0
    p[mask] = dL[mask]**2 / Z
    return p

def _gaussian_pos_pdf(ra0: float, dec0: float, sigma: float, nside: int) -> np.ndarray:
    """Circular Gaussian on the sphere (small-angle approx) centered at (ra0, dec0)."""
    npix = hp_utils.nside2npix(nside)
    pix = np.arange(npix)
    theta, phi = hp_utils.pix2ang(nside, pix)
    ra = phi
    dec = np.pi/2 - theta
    cosd = np.sin(dec0)*np.sin(dec) + np.cos(dec0)*np.cos(dec)*np.cos(ra - ra0)
    ang = np.arccos(np.clip(cosd, -1.0, 1.0))
    p = np.exp(-0.5*(ang/ max(sigma, 1e-9))**2)
    p = p / (p.sum() + 1e-300)
    return p

def spatial_overlap_Iw(pGW_Omega: np.ndarray,
                       pEM_Omega: np.ndarray,
                       prior_Omega: np.ndarray) -> float:
    """
    I_W = ∑_pix [ pGW(Ω_pix) * pEM(Ω_pix) / prior(Ω_pix) ] * dΩ
    All inputs normalized over pixels; dΩ = 4π/npix for HEALPix equal-area.
    """
    assert pGW_Omega.shape == pEM_Omega.shape == prior_Omega.shape
    npix = pGW_Omega.size
    domega = 4*np.pi/npix
    num = pGW_Omega * pEM_Omega
    den = np.clip(prior_Omega, 1e-300, np.inf)
    return float((num/den).sum() * domega)

def distance_overlap_ID(pGW_DL: np.ndarray, pEM_DL: np.ndarray, dDL: float) -> float:
    """
    I_D = ∫ pGW(DL) pEM(DL) / π(DL) dDL ; MVP uses flat prior π(DL)=const.
    Inputs are normalized PDFs sharing the same uniform DL grid (spacing dDL).
    """
    return float((pGW_DL * pEM_DL).sum() * dDL)

def distance_overlap_ID_lensed(pGW_DL: np.ndarray,
                               pEM_DL: np.ndarray,
                               dDL: float,
                               DL_grid: np.ndarray,
                               mu_grid: np.ndarray | None = None,
                               p_mu: np.ndarray | None = None) -> float:
    """
    Lensing-aware I_D: I_D^(lens) = ∫ dμ p(μ) ∫ dDL pGW(DL) pEM(DL*sqrt(μ)) dDL
    Implemented via discrete sum on (μ, DL). pEM is resampled as needed.
    """
    from scipy.interpolate import interp1d
    if mu_grid is None:
        mu_grid = np.concatenate([
            np.linspace(0.3, 0.9, 15),
            np.linspace(0.9, 3.0, 30),
            np.linspace(3.0, 10.0, 20)
        ])
    if p_mu is None:
        w = np.where(mu_grid < 1.0, mu_grid, mu_grid**-3)
        p_mu = w / (np.trapz(w, mu_grid) + 1e-300)

    f_pEM = interp1d(DL_grid, pEM_DL, kind='linear', bounds_error=False, fill_value=0.0)
    acc = 0.0
    for mu, wmu in zip(mu_grid, p_mu):
        pEM_mu = f_pEM(DL_grid * np.sqrt(mu))
        acc += wmu * (pGW_DL * pEM_mu).sum() * dDL
    return float(acc)

def joint_overlap_I3D(gw_density,
                      em,  # EMTransient
                      prior_Omega: np.ndarray,
                      DL_grid: np.ndarray) -> float:
    """
    I_3D ≈ ∑_pix ∫ dDL [ pGW(DL,Ω_pix) pEM(DL) pEM(Ω_pix) / prior(Ω_pix) ].
    EM sky is a Gaussian at (ra,dec) with width em.sigma_pos.
    """
    npix = hp_utils.nside2npix(gw_density.nside)
    pix = np.arange(npix)
    pEM_Omega = _gaussian_pos_pdf(em.ra, em.dec, max(em.sigma_pos, 1e-6), gw_density.nside)
    prior = prior_Omega / (prior_Omega.sum() + 1e-300)

    dDL = float(np.mean(np.diff(DL_grid)))
    acc = 0.0
    for p in pix:
        pGW_D = gw_density.pdf_D_omega(DL_grid, np.array([p])).squeeze()
        num = (pGW_D * em.p_DL) * pEM_Omega[p]
        den = np.clip(prior[p], 1e-300, np.inf)
        acc += (num/den).sum() * dDL
    domega = 4*np.pi/npix
    return float(acc * domega)

def IdL_at_host(ra_rad, dec_rad, z_host, cosmo, los_pdf_fn,
                dmin=DL_MIN_MPC, dmax=DL_MAX_MPC):
    """
    I_dL = p_LOS(dL_em) / p0(dL_em), evaluated at the EM host distance.

    Parameters
    ----------
    ra_rad, dec_rad : float (radians)
    z_host : float
    cosmo : astropy.cosmology.FlatLambdaCDM
    los_pdf_fn : callable (ra_rad, dec_rad, dL_mpc) -> p_LOS(dL | Ω) in 1/Mpc
        NOTE: los_pdf_fn must already include r^2 * Normal(mu, sigma) * distnorm(Ω).
    """
    dL_em = cosmo.luminosity_distance(z_host).to_value('Mpc')
    p_los = float(los_pdf_fn(ra_rad, dec_rad, dL_em))
    p0 = float(prior_dl2(dL_em, dmin, dmax))
    return p_los / (p0 if p0 > 1e-300 else 1e-300)

def IOmega_point(ra_rad, dec_rad, gw_prob, nside, nest=True):
    """
    I_Ω for a point host: gw sky probability density at (ra,dec) vs sky prior.
    gw_prob: HEALPix probability map normalized to 1 over the sphere.
    """
    theta = np.pi/2 - dec_rad
    phi = ra_rad
    ipix = hp_utils.ang2pix(nside, theta, phi, nest=nest)
    pix_area = hp_utils.nside2pixarea(nside)      # steradian
    p_density = gw_prob[ipix] / pix_area    # per steradian
    sky_prior = 1.0 / (4.0 * np.pi)
    return float(p_density / sky_prior)

def IOmega_maps(gw_prob, ext_prob, nside):
    """
    Map–map overlap with area weighting and proper normalization.
    Both inputs should be probabilities summing to 1 over the sphere.
    """
    area = hp_utils.nside2pixarea(nside)  # steradian
    num = float(np.sum(gw_prob * ext_prob) * area)                   # ∫ p_gw p_ext dΩ
    den = float((np.sum(gw_prob) * area) * (np.sum(ext_prob) * area))# (=1*1 if normalized)
    return num / (den if den > 1e-300 else 1e-300)

def temporal_factor(R_EM_per_day: float, delta_t_hours: float) -> float:
    """Return 1 / (R_EM * Δt). R_EM in day^-1, Δt in hours."""
    return 1.0 / (R_EM_per_day * max(delta_t_hours/24.0, 1e-9))

def posterior_odds(I_value: float, R_EM_per_day: float, delta_t_hours: float) -> float:
    """O_C/R ∝ I / (R_EM * Δt). Proportional odds (no absolute prior normalization)."""
    return float(I_value * temporal_factor(R_EM_per_day, delta_t_hours))
