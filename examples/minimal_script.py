# examples/minimal_script.py  — REAL calculation (3D FITS path)
import numpy as np
from gw_assoc.ingest import load_gw_3d_fits, load_transient
from gw_assoc.density import GWPosteriorDensity
from gw_assoc.fov import fullsky_prior, chime_fov_prior
from gw_assoc.stats import joint_overlap_I3D, distance_overlap_ID, distance_overlap_ID_lensed, posterior_odds
from gw_assoc.plots import plot_los_distance, render_odds_table

# --- EDIT THESE ---
GW_3D = "Bilby.offline1.multiorder.fits.gz"  # your file
EM_RA, EM_DEC = 255.72, 21.52
EM_TIME = 58630.5
EM_Z, EM_SIGMA_Z = 0.03136, 0.0005
R_EM_PER_DAY = 1.6
DELTA_T_HOURS = 26.0
USE_LENSING = False
# ------------------

# Build GW density from 3D FITS cube (p(D,Ω))
skymap = load_gw_3d_fits(GW_3D)          # NOTE: MVP loader expects DIST_MID/PROB HDUs
gw = GWPosteriorDensity.from_3d_fits(skymap)

# EM transient
em = load_transient(EM_RA, EM_DEC, t_em=EM_TIME, z=EM_Z, sigma_z=EM_SIGMA_Z)

prior_full = fullsky_prior(gw.nside)
prior_fov  = chime_fov_prior(gw.nside, ra_center=np.deg2rad(EM_RA), width_hours=5.0)

I3_full = joint_overlap_I3D(gw, em, prior_full, em.DL_grid)
I3_fov  = joint_overlap_I3D(gw, em, prior_fov,  em.DL_grid)

# LOS from the 3D map
pGW_LOS = gw.pdf_LOS(em.DL_grid, em.ra, em.dec)
dDL = float(np.mean(np.diff(em.DL_grid)))
ID = (distance_overlap_ID_lensed(pGW_LOS, em.p_DL, dDL, em.DL_grid)
      if USE_LENSING else
      distance_overlap_ID(pGW_LOS, em.p_DL, dDL))

O3_full = posterior_odds(I3_full, R_EM_PER_DAY, DELTA_T_HOURS)
O3_fov  = posterior_odds(I3_fov,  R_EM_PER_DAY, DELTA_T_HOURS)
O1D     = posterior_odds(ID,      R_EM_PER_DAY, DELTA_T_HOURS)

print("\n=== gw_assoc MVP results (real overlaps) ===")
render_odds_table([
    ("I_3D + Full-sky prior", O3_full),
    ("I_3D + FOV prior",      O3_fov),
    ("LOS distance only",     O1D),
])
