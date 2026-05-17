import matplotlib.pyplot as plt
import numpy as np

from ..utils import healpix as hp_utils

def plot_skymap(gw, transient, out_file="skymap.png"):
    '''
    Enhanced skymap plotting with proper healpy visualization
    '''
    fig = plt.figure(figsize=(12, 6))
    
    try:
        # Check if we have skymap data
        if hasattr(gw, 'skymap') and gw.skymap is not None:
            skymap_data = gw.skymap
        elif hasattr(gw, 'data') and gw.data is not None:
            skymap_data = gw.data
        else:
            # Load if needed
            if hasattr(gw, 'load_skymap'):
                gw.load_skymap()
                skymap_data = gw.skymap
            else:
                skymap_data = None
        
        if skymap_data is not None:
            # Plot with healpy (or raise informative error if unavailable)
            hp_utils.mollview(skymap_data, 
                       title="GW Skymap with EM Candidate",
                       unit="Probability",
                       fig=fig.number,
                       cmap='YlOrRd')
            
            # Add transient position
            if transient is not None and hasattr(transient, 'ra') and hasattr(transient, 'dec'):
                hp_utils.projscatter(transient.ra, transient.dec,
                              lonlat=True, marker='*', s=500, 
                              color='blue', edgecolor='white', linewidth=2,
                              label=f'EM Transient')
                
            # Add graticule
            hp_utils.graticule(dpar=30, dmer=30, alpha=0.3)
            
        else:
            # Fallback to simple plot
            ax = fig.add_subplot(111)
            ax.set_title("Skymap (data not loaded)")
            if transient is not None and hasattr(transient, 'ra') and hasattr(transient, 'dec'):
                ax.plot(transient.ra, transient.dec, '*', markersize=20, 
                       color='blue', label="EM Transient")
                ax.set_xlabel("RA [deg]")
                ax.set_ylabel("Dec [deg]")
                ax.legend()
                ax.grid(True, alpha=0.3)
            
    except Exception as e:
        print(f"Error in skymap plotting: {e}")
        # Basic fallback plot
        ax = fig.add_subplot(111)
        ax.set_title("Skymap (error in plotting)")
        if transient is not None and hasattr(transient, 'ra') and hasattr(transient, 'dec'):
            ax.plot(transient.ra, transient.dec, '*', markersize=20, color='blue')
            ax.set_xlabel("RA [deg]")
            ax.set_ylabel("Dec [deg]")
    
    plt.tight_layout()
    fig.savefig(out_file, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig