import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from astropy.cosmology import Planck15
import astropy.units as u

def plot_distance_posteriors(gw_event, em_transient, output_file="distance_posteriors.png"):
    '''Plot GW and EM distance posteriors and their overlap'''
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Distance range
    dL_range = np.linspace(1, 500, 1000)
    
    # GW distance posterior (if available)
    ax1 = axes[0]
    if hasattr(gw_event, 'distances') and gw_event.distances is not None:
        dL_gw_mean = gw_event.distances.get('distmean', 100)
        dL_gw_std = gw_event.distances.get('diststd', 50)
        
        # Plot GW posterior
        gw_posterior = stats.norm.pdf(dL_range, dL_gw_mean, dL_gw_std)
        ax1.fill_between(dL_range, gw_posterior, alpha=0.5, 
                        color='blue', label='GW Posterior')
        ax1.axvline(dL_gw_mean, color='blue', linestyle='--', 
                   label=f'GW: {dL_gw_mean:.1f}±{dL_gw_std:.1f} Mpc')
    
    # EM distance posterior (from redshift)
    if hasattr(em_transient, 'z') and em_transient.z is not None:
        dL_em = Planck15.luminosity_distance(em_transient.z).to(u.Mpc).value
        
        # Approximate uncertainty
        z_err = em_transient.z_err if hasattr(em_transient, 'z_err') and em_transient.z_err else 0.01
        dL_em_err = dL_em * (z_err / em_transient.z) if em_transient.z > 0 else 10
        
        em_posterior = stats.norm.pdf(dL_range, dL_em, dL_em_err)
        ax1.fill_between(dL_range, em_posterior, alpha=0.5,
                        color='red', label='EM Posterior')
        ax1.axvline(dL_em, color='red', linestyle='--',
                   label=f'EM: {dL_em:.1f}±{dL_em_err:.1f} Mpc')
    
    ax1.set_xlabel('Luminosity Distance (Mpc)', fontsize=12)
    ax1.set_ylabel('Probability Density', fontsize=12)
    ax1.set_title('Distance Posteriors', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Overlap integral visualization
    ax2 = axes[1]
    if (hasattr(gw_event, 'distances') and gw_event.distances is not None and 
        hasattr(em_transient, 'z') and em_transient.z is not None):
        # Calculate overlap
        overlap = gw_posterior * em_posterior
        overlap_integral = np.trapz(overlap, dL_range)
        
        ax2.fill_between(dL_range, overlap, alpha=0.7,
                        color='purple', label=f'Overlap: {overlap_integral:.3f}')
        ax2.set_xlabel('Luminosity Distance (Mpc)', fontsize=12)
        ax2.set_ylabel('Overlap Density', fontsize=12)
        ax2.set_title(f'Distance Overlap Integral: {overlap_integral:.3f}', fontsize=14)
        ax2.legend()
        ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return fig

def plot_temporal_distribution(gw_event, em_transient, output_file="temporal.png"):
    '''Plot temporal probability distributions'''
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Time array (in days from GW event)
    t_days = np.linspace(0, 30, 1000)
    t_seconds = t_days * 86400
    
    # Different EM counterpart models
    models = {
        'Kilonova': (1.0, 2.0, 'red'),
        'Short GRB': (0.001, 0.01, 'blue'),
        'Afterglow': (0.1, 5.0, 'green')
    }
    
    for model_name, (t_peak, t_width, color) in models.items():
        t_peak_sec = t_peak * 86400
        t_width_sec = t_width * 86400
        prob = np.exp(-0.5 * ((t_seconds - t_peak_sec) / t_width_sec)**2)
        ax.plot(t_days, prob / np.max(prob), label=model_name, 
               color=color, linewidth=2)
    
    # Mark actual detection time if available
    if (hasattr(em_transient, 'time') and em_transient.time is not None and
        hasattr(gw_event, 'event_time')):
        dt_days = (em_transient.time - gw_event.event_time) / 86400
        if 0 <= dt_days <= 30:
            ax.axvline(dt_days, color='black', linestyle='--', 
                      label=f'Detection: {dt_days:.2f} days')
    
    ax.set_xlabel('Days since merger', fontsize=12)
    ax.set_ylabel('Temporal Probability', fontsize=12)
    ax.set_title('Temporal Association Probability', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 10)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return fig