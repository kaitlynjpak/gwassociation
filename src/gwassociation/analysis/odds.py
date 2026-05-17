import numpy as np
from typing import Dict, Optional
from .spatial import SpatialOverlap
from .los import DistanceOverlap
from .temporal import TemporalOverlap
from .radial import RadialOverlap

def compute_posterior_odds(gw, transient, **kwargs):
    '''
    Compute posterior odds for GW-EM association
    
    O_posterior = O_prior × BF
    where BF = P(data|associated) / P(data|not associated)
         = I_Ω × I_D_L × I_t / P_chance
    '''
    
    # Initialize calculators
    spatial_calc = SpatialOverlap()
    distance_calc = DistanceOverlap()
    temporal_calc = TemporalOverlap()
    
    # Calculate overlap integrals
    I_omega = kwargs.get('I_omega')
    if I_omega is None:
        I_omega = spatial_calc.compute(gw, transient, 
                                       search_radius=kwargs.get('search_radius', 0.0))
    
    I_dl = kwargs.get('I_dl')
    if I_dl is None:
        I_dl = distance_calc.compute(gw, transient,
                                     H0_uncertainty=kwargs.get('H0_uncertainty', 7.0))
    
    I_t = kwargs.get('I_t')
    if I_t is None:
        I_t = temporal_calc.compute(gw, transient,
                                   model=kwargs.get('em_model', 'kilonova'))
    
    # Get priors and chance coincidence rate
    prior_odds = kwargs.get('prior_odds', 1.0)
    P_chance = kwargs.get('chance_coincidence_rate', 1e-4)
    
    # Calculate Bayes factor
    bayes_factor = (I_omega * I_dl * I_t) / P_chance
    
    # Posterior odds
    posterior_odds = prior_odds * bayes_factor
    
    # Log odds for numerical stability
    log_posterior_odds = np.log10(posterior_odds) if posterior_odds > 0 else -np.inf
    
    # Convert to probability
    confidence = 1 - 1/(1 + posterior_odds)
    
    return {
        'gw_file': gw.skymap_path if hasattr(gw, 'skymap_path') else None,
        'transient': {
            'ra': getattr(transient, 'ra', None),
            'dec': getattr(transient, 'dec', None),
            'z': getattr(transient, 'z', None),
            'time': getattr(transient, 'time', None),
        },
        'I_omega': I_omega,
        'I_dl': I_dl,
        'I_t': I_t,
        'bayes_factor': bayes_factor,
        'prior_odds': prior_odds,
        'posterior_odds': posterior_odds,
        'log_posterior_odds': log_posterior_odds,
        'P_chance': P_chance,
        'associated': posterior_odds > 1.0,
        'confidence': confidence,
        'note': 'Full Bayesian calculation'
    }


def compute_coincident_odds(primary_event, secondary_event, **kwargs):
    """
    Compute posterior odds for two skymaps (e.g., GW vs EM localization) following
    the coincidence framework described in "Coincident Detection Significance in
    Multimessenger Astronomy".
    """
    spatial_calc = SpatialOverlap()
    radial_calc = RadialOverlap()

    I_omega = kwargs.get('I_omega')
    if I_omega is None:
        I_omega = spatial_calc.compute_map_overlap(primary_event, secondary_event)

    I_dl = kwargs.get('I_dl')
    if I_dl is None:
        I_dl = radial_calc.compute_map_overlap(primary_event, secondary_event)

    I_t = kwargs.get('I_t', 1.0)

    prior_odds = kwargs.get('prior_odds', 1.0)
    P_chance = kwargs.get('chance_coincidence_rate', 1e-4)

    bayes_factor = (I_omega * I_dl * I_t) / P_chance
    posterior_odds = prior_odds * bayes_factor
    log_posterior_odds = np.log10(posterior_odds) if posterior_odds > 0 else -np.inf
    confidence = 1 - 1 / (1 + posterior_odds)

    return {
        'gw_file': getattr(primary_event, 'skymap_path', None),
        'secondary_gw_file': getattr(secondary_event, 'skymap_path', None),
        'I_omega': I_omega,
        'I_dl': I_dl,
        'I_t': I_t,
        'bayes_factor': bayes_factor,
        'prior_odds': prior_odds,
        'posterior_odds': posterior_odds,
        'log_posterior_odds': log_posterior_odds,
        'P_chance': P_chance,
        'associated': posterior_odds > 1.0,
        'confidence': confidence,
        'mode': 'skymap_coincidence',
        'note': 'Skymap coincidence Bayesian calculation'
    }