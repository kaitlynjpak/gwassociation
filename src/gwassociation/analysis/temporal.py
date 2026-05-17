"""Temporal compatibility models for GW--EM counterpart candidates.

The :class:`TemporalOverlap` calculator compares the elapsed time between a
primary event and a transient against simple phenomenological delay models for
kilonovae, short gamma-ray bursts, and afterglows.
"""

# temporal overlap integrals
import numpy as np

class TemporalOverlap:
    '''Calculate temporal overlap integral I_t between GW and EM times'''
    
    @staticmethod
    def compute(gw_event, em_transient, model: str = 'kilonova') -> float:
        '''
        Calculate temporal overlap integral I_t
        
        Parameters
        ----------
        gw_event : GWEvent
            Event-like object with a merger time.
        em_transient : Transient
            Candidate-like object with a detection time.
        model : str
            Counterpart model. Supported values are ``'kilonova'``, ``'grb'``,
            and ``'afterglow'``.

        Returns
        -------
        float
            Temporal overlap probability.
        '''
        event_time = getattr(em_transient, "time", None)
        if event_time is None:
            return 1.0  # No time constraint
        
        dt = event_time - gw_event.event_time  # Time delay
        
        if model == 'kilonova':
            # Kilonova light curve model (simplified)
            t_peak = 1.0 * 86400  # 1 day in seconds
            sigma_t = 2.0 * 86400  # 2 day width
            
            # Log-normal distribution for rise and exponential decay
            if dt > 0:
                I_t = np.exp(-0.5 * (np.log(dt/t_peak))**2 / (sigma_t/t_peak)**2)
            else:
                I_t = 0.0
                
        elif model == 'grb':
            # GRB prompt emission
            t_peak = 2.0  # seconds
            sigma_t = 5.0
            I_t = np.exp(-0.5 * (dt - t_peak)**2 / sigma_t**2)
            
        elif model == 'afterglow':
            # GRB afterglow
            t_peak = 1.0 * 86400  # 1 day
            if dt > 0:
                # Power-law decay
                I_t = (dt / t_peak) ** (-0.7) if dt > t_peak else 1.0
            else:
                I_t = 0.0
        else:
            I_t = 1.0
            
        return I_t