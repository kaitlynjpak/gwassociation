"""
Data ingestion utilities for GW-EM association
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from typing import Dict, List, Union, Optional

try:
    from astropy.table import Table
    from astropy.io import fits
    ASTROPY_AVAILABLE = True
except ImportError:
    ASTROPY_AVAILABLE = False

# --- Back-compat shim so __init__ can import this safely ---
def load_gw_posteriors(path, *args, **kwargs):
    """
    Back-compat loader. Delegates to an existing 3D GW skymap/ posterior loader.
    Prefer calling the concrete loader directly (e.g., load_gw_3d) in new code.
    """
    # Try the most likely concrete loaders you already have:
    if 'load_gw_3d' in globals():
        return load_gw_3d(path, *args, **kwargs)
    if 'load_skymap' in globals():
        return load_skymap(path, *args, **kwargs)
    raise NotImplementedError(
        "load_gw_posteriors() shim: no underlying loader found "
        "(expected load_gw_3d or load_skymap in ingest.py)."
    )

def ingest_transient_list(filename: str) -> List[Dict]:
    """
    Ingest a list of transient candidates from various formats
    
    Parameters
    ----------
    filename : str
        Path to file containing transient data
        
    Returns
    -------
    list
        List of dictionaries containing transient information
    """
    filename = str(filename)  # Convert Path objects to string
    
    if filename.endswith('.csv'):
        df = pd.read_csv(filename)
        # Replace NaN with None for JSON compatibility
        df = df.where(pd.notnull(df), None)
        return df.to_dict('records')
    
    elif filename.endswith('.json'):
        with open(filename, 'r') as f:
            return json.load(f)
    
    elif filename.endswith('.fits') and ASTROPY_AVAILABLE:
        table = Table.read(filename)
        return table.to_pandas().to_dict('records')
    
    elif filename.endswith('.txt') or filename.endswith('.dat'):
        # Try to read as whitespace-delimited file
        df = pd.read_csv(filename, delim_whitespace=True)
        df = df.where(pd.notnull(df), None)
        return df.to_dict('records')
    
    else:
        raise ValueError(f"Unsupported file format: {filename}")


def load_gw_posteriors(filename: str) -> Dict:
    """
    Load GW posterior samples or parameters
    
    Parameters
    ----------
    filename : str
        Path to file containing GW posteriors
        
    Returns
    -------
    dict
        Dictionary containing posterior samples or summary statistics
    """
    filename = str(filename)
    
    if filename.endswith('.json'):
        with open(filename, 'r') as f:
            return json.load(f)
    
    elif filename.endswith('.h5') or filename.endswith('.hdf5'):
        # HDF5 format often used for posterior samples
        try:
            import h5py
            posteriors = {}
            with h5py.File(filename, 'r') as f:
                # Try to extract common parameters
                if 'posterior_samples' in f:
                    group = f['posterior_samples']
                    for key in group.keys():
                        posteriors[key] = group[key][:]
                else:
                    # Fallback: get all datasets
                    for key in f.keys():
                        if isinstance(f[key], h5py.Dataset):
                            posteriors[key] = f[key][:]
            return posteriors
        except ImportError:
            raise ImportError("h5py required to read HDF5 posterior files")
    
    elif filename.endswith('.dat') or filename.endswith('.txt'):
        # ASCII format posteriors
        df = pd.read_csv(filename, delim_whitespace=True)
        return df.to_dict('list')
    
    elif filename.endswith('.fits') and ASTROPY_AVAILABLE:
        # FITS table format
        from astropy.table import Table
        table = Table.read(filename)
        return {col: table[col].data.tolist() for col in table.colnames}
    
    else:
        # Default: try to load as JSON
        try:
            with open(filename, 'r') as f:
                return json.load(f)
        except:
            raise ValueError(f"Cannot read posterior file: {filename}")


def ingest_gracedb_data(superevent_id: str) -> Dict:
    """
    Download and ingest data from GraceDB (placeholder implementation)
    
    Parameters
    ----------
    superevent_id : str
        GraceDB superevent ID (e.g., 'S190425z')
        
    Returns
    -------
    dict
        Dictionary containing event information and file URLs
    """
    # This is a placeholder implementation
    # In production, this would use the GraceDB API
    
    base_url = "https://gracedb.ligo.org/superevents"
    
    return {
        'superevent_id': superevent_id,
        'skymap_url': f'{base_url}/{superevent_id}/files/bayestar.fits',
        'classification': None,
        'event_time': None,
        'far': None,  # False alarm rate
        'instruments': None,
        'note': 'Placeholder - implement GraceDB API access for real data'
    }


def format_for_analysis(transient_dict: Dict) -> Dict:
    """
    Format transient dictionary for analysis, ensuring required fields
    
    Parameters
    ----------
    transient_dict : dict
        Raw transient data
        
    Returns
    -------
    dict
        Formatted dictionary ready for analysis
    """
    required = ['ra', 'dec']
    optional = ['z', 'z_err', 'time', 'time_err', 'magnitude', 'filter_band', 'name']
    
    formatted = {}
    
    # Check required fields
    for key in required:
        if key not in transient_dict:
            raise ValueError(f"Missing required field: {key}")
        formatted[key] = float(transient_dict[key]) if transient_dict[key] is not None else None
    
    # Add optional fields
    for key in optional:
        value = transient_dict.get(key, None)
        if value is not None:
            # Convert to appropriate type
            if key in ['z', 'z_err', 'time', 'time_err', 'magnitude']:
                formatted[key] = float(value) if value is not None else None
            else:
                formatted[key] = value
        else:
            formatted[key] = None
    
    return formatted


def load_skymap_metadata(skymap_file: str) -> Dict:
    """
    Load metadata from a GW skymap file
    
    Parameters
    ----------
    skymap_file : str
        Path to skymap FITS file
        
    Returns
    -------
    dict
        Metadata from skymap header
    """
    metadata = {}
    
    if ASTROPY_AVAILABLE:
        try:
            from astropy.io import fits
            with fits.open(skymap_file) as hdul:
                header = hdul[1].header  # Extension 1 usually has the skymap
                
                # Extract common metadata
                metadata['gps_time'] = header.get('GPS_TIME', None)
                metadata['distance'] = header.get('DISTMEAN', None)
                metadata['distance_std'] = header.get('DISTSTD', None)
                metadata['instruments'] = header.get('INSTRUME', None)
                metadata['creator'] = header.get('CREATOR', 'unknown')
                metadata['origin'] = header.get('ORIGIN', 'unknown')
                
        except Exception as e:
            print(f"Warning: Could not read skymap metadata: {e}")
    
    return metadata


def batch_ingest_transients(file_list: List[str]) -> List[Dict]:
    """
    Ingest transients from multiple files
    
    Parameters
    ----------
    file_list : list
        List of file paths
        
    Returns
    -------
    list
        Combined list of all transients
    """
    all_transients = []
    
    for filename in file_list:
        try:
            transients = ingest_transient_list(filename)
            all_transients.extend(transients)
            print(f"Loaded {len(transients)} transients from {filename}")
        except Exception as e:
            print(f"Warning: Could not load {filename}: {e}")
    
    return all_transients


def save_results(results: Dict, output_file: str, format: str = 'json'):
    """
    Save analysis results to file
    
    Parameters
    ----------
    results : dict
        Analysis results
    output_file : str
        Output file path
    format : str
        Output format ('json', 'csv', 'fits')
    """
    output_file = str(output_file)
    
    # Convert numpy types to Python native types for JSON serialization
    def convert_numpy(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.int32, np.int64)):
            return int(obj)
        elif isinstance(obj, dict):
            return {k: convert_numpy(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_numpy(item) for item in obj]
        else:
            return obj
    
    results = convert_numpy(results)
    
    if format == 'json' or output_file.endswith('.json'):
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
    
    elif format == 'csv' or output_file.endswith('.csv'):
        df = pd.DataFrame([results] if isinstance(results, dict) else results)
        df.to_csv(output_file, index=False)
    
    elif format == 'fits' and ASTROPY_AVAILABLE:
        from astropy.table import Table
        table = Table(results)
        table.write(output_file, format='fits', overwrite=True)
    
    else:
        raise ValueError(f"Unsupported output format: {format}")
    
    print(f"Results saved to: {output_file}")


# Make sure all functions are available for import
__all__ = [
    'ingest_transient_list',
    'load_gw_posteriors',
    'ingest_gracedb_data',
    'format_for_analysis',
    'load_skymap_metadata',
    'batch_ingest_transients',
    'save_results'
]