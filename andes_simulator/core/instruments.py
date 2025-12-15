"""
Instrument configurations for ANDES spectrograph bands.

Defines the fiber counts, detector sizes, HDF models, and other 
band-specific parameters for all ANDES spectral channels.
"""

from pathlib import Path
from typing import Dict, List, Any

# Instrument configurations for each spectral band
INSTRUMENTS = {
    # Near-infrared bands (75 fibers)
    'Y': {
        'n_fibers': 75,
        'detector_size': (4096, 4096),
        'pixel_size': 15,  # microns
        'hdf_models': {
            'default': 'ANDES_75fibre_Y',
            'with_fiber_eff': 'ANDES_Y01_wFiberEff'
        },
        'zemax_files': {
            'Y': 'HIRES_Y_21jan2022_sconf_noap.zmx'
        },
        'diffraction_orders': list(range(109, 127)),
        'fp_spectrum': 'FP_simulation_YJH_finesse_26.csv',
        'sampling': 2.1,
        'skip_fibers': [3, 4, 36, 37, 39, 40, 72, 73]
    },
    'J': {
        'n_fibers': 75,
        'detector_size': (4096, 4096),
        'pixel_size': 15,
        'hdf_models': {
            'default': 'ANDES_75fibre_J'
        },
        'zemax_files': {
            'J': 'HIRES_J_21jan2022_sconf_noap.zmx'
        },
        'diffraction_orders': list(range(90, 108)),
        'fp_spectrum': 'FP_simulation_YJH_finesse_26.csv',
        'sampling': 2.1,
        'skip_fibers': [3, 4, 36, 37, 39, 40, 72, 73]
    },
    'H': {
        'n_fibers': 75,
        'detector_size': (4096, 4096),
        'pixel_size': 15,
        'hdf_models': {
            'default': 'ANDES_75fibre_H'
        },
        'zemax_files': {
            'H': 'HIRES_H_21jan2022_sconf.zmx'
        },
        'diffraction_orders': list(range(68, 83)),
        'fp_spectrum': 'FP_simulation_YJH_finesse_26.csv',
        'sampling': 2.1,
        'skip_fibers': [3, 4, 36, 37, 39, 40, 72, 73]
    },
    
    # Optical bands (66 fibers)
    'R': {
        'n_fibers': 66,
        'detector_size': (9216, 9232),
        'pixel_size': 15,
        'hdf_models': {
            'default': 'ANDES_123_R3',
            'thermal_variants': [
                'ANDES_full_F18A33_win_jmr_MC_T0019_Rband_p0',
                'Andes_full_F18A33_win_jmr_MC_T0108_Rband_P0_cfg1'
            ]
        },
        'fp_spectrum': 'FP_simulation_RIZ_finesse_23.csv',
        'sampling': 4.0,
        'skip_fibers': [32, 35]
    },
    'IZ': {
        'n_fibers': 66,
        'detector_size': (9216, 9232),
        'pixel_size': 15,
        'hdf_models': {
            'default': 'ANDES_123_IZ3',
            'thermal_variants': [
                'Andes_F18A33_VM246aa_win_jmr9_MC_T0045_IZband_P0_cf1',
                'Andes_full_F18A33_win_jmr_MC_T0028_IZband_P0'
            ]
        },
        'fp_spectrum': 'FP_simulation_RIZ_finesse_23.csv',
        'sampling': 4.0,
        'skip_fibers': [32, 35]
    },
    'U': {
        'n_fibers': 66,
        'detector_size': (9216, 9232),
        'pixel_size': 15,
        'hdf_models': {
            'default': 'ANDES_123_U3'
        },
        'fp_spectrum': 'FP_simulation_UBV_finesse_23.csv',
        'sampling': 4.0,
        'skip_fibers': [32, 35]
    },
    'B': {
        'n_fibers': 66,
        'detector_size': (9216, 9232),
        'pixel_size': 15,
        'hdf_models': {
            'default': 'ANDES_123_B3'
        },
        'fp_spectrum': 'FP_simulation_UBV_finesse_26.csv',
        'sampling': 4.0,
        'skip_fibers': [32, 35]
    },
    'V': {
        'n_fibers': 66,
        'detector_size': (9216, 9232),
        'pixel_size': 15,
        'hdf_models': {
            'default': 'ANDES_123_V3'
        },
        'fp_spectrum': 'FP_simulation_UBV_finesse_26.csv',
        'sampling': 4.0,
        'skip_fibers': [32, 35]
    }
}

# Telescope configuration (ELT)
TELESCOPE = {
    'primary_diameter': 39.3,  # meters
    'central_obstruction': 4.09  # meters
}

# Standard fiber sizes and configurations
FIBER_CONFIG = {
    'UBVRIZ': {
        'n_fibers': 66,
        'slits': {
            'slitA': list(range(1, 32)),
            'slitB': list(range(35, 66)),
            'cal_fibers': [33, 34]
        }
    },
    'YJH_SL': {
        'fiber_size': 474,  # micrometers
        'n_fibers': 75,
        'slits': {
            'slitA': list(range(4, 35)),
            'slitB': list(range(42, 73)),
            'cal_fibers': [1,37,38,39,75]
        }
    },
    'YJH_IFU': {
        'fiber_size': 474,  # micrometers
        'n_fibers': 75,
        'slits': {
            'ring0': [3],
            'ring1': [6,8,10,12,14,16],
            'ring2': list(range(18, 18+12)),
            'ring3': list(range(31, 31+18)),
            'ring4': list(range(50, 50+24)),
            'cal_fibers': [1,75]
        }
    }
}


def get_instrument_config(band: str) -> Dict[str, Any]:
    """
    Get instrument configuration for a specific band.
    
    Parameters
    ----------
    band : str
        Spectral band name (Y, J, H, R, IZ, U, B, V)
        
    Returns
    -------
    dict
        Complete instrument configuration including telescope and fiber setup
    """
    if band not in INSTRUMENTS:
        raise ValueError(f"Unknown band '{band}'. Available bands: {list(INSTRUMENTS.keys())}")
    
    config = INSTRUMENTS[band].copy()
    config['telescope'] = TELESCOPE.copy()
    
    # Add appropriate fiber configuration
    if band in ['Y', 'J', 'H']:
        config['fiber_config'] = FIBER_CONFIG['YJH_SL'].copy()
    else:
        config['fiber_config'] = FIBER_CONFIG['UBVRIZ'].copy()
    
    return config


def get_all_bands() -> List[str]:
    """Get list of all available spectral bands."""
    return list(INSTRUMENTS.keys())


def get_nir_bands() -> List[str]:
    """Get list of near-infrared bands (YJH)."""
    return ['Y', 'J', 'H']


def get_optical_bands() -> List[str]:
    """Get list of optical bands (RIUZV).""" 
    return ['R', 'IZ', 'U', 'B', 'V']


def get_hdf_model_path(band: str, model_type: str = 'default', project_root: Path = None) -> Path:
    """
    Get the full path to an HDF model file.
    
    Parameters
    ----------
    band : str
        Spectral band name
    model_type : str
        Type of model ('default', 'with_fiber_eff', or specific thermal variant)
    project_root : Path, optional
        Project root directory. If None, uses current file location.
        
    Returns
    -------
    Path
        Full path to HDF model file
    """
    if project_root is None:
        project_root = Path(__file__).parent.parent.parent
    
    config = get_instrument_config(band)
    
    if model_type == 'default':
        model_name = config['hdf_models']['default']
    elif model_type in config['hdf_models']:
        if isinstance(config['hdf_models'][model_type], list):
            # For thermal variants, return the first one as default
            model_name = config['hdf_models'][model_type][0]
        else:
            model_name = config['hdf_models'][model_type]
    else:
        # Assume it's a specific model name
        model_name = model_type
    
    return project_root / 'HDF' / f'{model_name}.hdf'


def get_sed_path(band: str, project_root: Path = None) -> Path:
    """
    Get the path to the Fabry-Perot SED file for a band.
    
    Parameters
    ----------
    band : str
        Spectral band name
    project_root : Path, optional
        Project root directory
        
    Returns
    -------
    Path
        Full path to SED file
    """
    if project_root is None:
        project_root = Path(__file__).parent.parent.parent
    
    config = get_instrument_config(band)
    return project_root / 'SED' / config['fp_spectrum']