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
            'default': 'ANDES_U_v88'
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
            'default': 'ANDES_B_v88'
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
            'default': 'ANDES_V_v88'
        },
        'fp_spectrum': 'FP_simulation_UBV_finesse_26.csv',
        'sampling': 4.0,
        'skip_fibers': [32, 35]
    }
}

# Cache for wavelength ranges read from HDF files
_wavelength_range_cache: Dict[str, tuple] = {}

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


def _read_wavelength_range_from_hdf(hdf_path: Path) -> tuple:
    """
    Read wavelength range from an HDF model file.

    Returns
    -------
    tuple
        (wl_min, wl_max) in nm
    """
    import h5py

    with h5py.File(hdf_path, 'r') as f:
        ccd = f['CCD_1']
        fiber_key = next(k for k in ccd.keys() if k.startswith('fiber_'))
        fiber = ccd[fiber_key]

        wavelengths = []
        for key in fiber.keys():
            if key.startswith('psf_order'):
                psf_grp = fiber[key]
                for wl_key in psf_grp.keys():
                    if wl_key.startswith('wavelength_'):
                        wl = float(wl_key.replace('wavelength_', ''))
                        wavelengths.append(wl)

    if not wavelengths:
        raise ValueError(f"No wavelength data found in {hdf_path}")

    # Convert micrometers to nm
    return min(wavelengths) * 1000, max(wavelengths) * 1000


def get_band_wavelength_range(band: str, project_root: Path = None) -> tuple:
    """
    Get wavelength range for a band from its default HDF model.

    Results are cached to avoid repeated file reads.

    Parameters
    ----------
    band : str
        Spectral band name (U, B, V, R, IZ, Y, J, H)
    project_root : Path, optional
        Project root directory

    Returns
    -------
    tuple
        (wl_min, wl_max) in nm
    """
    if band in _wavelength_range_cache:
        return _wavelength_range_cache[band]

    hdf_path = get_hdf_model_path(band, 'default', project_root)
    wl_range = _read_wavelength_range_from_hdf(hdf_path)
    _wavelength_range_cache[band] = wl_range
    return wl_range


def get_all_band_wavelength_ranges(project_root: Path = None) -> Dict[str, tuple]:
    """
    Get wavelength ranges for all bands.

    Returns
    -------
    dict
        Mapping of band name to (wl_min, wl_max) in nm
    """
    return {band: get_band_wavelength_range(band, project_root)
            for band in INSTRUMENTS.keys()}


def infer_band_from_hdf(hdf_path: Path, project_root: Path = None) -> str:
    """
    Infer spectral band from HDF file by reading wavelength coverage.

    Parameters
    ----------
    hdf_path : Path
        Path to HDF model file
    project_root : Path, optional
        Project root directory (for reading default HDF models)

    Returns
    -------
    str
        Inferred band name

    Raises
    ------
    ValueError
        If band cannot be determined from wavelength range
    """
    wl_min, wl_max = _read_wavelength_range_from_hdf(hdf_path)
    wl_center = (wl_min + wl_max) / 2

    for band in INSTRUMENTS.keys():
        low, high = get_band_wavelength_range(band, project_root)
        if low <= wl_center <= high:
            return band

    raise ValueError(f"Cannot infer band from wavelength range {wl_min:.0f}-{wl_max:.0f}nm")


def infer_band_from_wavelengths(wl_min: float = None, wl_max: float = None,
                                project_root: Path = None) -> str:
    """
    Infer spectral band from wavelength limits.

    Parameters
    ----------
    wl_min : float, optional
        Minimum wavelength in nm
    wl_max : float, optional
        Maximum wavelength in nm
    project_root : Path, optional
        Project root directory

    Returns
    -------
    str
        Inferred band name

    Raises
    ------
    ValueError
        If band cannot be uniquely determined from wavelength range
    """
    if wl_min is None and wl_max is None:
        raise ValueError("At least one of wl_min or wl_max must be provided")

    matching_bands = []
    for band in INSTRUMENTS.keys():
        low, high = get_band_wavelength_range(band, project_root)
        req_min = wl_min if wl_min is not None else low
        req_max = wl_max if wl_max is not None else high
        if req_min <= high and req_max >= low:
            matching_bands.append(band)

    if len(matching_bands) == 0:
        raise ValueError(
            f"Wavelength range {wl_min}-{wl_max}nm does not match any band")
    if len(matching_bands) > 1:
        raise ValueError(
            f"Wavelength range {wl_min}-{wl_max}nm is ambiguous, "
            f"could match: {', '.join(matching_bands)}. Use --band to specify.")

    return matching_bands[0]


def validate_wavelength_range(band: str, wl_min: float = None, wl_max: float = None,
                              project_root: Path = None) -> None:
    """
    Validate that wavelength limits fall within the band's range.

    Parameters
    ----------
    band : str
        Spectral band name
    wl_min : float, optional
        Minimum wavelength in nm
    wl_max : float, optional
        Maximum wavelength in nm
    project_root : Path, optional
        Project root directory

    Raises
    ------
    ValueError
        If wavelength limits are outside the band's range
    """
    if wl_min is None and wl_max is None:
        return

    if band not in INSTRUMENTS:
        raise ValueError(f"Unknown band '{band}'")

    band_min, band_max = get_band_wavelength_range(band, project_root)

    if wl_min is not None and wl_min < band_min:
        raise ValueError(
            f"wl_min={wl_min}nm is below {band}-band range ({band_min:.1f}-{band_max:.1f}nm)")
    if wl_max is not None and wl_max > band_max:
        raise ValueError(
            f"wl_max={wl_max}nm is above {band}-band range ({band_min:.1f}-{band_max:.1f}nm)")
    if wl_min is not None and wl_max is not None and wl_min >= wl_max:
        raise ValueError(f"wl_min={wl_min}nm must be less than wl_max={wl_max}nm")


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