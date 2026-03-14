"""
Instrument configuration registry.

Merges band definitions from instrument-specific modules (andes, mosaic)
into a single namespace. All downstream code imports from here.
"""

from pathlib import Path
from typing import Dict, List, Any

from . import andes as _andes
from . import mosaic as _mosaic

_INSTRUMENT_MODULES = [_andes, _mosaic]

# Merged dicts — band names must be unique across instruments
INSTRUMENTS = {**_andes.INSTRUMENTS, **_mosaic.INSTRUMENTS}
BAND_ORDER_ESTIMATES = {**_andes.BAND_ORDER_ESTIMATES, **_mosaic.BAND_ORDER_ESTIMATES}
DEFAULT_SCALING = {**_andes.DEFAULT_SCALING, **_mosaic.DEFAULT_SCALING}
DEFAULT_FINESSE = {**_andes.DEFAULT_FINESSE, **_mosaic.DEFAULT_FINESSE}
FIBER_CONFIG = {**_andes.FIBER_CONFIG, **_mosaic.FIBER_CONFIG}

# Both on ELT; use ANDES telescope as canonical (identical values)
TELESCOPE = _andes.TELESCOPE.copy()

# Cache for wavelength ranges read from HDF files
_wavelength_range_cache: Dict[str, tuple] = {}


def get_instrument_name(band: str) -> str:
    """Get the instrument name for a band."""
    for mod in _INSTRUMENT_MODULES:
        if band in mod.INSTRUMENTS:
            return mod.INSTRUMENT_NAME
    raise ValueError(f"Unknown band '{band}'. Available: {list(INSTRUMENTS.keys())}")


def get_instrument_config(band: str) -> Dict[str, Any]:
    """Get instrument configuration for a specific band."""
    if band not in INSTRUMENTS:
        raise ValueError(f"Unknown band '{band}'. Available: {list(INSTRUMENTS.keys())}")

    config = INSTRUMENTS[band].copy()
    config['telescope'] = TELESCOPE.copy()

    fc_key = config.get('fiber_config_key')
    if fc_key and fc_key in FIBER_CONFIG:
        config['fiber_config'] = FIBER_CONFIG[fc_key].copy()

    ifu_key = config.get('ifu_config_key')
    if ifu_key and ifu_key in FIBER_CONFIG:
        config['ifu_config'] = FIBER_CONFIG[ifu_key].copy()

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
    """Get the full path to an HDF model file."""
    if project_root is None:
        project_root = Path(__file__).parent.parent.parent

    config = get_instrument_config(band)

    if model_type == 'default':
        model_name = config['hdf_models']['default']
    elif model_type in config['hdf_models']:
        if isinstance(config['hdf_models'][model_type], list):
            model_name = config['hdf_models'][model_type][0]
        else:
            model_name = config['hdf_models'][model_type]
    else:
        model_name = model_type

    return project_root / 'HDF' / f'{model_name}.hdf'


def _read_wavelength_range_from_hdf(hdf_path: Path) -> tuple:
    """Read wavelength range from an HDF model file. Returns (wl_min, wl_max) in nm."""
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

    return min(wavelengths) * 1000, max(wavelengths) * 1000


def get_band_wavelength_range(band: str, project_root: Path = None) -> tuple:
    """Get wavelength range for a band from its default HDF model. Returns (wl_min, wl_max) in nm."""
    if band in _wavelength_range_cache:
        return _wavelength_range_cache[band]

    hdf_path = get_hdf_model_path(band, 'default', project_root)
    wl_range = _read_wavelength_range_from_hdf(hdf_path)
    _wavelength_range_cache[band] = wl_range
    return wl_range


def get_all_band_wavelength_ranges(project_root: Path = None) -> Dict[str, tuple]:
    """Get wavelength ranges for all bands."""
    return {band: get_band_wavelength_range(band, project_root)
            for band in INSTRUMENTS.keys()}


def infer_band_from_hdf(hdf_path: Path, project_root: Path = None) -> str:
    """Infer spectral band from HDF file by reading wavelength coverage."""
    wl_min, wl_max = _read_wavelength_range_from_hdf(hdf_path)
    wl_center = (wl_min + wl_max) / 2

    for band in INSTRUMENTS.keys():
        low, high = get_band_wavelength_range(band, project_root)
        if low <= wl_center <= high:
            return band

    raise ValueError(f"Cannot infer band from wavelength range {wl_min:.0f}-{wl_max:.0f}nm")


def infer_band_from_wavelengths(wl_min: float = None, wl_max: float = None,
                                project_root: Path = None) -> str:
    """Infer spectral band from wavelength limits."""
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
    """Validate that wavelength limits fall within the band's range."""
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
