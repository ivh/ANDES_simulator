"""MOSAIC instrument configuration."""

INSTRUMENT_NAME = 'MOSAIC'

MOSAIC_BANDS = ['LR-blue', 'LR-red', 'LR-J', 'LR-H']

SUBSLIT_CHOICES = ['all', 'even', 'odd']

# Single VPH order for all bands
BAND_ORDER_ESTIMATES = {
    'LR-blue': 1,
    'LR-red': 1,
    'LR-J': 1,
    'LR-H': 1,
}

DEFAULT_SCALING = {
    'LR-blue': 1e3,
    'LR-red': 1e3,
    'LR-J': 1e3,
    'LR-H': 1e3,
}

DEFAULT_FINESSE = {
    'LR-blue': 23,
    'LR-red': 23,
    'LR-J': 23,
    'LR-H': 23,
}

TELESCOPE = {
    'primary_diameter': 39.3,
    'central_obstruction': 4.09,
}

_VIS_COMMON = {
    'instrument_name': 'MOSAIC',
    'n_fibers': 980,
    'detector_size': (12288, 12288),
    'pixel_size': 15,
    'diffraction_orders': [1],
    'sampling': 5.8,
    'skip_fibers': [],
    'fiber_config_key': 'MOSAIC_VIS',
    'bundle_size': 7,
    'n_bundles': 140,
}

_NIR_COMMON = {
    'instrument_name': 'MOSAIC',
    'n_fibers': 630,
    'detector_size': (4096, 4096),
    'pixel_size': 15,
    'diffraction_orders': [1],
    'sampling': 4.0,
    'skip_fibers': [],
    'fiber_config_key': 'MOSAIC_NIR',
    'bundle_size': 7,
    'n_bundles': 90,
}

INSTRUMENTS = {
    'LR-blue': {
        **_VIS_COMMON,
        'wavelength_range': (390, 625),
        'hdf_models': {'default': 'MOSAIC_VIS_LR_Blue'},
    },
    'LR-red': {
        **_VIS_COMMON,
        'wavelength_range': (595, 952),
        'hdf_models': {'default': 'MOSAIC_VIS_LR_Red'},
    },
    'LR-J': {
        **_NIR_COMMON,
        'wavelength_range': (950, 1340),
        'hdf_models': {'default': 'MOSAIC_NIR_LR_J'},
    },
    'LR-H': {
        **_NIR_COMMON,
        'wavelength_range': (1430, 1800),
        'hdf_models': {'default': 'MOSAIC_NIR_LR_H'},
    },
}

FIBER_CONFIG = {
    'MOSAIC_VIS': {
        'n_fibers': 980,
        'slits': {},
    },
    'MOSAIC_NIR': {
        'n_fibers': 630,
        'slits': {},
    },
}
