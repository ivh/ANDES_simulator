"""MOSAIC instrument configuration."""

INSTRUMENT_NAME = 'MOSAIC'

MOSAIC_BANDS = ['VIS']

SUBSLIT_CHOICES = ['all', 'even', 'odd']

# Single VPH order
BAND_ORDER_ESTIMATES = {
    'VIS': 1,
}

DEFAULT_SCALING = {
    'VIS': 1e3,
}

DEFAULT_FINESSE = {
    'VIS': 23,
}

TELESCOPE = {
    'primary_diameter': 39.3,
    'central_obstruction': 4.09,
}

INSTRUMENTS = {
    'VIS': {
        'instrument_name': 'MOSAIC',
        'n_fibers': 980,
        'detector_size': (12288, 12288),
        'pixel_size': 15,
        'hdf_models': {
            'default': 'MOSAIC_VIS_LR_Blue_full',
            'test_2bundle': 'MOSAIC_VIS_LR_Blue',
        },
        'diffraction_orders': [1],
        'sampling': 5.8,
        'skip_fibers': [],
        'fiber_config_key': 'MOSAIC_VIS',
        'bundle_size': 7,
        'n_bundles': 140,
    },
}

FIBER_CONFIG = {
    'MOSAIC_VIS': {
        'n_fibers': 980,
        'slits': {},
    },
}
