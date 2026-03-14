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
        'n_fibers': 75,
        'detector_size': (4096, 4096),
        'pixel_size': 15,
        'hdf_models': {
            'default': 'MOSAIC_VIS',
        },
        'diffraction_orders': [1],
        'sampling': 4.0,
        'skip_fibers': [],
        'fiber_config_key': 'MOSAIC_VIS',
    },
}

FIBER_CONFIG = {
    'MOSAIC_VIS': {
        'n_fibers': 75,
        'slits': {},
    },
}
