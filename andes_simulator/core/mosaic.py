"""MOSAIC instrument configuration."""

INSTRUMENT_NAME = 'MOSAIC'

MOSAIC_BANDS = ['LR-blue', 'LR-red', 'LR-J', 'LR-H',
                'HR-B1', 'HR-R1', 'HR-B2', 'HR-R2', 'HR-H']

SUBSLIT_CHOICES = ['all', 'even', 'odd']

# Single VPH order for all bands
BAND_ORDER_ESTIMATES = {
    'LR-blue': 1, 'LR-red': 1, 'LR-J': 1, 'LR-H': 1,
    'HR-B1': 1, 'HR-R1': 1, 'HR-B2': 1, 'HR-R2': 1, 'HR-H': 1,
}

DEFAULT_SCALING = {
    'LR-blue': 1e3, 'LR-red': 1e3, 'LR-J': 1e3, 'LR-H': 1e3,
    'HR-B1': 1e3, 'HR-R1': 1e3, 'HR-B2': 1e3, 'HR-R2': 1e3, 'HR-H': 1e3,
}

DEFAULT_FINESSE = {
    'LR-blue': 23, 'LR-red': 23, 'LR-J': 23, 'LR-H': 23,
    'HR-B1': 23, 'HR-R1': 23, 'HR-B2': 23, 'HR-R2': 23, 'HR-H': 23,
}

TELESCOPE = {
    'primary_diameter': 39.3,
    'central_obstruction': 4.09,
}

# VIS LR: 140 bundles x 7 fibers = 980, 177um pitch
_VIS_LR = {
    'instrument_name': 'MOSAIC',
    'n_fibers': 980,
    'detector_size': (12288, 12288),
    'pixel_size': 15,
    'diffraction_orders': [1],
    'sampling': 5.8,
    'skip_fibers': [],
    'fiber_config_key': 'MOSAIC_VIS_LR',
    'bundle_size': 7,
    'n_bundles': 140,
}

# VIS HR: 60 bundles x 19 fibers = 1140, 152um pitch
_VIS_HR = {
    'instrument_name': 'MOSAIC',
    'n_fibers': 1140,
    'detector_size': (12288, 12288),
    'pixel_size': 15,
    'diffraction_orders': [1],
    'sampling': 4.7,
    'skip_fibers': [],
    'fiber_config_key': 'MOSAIC_VIS_HR',
    'bundle_size': 19,
    'n_bundles': 60,
}

# NIR: 90 bundles x 7 fibers = 630, 203um pitch
_NIR = {
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
    # VIS LR
    'LR-blue': {
        **_VIS_LR,
        'wavelength_range': (390, 625),
        'hdf_models': {'default': 'MOSAIC_VIS_LR_Blue'},
    },
    'LR-red': {
        **_VIS_LR,
        'wavelength_range': (595, 952),
        'hdf_models': {'default': 'MOSAIC_VIS_LR_Red'},
    },
    # VIS HR
    'HR-B1': {
        **_VIS_HR,
        'wavelength_range': (505, 580),
        'hdf_models': {'default': 'MOSAIC_VIS_HR_B1'},
    },
    'HR-R1': {
        **_VIS_HR,
        'wavelength_range': (610, 680),
        'hdf_models': {'default': 'MOSAIC_VIS_HR_R1'},
    },
    'HR-B2': {
        **_VIS_HR,
        'wavelength_range': (393, 458),
        'hdf_models': {'default': 'MOSAIC_VIS_HR_B2'},
    },
    'HR-R2': {
        **_VIS_HR,
        'wavelength_range': (765, 890),
        'hdf_models': {'default': 'MOSAIC_VIS_HR_R2'},
    },
    # NIR LR
    'LR-J': {
        **_NIR,
        'wavelength_range': (950, 1340),
        'hdf_models': {'default': 'MOSAIC_NIR_LR_J'},
    },
    'LR-H': {
        **_NIR,
        'wavelength_range': (1430, 1800),
        'hdf_models': {'default': 'MOSAIC_NIR_LR_H'},
    },
    # NIR HR
    'HR-H': {
        **_NIR,
        'wavelength_range': (1520, 1620),
        'hdf_models': {'default': 'MOSAIC_NIR_HR'},
    },
}

FIBER_CONFIG = {
    'MOSAIC_VIS_LR': {
        'n_fibers': 980,
        'slits': {},
    },
    'MOSAIC_VIS_HR': {
        'n_fibers': 1140,
        'slits': {},
    },
    'MOSAIC_NIR': {
        'n_fibers': 630,
        'slits': {},
    },
}
