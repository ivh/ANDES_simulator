"""MOSAIC instrument configuration."""

INSTRUMENT_NAME = 'MOSAIC'

# Band names follow the E2E `ESO INS MODE` convention: <config>_<resolution>,
# e.g. B_LR (blue low-res), B1_HR (blue high-res setting 1), J_LR, H_HR.
MOSAIC_BANDS = ['B_LR', 'R_LR', 'J_LR', 'H_LR',
                'B1_HR', 'R1_HR', 'B2_HR', 'R2_HR', 'H_HR']

SUBSLIT_CHOICES = ['all', 'even', 'odd']

# Single VPH order for all bands
BAND_ORDER_ESTIMATES = {
    'B_LR': 1, 'R_LR': 1, 'J_LR': 1, 'H_LR': 1,
    'B1_HR': 1, 'R1_HR': 1, 'B2_HR': 1, 'R2_HR': 1, 'H_HR': 1,
}

DEFAULT_SCALING = {
    'B_LR': 1e3, 'R_LR': 1e3, 'J_LR': 1e3, 'H_LR': 1e3,
    'B1_HR': 1e3, 'R1_HR': 1e3, 'B2_HR': 1e3, 'R2_HR': 1e3, 'H_HR': 1e3,
}

DEFAULT_FINESSE = {
    'B_LR': 23, 'R_LR': 23, 'J_LR': 23, 'H_LR': 23,
    'B1_HR': 23, 'R1_HR': 23, 'B2_HR': 23, 'R2_HR': 23, 'H_HR': 23,
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
    'B_LR': {
        **_VIS_LR,
        'wavelength_range': (390, 625),
        'hdf_models': {'default': 'MOSAIC_VIS_LR_Blue'},
    },
    'R_LR': {
        **_VIS_LR,
        'wavelength_range': (595, 952),
        'hdf_models': {'default': 'MOSAIC_VIS_LR_Red'},
    },
    # VIS HR
    'B1_HR': {
        **_VIS_HR,
        'wavelength_range': (505, 580),
        'hdf_models': {'default': 'MOSAIC_VIS_HR_B1'},
    },
    'R1_HR': {
        **_VIS_HR,
        'wavelength_range': (610, 680),
        'hdf_models': {'default': 'MOSAIC_VIS_HR_R1'},
    },
    'B2_HR': {
        **_VIS_HR,
        'wavelength_range': (393, 458),
        'hdf_models': {'default': 'MOSAIC_VIS_HR_B2'},
    },
    'R2_HR': {
        **_VIS_HR,
        'wavelength_range': (765, 890),
        'hdf_models': {'default': 'MOSAIC_VIS_HR_R2'},
    },
    # NIR LR
    'J_LR': {
        **_NIR,
        'wavelength_range': (950, 1340),
        'hdf_models': {'default': 'MOSAIC_NIR_LR_J'},
    },
    'H_LR': {
        **_NIR,
        'wavelength_range': (1430, 1800),
        'hdf_models': {'default': 'MOSAIC_NIR_LR_H'},
    },
    # NIR HR
    'H_HR': {
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
