"""ANDES instrument configuration."""

INSTRUMENT_NAME = 'ANDES'

ANDES_BANDS = ['U', 'B', 'V', 'R', 'IZ', 'Y', 'J', 'H']

SUBSLIT_CHOICES = ['all', 'even', 'odd', 'slitA', 'slitB', 'cal_sl', 'cal_ifu',
                   'ifu', 'ring0', 'ring1', 'ring2', 'ring3', 'ring4']

BAND_ORDER_ESTIMATES = {
    'U': 130,
    'B': 100,
    'V': 85,
    'R': 90,
    'IZ': 55,
    'Y': 118,
    'J': 99,
    'H': 76,
}

DEFAULT_SCALING = {
    'U': 2.1e5, 'B': 1e5, 'V': 1e5, 'R': 8e4,
    'IZ': 9e4, 'Y': 1.7e4, 'J': 1.4e4, 'H': 1.1e4,
}

DEFAULT_FINESSE = {
    'U': 23, 'R': 23, 'IZ': 23,
    'B': 26, 'V': 26, 'Y': 26, 'J': 26, 'H': 26,
}

TELESCOPE = {
    'primary_diameter': 39.3,
    'central_obstruction': 4.09,
}

INSTRUMENTS = {
    'Y': {
        'instrument_name': 'ANDES',
        'n_fibers': 75,
        'detector_size': (4096, 4096),
        'pixel_size': 15,
        'hdf_models': {
            'default': 'ANDES_75fibre_Y',
            'with_fiber_eff': 'ANDES_Y01_wFiberEff'
        },
        'zemax_files': {
            'Y': 'HIRES_Y_21jan2022_sconf_noap.zmx'
        },
        'diffraction_orders': list(range(109, 127)),
        'sampling': 2.1,
        'skip_fibers': [3, 4, 36, 37, 39, 40, 72, 73],
        'fiber_config_key': 'YJH_SL',
        'ifu_config_key': 'YJH_IFU',
    },
    'J': {
        'instrument_name': 'ANDES',
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
        'sampling': 2.1,
        'skip_fibers': [3, 4, 36, 37, 39, 40, 72, 73],
        'fiber_config_key': 'YJH_SL',
        'ifu_config_key': 'YJH_IFU',
    },
    'H': {
        'instrument_name': 'ANDES',
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
        'sampling': 2.1,
        'skip_fibers': [3, 4, 36, 37, 39, 40, 72, 73],
        'fiber_config_key': 'YJH_SL',
        'ifu_config_key': 'YJH_IFU',
    },
    'R': {
        'instrument_name': 'ANDES',
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
        'sampling': 4.0,
        'skip_fibers': [32, 35],
        'fiber_config_key': 'UBVRIZ',
    },
    'IZ': {
        'instrument_name': 'ANDES',
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
        'sampling': 4.0,
        'skip_fibers': [32, 35],
        'fiber_config_key': 'UBVRIZ',
    },
    'U': {
        'instrument_name': 'ANDES',
        'n_fibers': 66,
        'detector_size': (9216, 9232),
        'pixel_size': 15,
        'hdf_models': {
            'default': 'ANDES_U_v88'
        },
        'sampling': 4.0,
        'skip_fibers': [32, 35],
        'fiber_config_key': 'UBVRIZ',
    },
    'B': {
        'instrument_name': 'ANDES',
        'n_fibers': 66,
        'detector_size': (9216, 9232),
        'pixel_size': 15,
        'hdf_models': {
            'default': 'ANDES_B_v88'
        },
        'sampling': 4.0,
        'skip_fibers': [32, 35],
        'fiber_config_key': 'UBVRIZ',
    },
    'V': {
        'instrument_name': 'ANDES',
        'n_fibers': 66,
        'detector_size': (9216, 9232),
        'pixel_size': 15,
        'hdf_models': {
            'default': 'ANDES_V_v88'
        },
        'sampling': 4.0,
        'skip_fibers': [32, 35],
        'fiber_config_key': 'UBVRIZ',
    },
}

FIBER_CONFIG = {
    'UBVRIZ': {
        'n_fibers': 66,
        'slits': {
            'slitA': list(range(1, 32)),
            'slitB': list(range(36, 67)),
            'cal_fibers': [33, 34]
        }
    },
    'YJH_SL': {
        'fiber_size': 474,
        'n_fibers': 75,
        'slits': {
            'slitA': list(range(4, 35)),
            'slitB': list(range(42, 73)),
            'cal_fibers': [1, 37, 38, 39, 75]
        }
    },
    'YJH_IFU': {
        'fiber_size': 474,
        'n_fibers': 75,
        'slits': {
            'ring0': [3],
            'ring1': [6, 8, 10, 12, 14, 16],
            'ring2': list(range(18, 18 + 12)),
            'ring3': list(range(31, 31 + 18)),
            'ring4': list(range(50, 50 + 24)),
            'cal_fibers': [1, 75]
        }
    }
}
