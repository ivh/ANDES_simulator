"""Tests for AndesSimulator helper methods that don't require PyEchelle.

_parse_fib_eff handles a user-facing string format and is easy to break.
_update_fits_headers writes simulation metadata that downstream tools
(combine, analysis scripts) rely on for reproducibility.
"""

import numpy as np
import pytest
from astropy.io import fits

from andes_simulator.core.config import (
    SimulationConfig, SourceConfig, FiberConfig
)
from andes_simulator.core.simulator import AndesSimulator


def _make_sim(**overrides) -> AndesSimulator:
    base = dict(
        simulation_type='flat_field',
        band='R',
        source=SourceConfig(type='constant'),
        fibers=FiberConfig(mode='all'),
    )
    base.update(overrides)
    return AndesSimulator(SimulationConfig(**base))


# ---------- _parse_fib_eff ----------

def test_parse_fib_eff_range():
    sim = _make_sim()
    assert sim._parse_fib_eff('0.7-0.9') == (0.7, 0.9)


def test_parse_fib_eff_scalar_returns_equal_min_max():
    sim = _make_sim()
    assert sim._parse_fib_eff('0.85') == (0.85, 0.85)


def test_parse_fib_eff_one():
    sim = _make_sim()
    assert sim._parse_fib_eff('1.0') == (1.0, 1.0)


def test_parse_fib_eff_invalid_format_three_parts_rejected():
    sim = _make_sim()
    with pytest.raises(ValueError, match="Invalid fib_eff"):
        sim._parse_fib_eff('0.7-0.8-0.9')


def test_parse_fib_eff_non_numeric_rejected():
    sim = _make_sim()
    with pytest.raises(ValueError):
        sim._parse_fib_eff('high')


# ---------- _update_fits_headers ----------

@pytest.fixture
def fits_file(tmp_path):
    """Write a small empty FITS file ready for header updates."""
    path = tmp_path / "out.fits"
    fits.PrimaryHDU(np.zeros((4, 4), dtype=np.uint16)).writeto(path)
    return path


def test_update_fits_headers_writes_core_metadata(fits_file, tmp_path):
    sim = _make_sim(exposure_time=42.0)
    sim._update_fits_headers(fits_file, hdf_path=tmp_path / "ANDES_R.hdf")
    with fits.open(fits_file) as hdul:
        hdr = hdul[0].header
        assert hdr['INSTRUME'] == 'ANDES'
        assert hdr['BAND'] == 'R'
        assert hdr['SIMTYPE'] == 'flat_field'
        assert hdr['EXPTIME'] == 42.0
        assert hdr['HDFMODEL'] == 'ANDES_R.hdf'
        assert hdr['DATE-OBS']  # set


def test_update_fits_headers_source_fields(fits_file):
    sim = _make_sim(source=SourceConfig(type='constant', flux=0.5,
                                        scaling_factor=1e5))
    sim._update_fits_headers(fits_file)
    with fits.open(fits_file) as hdul:
        hdr = hdul[0].header
        assert hdr['SRCTYPE'] == 'constant'
        assert hdr['SRCFLUX'] == 0.5
        assert hdr['SRCSCALE'] == 1e5


def test_update_fits_headers_fiber_mode(fits_file):
    sim = _make_sim(fibers=FiberConfig(mode='slitA'))
    sim._update_fits_headers(fits_file)
    with fits.open(fits_file) as hdul:
        assert hdul[0].header['FIBMODE'] == 'slitA'


def test_update_fits_headers_velocity_shift_scalar(fits_file):
    sim = _make_sim(velocity_shift=250.0)
    sim._update_fits_headers(fits_file)
    with fits.open(fits_file) as hdul:
        assert hdul[0].header['VSHIFT'] == 250.0


def test_update_fits_headers_velocity_shift_per_fiber_list(fits_file):
    sim = _make_sim(velocity_shift=[100.0, 200.0, 300.0])
    sim._update_fits_headers(fits_file)
    with fits.open(fits_file) as hdul:
        # Per-fiber lists become the literal string 'per-fiber'
        assert hdul[0].header['VSHIFT'] == 'per-fiber'


def test_update_fits_headers_no_velocity_shift_omits_key(fits_file):
    sim = _make_sim()
    sim._update_fits_headers(fits_file)
    with fits.open(fits_file) as hdul:
        assert 'VSHIFT' not in hdul[0].header


def test_update_fits_headers_x_shift(fits_file):
    sim = _make_sim(x_shift=0.5)
    sim._update_fits_headers(fits_file)
    with fits.open(fits_file) as hdul:
        assert hdul[0].header['XSHIFT'] == 0.5


def test_update_fits_headers_wavelength_range(fits_file):
    sim = _make_sim(wl_min=650.0, wl_max=700.0)
    sim._update_fits_headers(fits_file)
    with fits.open(fits_file) as hdul:
        hdr = hdul[0].header
        assert hdr['WL_MIN'] == 650.0
        assert hdr['WL_MAX'] == 700.0


def test_update_fits_headers_fib_eff(fits_file):
    sim = _make_sim(fib_eff='0.85')
    sim._update_fits_headers(fits_file)
    with fits.open(fits_file) as hdul:
        assert hdul[0].header['FIBEFF'] == '0.85'


def test_update_fits_headers_e2e_marker(fits_file):
    sim = _make_sim()
    sim._update_fits_headers(fits_file)
    with fits.open(fits_file) as hdul:
        # HIERARCH key for the E2E_SIM marker
        assert hdul[0].header['E2E_SIM'] is True
