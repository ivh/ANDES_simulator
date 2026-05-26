"""Tests for FabryPerotSource — Airy transmission physics.

The gap thickness auto-compute and the peak wavelength formula
lambda_k = 2*n*d / k are the load-bearing equations. Regressions here
would shift every FP calibration line and break wavelength solution.
"""

import numpy as np
import pytest

from andes_simulator.sources.fabry_perot import FabryPerotSource
from andes_simulator.core.instruments import (
    get_band_wavelength_range, BAND_ORDER_ESTIMATES,
)


# ---------- gap auto-compute ----------

def test_auto_gap_matches_formula_R_band():
    # d = N * m * lambda_c / (2*n)   with d in mm, lambda in nm * 1e-6
    fp = FabryPerotSource(band="R", finesse=23, lines_per_order=100,
                          scaling_factor=1e6, n_refr=1.0)
    wl_min, wl_max = get_band_wavelength_range("R")
    wl_center = (wl_min + wl_max) / 2
    m = BAND_ORDER_ESTIMATES["R"]
    expected = (100 * m * wl_center * 1e-6) / 2.0
    assert fp.gap_mm == pytest.approx(expected)


def test_explicit_gap_overrides_auto():
    fp = FabryPerotSource(band="R", finesse=23, gap_mm=0.5, scaling_factor=1e6)
    assert fp.gap_mm == 0.5


def test_gap_scales_with_lines_per_order():
    fp_low = FabryPerotSource(band="R", finesse=23, lines_per_order=50, scaling_factor=1e6)
    fp_high = FabryPerotSource(band="R", finesse=23, lines_per_order=100, scaling_factor=1e6)
    # More lines per order -> larger gap
    assert fp_high.gap_mm == pytest.approx(2 * fp_low.gap_mm)


# ---------- peak wavelengths ----------

def test_peak_wavelengths_satisfy_2nd_over_k_formula():
    fp = FabryPerotSource(band="R", finesse=23, gap_mm=3.0,
                          scaling_factor=1e6, n_refr=1.0)
    peaks = fp.get_peak_wavelengths()
    gap_nm = 3.0 * 1e6
    # Reconstruct k from each peak: k = 2*n*d / lambda
    k_values = np.round(2.0 * gap_nm / peaks).astype(int)
    # Each peak should equal 2*n*d / k exactly (within float)
    reconstructed = 2.0 * gap_nm / k_values
    np.testing.assert_allclose(reconstructed, peaks, rtol=1e-12)


def test_peak_wavelengths_within_band():
    fp = FabryPerotSource(band="R", finesse=23, scaling_factor=1e6)
    wl_min, wl_max = get_band_wavelength_range("R")
    peaks = fp.get_peak_wavelengths()
    assert peaks.min() >= wl_min
    assert peaks.max() <= wl_max


def test_peak_wavelengths_clip_respects_user_range():
    fp = FabryPerotSource(band="R", finesse=23, scaling_factor=1e6)
    peaks = fp.get_peak_wavelengths(wl_min=650.0, wl_max=700.0)
    assert peaks.min() >= 650.0
    assert peaks.max() <= 700.0


def test_peak_count_scales_with_lines_per_order():
    fp_low = FabryPerotSource(band="R", finesse=23, lines_per_order=50, scaling_factor=1e6)
    fp_high = FabryPerotSource(band="R", finesse=23, lines_per_order=100, scaling_factor=1e6)
    # ~2x more lines per order should give ~2x more peaks in the band
    n_low = len(fp_low.get_peak_wavelengths())
    n_high = len(fp_high.get_peak_wavelengths())
    assert 1.8 < n_high / n_low < 2.2


# ---------- spectrum info ----------

def test_get_spectrum_info_basic_fields():
    fp = FabryPerotSource(band="R", finesse=23, scaling_factor=1e6)
    info = fp.get_spectrum_info()
    assert info["band"] == "R"
    assert info["finesse"] == 23
    assert info["n_refr"] == 1.0
    assert info["gap_mm"] == fp.gap_mm
    assert info["n_samples"] > 0


def test_fsr_and_fwhm_are_consistent_with_finesse():
    # finesse F = FSR / FWHM, so FWHM = FSR / F (by construction)
    fp = FabryPerotSource(band="R", finesse=23, scaling_factor=1e6)
    info = fp.get_spectrum_info()
    assert info["fwhm_nm_at_center"] == pytest.approx(info["fsr_nm_at_center"] / 23)


# ---------- default finesse from band ----------

def test_default_finesse_used_when_not_specified():
    from andes_simulator.core.instruments import DEFAULT_FINESSE
    fp = FabryPerotSource(band="R", scaling_factor=1e6)
    assert fp.finesse == DEFAULT_FINESSE["R"]


# ---------- error paths ----------

def test_unknown_band_raises():
    with pytest.raises((KeyError, ValueError)):
        FabryPerotSource(band="NOT_A_BAND", scaling_factor=1e6)


def test_unknown_band_with_explicit_finesse_raises_value_error():
    # Passing finesse skips DEFAULT_FINESSE lookup, so the band check fires
    with pytest.raises(ValueError, match="Unknown band"):
        FabryPerotSource(band="NOT_A_BAND", finesse=23, scaling_factor=1e6)


# ---------- single fiber selection ----------

def test_get_single_fiber_sources_length_matches_n_fibers():
    fp = FabryPerotSource(band="R", finesse=23, scaling_factor=1e6)
    sources = fp.get_single_fiber_sources(fiber_num=5, n_fibers=10)
    assert len(sources) == 10


def test_get_single_fiber_out_of_range_raises():
    fp = FabryPerotSource(band="R", finesse=23, scaling_factor=1e6)
    with pytest.raises(ValueError, match="out of range"):
        fp.get_single_fiber_sources(fiber_num=11, n_fibers=10)


def test_get_all_fibers_respects_skip_fibers():
    from pyechelle.sources import ConstantPhotonFlux
    fp = FabryPerotSource(band="R", finesse=23, scaling_factor=1e6)
    sources = fp.get_all_fibers_sources(n_fibers=10, skip_fibers=[4])
    assert isinstance(sources[3], ConstantPhotonFlux)  # fiber 4 dark
