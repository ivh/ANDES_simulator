"""Tests for LFCSource — laser frequency comb physics.

LFC lines are equidistant in velocity (log-uniform in wavelength). The
velocity spacing is derived from band order and a target lines-per-order;
mis-deriving it would shift every calibration line. We pin the formula,
the log-uniform property, and clip behavior.
"""

import numpy as np
import pytest

from andes_simulator.sources.lfc import LFCSource


# ---------- velocity spacing physics ----------

def test_velocity_spacing_formula_for_R_band():
    # delta_v/c = 1/(N * m); R-band m_center=90, N=100 default
    lfc = LFCSource(band="R", lines_per_order=100)
    c = 299792458.0
    expected = c / (100 * 90)
    assert lfc._calculate_velocity_spacing() == pytest.approx(expected)


def test_velocity_spacing_scales_inversely_with_lines_per_order():
    fast = LFCSource(band="R", lines_per_order=50)
    slow = LFCSource(band="R", lines_per_order=100)
    # halving lines_per_order doubles the spacing
    assert fast._calculate_velocity_spacing() == pytest.approx(
        2 * slow._calculate_velocity_spacing())


def test_velocity_spacing_around_33km_per_second_for_R():
    # Documented behavior: ~33 km/s for R-band
    lfc = LFCSource(band="R", lines_per_order=100)
    delta_v_km_s = lfc._calculate_velocity_spacing() / 1000
    assert 30 < delta_v_km_s < 36


# ---------- wavelength generation ----------

def test_generated_wavelengths_are_log_uniform():
    lfc = LFCSource(band="R", lines_per_order=100)
    wl, _ = lfc._generate_lfc_wavelengths()
    ratios = wl[1:] / wl[:-1]
    # All ratios should be equal -> std/mean is essentially zero
    assert np.std(ratios) / np.mean(ratios) < 1e-10


def test_generated_wavelengths_cover_band_range():
    from andes_simulator.core.instruments import get_band_wavelength_range
    lfc = LFCSource(band="R", lines_per_order=100)
    wl, _ = lfc._generate_lfc_wavelengths()
    band_min, band_max = get_band_wavelength_range("R")
    assert wl[0] == pytest.approx(band_min)
    assert wl[-1] <= band_max


def test_generated_flux_is_constant_per_line():
    lfc = LFCSource(band="R", flux_per_line=5e4)
    _, fluxes = lfc._generate_lfc_wavelengths()
    assert np.all(fluxes == 5e4)


def test_get_line_info_reports_count_and_spacing():
    lfc = LFCSource(band="R", lines_per_order=100, flux_per_line=1e5)
    info = lfc.get_line_info()
    assert info["band"] == "R"
    assert info["target_lines_per_order"] == 100
    assert info["flux_per_line"] == 1e5
    assert info["n_lines"] > 0
    assert info["velocity_spacing_m_s"] > 0


# ---------- clip behavior ----------

def test_wavelength_clip_restricts_output():
    lfc = LFCSource(band="R", lines_per_order=100)
    wl, _ = lfc._generate_lfc_wavelengths(wl_min=650.0, wl_max=700.0)
    assert wl[0] >= 650.0
    assert wl[-1] <= 700.0


def test_clip_returns_empty_array_when_range_outside_band():
    lfc = LFCSource(band="R", lines_per_order=100)
    # Note: clip is applied after full-band generation, so out-of-band returns []
    wl, _ = lfc._generate_lfc_wavelengths(wl_min=10.0, wl_max=20.0)
    assert wl.size == 0


def test_create_lfc_source_raises_on_empty_clip():
    lfc = LFCSource(band="R", lines_per_order=100)
    with pytest.raises(ValueError, match="no lines within"):
        lfc._create_lfc_source(wl_min=10.0, wl_max=20.0)


# ---------- bad band ----------

def test_unknown_band_raises():
    with pytest.raises(ValueError, match="Unknown band"):
        LFCSource(band="NOT_A_BAND")


# ---------- single fiber selection ----------

def test_get_single_fiber_sources_correct_length():
    lfc = LFCSource(band="R", flux_per_line=1e5)
    sources = lfc.get_single_fiber_sources(fiber_num=5, n_fibers=10)
    assert len(sources) == 10


def test_get_single_fiber_rejects_out_of_range():
    lfc = LFCSource(band="R", flux_per_line=1e5)
    with pytest.raises(ValueError, match="out of range"):
        lfc.get_single_fiber_sources(fiber_num=11, n_fibers=10)


def test_get_all_fibers_skipped_have_dark_source():
    from pyechelle.sources import ConstantPhotonFlux
    lfc = LFCSource(band="R", flux_per_line=1e5)
    sources = lfc.get_all_fibers_sources(n_fibers=10, skip_fibers=[3, 7])
    assert isinstance(sources[2], ConstantPhotonFlux)  # fiber 3
    assert isinstance(sources[6], ConstantPhotonFlux)  # fiber 7
