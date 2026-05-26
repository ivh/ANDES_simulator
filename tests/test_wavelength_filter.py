"""Tests for WavelengthFilteredSource — masking input wavelengths to a band.

This wrapper sits between PyEchelle's raytracer and the source spectrum.
The unit auto-detection (samples < 10 -> microns) is fragile in principle
but is what PyEchelle's bare-float calls require. We pin both code paths.
"""

import numpy as np
import pytest

from andes_simulator.sources.wavelength_filter import WavelengthFilteredSource


class _StubSource:
    """Minimal stub mimicking PyEchelle source interface."""
    list_like = False
    name = "stub"

    def __init__(self, value=10.0):
        self.value = value

    def get_counts(self, wl, integration_time):
        return np.full(np.asarray(wl).shape, self.value, dtype=float)


# ---------- masking with bounds ----------

def test_mask_zeros_below_wl_min():
    w = WavelengthFilteredSource(_StubSource(10.0), wl_min=400.0)
    wl = np.array([350.0, 400.0, 450.0])
    out = w.get_counts(wl, 1)
    assert out[0] == 0  # below min
    assert out[1] == 10  # at min (inclusive)
    assert out[2] == 10


def test_mask_zeros_above_wl_max():
    w = WavelengthFilteredSource(_StubSource(10.0), wl_max=500.0)
    wl = np.array([450.0, 500.0, 550.0])
    out = w.get_counts(wl, 1)
    assert out[0] == 10
    assert out[1] == 10  # at max (inclusive)
    assert out[2] == 0   # above max


def test_mask_with_both_bounds():
    w = WavelengthFilteredSource(_StubSource(10.0), wl_min=400.0, wl_max=500.0)
    wl = np.array([350.0, 400.0, 450.0, 500.0, 550.0])
    out = w.get_counts(wl, 1)
    np.testing.assert_array_equal(out, [0, 10, 10, 10, 0])


def test_no_bounds_passes_through_unchanged():
    w = WavelengthFilteredSource(_StubSource(10.0))
    wl = np.array([100.0, 1000.0])
    out = w.get_counts(wl, 1)
    np.testing.assert_array_equal(out, [10, 10])


# ---------- unit auto-detection ----------

def test_micron_input_auto_detected_and_converted_to_nm_for_comparison():
    # Bounds are in nm; raytracing passes microns (sample < 10).
    w = WavelengthFilteredSource(_StubSource(10.0), wl_min=400.0, wl_max=500.0)
    wl_um = np.array([0.35, 0.40, 0.45, 0.50, 0.55])
    out = w.get_counts(wl_um, 1)
    np.testing.assert_array_equal(out, [0, 10, 10, 10, 0])


def test_nm_input_when_sample_above_10():
    # First sample is large (>10) -> treated as nm already
    w = WavelengthFilteredSource(_StubSource(10.0), wl_min=400.0, wl_max=500.0)
    wl_nm = np.array([350.0, 400.0, 500.0, 550.0])
    out = w.get_counts(wl_nm, 1)
    np.testing.assert_array_equal(out, [0, 10, 10, 0])


# ---------- attribute propagation ----------

def test_list_like_attribute_propagated():
    src = _StubSource()
    src.list_like = True
    w = WavelengthFilteredSource(src)
    assert w.list_like is True


def test_name_attribute_propagated():
    src = _StubSource()
    src.name = "my_source"
    w = WavelengthFilteredSource(src)
    assert w.name == "my_source"


def test_default_list_like_when_missing():
    class Bare:
        def get_counts(self, wl, t):
            return np.zeros_like(wl)
    w = WavelengthFilteredSource(Bare())
    assert w.list_like is False
    assert w.name == "filtered_source"


# ---------- get_counts_rv (radial velocity variant) ----------

def test_get_counts_rv_delegates_when_supported():
    class WithRV(_StubSource):
        def get_counts_rv(self, wl, integration_time, rv=0):
            return np.full(np.asarray(wl).shape, 99.0)
    w = WavelengthFilteredSource(WithRV(), wl_min=400.0)
    out = w.get_counts_rv(np.array([350.0, 450.0]), 1, rv=100)
    assert out[0] == 0    # masked
    assert out[1] == 99   # from RV-aware source


def test_get_counts_rv_falls_back_to_get_counts():
    w = WavelengthFilteredSource(_StubSource(10.0), wl_min=400.0)
    out = w.get_counts_rv(np.array([350.0, 450.0]), 1, rv=0)
    assert out[0] == 0
    assert out[1] == 10
