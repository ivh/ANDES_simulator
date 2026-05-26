"""Tests for the instrument registry (core/instruments.py).

Covers band lookup, wavelength validation, and the band-inference helpers
that the CLI uses to skip --band when --hdf or --wl-min/max are given.
"""

import pytest

from andes_simulator.core.instruments import (
    get_instrument_config,
    get_instrument_name,
    get_all_bands,
    get_nir_bands,
    get_optical_bands,
    get_band_wavelength_range,
    infer_band_from_wavelengths,
    validate_wavelength_range,
)


# ---------- get_instrument_config ----------

def test_R_band_config_has_expected_fields():
    cfg = get_instrument_config("R")
    assert cfg["n_fibers"] == 66
    assert cfg["detector_size"] == (9216, 9232)
    assert cfg["pixel_size"] == 15
    assert "telescope" in cfg
    assert "fiber_config" in cfg


def test_Y_band_config_has_ifu():
    cfg = get_instrument_config("Y")
    assert cfg["n_fibers"] == 75
    assert "fiber_config" in cfg
    assert "ifu_config" in cfg


def test_R_band_config_has_no_ifu():
    cfg = get_instrument_config("R")
    assert "ifu_config" not in cfg


def test_mosaic_LR_blue_config():
    cfg = get_instrument_config("LR-blue")
    assert cfg["n_fibers"] == 980
    assert cfg["bundle_size"] == 7
    assert cfg["n_bundles"] == 140


def test_mosaic_HR_B1_has_19_fibers_per_bundle():
    cfg = get_instrument_config("HR-B1")
    assert cfg["n_fibers"] == 1140
    assert cfg["bundle_size"] == 19
    assert cfg["n_bundles"] == 60


def test_unknown_band_raises_value_error():
    with pytest.raises(ValueError, match="Unknown band"):
        get_instrument_config("NOT_A_BAND")


def test_instrument_config_is_copy_not_reference():
    cfg1 = get_instrument_config("R")
    cfg1["n_fibers"] = 999  # mutate the copy
    cfg2 = get_instrument_config("R")
    assert cfg2["n_fibers"] == 66  # original untouched


# ---------- get_instrument_name ----------

def test_andes_bands_resolve_to_ANDES():
    for band in ["U", "B", "V", "R", "IZ", "Y", "J", "H"]:
        assert get_instrument_name(band) == "ANDES"


def test_mosaic_bands_resolve_to_MOSAIC():
    for band in ["LR-blue", "LR-red", "HR-B1", "HR-H"]:
        assert get_instrument_name(band) == "MOSAIC"


def test_get_instrument_name_unknown_band():
    with pytest.raises(ValueError, match="Unknown band"):
        get_instrument_name("NOT_A_BAND")


# ---------- band groupings ----------

def test_all_bands_includes_andes_and_mosaic():
    bands = get_all_bands()
    assert "R" in bands
    assert "LR-blue" in bands


def test_nir_bands_are_YJH():
    assert get_nir_bands() == ["Y", "J", "H"]


def test_optical_bands_includes_R():
    assert "R" in get_optical_bands()


# ---------- get_band_wavelength_range ----------

def test_R_band_wavelength_range_reasonable():
    wl_min, wl_max = get_band_wavelength_range("R")
    # R-band roughly 620-770 nm
    assert 600 < wl_min < 650
    assert 700 < wl_max < 800


def test_mosaic_LR_blue_uses_config_range():
    # LR-blue has wavelength_range in config (no HDF lookup needed)
    wl_min, wl_max = get_band_wavelength_range("LR-blue")
    assert wl_min == 390
    assert wl_max == 625


def test_wavelength_range_is_cached():
    # Second call uses internal cache; verify it returns same value
    r1 = get_band_wavelength_range("R")
    r2 = get_band_wavelength_range("R")
    assert r1 == r2


# ---------- validate_wavelength_range ----------

def test_validate_accepts_in_band_range():
    # Should not raise
    validate_wavelength_range("R", wl_min=650.0, wl_max=700.0)


def test_validate_accepts_none():
    validate_wavelength_range("R", wl_min=None, wl_max=None)


def test_validate_rejects_wl_min_below_band():
    with pytest.raises(ValueError, match="below R-band"):
        validate_wavelength_range("R", wl_min=100.0)


def test_validate_rejects_wl_max_above_band():
    with pytest.raises(ValueError, match="above R-band"):
        validate_wavelength_range("R", wl_max=2000.0)


def test_validate_rejects_inverted_range():
    with pytest.raises(ValueError, match="must be less than"):
        validate_wavelength_range("R", wl_min=700.0, wl_max=650.0)


def test_validate_rejects_unknown_band():
    with pytest.raises(ValueError, match="Unknown band"):
        validate_wavelength_range("NOT_A_BAND", wl_min=500.0)


# ---------- infer_band_from_wavelengths ----------

def test_infer_unique_band_within_restriction():
    # Restricting to ANDES NIR removes MOSAIC LR-J ambiguity
    assert infer_band_from_wavelengths(1000, 1100, restrict_to=["Y", "J", "H"]) == "Y"


def test_infer_ambiguous_range_raises_with_candidates():
    # ~1000nm overlaps Y and LR-J across the two instruments
    with pytest.raises(ValueError, match="ambiguous"):
        infer_band_from_wavelengths(1000, 1050)


def test_infer_no_match_raises():
    with pytest.raises(ValueError, match="does not match"):
        infer_band_from_wavelengths(50, 60)


def test_infer_requires_at_least_one_bound():
    with pytest.raises(ValueError, match="At least one"):
        infer_band_from_wavelengths(wl_min=None, wl_max=None)


def test_infer_only_wl_min_supplied():
    # When only wl_min given, requested range is [wl_min, band_max]
    band = infer_band_from_wavelengths(wl_min=650.0, restrict_to=["R"])
    assert band == "R"
