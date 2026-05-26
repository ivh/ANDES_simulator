"""Tests for FlatFieldSource — fiber illumination patterns for flat-field calibration.

Each method returns a per-fiber list of PyEchelle source objects. We verify
which fibers are illuminated by comparing against the source's own
illuminated_source/dark_source instances, so a swap of the two (the classic
off-by-one regression here) would be caught immediately.
"""

import pytest

from andes_simulator.sources.flat_field import FlatFieldSource


def _illuminated_mask(sources, ff: FlatFieldSource):
    """Return list of bools: True where source is the illuminated one."""
    return [src is ff.illuminated_source for src in sources]


# ---------- single fiber ----------

def test_single_fiber_illuminates_only_that_fiber():
    ff = FlatFieldSource()
    sources = ff.get_single_fiber_sources(fiber_num=5, n_fibers=10)
    mask = _illuminated_mask(sources, ff)
    assert sum(mask) == 1
    assert mask[4] is True  # 1-based -> index 4


def test_single_fiber_first_fiber():
    ff = FlatFieldSource()
    sources = ff.get_single_fiber_sources(fiber_num=1, n_fibers=10)
    mask = _illuminated_mask(sources, ff)
    assert mask[0] is True
    assert sum(mask) == 1


def test_single_fiber_last_fiber():
    ff = FlatFieldSource()
    sources = ff.get_single_fiber_sources(fiber_num=10, n_fibers=10)
    mask = _illuminated_mask(sources, ff)
    assert mask[-1] is True
    assert sum(mask) == 1


def test_single_fiber_out_of_range_raises():
    ff = FlatFieldSource()
    with pytest.raises(ValueError, match="out of range"):
        ff.get_single_fiber_sources(fiber_num=11, n_fibers=10)


def test_single_fiber_zero_raises():
    ff = FlatFieldSource()
    with pytest.raises(ValueError, match="out of range"):
        ff.get_single_fiber_sources(fiber_num=0, n_fibers=10)


# ---------- all fibers ----------

def test_all_fibers_illuminates_every_fiber():
    ff = FlatFieldSource()
    sources = ff.get_all_fibers_sources(n_fibers=66)
    assert all(src is ff.illuminated_source for src in sources)
    assert len(sources) == 66


def test_all_fibers_with_skip_darkens_specified():
    ff = FlatFieldSource()
    sources = ff.get_all_fibers_sources(n_fibers=66, skip_fibers=[32, 35])
    mask = _illuminated_mask(sources, ff)
    # Indices 31 and 34 should be dark
    assert mask[31] is False
    assert mask[34] is False
    # All others illuminated
    assert sum(mask) == 64


def test_all_fibers_skip_out_of_range_is_ignored():
    ff = FlatFieldSource()
    # 999 is outside [1, 66] and should silently not affect the pattern
    sources = ff.get_all_fibers_sources(n_fibers=66, skip_fibers=[999])
    assert sum(_illuminated_mask(sources, ff)) == 66


# ---------- even/odd ----------

def test_even_illuminates_only_even_fibers():
    ff = FlatFieldSource()
    sources = ff.get_even_odd_sources(n_fibers=10, illuminate_even=True)
    mask = _illuminated_mask(sources, ff)
    # 1-based: fibers 2,4,6,8,10 illuminated => indices 1,3,5,7,9
    assert mask == [False, True, False, True, False, True, False, True, False, True]


def test_odd_illuminates_only_odd_fibers():
    ff = FlatFieldSource()
    sources = ff.get_even_odd_sources(n_fibers=10, illuminate_even=False)
    mask = _illuminated_mask(sources, ff)
    assert mask == [True, False, True, False, True, False, True, False, True, False]


# ---------- pseudo-slits (optical bands, 66 fibers) ----------

def test_first_slit_illuminates_first_31_fibers():
    ff = FlatFieldSource()
    sources = ff.get_first_slit_sources(n_fibers=66)
    mask = _illuminated_mask(sources, ff)
    assert sum(mask) == 31
    assert all(mask[i] for i in range(31))
    assert not any(mask[i] for i in range(31, 66))


def test_first_slit_rejects_non_optical_band_size():
    ff = FlatFieldSource()
    with pytest.raises(ValueError, match="66 fibers"):
        ff.get_first_slit_sources(n_fibers=75)


def test_second_slit_illuminates_fibers_35_to_66():
    ff = FlatFieldSource()
    sources = ff.get_second_slit_sources(n_fibers=66)
    mask = _illuminated_mask(sources, ff)
    # 1-based fibers 35..66 => indices 34..65 illuminated, 32 fibers
    assert sum(mask) == 32
    assert all(mask[i] for i in range(34, 66))
    assert not any(mask[i] for i in range(34))


def test_second_slit_rejects_non_optical_band_size():
    ff = FlatFieldSource()
    with pytest.raises(ValueError, match="66 fibers"):
        ff.get_second_slit_sources(n_fibers=75)


# ---------- calibration ----------

def test_calibration_illuminates_fibers_33_and_34():
    ff = FlatFieldSource()
    sources = ff.get_calibration_sources(n_fibers=66)
    mask = _illuminated_mask(sources, ff)
    assert mask[32] is True  # fiber 33
    assert mask[33] is True  # fiber 34
    assert sum(mask) == 2


def test_calibration_rejects_non_optical_band_size():
    ff = FlatFieldSource()
    with pytest.raises(ValueError, match="66 fibers"):
        ff.get_calibration_sources(n_fibers=75)


# ---------- custom ----------

def test_custom_illuminates_specified_fibers():
    ff = FlatFieldSource()
    sources = ff.get_custom_sources(illuminated_fibers=[2, 4, 6], n_fibers=10)
    mask = _illuminated_mask(sources, ff)
    assert mask == [False, True, False, True, False, True, False, False, False, False]


def test_custom_rejects_out_of_range_fiber():
    ff = FlatFieldSource()
    with pytest.raises(ValueError, match="out of range"):
        ff.get_custom_sources(illuminated_fibers=[1, 11], n_fibers=10)


# ---------- flux level wiring ----------

def test_flux_level_is_stored_on_illuminated_source():
    ff = FlatFieldSource(flux_level=0.005)
    # The illuminated source carries the flux value; the dark source carries 0.
    # We don't peek into PyEchelle internals but verify the instance is reused.
    sources = ff.get_all_fibers_sources(n_fibers=3)
    assert all(src is ff.illuminated_source for src in sources)
    assert ff.flux_level == 0.005
