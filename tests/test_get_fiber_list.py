"""Tests for SimulationConfig.get_fiber_list — the fiber-selection heart of every run.

These map every fiber mode to the concrete fiber list that flows into PyEchelle.
A regression here silently illuminates the wrong fibers, so the modes are pinned
to the exact lists encoded in core/andes.py and core/mosaic.py.
"""

import pytest

from andes_simulator.core.config import SimulationConfig, FiberConfig, SourceConfig


def _make(band, mode, fibers="all", skip_fibers=None, source_type="constant", filepath=None):
    """Build a minimal SimulationConfig for fiber-list testing."""
    return SimulationConfig(
        simulation_type="flat_field",
        band=band,
        source=SourceConfig(type=source_type, filepath=filepath),
        fibers=FiberConfig(mode=mode, fibers=fibers, skip_fibers=skip_fibers or []),
    )


# ---------- "all" mode ----------

def test_all_mode_R_band_returns_1_to_66():
    cfg = _make("R", "all")
    assert cfg.get_fiber_list() == list(range(1, 67))


def test_all_mode_Y_band_returns_1_to_75():
    cfg = _make("Y", "all")
    assert cfg.get_fiber_list() == list(range(1, 76))


def test_all_mode_LR_blue_returns_1_to_980():
    cfg = _make("B_LR", "all")
    fibers = cfg.get_fiber_list()
    assert fibers[0] == 1
    assert fibers[-1] == 980
    assert len(fibers) == 980


def test_all_mode_respects_skip_fibers_from_config():
    cfg = _make("R", "all", skip_fibers=[10, 20, 30])
    fibers = cfg.get_fiber_list()
    assert 10 not in fibers
    assert 20 not in fibers
    assert 30 not in fibers
    assert len(fibers) == 63


# ---------- "single" mode ----------

def test_single_mode_returns_one_fiber():
    cfg = _make("R", "single", fibers=[21])
    assert cfg.get_fiber_list() == [21]


def test_single_mode_with_int_fiber():
    cfg = _make("R", "single", fibers=21)
    assert cfg.get_fiber_list() == [21]


def test_single_mode_rejects_string_fiber():
    with pytest.raises(ValueError, match="Specific fiber"):
        _make("R", "single", fibers="all")


# ---------- "even" / "odd" ----------

def test_even_mode_R_band():
    cfg = _make("R", "even")
    fibers = cfg.get_fiber_list()
    assert fibers[0] == 2
    assert fibers[-1] == 66
    assert all(f % 2 == 0 for f in fibers)
    assert len(fibers) == 33


def test_odd_mode_R_band():
    cfg = _make("R", "odd")
    fibers = cfg.get_fiber_list()
    assert fibers[0] == 1
    assert fibers[-1] == 65
    assert all(f % 2 == 1 for f in fibers)
    assert len(fibers) == 33


def test_even_mode_Y_band():
    cfg = _make("Y", "even")
    fibers = cfg.get_fiber_list()
    assert all(f % 2 == 0 for f in fibers)
    assert len(fibers) == 37  # 75 fibers -> 37 even


# ---------- "slitA" / "slitB" / "cal_sl" ----------

def test_slitA_R_band_returns_fibers_1_through_31():
    cfg = _make("R", "slitA")
    assert cfg.get_fiber_list() == list(range(1, 32))


def test_slitB_R_band_returns_fibers_36_through_66():
    cfg = _make("R", "slitB")
    assert cfg.get_fiber_list() == list(range(36, 67))


def test_cal_sl_R_band_returns_calibration_fibers():
    cfg = _make("R", "cal_sl")
    assert cfg.get_fiber_list() == [33, 34]


def test_cal_sl_Y_band_returns_YJH_calibration_fibers():
    cfg = _make("Y", "cal_sl")
    assert cfg.get_fiber_list() == [1, 37, 38, 39, 75]


def test_slitA_Y_band_uses_YJH_layout():
    cfg = _make("Y", "slitA")
    assert cfg.get_fiber_list() == list(range(4, 35))


def test_slitB_Y_band_uses_YJH_layout():
    cfg = _make("Y", "slitB")
    assert cfg.get_fiber_list() == list(range(42, 73))


# ---------- IFU modes (YJH only) ----------

def test_ring0_Y_band():
    cfg = _make("Y", "ring0")
    assert cfg.get_fiber_list() == [3]


def test_ring1_Y_band():
    cfg = _make("Y", "ring1")
    assert cfg.get_fiber_list() == [6, 8, 10, 12, 14, 16]


def test_ring2_Y_band_has_12_fibers():
    cfg = _make("Y", "ring2")
    fibers = cfg.get_fiber_list()
    assert fibers == list(range(18, 30))


def test_ring3_Y_band_has_18_fibers():
    cfg = _make("Y", "ring3")
    fibers = cfg.get_fiber_list()
    assert fibers == list(range(31, 49))


def test_ring4_Y_band_has_24_fibers():
    cfg = _make("Y", "ring4")
    fibers = cfg.get_fiber_list()
    assert fibers == list(range(50, 74))


def test_ifu_Y_band_concatenates_all_rings():
    cfg = _make("Y", "ifu")
    fibers = cfg.get_fiber_list()
    # ring0 (1) + ring1 (6) + ring2 (12) + ring3 (18) + ring4 (24) = 61
    assert len(fibers) == 61
    assert fibers[0] == 3  # ring0 first


def test_cal_ifu_Y_band():
    cfg = _make("Y", "cal_ifu")
    assert cfg.get_fiber_list() == [1, 75]


def test_ifu_modes_unavailable_on_optical_bands():
    cfg = _make("R", "ring0")
    with pytest.raises(ValueError, match="not available"):
        cfg.get_fiber_list()


def test_cal_ifu_unavailable_on_optical_bands():
    cfg = _make("R", "cal_ifu")
    with pytest.raises(ValueError, match="not available"):
        cfg.get_fiber_list()


# ---------- "custom" mode ----------

def test_custom_mode_returns_explicit_list():
    cfg = _make("R", "custom", fibers=[5, 10, 15])
    assert cfg.get_fiber_list() == [5, 10, 15]


def test_custom_mode_requires_list():
    cfg = _make("R", "custom", fibers="all")
    with pytest.raises(ValueError, match="explicit fiber list"):
        cfg.get_fiber_list()


# ---------- skip_fibers applied after mode resolution ----------

def test_skip_fibers_filters_slit_selection():
    cfg = _make("R", "slitA", skip_fibers=[5, 6])
    fibers = cfg.get_fiber_list()
    assert 5 not in fibers and 6 not in fibers
    assert len(fibers) == 29  # 31 - 2


def test_skip_fibers_filters_ifu_selection():
    cfg = _make("Y", "ring1", skip_fibers=[10])
    assert cfg.get_fiber_list() == [6, 8, 12, 14, 16]
