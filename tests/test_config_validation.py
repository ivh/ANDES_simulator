"""Tests for SimulationConfig.validate — guard against silently bad configs.

The dataclass validates in __post_init__; these tests assert that mistyped
bands, modes, source types, and out-of-range wavelengths fail loudly at
construction rather than mid-simulation.
"""

import pytest

from andes_simulator.core.config import (
    SimulationConfig, FiberConfig, SourceConfig
)


def _make(**kwargs):
    """Build a baseline-valid SimulationConfig with optional overrides."""
    base = dict(
        simulation_type="flat_field",
        band="R",
        source=SourceConfig(type="constant"),
        fibers=FiberConfig(mode="all"),
    )
    base.update(kwargs)
    return SimulationConfig(**base)


# ---------- band ----------

def test_unknown_band_rejected():
    with pytest.raises(ValueError, match="Invalid band"):
        _make(band="ZZZ")


def test_all_known_bands_accepted():
    for band in ["U", "B", "V", "R", "IZ", "Y", "J", "H",
                 "LR-blue", "LR-red", "LR-J", "LR-H",
                 "HR-B1", "HR-R1", "HR-B2", "HR-R2", "HR-H"]:
        _make(band=band)  # should not raise


# ---------- simulation_type ----------

def test_unknown_simulation_type_rejected():
    with pytest.raises(ValueError, match="Invalid simulation_type"):
        _make(simulation_type="bogus")


@pytest.mark.parametrize("sim_type", ["flat_field", "fabry_perot", "spectrum", "lfc", "hdf_generation"])
def test_known_simulation_types_accepted(sim_type):
    # spectrum needs csv filepath
    if sim_type == "spectrum":
        _make(simulation_type=sim_type, source=SourceConfig(type="csv", filepath="x.csv"))
    else:
        _make(simulation_type=sim_type)


# ---------- fiber mode ----------

def test_unknown_fiber_mode_rejected():
    with pytest.raises(ValueError, match="Invalid fiber mode"):
        _make(fibers=FiberConfig(mode="bogus"))


# ---------- source type ----------

def test_unknown_source_type_rejected():
    with pytest.raises(ValueError, match="Invalid source type"):
        _make(source=SourceConfig(type="bogus"))


def test_csv_spectrum_requires_filepath():
    with pytest.raises(ValueError, match="CSV filepath required"):
        _make(simulation_type="spectrum", source=SourceConfig(type="csv", filepath=None))


def test_single_mode_requires_explicit_fiber():
    with pytest.raises(ValueError, match="Specific fiber"):
        _make(fibers=FiberConfig(mode="single", fibers="all"))


# ---------- wavelength range ----------

def test_wl_min_below_band_range_rejected():
    # R-band lower edge is ~621 nm
    with pytest.raises(ValueError, match="below R-band range"):
        _make(wl_min=100.0)


def test_wl_max_above_band_range_rejected():
    # R-band upper edge is ~764 nm
    with pytest.raises(ValueError, match="above R-band range"):
        _make(wl_max=2000.0)


def test_wl_min_greater_than_wl_max_rejected():
    with pytest.raises(ValueError, match="must be less than"):
        _make(wl_min=700.0, wl_max=650.0)


def test_wl_range_within_band_accepted():
    cfg = _make(wl_min=650.0, wl_max=700.0)
    assert cfg.wl_min == 650.0
    assert cfg.wl_max == 700.0


def test_wl_min_only_accepted():
    cfg = _make(wl_min=650.0)
    assert cfg.wl_max is None


# ---------- derived properties ----------

def test_n_fibers_property_R_band():
    cfg = _make(band="R")
    assert cfg.n_fibers == 66


def test_n_fibers_property_Y_band():
    cfg = _make(band="Y")
    assert cfg.n_fibers == 75


def test_detector_size_property_R_band():
    cfg = _make(band="R")
    assert cfg.detector_size == (9216, 9232)
