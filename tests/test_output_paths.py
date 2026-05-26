"""Tests for SimulationConfig.get_output_path — filename naming conventions.

Downstream post-processing tools (combine, psf-process) parse filenames with
patterns like "R_FP_fiber{fib:02d}_*.fits". A regression in the naming would
silently break those wildcards, so each simulation type's filename shape is
pinned here.
"""

from andes_simulator.core.config import (
    SimulationConfig, FiberConfig, SourceConfig, OutputConfig
)


def _make(**kwargs):
    base = dict(
        band="R",
        source=SourceConfig(type="constant"),
        fibers=FiberConfig(mode="all"),
        output=OutputConfig(directory="."),
        exposure_time=30.0,
    )
    base.update(kwargs)
    return SimulationConfig(**base)


# ---------- flat field ----------

def test_flat_field_single_fiber_filename():
    cfg = _make(simulation_type="flat_field",
                fibers=FiberConfig(mode="single", fibers=[21]))
    assert cfg.get_output_path(fiber_num=21).name == "R_FF_fiber21_30s.fits"


def test_flat_field_single_fiber_zero_padded():
    cfg = _make(simulation_type="flat_field",
                fibers=FiberConfig(mode="single", fibers=[7]))
    assert cfg.get_output_path(fiber_num=7).name == "R_FF_fiber07_30s.fits"


def test_flat_field_all_fibers_uses_mode_in_name():
    cfg = _make(simulation_type="flat_field", fibers=FiberConfig(mode="all"))
    assert cfg.get_output_path().name == "R_FF_all_30s.fits"


def test_flat_field_even_mode_name():
    cfg = _make(simulation_type="flat_field", fibers=FiberConfig(mode="even"))
    assert cfg.get_output_path().name == "R_FF_even_30s.fits"


# ---------- Fabry-Perot ----------

def test_fp_single_fiber_filename():
    cfg = _make(simulation_type="fabry_perot",
                source=SourceConfig(type="fabry_perot"),
                fibers=FiberConfig(mode="single", fibers=[21]))
    assert cfg.get_output_path(fiber_num=21).name == "R_FP_fiber21.fits"


def test_fp_all_fibers_filename():
    cfg = _make(simulation_type="fabry_perot",
                source=SourceConfig(type="fabry_perot"),
                fibers=FiberConfig(mode="all"))
    assert cfg.get_output_path().name == "R_FP_all_30s.fits"


def test_fp_fiber_with_shift_suffix():
    cfg = _make(simulation_type="fabry_perot",
                source=SourceConfig(type="fabry_perot"),
                fibers=FiberConfig(mode="single", fibers=[21]))
    assert cfg.get_output_path(fiber_num=21, suffix="0p5").name == "R_FP_fiber21_shift0p5.fits"


# ---------- LFC ----------

def test_lfc_single_fiber_filename():
    cfg = _make(simulation_type="lfc",
                source=SourceConfig(type="lfc"),
                fibers=FiberConfig(mode="single", fibers=[21]),
                band="Y")
    assert cfg.get_output_path(fiber_num=21).name == "Y_LFC_fiber21.fits"


def test_lfc_all_fibers_filename():
    cfg = _make(simulation_type="lfc",
                source=SourceConfig(type="lfc"),
                fibers=FiberConfig(mode="all"),
                band="Y")
    assert cfg.get_output_path().name == "Y_LFC_all_30s.fits"


# ---------- spectrum ----------

def test_spectrum_single_fiber_filename():
    cfg = _make(simulation_type="spectrum",
                source=SourceConfig(type="csv", filepath="star.csv"),
                fibers=FiberConfig(mode="single", fibers=[33]),
                band="J")
    assert cfg.get_output_path(fiber_num=33).name == "J_spectrum_fiber33.fits"


def test_spectrum_all_fibers_filename():
    cfg = _make(simulation_type="spectrum",
                source=SourceConfig(type="csv", filepath="star.csv"),
                fibers=FiberConfig(mode="all"),
                band="J")
    assert cfg.get_output_path().name == "J_spectrum_30s.fits"


# ---------- wavelength suffix ----------

def test_wl_suffix_integer_bounds():
    cfg = _make(simulation_type="flat_field",
                fibers=FiberConfig(mode="single", fibers=[21]),
                wl_min=625.0, wl_max=720.0)
    assert cfg.get_output_path(fiber_num=21).name == "R_FF_fiber21_30s_wl625-720.fits"


def test_wl_suffix_decimal_bounds_use_p_separator():
    cfg = _make(simulation_type="flat_field",
                fibers=FiberConfig(mode="single", fibers=[21]),
                wl_min=625.0, wl_max=720.5)
    assert cfg.get_output_path(fiber_num=21).name == "R_FF_fiber21_30s_wl625-720p5.fits"


def test_wl_suffix_only_min_uses_x_placeholder():
    cfg = _make(simulation_type="flat_field",
                fibers=FiberConfig(mode="single", fibers=[21]),
                wl_min=625.0)
    assert cfg.get_output_path(fiber_num=21).name == "R_FF_fiber21_30s_wl625-x.fits"


def test_wl_suffix_only_max_uses_x_placeholder():
    cfg = _make(simulation_type="flat_field",
                fibers=FiberConfig(mode="single", fibers=[21]),
                wl_max=720.0)
    assert cfg.get_output_path(fiber_num=21).name == "R_FF_fiber21_30s_wlx-720.fits"


def test_wl_suffix_strips_trailing_zeros():
    cfg = _make(simulation_type="flat_field",
                fibers=FiberConfig(mode="single", fibers=[21]),
                wl_min=625.100, wl_max=720.000)
    # 625.100 -> 625p1, 720.000 -> 720
    assert "wl625p1-720" in cfg.get_output_path(fiber_num=21).name


# ---------- explicit filename override ----------

def test_explicit_filename_overrides_template():
    cfg = _make(simulation_type="flat_field",
                fibers=FiberConfig(mode="single", fibers=[21]),
                output=OutputConfig(directory=".", filename="custom_name.fits"))
    assert cfg.get_output_path(fiber_num=21).name == "custom_name.fits"
