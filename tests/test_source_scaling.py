"""Tests for resolve_source_scaling — flux/scaling resolution from CLI options.

This function controls the final ph/s reaching PyEchelle. Each source type has
distinct semantics (LFC and FP divide their scaling to compensate for
line/peak concentration), and CSV uses a separate flag to opt into file-header
scaling. Mis-routing the branches changes the brightness of every simulation.
"""

import pytest

from andes_simulator.cli.utils import resolve_source_scaling
from andes_simulator.core.instruments import DEFAULT_SCALING


# ---------- CSV branch ----------

def test_csv_no_user_scaling_returns_flux_and_enables_file_scaling():
    value, use_file = resolve_source_scaling("csv", "R", flux=2.0, user_scaling=None)
    assert value == 2.0
    assert use_file is True


def test_csv_with_user_scaling_multiplies_and_disables_file_scaling():
    value, use_file = resolve_source_scaling("csv", "R", flux=2.0, user_scaling=10.0)
    assert value == 20.0
    assert use_file is False


def test_csv_flux_equal_one_passes_through_when_no_scaling():
    value, use_file = resolve_source_scaling("csv", "Y", flux=1.0, user_scaling=None)
    assert value == 1.0
    assert use_file is True


# ---------- constant (flat-field) branch ----------

def test_constant_uses_band_default_when_no_user_scaling():
    value, flag = resolve_source_scaling("constant", "R", flux=1.0, user_scaling=None)
    assert value == DEFAULT_SCALING["R"]
    assert flag is None  # only csv exposes the file-scaling flag


def test_constant_flux_multiplies_band_default():
    value, _ = resolve_source_scaling("constant", "R", flux=2.0, user_scaling=None)
    assert value == 2.0 * DEFAULT_SCALING["R"]


def test_constant_user_scaling_replaces_band_default():
    value, _ = resolve_source_scaling("constant", "R", flux=1.0, user_scaling=5e5)
    assert value == 5e5


def test_constant_user_scaling_and_flux_multiply():
    value, _ = resolve_source_scaling("constant", "R", flux=3.0, user_scaling=1e4)
    assert value == 3.0 * 1e4


def test_constant_band_default_varies_with_band():
    v_r, _ = resolve_source_scaling("constant", "R", flux=1.0, user_scaling=None)
    v_h, _ = resolve_source_scaling("constant", "H", flux=1.0, user_scaling=None)
    assert v_r == DEFAULT_SCALING["R"]
    assert v_h == DEFAULT_SCALING["H"]
    assert v_r != v_h  # bands have different defaults


# ---------- LFC branch (divides by 20) ----------

def test_lfc_default_divides_band_default_by_20():
    value, _ = resolve_source_scaling("lfc", "R", flux=1.0, user_scaling=None)
    assert value == pytest.approx(DEFAULT_SCALING["R"] / 20)


def test_lfc_flux_propagates_through_divisor():
    value, _ = resolve_source_scaling("lfc", "R", flux=5.0, user_scaling=None)
    assert value == pytest.approx(5.0 * DEFAULT_SCALING["R"] / 20)


def test_lfc_user_scaling_replaces_default_then_divides():
    value, _ = resolve_source_scaling("lfc", "R", flux=1.0, user_scaling=1e6)
    assert value == pytest.approx(1e6 / 20)


# ---------- Fabry-Perot branch (divides by 100) ----------

def test_fp_default_divides_band_default_by_100():
    value, _ = resolve_source_scaling("fabry_perot", "R", flux=1.0, user_scaling=None)
    assert value == pytest.approx(DEFAULT_SCALING["R"] / 100)


def test_fp_flux_propagates_through_divisor():
    value, _ = resolve_source_scaling("fabry_perot", "R", flux=2.0, user_scaling=None)
    assert value == pytest.approx(2.0 * DEFAULT_SCALING["R"] / 100)


def test_fp_user_scaling_replaces_default_then_divides():
    value, _ = resolve_source_scaling("fabry_perot", "R", flux=1.0, user_scaling=1e6)
    assert value == pytest.approx(1e6 / 100)


# ---------- Unknown band falls back to a sane default ----------

def test_unknown_band_constant_uses_default_fallback():
    # falls back to DEFAULT_SCALING.get(band, 1e5)
    value, _ = resolve_source_scaling("constant", "NOT_A_BAND", flux=1.0, user_scaling=None)
    assert value == 1e5
