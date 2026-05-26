"""Tests for SourceFactory utility helpers.

These convert flux units and wavelengths into PyEchelle's expected
formats and parse CSV scaling headers. Each is a single small function
that all source-loading downstream relies on, so regressions cascade.
"""

import textwrap
from pathlib import Path

import pytest

from andes_simulator.core.sources import SourceFactory


@pytest.fixture
def factory(tmp_path):
    return SourceFactory(project_root=tmp_path)


# ---------- _to_nanometers ----------

def test_to_nanometers_passthrough_for_nm(factory):
    assert factory._to_nanometers(500.0, "nm") == 500.0


def test_to_nanometers_from_AA(factory):
    assert factory._to_nanometers(5000.0, "AA") == pytest.approx(500.0)


def test_to_nanometers_from_angstrom_synonym(factory):
    assert factory._to_nanometers(5000.0, "angstrom") == pytest.approx(500.0)


def test_to_nanometers_from_micron(factory):
    assert factory._to_nanometers(0.5, "um") == pytest.approx(500.0)


def test_to_nanometers_from_micron_synonyms(factory):
    assert factory._to_nanometers(0.5, "microns") == pytest.approx(500.0)
    assert factory._to_nanometers(0.5, "micron") == pytest.approx(500.0)


def test_to_nanometers_unknown_unit_assumes_nm(factory):
    # The factory logs a warning and returns the value unchanged
    assert factory._to_nanometers(500.0, "bogus_unit") == 500.0


def test_to_nanometers_none_unit_passthrough(factory):
    assert factory._to_nanometers(500.0, None) == 500.0


# ---------- _convert_flux_units ----------

def test_convert_flux_ph_per_s_per_AA_passthrough(factory):
    assert factory._convert_flux_units(10.0, "ph/s/AA") == 10.0


def test_convert_flux_ph_per_s_per_nm_divides_by_10(factory):
    # 1 nm = 10 AA, so ph/s/nm -> ph/s/AA divides by 10
    assert factory._convert_flux_units(10.0, "ph/s/nm") == pytest.approx(1.0)


def test_convert_flux_ph_per_s_divides_by_default_bandwidth(factory):
    # Default bandwidth is 100 AA
    assert factory._convert_flux_units(100.0, "ph/s") == pytest.approx(1.0)


def test_convert_flux_unknown_unit_raises(factory):
    with pytest.raises(ValueError, match="Unknown flux_unit"):
        factory._convert_flux_units(10.0, "bogus")


# ---------- _parse_csv_scaling ----------

def test_parse_csv_scaling_with_header_comment(tmp_path):
    csv_path = tmp_path / "spec.csv"
    csv_path.write_text(textwrap.dedent("""\
        # scaling: 1.5e6
        500.0, 100.0
        510.0, 110.0
    """))
    assert SourceFactory._parse_csv_scaling(csv_path) == 1.5e6


def test_parse_csv_scaling_returns_none_without_header(tmp_path):
    csv_path = tmp_path / "spec.csv"
    csv_path.write_text("500.0, 100.0\n510.0, 110.0\n")
    assert SourceFactory._parse_csv_scaling(csv_path) is None


def test_parse_csv_scaling_with_multiple_comments(tmp_path):
    csv_path = tmp_path / "spec.csv"
    csv_path.write_text(textwrap.dedent("""\
        # comment line
        # scaling: 2.0
        # another comment
        500.0, 100.0
    """))
    assert SourceFactory._parse_csv_scaling(csv_path) == 2.0


def test_parse_csv_scaling_case_insensitive(tmp_path):
    csv_path = tmp_path / "spec.csv"
    csv_path.write_text("# Scaling: 3.5\n500.0, 100.0\n")
    assert SourceFactory._parse_csv_scaling(csv_path) == 3.5


def test_parse_csv_scaling_invalid_value_returns_none(tmp_path):
    csv_path = tmp_path / "spec.csv"
    csv_path.write_text("# scaling: not_a_number\n500.0, 100.0\n")
    assert SourceFactory._parse_csv_scaling(csv_path) is None


def test_parse_csv_scaling_stops_at_first_data_line(tmp_path):
    # The header comment must come before any data
    csv_path = tmp_path / "spec.csv"
    csv_path.write_text("500.0, 100.0\n# scaling: 99\n")
    assert SourceFactory._parse_csv_scaling(csv_path) is None


# ---------- create_dark_source ----------

def test_create_dark_source_returns_zero_flux(factory):
    from pyechelle.sources import ConstantPhotonFlux
    src = factory.create_dark_source()
    assert isinstance(src, ConstantPhotonFlux)


# ---------- _scale_csv_flux ----------

def test_scale_csv_flux_multiplies_flux_column(factory, tmp_path):
    csv_path = tmp_path / "in.csv"
    csv_path.write_text("500.0,100.0\n510.0,200.0\n")
    out_path = factory._scale_csv_flux(csv_path, 2.5)
    lines = Path(out_path).read_text().strip().split("\n")
    assert lines[0].split(",") == ["500.0", "250.0"]
    assert lines[1].split(",") == ["510.0", "500.0"]


def test_scale_csv_flux_is_cached(factory, tmp_path):
    csv_path = tmp_path / "in.csv"
    csv_path.write_text("500.0,100.0\n")
    p1 = factory._scale_csv_flux(csv_path, 2.0)
    p2 = factory._scale_csv_flux(csv_path, 2.0)
    assert p1 == p2  # cache hit returns same path
