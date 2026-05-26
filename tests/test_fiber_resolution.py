"""Tests for resolve_negative_fiber (Python-style negative --fiber indexing)."""

import click
import pytest

from andes_simulator.cli.main import resolve_negative_fiber
from andes_simulator.core.instruments import get_instrument_config


# Bands with known fiber counts, sampled across ANDES and MOSAIC.
# (band, n_fibers)
BAND_CASES = [
    ('R', 66),       # ANDES UBVRIZ
    ('Y', 75),       # ANDES YJH
    ('LR-blue', 980),  # MOSAIC VIS LR
    ('HR-B1', 1140),   # MOSAIC VIS HR
    ('HR-H', 630),     # MOSAIC NIR
]


@pytest.mark.parametrize('band, n_fibers', BAND_CASES)
def test_config_n_fibers_matches_expected(band, n_fibers):
    """Guardrail: if instrument configs change n_fibers, this suite needs updating."""
    assert get_instrument_config(band)['n_fibers'] == n_fibers


@pytest.mark.parametrize('band, n_fibers', BAND_CASES)
def test_negative_one_maps_to_last_fiber(band, n_fibers):
    assert resolve_negative_fiber(-1, band) == n_fibers


@pytest.mark.parametrize('band, n_fibers', BAND_CASES)
def test_negative_two_maps_to_second_last(band, n_fibers):
    assert resolve_negative_fiber(-2, band) == n_fibers - 1


@pytest.mark.parametrize('band, n_fibers', BAND_CASES)
def test_most_negative_valid_maps_to_first_fiber(band, n_fibers):
    assert resolve_negative_fiber(-n_fibers, band) == 1


@pytest.mark.parametrize('band, n_fibers', BAND_CASES)
def test_out_of_range_negative_raises(band, n_fibers):
    with pytest.raises(click.UsageError, match='out of range'):
        resolve_negative_fiber(-(n_fibers + 1), band)


@pytest.mark.parametrize('value', [1, 21, 66, 100])
def test_positive_int_passthrough(value):
    assert resolve_negative_fiber(value, 'R') == value


def test_zero_passthrough():
    # Zero isn't a valid fiber but resolution shouldn't touch it; downstream code rejects it.
    assert resolve_negative_fiber(0, 'R') == 0


@pytest.mark.parametrize('value', ['all', 'slitA', 'cal_sl', 'ifu', 'bundle:1', 'bundle:1-10'])
def test_string_modes_passthrough(value):
    assert resolve_negative_fiber(value, 'R') == value


def test_emits_echo_on_resolution(capsys):
    resolve_negative_fiber(-3, 'R')
    captured = capsys.readouterr()
    assert 'Fiber -3 -> 64' in captured.out
    assert 'R-band' in captured.out


def test_no_echo_for_positive(capsys):
    resolve_negative_fiber(21, 'R')
    assert capsys.readouterr().out == ''


def test_no_echo_for_string_mode(capsys):
    resolve_negative_fiber('all', 'R')
    assert capsys.readouterr().out == ''


def test_unknown_band_raises():
    with pytest.raises(ValueError, match='Unknown band'):
        resolve_negative_fiber(-1, 'NOT_A_BAND')


def test_error_message_includes_valid_range():
    with pytest.raises(click.UsageError) as exc:
        resolve_negative_fiber(-67, 'R')
    assert '-66..-1' in str(exc.value)
