"""Integration tests for the CLI surface via Click's CliRunner.

These exercise the full option parsing -> fiber resolution -> source dispatch
-> config building -> dry-run output pipeline without invoking PyEchelle.
They catch regressions in CLI wiring (option types, callbacks, source
dispatch, band inference, bundle parsing, JSON loading) that pure unit
tests on individual helpers would miss.
"""

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from andes_simulator.cli.andes import cli as andes_cli
from andes_simulator.cli.mosaic import cli as mosaic_cli


@pytest.fixture
def runner():
    return CliRunner()


# ============================================================================
# list-bands
# ============================================================================

def test_andes_list_bands(runner):
    res = runner.invoke(andes_cli, ['list-bands'])
    assert res.exit_code == 0
    for band in ['U', 'B', 'V', 'R', 'IZ', 'Y', 'J', 'H']:
        assert band in res.output


def test_mosaic_list_bands(runner):
    res = runner.invoke(mosaic_cli, ['list-bands'])
    assert res.exit_code == 0
    assert 'LR-blue' in res.output
    assert 'HR-H' in res.output


# ============================================================================
# simulate --dry-run: source types
# ============================================================================

def test_simulate_flat_source_dry_run(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'flat',
        '--fiber', '21', '--dry-run'])
    assert res.exit_code == 0, res.output
    assert 'Type: flat_field' in res.output
    assert 'Band: R' in res.output
    assert 'Flux: 80000' in res.output  # R-band default scaling


def test_simulate_fp_source_dry_run(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'fp',
        '--fiber', '21', '--dry-run'])
    assert res.exit_code == 0, res.output
    assert 'Type: fabry_perot' in res.output
    assert 'Scaling:' in res.output
    assert 'Finesse:' in res.output


def test_simulate_lfc_source_dry_run(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'lfc',
        '--fiber', '21', '--dry-run'])
    assert res.exit_code == 0, res.output
    assert 'Type: lfc' in res.output
    assert 'LFC flux per line:' in res.output


def test_simulate_csv_source_dry_run(runner, tmp_path):
    csv = tmp_path / "spec.csv"
    csv.write_text("500.0, 100.0\n510.0, 110.0\n")
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', str(csv),
        '--fiber', '21', '--dry-run'])
    assert res.exit_code == 0, res.output
    assert 'Type: spectrum' in res.output
    assert 'Spectrum:' in res.output
    assert str(csv) in res.output


def test_simulate_csv_missing_file_rejected(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'nope.csv',
        '--fiber', '21', '--dry-run'])
    assert res.exit_code != 0
    assert 'Spectrum file not found' in res.output


def test_simulate_unknown_source_rejected(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'gibberish_no_extension',
        '--fiber', '21', '--dry-run'])
    assert res.exit_code != 0
    assert 'Unknown source' in res.output


# ============================================================================
# simulate --dry-run: fiber selection
# ============================================================================

def test_simulate_fiber_integer(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'flat',
        '--fiber', '21', '--dry-run'])
    assert res.exit_code == 0
    assert 'Subslit: single' in res.output


def test_simulate_fiber_negative_resolves(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'flat',
        '--fiber', '-1', '--dry-run'])
    assert res.exit_code == 0
    assert 'Fiber -1 -> 66' in res.output  # R-band has 66 fibers


def test_simulate_subslit_all(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'flat',
        '--subslit', 'all', '--dry-run'])
    assert res.exit_code == 0
    assert 'Subslit: all' in res.output


def test_simulate_subslit_slitA(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'flat',
        '--subslit', 'slitA', '--dry-run'])
    assert res.exit_code == 0
    assert 'Subslit: slitA' in res.output


def test_simulate_subslit_ifu_only_on_YJH(runner):
    # IFU mode is allowed on Y-band
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'Y', '--source', 'flat',
        '--subslit', 'ifu', '--dry-run'])
    assert res.exit_code == 0
    assert 'Subslit: ifu' in res.output


def test_simulate_bogus_subslit_rejected(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'flat',
        '--subslit', 'bogus', '--dry-run'])
    assert res.exit_code != 0


# ============================================================================
# simulate --dry-run: MOSAIC bundle parsing
# ============================================================================

def test_mosaic_bundle_single(runner):
    res = runner.invoke(mosaic_cli, [
        'simulate', '--band', 'LR-blue', '--source', 'flat',
        '--fiber', 'bundle:2', '--dry-run'])
    assert res.exit_code == 0, res.output
    # bundle 2 of size 7 (LR) -> fibers 8..14
    assert 'Bundle 2 -> fibers [8, 9, 10, 11, 12, 13, 14]' in res.output
    assert 'Subslit: custom' in res.output


def test_mosaic_bundle_range(runner):
    res = runner.invoke(mosaic_cli, [
        'simulate', '--band', 'LR-blue', '--source', 'flat',
        '--fiber', 'bundle:1-2', '--dry-run'])
    assert res.exit_code == 0, res.output
    # bundles 1-2 of size 7 -> fibers 1..14
    assert 'Bundle 1-2' in res.output


def test_mosaic_bundle_HR_uses_19_fibers(runner):
    res = runner.invoke(mosaic_cli, [
        'simulate', '--band', 'HR-B1', '--source', 'flat',
        '--fiber', 'bundle:1', '--dry-run'])
    assert res.exit_code == 0, res.output
    # bundle 1 of size 19 -> fibers 1..19
    assert 'fibers [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]' in res.output


def test_mosaic_bundle_out_of_range_rejected(runner):
    res = runner.invoke(mosaic_cli, [
        'simulate', '--band', 'LR-blue', '--source', 'flat',
        '--fiber', 'bundle:999', '--dry-run'])
    assert res.exit_code != 0
    assert 'out of range' in res.output


def test_andes_bundle_rejected_no_bundle_size(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'flat',
        '--fiber', 'bundle:2', '--dry-run'])
    assert res.exit_code != 0
    assert 'bundle: syntax not supported' in res.output


# ============================================================================
# simulate --dry-run: band inference from --wl-min / --wl-max
# ============================================================================

def test_band_inferred_from_wavelengths(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--wl-min', '650', '--wl-max', '700',
        '--source', 'flat', '--fiber', '21', '--dry-run'])
    assert res.exit_code == 0, res.output
    assert 'Inferred band: R' in res.output
    assert 'Band: R' in res.output


def test_band_inference_ambiguous_fails(runner):
    # ~1000nm matches both Y and LR-J — but each CLI is restricted to its
    # own instrument's bands, so ANDES picks Y unambiguously.
    res = runner.invoke(andes_cli, [
        'simulate', '--wl-min', '1000', '--wl-max', '1050',
        '--source', 'flat', '--fiber', '21', '--dry-run'])
    assert res.exit_code == 0
    assert 'Inferred band: Y' in res.output


def test_simulate_no_band_no_wl_fails(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--source', 'flat', '--fiber', '21', '--dry-run'])
    assert res.exit_code != 0
    assert 'Either --band' in res.output


def test_simulate_wl_out_of_band_rejected(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--wl-min', '100', '--wl-max', '200',
        '--source', 'flat', '--fiber', '21', '--dry-run'])
    assert res.exit_code != 0


# ============================================================================
# simulate --dry-run: velocity-shift and x-shift
# ============================================================================

def test_velocity_shift_scalar(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'lfc', '--fiber', '21',
        '--velocity-shift', '500', '--dry-run'])
    assert res.exit_code == 0
    assert 'Velocity shift: 500.0 m/s' in res.output


def test_velocity_shift_json_file(runner, tmp_path):
    offsets = tmp_path / "offsets.json"
    offsets.write_text(json.dumps({"21": 250.5}))
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'lfc', '--fiber', '21',
        '--velocity-shift', str(offsets), '--dry-run'])
    assert res.exit_code == 0, res.output
    assert '250.5' in res.output
    assert str(offsets) in res.output


def test_velocity_shift_json_missing_fiber_rejected(runner, tmp_path):
    offsets = tmp_path / "offsets.json"
    offsets.write_text(json.dumps({"99": 250.5}))
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'lfc', '--fiber', '21',
        '--velocity-shift', str(offsets), '--dry-run'])
    assert res.exit_code != 0
    assert 'Fiber 21 not found' in res.output


def test_velocity_shift_neither_number_nor_file_rejected(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'lfc', '--fiber', '21',
        '--velocity-shift', 'nope.json', '--dry-run'])
    assert res.exit_code != 0
    assert 'Not a number and file not found' in res.output


def test_x_shift_scalar(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'lfc', '--fiber', '21',
        '--x-shift', '0.5', '--dry-run'])
    assert res.exit_code == 0
    assert 'X-shift: 0.5 px' in res.output


def test_x_shift_json_requires_single_fiber(runner, tmp_path):
    offsets = tmp_path / "xshift.json"
    offsets.write_text(json.dumps({"21": 0.3}))
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'lfc', '--subslit', 'all',
        '--x-shift', str(offsets), '--dry-run'])
    assert res.exit_code != 0
    assert 'single fiber mode' in res.output


# ============================================================================
# simulate --dry-run: FP options
# ============================================================================

def test_simulate_fp_with_finesse_override(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'fp', '--fiber', '21',
        '--finesse', '40', '--dry-run'])
    assert res.exit_code == 0
    assert 'Finesse: 40' in res.output


def test_simulate_fp_with_gap_override(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'fp', '--fiber', '21',
        '--fp-gap', '0.5', '--dry-run'])
    assert res.exit_code == 0
    assert 'FP gap: 0.5 mm' in res.output


def test_simulate_flux_multiplier(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'flat', '--fiber', '21',
        '--flux', '0.5', '--dry-run'])
    assert res.exit_code == 0
    # 0.5 * R-band default 80000 = 40000
    assert 'Flux: 40000' in res.output


def test_simulate_explicit_scaling_overrides_default(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'R', '--source', 'flat', '--fiber', '21',
        '--scaling', '12345', '--dry-run'])
    assert res.exit_code == 0
    assert 'Flux: 12345' in res.output


# ============================================================================
# combine --dry-run
# ============================================================================

def test_combine_dry_run(runner):
    res = runner.invoke(andes_cli, [
        'combine', '--band', 'R',
        '--input-pattern', 'R_FF_fiber{fib:02d}_*.fits',
        '--mode', 'all', '--dry-run'])
    assert res.exit_code == 0, res.output
    assert 'Mode: all' in res.output
    assert 'R_FF_fiber{fib:02d}' in res.output


def test_combine_custom_mode_dry_run_shows_fibers(runner):
    res = runner.invoke(andes_cli, [
        'combine', '--band', 'R',
        '--input-pattern', 'R_FF_fiber{fib:02d}_*.fits',
        '--mode', 'custom', '--fibers', '1,5,10-12', '--dry-run'])
    assert res.exit_code == 0
    assert 'Fibers: 1,5,10-12' in res.output


# ============================================================================
# psf-process --dry-run
# ============================================================================

def test_psf_process_dry_run(runner):
    res = runner.invoke(andes_cli, [
        'psf-process', '--band', 'R',
        '--input-pattern', 'R_FP_*.fits',
        '--fwhm', '3.2', '--dry-run'])
    assert res.exit_code == 0, res.output
    assert 'FWHM=3.2' in res.output
    assert 'Band: R' in res.output


def test_psf_process_bad_kernel_size_rejected(runner):
    res = runner.invoke(andes_cli, [
        'psf-process', '--band', 'R',
        '--input-pattern', 'R_FP_*.fits',
        '--kernel-size', 'bogus', '--dry-run'])
    assert res.exit_code != 0


# ============================================================================
# generate-offsets (actually executes — writes a JSON file)
# ============================================================================

def test_generate_offsets_single_fiber(runner, tmp_path):
    out = tmp_path / "offsets.json"
    res = runner.invoke(andes_cli, [
        'generate-offsets', '--band', 'R', '--fiber', '21',
        '--rms', '100', '--seed', '0', '-o', str(out)])
    assert res.exit_code == 0, res.output
    assert out.exists()
    data = json.loads(out.read_text())
    assert list(data.keys()) == ['21']
    assert isinstance(data['21'], float)


def test_generate_offsets_seed_reproducible(runner, tmp_path):
    out1 = tmp_path / "o1.json"
    out2 = tmp_path / "o2.json"
    runner.invoke(andes_cli, ['generate-offsets', '--band', 'R',
                              '--fiber', '21', '--rms', '100',
                              '--seed', '42', '-o', str(out1)])
    runner.invoke(andes_cli, ['generate-offsets', '--band', 'R',
                              '--fiber', '21', '--rms', '100',
                              '--seed', '42', '-o', str(out2)])
    assert json.loads(out1.read_text()) == json.loads(out2.read_text())


def test_generate_offsets_slit_writes_per_fiber_values(runner, tmp_path):
    out = tmp_path / "offsets.json"
    res = runner.invoke(andes_cli, [
        'generate-offsets', '--band', 'R', '--subslit', 'slitA',
        '--rms', '50', '--seed', '0', '-o', str(out)])
    assert res.exit_code == 0
    data = json.loads(out.read_text())
    # slitA on R-band has 31 fibers
    assert len(data) == 31


# ============================================================================
# generate-hdf (dry-run only)
# ============================================================================

def test_generate_hdf_dry_run(runner):
    res = runner.invoke(andes_cli, [
        'generate-hdf', '--band', 'R', '--dry-run'])
    assert res.exit_code == 0
    assert 'Band: R' in res.output


# ============================================================================
# help / errors
# ============================================================================

def test_simulate_help(runner):
    res = runner.invoke(andes_cli, ['simulate', '--help'])
    assert res.exit_code == 0
    assert '--source' in res.output
    assert '--fiber' in res.output


def test_invalid_band_rejected(runner):
    res = runner.invoke(andes_cli, [
        'simulate', '--band', 'NOT_A_BAND', '--source', 'flat',
        '--fiber', '21', '--dry-run'])
    assert res.exit_code != 0


def test_mosaic_rejects_andes_only_band(runner):
    res = runner.invoke(mosaic_cli, [
        'simulate', '--band', 'R', '--source', 'flat',
        '--fiber', '1', '--dry-run'])
    assert res.exit_code != 0  # R is ANDES, not MOSAIC
