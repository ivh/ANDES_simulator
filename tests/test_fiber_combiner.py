"""Tests for FiberCombiner using synthetic FITS fixtures.

We write tiny FITS files containing per-fiber arrays of known content
(typically value = fiber number), then verify combination modes produce
the expected pixel sums. Skip-fibers, header propagation, and per-mode
arithmetic are all checked without involving PyEchelle.
"""

from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from andes_simulator.postprocess.combine import FiberCombiner


# Small detector for fast tests. Note: FiberCombiner reads the band config
# for n_fibers and detector_size, so we patch the detector_size after init
# to keep the arrays tiny.
SMALL_SHAPE = (8, 8)


def _write_fiber_fits(path: Path, value: int, headers: dict = None):
    """Write a small FITS file filled with a single value."""
    data = np.full(SMALL_SHAPE, value, dtype=np.uint16)
    hdu = fits.PrimaryHDU(data)
    if headers:
        for k, v in headers.items():
            hdu.header[k] = v
    hdu.writeto(path, overwrite=True)


@pytest.fixture
def combiner_R(tmp_path):
    """FiberCombiner for R-band, with tiny detector size and tmp input dir."""
    c = FiberCombiner(band='R', input_dir=tmp_path)
    c.detector_size = SMALL_SHAPE  # patch to keep fixtures small
    return c


@pytest.fixture
def fiber_files_R(tmp_path):
    """Write FITS files for all R-band fibers (1..66); value = fiber number."""
    for fib in range(1, 67):
        path = tmp_path / f"R_FF_fiber{fib:02d}_30s.fits"
        _write_fiber_fits(path, value=fib,
                          headers={'HDFMODEL': 'ANDES_123_R3',
                                   'EXPTIME': 30.0})
    return tmp_path


# ---------- load_fiber_data ----------

def test_load_fiber_data_returns_array(combiner_R, fiber_files_R):
    arr = combiner_R.load_fiber_data(21, "R_FF_fiber{fib:02d}_*.fits")
    assert arr.shape == SMALL_SHAPE
    assert np.all(arr == 21)


def test_load_fiber_data_skip_dark_returns_zeros(combiner_R, fiber_files_R):
    # R-band has skip_fibers=[32, 35]
    arr = combiner_R.load_fiber_data(32, "R_FF_fiber{fib:02d}_*.fits",
                                     skip_dark=True)
    assert np.all(arr == 0)


def test_load_fiber_data_skip_dark_false_loads_file(combiner_R, fiber_files_R):
    arr = combiner_R.load_fiber_data(32, "R_FF_fiber{fib:02d}_*.fits",
                                     skip_dark=False)
    assert np.all(arr == 32)


def test_load_fiber_data_missing_file_returns_zeros(combiner_R, tmp_path):
    # No files written; expect zeros (not None — the function fills with 0)
    arr = combiner_R.load_fiber_data(21, "R_FF_fiber{fib:02d}_*.fits")
    assert np.all(arr == 0)


# ---------- combine_all_fibers ----------

def test_combine_all_fibers_sums_all_except_skipped(combiner_R, fiber_files_R):
    combined = combiner_R.combine_all_fibers("R_FF_fiber{fib:02d}_*.fits")
    # Sum of 1..66 minus the skipped fibers 32 and 35
    expected = sum(range(1, 67)) - 32 - 35
    assert combined[0, 0] == expected


def test_combine_all_fibers_skip_dark_false_includes_skipped(combiner_R, fiber_files_R):
    combined = combiner_R.combine_all_fibers(
        "R_FF_fiber{fib:02d}_*.fits", skip_dark=False)
    expected = sum(range(1, 67))
    assert combined[0, 0] == expected


def test_combine_all_fibers_returns_correct_shape(combiner_R, fiber_files_R):
    combined = combiner_R.combine_all_fibers("R_FF_fiber{fib:02d}_*.fits")
    assert combined.shape == SMALL_SHAPE


# ---------- combine_fiber_subset ----------

def test_combine_fiber_subset(combiner_R, fiber_files_R):
    combined = combiner_R.combine_fiber_subset(
        [1, 2, 3], "R_FF_fiber{fib:02d}_*.fits")
    assert combined[0, 0] == 1 + 2 + 3


def test_combine_fiber_subset_skipped_fibers_pass_through_load_logic(
        combiner_R, fiber_files_R):
    # combine_fiber_subset calls load_fiber_data with default skip_dark=True,
    # so skipped fibers in the subset contribute zero.
    combined = combiner_R.combine_fiber_subset(
        [1, 32, 35], "R_FF_fiber{fib:02d}_*.fits")
    assert combined[0, 0] == 1  # 32 and 35 are zero


# ---------- combine_even_odd ----------

def test_combine_even_odd_partitions_correctly(combiner_R, fiber_files_R):
    results = combiner_R.combine_even_odd_fibers("R_FF_fiber{fib:02d}_*.fits")
    assert set(results.keys()) == {'even', 'odd'}
    # Even fibers 2..66 sum (skipped fiber 32 is even, returns zero)
    expected_even = sum(f for f in range(2, 67, 2) if f != 32)
    expected_odd = sum(f for f in range(1, 67, 2) if f != 35)
    assert results['even'][0, 0] == expected_even
    assert results['odd'][0, 0] == expected_odd


# ---------- combine_pseudo_slits ----------

def test_combine_pseudo_slits_R_band(combiner_R, fiber_files_R):
    results = combiner_R.combine_pseudo_slits("R_FF_fiber{fib:02d}_*.fits")
    # R-band has slitA=1..31, slitB=36..66, cal_fibers=[33, 34]
    assert results['first_slit'][0, 0] == sum(range(1, 32))
    assert results['second_slit'][0, 0] == sum(range(36, 67))
    assert results['calibration'][0, 0] == 33 + 34


def test_combine_pseudo_slits_rejects_YJH_band(tmp_path):
    combiner = FiberCombiner(band='Y', input_dir=tmp_path)
    with pytest.raises(ValueError, match="not applicable to YJH"):
        combiner.combine_pseudo_slits("Y_FF_*.fits")


# ---------- save_combined_image ----------

def test_save_combined_image_writes_with_metadata(combiner_R, fiber_files_R, tmp_path):
    combined = combiner_R.combine_all_fibers("R_FF_fiber{fib:02d}_*.fits")
    out_path = tmp_path / "out.fits"
    combiner_R.save_combined_image(combined, out_path, {'mode': 'all'})

    with fits.open(out_path) as hdul:
        hdr = hdul[0].header
        assert hdr['BAND'] == 'R'
        assert hdr['INSTRUME'] == 'ANDES'
        assert hdr['NFIBERS'] == 66
        assert hdr['COMBINED'] is True
        assert hdr['MODE'] == 'all'
        # skip_fibers header
        assert '32' in str(hdr['SKIPFIB'])
        # Propagated header from source files
        assert hdr['HDFMODEL'] == 'ANDES_123_R3'


def test_save_combined_image_propagates_uniform_headers(combiner_R, fiber_files_R, tmp_path):
    combined = combiner_R.combine_all_fibers("R_FF_fiber{fib:02d}_*.fits")
    out_path = tmp_path / "out.fits"
    combiner_R.save_combined_image(combined, out_path)
    with fits.open(out_path) as hdul:
        # EXPTIME was 30.0 for all fibers -> uniform, propagated as single key
        assert hdul[0].header['EXPTIME'] == 30.0


def test_save_combined_image_creates_per_fiber_keys_for_varying_headers(combiner_R, tmp_path):
    """When a header differs per fiber, it should be propagated as VSHIFTNN."""
    # Custom files where each fiber has a different VSHIFT
    for fib in range(1, 5):
        path = tmp_path / f"R_FF_fiber{fib:02d}_30s.fits"
        _write_fiber_fits(path, value=fib,
                          headers={'VSHIFT': fib * 100.0})

    # Patch n_fibers down so we don't need to write all 66
    combiner_R.n_fibers = 4
    combined = combiner_R.combine_all_fibers("R_FF_fiber{fib:02d}_*.fits")
    out_path = tmp_path / "out.fits"
    combiner_R.save_combined_image(combined, out_path)

    with fits.open(out_path) as hdul:
        hdr = hdul[0].header
        # VSHIFT differed per fiber -> emitted as VSHIF01, VSHIF02, ...
        assert 'VSHIF01' in hdr
        assert hdr['VSHIF01'] == 100.0
        assert hdr['VSHIF04'] == 400.0


# ---------- analyze_fiber_contributions ----------

def test_analyze_fiber_contributions(combiner_R, fiber_files_R):
    analysis = combiner_R.analyze_fiber_contributions(
        "R_FF_fiber{fib:02d}_*.fits")
    assert analysis['band'] == 'R'
    assert analysis['n_active_fibers'] == 64  # 66 - 2 skipped
    assert analysis['n_skipped_fibers'] == 2

    # Each fiber's total_flux = fiber_number * pixel_count
    stats = analysis['fiber_stats']
    assert stats[1]['total_flux'] == 1 * np.prod(SMALL_SHAPE)
    assert stats[21]['total_flux'] == 21 * np.prod(SMALL_SHAPE)
    assert stats[32]['skipped'] is True


def test_analyze_fiber_contributions_fractions_sum_to_one(combiner_R, fiber_files_R):
    analysis = combiner_R.analyze_fiber_contributions(
        "R_FF_fiber{fib:02d}_*.fits")
    fractions = [s['flux_fraction'] for s in analysis['fiber_stats'].values()
                 if not s['skipped']]
    assert sum(fractions) == pytest.approx(1.0)


# ---------- create_combination_report ----------

def test_create_combination_report(combiner_R, fiber_files_R, tmp_path):
    report_path = combiner_R.create_combination_report(
        "R_FF_fiber{fib:02d}_*.fits", output_dir=tmp_path)
    assert report_path.exists()
    content = report_path.read_text()
    assert 'R-band' in content
    assert 'Total fibers: 66' in content
    assert 'Active fibers: 64' in content
    assert 'SKIPPED' in content  # for fibers 32 and 35
