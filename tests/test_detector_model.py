"""Tests for the detector layer (raw/detector.py)."""

import numpy as np
import pytest

from andes_simulator.raw.detector import DetectorModel, spectrograph_for_band


def test_spectrograph_lookup():
    assert spectrograph_for_band("R")[0] == "RIZ"
    assert spectrograph_for_band("Y")[0] == "YJH"
    with pytest.raises(ValueError):
        spectrograph_for_band("Y_iq15")


def test_bias_frame_statistics():
    det = DetectorModel("Y")
    frame = det.apply(None, exptime=0.0, rng=np.random.default_rng(1))
    assert frame.dtype == np.uint16
    assert frame.shape == (4096, 4096)
    assert abs(float(frame.mean()) - det.params["bias_adu"]) < 2.0
    ron_adu = det.ron_e / det.gain_e_adu
    assert float(frame.std()) == pytest.approx(ron_adu, rel=0.15)


def test_dark_frame_has_hot_pixels():
    det = DetectorModel("Y")
    bias = det.apply(None, exptime=0.0, rng=np.random.default_rng(2)).astype(float)
    dark = det.apply(None, exptime=1800.0, rng=np.random.default_rng(3)).astype(float)
    diff_mean = dark.mean() - bias.mean()
    expected = det.params["dark_e_s"] * 1800.0 / det.gain_e_adu
    assert diff_mean == pytest.approx(expected, rel=0.3)
    # hot pixels stand far above the normal dark level
    assert (dark > bias.mean() + 50 * expected).sum() > 100


def test_noise_differs_but_cosmetics_static():
    det = DetectorModel("Y")
    f1 = det.apply(None, exptime=1800.0, rng=np.random.default_rng(10)).astype(float)
    f2 = det.apply(None, exptime=1800.0, rng=np.random.default_rng(11)).astype(float)
    assert not np.array_equal(f1, f2)
    # hot pixels (strong outliers) must sit at identical positions
    hot1 = f1 > f1.mean() + 10 * f1.std()
    hot2 = f2 > f2.mean() + 10 * f2.std()
    n_common = (hot1 & hot2).sum()
    assert n_common > 0.9 * max(hot1.sum(), hot2.sum())


def test_reproducible_with_same_seed():
    det = DetectorModel("Y")
    f1 = det.apply(None, exptime=0.0, rng=np.random.default_rng(42))
    f2 = det.apply(None, exptime=0.0, rng=np.random.default_rng(42))
    assert np.array_equal(f1, f2)


def test_led_flat_level_and_prnu():
    det = DetectorModel("Y")
    exp = det.led_expectation(exptime=10.0)
    frame = det.apply(exp, exptime=10.0, rng=np.random.default_rng(5)).astype(float)
    signal_adu = 800.0 * 10.0 / det.gain_e_adu
    assert frame.mean() - det.params["bias_adu"] == pytest.approx(signal_adu, rel=0.05)
    # PRNU dominates over shot noise at this level: relative spread per pixel
    sig = frame - det.params["bias_adu"]
    assert sig.std() / sig.mean() > det.params["prnu_rms"] * 0.8


def test_ccd_binning_sums_charge():
    det = DetectorModel("R", readout="slow")
    exp = det.led_expectation(exptime=5.0)
    rng = np.random.default_rng(7)
    unbinned = det.apply(exp, exptime=5.0, rng=np.random.default_rng(7))
    binned = det.apply(exp, exptime=5.0, rng=rng, binx=2, biny=2)
    assert binned.shape == (unbinned.shape[0] // 2, unbinned.shape[1] // 2)
    sig_unb = unbinned.astype(float).mean() - det.params["bias_adu"]
    sig_bin = binned.astype(float).mean() - det.params["bias_adu"]
    assert sig_bin == pytest.approx(4 * sig_unb, rel=0.05)


def test_dead_columns_present_on_ccd():
    det = DetectorModel("R")
    assert len(det.bad_columns) == det.params["n_bad_columns"]
    frame = det.apply(det.led_expectation(5.0), exptime=5.0,
                      rng=np.random.default_rng(8)).astype(float)
    for col in det.bad_columns:
        assert frame[:, col].mean() < det.params["bias_adu"] + 5


def test_saturation_clipping():
    det = DetectorModel("Y")
    exp = np.full(det.shape, 1e6, dtype=np.float32)
    frame = det.apply(exp, exptime=1.0, rng=np.random.default_rng(9))
    assert frame.max() == 65535
