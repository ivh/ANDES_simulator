"""Tests for RawFrameBuilder and the night planner (detector-only paths)."""

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from andes_simulator.raw.builder import RawFrameBuilder
from andes_simulator.raw.cache import SimCache
from andes_simulator.raw.night import plan_night, run_night, _runs

PROJECT_ROOT = Path(__file__).parent.parent


@pytest.fixture
def builder(tmp_path):
    cache = SimCache(tmp_path / "cache", PROJECT_ROOT)
    return RawFrameBuilder(PROJECT_ROOT, cache, tmp_path / "out", seed=1)


def test_bias_mef_structure(builder):
    t0 = datetime(2026, 7, 6, 10, 0, tzinfo=timezone.utc)
    paths = builder.build(arm="YJH", dpr_type="BIAS", exptime=0.0, nexp=2,
                          tpl_start=t0, tpl_id="ANDES_gen_cal_bias")
    # timestamped names: DATE-OBS with colons replaced, 60s readout overhead
    assert [p.name for p in paths] == [
        "ANDES_YJH_2026-07-06T10_00_00.000.fits",
        "ANDES_YJH_2026-07-06T10_01_00.000.fits",
    ]

    with fits.open(paths[0]) as hdul:
        assert hdul[0].data is None
        assert [h.name for h in hdul[1:]] == ["Y", "J", "H"]
        hdr = hdul[0].header
        assert hdr["INSTRUME"] == "ANDES"
        assert hdr["HIERARCH ESO SEQ ARM"] == "YJH"
        assert hdr["HIERARCH ESO DPR CATG"] == "CALIB"
        assert hdr["HIERARCH ESO DPR TYPE"] == "BIAS"
        assert hdr["HIERARCH ESO DPR TECH"] == "IMAGE"
        assert hdr["HIERARCH ESO TPL START"] == "2026-07-06T10:00:00"
        assert hdr["HIERARCH ESO TPL NEXP"] == 2
        assert hdr["HIERARCH ESO TPL EXPNO"] == 1
        assert hdul["Y"].data.dtype == np.uint16

    with fits.open(paths[1]) as hdul:
        assert hdul[0].header["HIERARCH ESO TPL EXPNO"] == 2
        # same tpl.start, later mjd-obs
        assert hdul[0].header["HIERARCH ESO TPL START"] == "2026-07-06T10:00:00"


def test_led_flat_is_technical_and_illuminated(builder):
    paths = builder.build(arm="YJH", dpr_type="FLAT,LAMP", exptime=10.0)
    with fits.open(paths[0]) as hdul:
        hdr = hdul[0].header
        assert hdr["HIERARCH ESO DPR CATG"] == "TECHNICAL"
        assert hdr["HIERARCH ESO DPR TECH"] == "IMAGE"
        data = hdul["Y"].data.astype(float)
        assert data.mean() > 3000  # well above bias level


def test_seed_reproducibility(tmp_path):
    def gen(sub):
        cache = SimCache(tmp_path / "cache", PROJECT_ROOT)
        b = RawFrameBuilder(PROJECT_ROOT, cache, tmp_path / sub, seed=99)
        return b.build(arm="YJH", dpr_type="BIAS", exptime=0.0)[0]

    with fits.open(gen("a")) as h1, fits.open(gen("b")) as h2:
        assert np.array_equal(h1["Y"].data, h2["Y"].data)


def test_bands_subset_flags_extensions(builder):
    paths = builder.build(arm="YJH", dpr_type="FLAT,LAMP", exptime=5.0,
                          bands=["Y"])
    with fits.open(paths[0]) as hdul:
        assert hdul["Y"].header["HIERARCH ESO SIM SIMULATED"]
        assert not hdul["J"].header["HIERARCH ESO SIM SIMULATED"]


def test_headers_only_spectral_frame(tmp_path):
    cache = SimCache(tmp_path / "cache", PROJECT_ROOT)
    b = RawFrameBuilder(PROJECT_ROOT, cache, tmp_path / "out", seed=4,
                        headers_only=True)
    paths = b.build(arm="RIZ", dpr_type="WAVE,HCL,FP", calfib="OFF",
                    exptime=120.0, readout="fast")
    with fits.open(paths[0]) as hdul:
        hdr = hdul[0].header
        assert hdr["HIERARCH ESO DPR TYPE"] == "WAVE,HCL,FP"
        assert hdr["HIERARCH ESO DPR TECH"] == "ECHELLE,FIBER"
        assert hdr["HIERARCH ESO INS CALFIB"] == "OFF"
        assert hdr["HIERARCH ESO DET READOUT"] == "fast"
        assert hdul["R"].data.shape == (2, 2)
        assert not hdul["R"].header["HIERARCH ESO SIM SIMULATED"]
    # no pyechelle run happened: cache stayed empty
    assert not list((tmp_path / "cache").glob("*.fits"))


def test_calfib_keyword_written(builder):
    paths = builder.build(arm="YJH", dpr_type="FLAT,LAMP", exptime=5.0)
    with fits.open(paths[0]) as hdul:
        # detector-only LED flat carries no calfib keyword
        assert "HIERARCH ESO INS CALFIB" not in hdul[0].header


MINI_PLAN = {
    "setups": {
        "arms": {
            "RIZ": {"detectors": ["R", "IZ"], "detector_type": "CCD", "has_bias": True},
            "YJH": {"detectors": ["Y", "J", "H"], "detector_type": "HAWAII4RG", "has_bias": False},
        },
        "vis_configs": [
            {"name": "1x1_fast", "det.binx": 1, "det.biny": 1, "det.readout": "fast", "set": "daily"},
        ],
    },
    "procedures": [
        {"cp": "C-B", "template": "tpl_bias", "set": "detector",
         "applies_to": {"arms": ["RIZ"], "configs": "all_vis"},
         "dpr": {"catg": "CALIB", "tech": "IMAGE"},
         "exposures": [{"type": "BIAS", "n": 3, "exptime_s": 0}]},
        {"cp": "C-wl", "template": "tpl_wave", "set": "daily",
         "applies_to": {"arms": ["RIZ", "YJH"], "modes": ["SL-UNI"]},
         "dpr": {"catg": "CALIB", "tech": "ECHELLE,FIBER"},
         "exposures": [{"type": "WAVE,FP,OFF", "n": 2, "exptime_s": 60,
                        "keywords": {"ins.calfib": "FP"}},
                       {"type": "WAVE,OFF,FP", "n": 2, "exptime_s": 60,
                        "keywords": {"ins.calfib": "FP"}}]},
        {"cp": "ref", "reference": ["C-B"], "applies_to": {"arms": ["RIZ"]}},
    ],
}


def test_plan_night_filters_and_grouping():
    planned = plan_night(MINI_PLAN, arms=["RIZ"], sets=["detector"])
    assert len(planned) == 1
    assert planned[0].dpr_type == "BIAS"
    assert planned[0].n == 3
    assert planned[0].readout == "fast"

    planned = plan_night(MINI_PLAN, arms=["RIZ", "YJH"])
    wave = [p for p in planned if p.dpr_type.startswith("WAVE")]
    assert len(wave) == 4  # 2 exposure specs x 2 arms
    assert all(p.calfib == "FP" for p in wave)
    # wave exposures of one arm share the tpl.start group
    riz_groups = {p.group for p in wave if p.arm == "RIZ"}
    assert len(riz_groups) == 1
    # NIR arm gets no binning/readout keywords
    yjh_wave = [p for p in wave if p.arm == "YJH"][0]
    assert yjh_wave.readout is None


def test_run_night_shares_tpl_start_within_group(tmp_path):
    cache = SimCache(tmp_path / "cache", PROJECT_ROOT)
    b = RawFrameBuilder(PROJECT_ROOT, cache, tmp_path / "out", seed=3)
    planned = plan_night(MINI_PLAN, arms=["RIZ"], sets=["detector"])
    t0 = datetime(2026, 7, 6, 10, 0, tzinfo=timezone.utc)
    written = run_night(Path("plan.yaml"), tmp_path / "out", b, planned, start=t0)
    assert len(written) == 3
    starts, mjds, expnos = set(), [], []
    for p in written:
        with fits.open(p) as hdul:
            starts.add(hdul[0].header["HIERARCH ESO TPL START"])
            mjds.append(hdul[0].header["MJD-OBS"])
            expnos.append(hdul[0].header["HIERARCH ESO TPL EXPNO"])
    assert len(starts) == 1
    assert expnos == [1, 2, 3]
    assert mjds == sorted(mjds) and len(set(mjds)) == 3


def test_runs_collapse():
    assert _runs([60.0, 60.0, 120.0]) == [(60.0, 2), (120.0, 1)]
    assert _runs([1.0, 2.0, 1.0]) == [(1.0, 1), (2.0, 1), (1.0, 1)]
