"""Tests for the DPR.TYPE grammar parser (raw/dpr.py)."""

from pathlib import Path

import pytest

from andes_simulator.raw.dpr import (
    parse_dpr, subslit_fibers, mask_fibers, source_spec_for_token,
    peak_target_e,
)

PROJECT_ROOT = Path(__file__).parent.parent


def test_sl_wave_frame():
    spec = parse_dpr("WAVE,HCL,FP,FP", band="R")
    assert spec.kind == "WAVE"
    assert spec.catg == "CALIB"
    assert spec.tech == "ECHELLE,FIBER"
    assert spec.mode == "SL-UNI"
    assert not spec.detector_only
    assert [s.name for s in spec.slots] == ["A", "C", "B"]
    assert [s.token for s in spec.slots] == ["HCL", "FP", "FP"]
    assert spec.slots[0].fibers == list(range(1, 32))
    assert spec.slots[1].fibers == [33, 34]
    assert spec.slots[2].fibers == list(range(36, 67))


def test_dark_slots_omitted():
    spec = parse_dpr("WAVE,FP,FP,OFF", band="R")
    assert [s.name for s in spec.slots] == ["A", "C"]


def test_bias_dark_detector_only():
    for dpr_type in ("BIAS", "DARK"):
        spec = parse_dpr(dpr_type)
        assert spec.detector_only
        assert not spec.led
        assert spec.catg == "CALIB"
        assert spec.tech == "IMAGE"
        assert spec.slots == []


def test_led_flat():
    spec = parse_dpr("FLAT,LAMP")
    assert spec.detector_only
    assert spec.led
    assert spec.catg == "TECHNICAL"
    assert spec.tech == "IMAGE"


def test_ifu_flat_not_led():
    spec = parse_dpr("FLAT,LAMP", band="Y", mode="IFU-AO")
    assert not spec.detector_only
    assert spec.tech == "ECHELLE,IFU"
    assert len(spec.slots) == 1
    assert spec.slots[0].subslit == "ifu"


def test_ifu_wave_two_slots():
    spec = parse_dpr("WAVE,HCL,FP", band="Y", mode="IFU-AO")
    assert [s.subslit for s in spec.slots] == ["ifu", "cal_ifu"]
    assert [s.token for s in spec.slots] == ["HCL", "FP"]


def test_science_frame_no_kind():
    spec = parse_dpr("OBJECT,FP,SKY", band="R")
    assert spec.kind is None
    assert spec.catg == "SCIENCE"
    assert [s.token for s in spec.slots] == ["OBJECT", "FP", "SKY"]


def test_std_is_calib():
    spec = parse_dpr("STD,FLUX,OFF,SKY", band="R")
    assert spec.kind == "STD"
    assert spec.catg == "CALIB"
    assert [s.name for s in spec.slots] == ["A", "B"]


def test_mask_intersection():
    spec = parse_dpr("SLIT,FP,FP,FP", band="R", ins_mask="M1")
    all_masked = set(mask_fibers("R", "M1"))
    for slot in spec.slots:
        assert slot.fibers
        assert set(slot.fibers) <= all_masked
    slit_a = spec.slots[0].fibers
    assert slit_a == [f for f in range(1, 32) if (f - 1) % 3 == 0]


def test_mask_patterns_disjoint_and_complete():
    m1 = set(mask_fibers("R", "M1"))
    m2 = set(mask_fibers("R", "M2"))
    m3 = set(mask_fibers("R", "M3"))
    assert not (m1 & m2) and not (m1 & m3) and not (m2 & m3)
    assert m1 | m2 | m3 == set(range(1, 67))


def test_two_tokens_without_ifu_mode_rejected():
    with pytest.raises(ValueError):
        parse_dpr("WAVE,HCL,FP", band="R")


def test_unknown_token_rejected():
    with pytest.raises(ValueError):
        parse_dpr("WAVE,XYZ,FP,FP", band="R")


def test_source_token_mapping():
    assert source_spec_for_token("LAMP", "R", PROJECT_ROOT) == {"type": "constant"}
    assert source_spec_for_token("FP", "R", PROJECT_ROOT) == {"type": "fabry_perot"}
    assert source_spec_for_token("HCL", "R", PROJECT_ROOT) == {"type": "hcl"}

    sky_r = source_spec_for_token("SKY", "R", PROJECT_ROOT)
    assert sky_r == {"type": "csv", "filepath": "SED/sky_emission_R.csv"}

    sky_y = source_spec_for_token("SKY", "Y", PROJECT_ROOT)
    assert sky_y == {"type": "csv", "filepath": "SED/sky_emission_YK.csv"}

    # no sky spectrum covers B-band: faint continuum fallback
    assert source_spec_for_token("SKY", "B", PROJECT_ROOT) == {"type": "constant"}

    obj = source_spec_for_token("OBJECT", "J", PROJECT_ROOT)
    assert obj["type"] == "csv"


def test_peak_targets():
    assert peak_target_e("EFF", "SKY") > peak_target_e("OBJECT", "SKY")
    assert peak_target_e("WAVE", "FP") == 30000.0


def test_subslit_fibers_ifu():
    fibers = subslit_fibers("Y", "ifu")
    assert 3 in fibers and len(fibers) > 50
    assert subslit_fibers("Y", "cal_ifu") == [1, 75]
