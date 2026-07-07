"""Tests for the DPR.TYPE grammar parser (raw/dpr.py).

Grammar: Templates Manual v2.0 — two slots (A, B), calibration fibre C
via the calfib parameter / ins.calfib keyword.
"""

from pathlib import Path

import pytest

from andes_simulator.raw.dpr import (
    parse_dpr, subslit_fibers, mask_fibers, source_spec_for_token,
    peak_target_e,
)

PROJECT_ROOT = Path(__file__).parent.parent


def test_sl_wave_frame_with_calfib():
    spec = parse_dpr("WAVE,HCL,FP", band="R", calfib="FP")
    assert spec.kind == "WAVE"
    assert spec.catg == "CALIB"
    assert spec.tech == "ECHELLE,FIBER"
    assert spec.mode == "SL-UNI"
    assert spec.calfib == "FP"
    assert not spec.detector_only
    assert [s.name for s in spec.slots] == ["A", "B", "C"]
    assert [s.token for s in spec.slots] == ["HCL", "FP", "FP"]
    assert spec.slots[0].fibers == list(range(1, 32))
    assert spec.slots[1].fibers == list(range(36, 67))
    assert spec.slots[2].fibers == [33, 34]


def test_calfib_off_or_absent_gives_no_c_slot():
    for calfib in (None, "OFF"):
        spec = parse_dpr("WAVE,FP,FP", band="R", calfib=calfib)
        assert [s.name for s in spec.slots] == ["A", "B"]


def test_dark_slots_omitted():
    spec = parse_dpr("WAVE,FP,OFF", band="R", calfib="FP")
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


def test_ifu_wave_with_calfib():
    spec = parse_dpr("WAVE,HCL", band="Y", mode="IFU-AO", calfib="FP")
    assert [s.subslit for s in spec.slots] == ["ifu", "cal_ifu"]
    assert [s.token for s in spec.slots] == ["HCL", "FP"]


def test_science_frame_no_kind():
    spec = parse_dpr("OBJECT,SKY", band="R", calfib="FP")
    assert spec.kind is None
    assert spec.catg == "SCIENCE"
    assert [s.token for s in spec.slots] == ["OBJECT", "SKY", "FP"]


def test_tc_science_wave_token():
    spec = parse_dpr("OBJECT,WAVE", band="R")
    assert spec.catg == "SCIENCE"
    assert [s.token for s in spec.slots] == ["OBJECT", "WAVE"]
    assert source_spec_for_token("WAVE", "R", PROJECT_ROOT) == {"type": "fabry_perot"}


def test_std_is_calib():
    spec = parse_dpr("STD,FLUX,SKY", band="R")
    assert spec.kind == "STD"
    assert spec.catg == "CALIB"
    assert [s.name for s in spec.slots] == ["A", "B"]


def test_slitmask_with_mask_intersection():
    spec = parse_dpr("SLITMASK,FP,OFF", band="R", ins_mask="M1", calfib="FP")
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


def test_old_three_slot_grammar_rejected():
    with pytest.raises(ValueError):
        parse_dpr("WAVE,HCL,FP,FP", band="R")


def test_single_token_without_ifu_mode_rejected():
    with pytest.raises(ValueError):
        parse_dpr("WAVE,HCL", band="R")


def test_unknown_token_rejected():
    with pytest.raises(ValueError):
        parse_dpr("WAVE,XYZ,FP", band="R")


def test_unknown_calfib_rejected():
    with pytest.raises(ValueError):
        parse_dpr("WAVE,FP,FP", band="R", calfib="XYZ")


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
    # twilight sky flat is bright, science sky slot faint
    assert peak_target_e("FLAT", "SKY") > peak_target_e(None, "SKY")
    assert peak_target_e("WAVE", "FP") == 30000.0


def test_subslit_fibers_ifu():
    fibers = subslit_fibers("Y", "ifu")
    assert 3 in fibers and len(fibers) > 50
    assert subslit_fibers("Y", "cal_ifu") == [1, 75]
