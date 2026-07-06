"""Tests for the HCL (ThAr) source."""

from pathlib import Path

import numpy as np
import pytest

from andes_simulator.sources.hcl import HCLSource, HCLLineListUnavailable

PROJECT_ROOT = Path(__file__).parent.parent


def test_line_list_loads_and_clips():
    hcl = HCLSource(band="R", flux_scale=1e5, project_root=PROJECT_ROOT)
    wl, flux = hcl.get_lines()
    assert len(wl) > 500
    assert flux.max() == pytest.approx(1e5)
    assert np.all(np.diff(wl) >= 0)

    wl_clip, _ = hcl.get_lines(wl_min=650.0, wl_max=660.0)
    assert len(wl_clip) < len(wl)
    assert wl_clip.min() >= 650.0 and wl_clip.max() <= 660.0


def test_nir_band_has_lines():
    hcl = HCLSource(band="H", flux_scale=1e5, project_root=PROJECT_ROOT)
    wl, _ = hcl.get_lines()
    # J/H are beyond the ThAr atlas; the NIST list must cover them
    assert len(wl) > 200


def test_create_source_object():
    hcl = HCLSource(band="Y", flux_scale=2e4, project_root=PROJECT_ROOT)
    src = hcl._create_hcl_source(wl_min=1000.0, wl_max=1010.0)
    assert src is not None
    assert getattr(src, "list_like", True)


def test_missing_linelist_raises_or_falls_back(tmp_path):
    # a project root without SED/: everything unavailable
    hcl = HCLSource(band="R", flux_scale=1e5, project_root=tmp_path)
    with pytest.raises(HCLLineListUnavailable):
        hcl._load_line_list()


def test_factory_integration():
    from andes_simulator.core.sources import SourceFactory
    from andes_simulator.core.config import SourceConfig

    factory = SourceFactory(PROJECT_ROOT)
    cfg = SourceConfig(type="hcl", scaling_factor=1e5)
    src = factory.create_source(cfg, band="R", wl_min=650.0, wl_max=655.0)
    assert src is not None
