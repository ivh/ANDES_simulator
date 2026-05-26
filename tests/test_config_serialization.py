"""Tests for SimulationConfig serialization (from_dict, from_yaml, to_yaml).

YAML round-trip is the contract for the `run-config` command path; if the
nested dataclass handling drifts, every YAML-driven simulation breaks.

Also smoke-tests create_template_configs, which the gaps audit flagged as
likely broken — it currently uses mode='even_odd' and type='zemax' which
neither pass validate(). Test is xfail until the function is fixed.
"""

from pathlib import Path

import pytest
import yaml

from andes_simulator.core.config import (
    SimulationConfig, SourceConfig, FiberConfig, OutputConfig, PSFConfig,
    create_template_configs,
)


# ---------- from_dict ----------

def test_from_dict_minimal():
    cfg = SimulationConfig.from_dict({
        'simulation_type': 'flat_field',
        'band': 'R',
        'source': {'type': 'constant', 'flux': 0.5},
        'fibers': {'mode': 'single', 'fibers': [21]},
    })
    assert cfg.band == 'R'
    assert cfg.source.type == 'constant'
    assert cfg.source.flux == 0.5
    assert cfg.fibers.mode == 'single'
    assert cfg.fibers.fibers == [21]


def test_from_dict_with_nested_output():
    cfg = SimulationConfig.from_dict({
        'simulation_type': 'flat_field',
        'band': 'R',
        'source': {'type': 'constant'},
        'fibers': {'mode': 'all'},
        'output': {'directory': '/tmp/out', 'filename': 'custom.fits'},
    })
    assert cfg.output.directory == '/tmp/out'
    assert cfg.output.filename == 'custom.fits'


def test_from_dict_with_psf_config():
    cfg = SimulationConfig.from_dict({
        'simulation_type': 'flat_field',
        'band': 'R',
        'source': {'type': 'constant'},
        'fibers': {'mode': 'all'},
        'psf': {'enabled': True, 'sigma': 2.5},
    })
    assert cfg.psf.enabled is True
    assert cfg.psf.sigma == 2.5


def test_from_dict_propagates_validation_errors():
    with pytest.raises(ValueError, match="Invalid band"):
        SimulationConfig.from_dict({
            'simulation_type': 'flat_field',
            'band': 'NOT_A_BAND',
            'source': {'type': 'constant'},
            'fibers': {'mode': 'all'},
        })


# ---------- from_yaml / to_yaml round-trip ----------

def test_to_yaml_writes_file(tmp_path):
    cfg = SimulationConfig(
        simulation_type='flat_field',
        band='R',
        source=SourceConfig(type='constant', flux=0.5),
        fibers=FiberConfig(mode='single', fibers=[21]),
        output=OutputConfig(directory='/tmp/out'),
    )
    yaml_path = tmp_path / "cfg.yaml"
    cfg.to_yaml(yaml_path)
    assert yaml_path.exists()

    # Loaded with UnsafeLoader to bypass the !!python/tuple bug (see xfail
    # round-trip test below); confirms to_yaml writes the expected content.
    data = yaml.load(yaml_path.read_text(), Loader=yaml.UnsafeLoader)
    assert data['band'] == 'R'
    assert data['source']['type'] == 'constant'
    assert data['fibers']['mode'] == 'single'


@pytest.mark.xfail(reason=(
    "to_yaml emits PSFConfig.kernel_size as !!python/tuple, which "
    "from_yaml's safe_load rejects. Fix: convert tuple to list in to_yaml, "
    "or use yaml.safe_dump."
), strict=True)
def test_yaml_round_trip_preserves_core_fields(tmp_path):
    original = SimulationConfig(
        simulation_type='flat_field',
        band='R',
        exposure_time=30.0,
        source=SourceConfig(type='constant', flux=0.5, scaling_factor=1e5),
        fibers=FiberConfig(mode='single', fibers=[21]),
        output=OutputConfig(directory='/tmp/out'),
        wl_min=650.0,
        wl_max=700.0,
    )
    yaml_path = tmp_path / "cfg.yaml"
    original.to_yaml(yaml_path)
    loaded = SimulationConfig.from_yaml(yaml_path)

    assert loaded.band == original.band
    assert loaded.simulation_type == original.simulation_type
    assert loaded.source.type == original.source.type


def test_yaml_from_hand_written_file(tmp_path):
    """Loading a hand-written YAML (the common run-config use case) works."""
    yaml_path = tmp_path / "cfg.yaml"
    yaml_path.write_text("""
simulation_type: flat_field
band: R
exposure_time: 30.0
source:
  type: constant
  flux: 0.5
fibers:
  mode: single
  fibers: [21]
output:
  directory: /tmp/out
""")
    cfg = SimulationConfig.from_yaml(yaml_path)
    assert cfg.band == 'R'
    assert cfg.exposure_time == 30.0
    assert cfg.source.flux == 0.5
    assert cfg.fibers.fibers == [21]


def test_from_yaml_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        SimulationConfig.from_yaml(tmp_path / "nope.yaml")


# ---------- create_template_configs ----------

@pytest.mark.xfail(reason=(
    "create_template_configs uses mode='even_odd' and type='zemax', "
    "neither of which pass SimulationConfig.validate(). Bug found by "
    "tests; either remove the function or fix the templates."
), strict=True)
def test_create_template_configs_smoke(tmp_path):
    """Templates must be loadable by from_yaml after creation."""
    create_template_configs(tmp_path)
    for yaml_file in tmp_path.glob("*.yaml"):
        SimulationConfig.from_yaml(yaml_file)
