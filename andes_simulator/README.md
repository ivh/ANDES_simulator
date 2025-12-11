# ANDES Simulation Framework

A unified framework for running ANDES E2E spectrograph simulations including flat field calibrations, Fabry-Perot wavelength calibrations, stellar observations, and post-processing.

## Features

- **Unified Interface**: Single entry point for all simulation types
- **Configuration Management**: YAML-based configuration with validation  
- **Batch Processing**: Parallel execution of multiple simulations
- **Portable Paths**: Works from any directory without hardcoded paths
- **Post-Processing**: Built-in PSF convolution and fiber combination tools
- **Extensible**: Easy to add new simulation types and features

## Installation

### Using uv (Recommended)

```bash
# Install with basic dependencies
uv sync

# Install with all optional dependencies
uv sync --extra all

# Install for development
uv sync --extra development
```

### Using pip

```bash
pip install -e .
pip install -e ".[all]"  # With all optional dependencies
```

## Quick Start

### Command Line Interface

```bash
# List available spectral bands
andes-sim list-bands

# Generate flat field calibration
andes-sim flat-field --band Y --mode single --fiber 1

# Generate Fabry-Perot calibration  
andes-sim fabry-perot --band Y --mode all

# Run from configuration file
andes-sim run-config --config examples/flat_field_single_fiber.yaml

# Create template configurations
andes-sim create-templates

# Process with PSF convolution
andes-sim psf-process --band Y --input-pattern "Y_FP_fiber{fib:02d}_*.fits" --fwhm 3.2

# Combine fiber outputs
andes-sim combine --band Y --input-pattern "Y_FP_fiber{fib:02d}_*.fits" --mode all
```

### Python API

```python
from andes_simulator import AndesSimulator, SimulationConfig

# Quick flat field simulation
sim = AndesSimulator.quick_flat_field('Y', 'single', fiber=1)
result = sim.run_simulation()

# From configuration file
config = SimulationConfig.from_yaml('my_config.yaml')
sim = AndesSimulator(config)
result = sim.run_simulation()

# Batch processing
from andes_simulator.scripts.batch_runner import BatchRunner
runner = BatchRunner()
results = runner.quick_fiber_sweep('Y', 'fabry_perot')
```

## Simulation Types

### Flat Field Calibrations
- Single fiber illumination
- Even/odd fiber patterns  
- Pseudo-slit illumination
- Calibration fiber patterns

### Fabry-Perot Wavelength Calibrations
- All fiber illumination
- Single fiber with velocity shifts
- Thermal variation studies

### Stellar Spectrum Observations
- Custom CSV spectrum inputs
- Single fiber observations
- Flux scaling control

### HDF Model Generation
- ZEMAX integration
- Automated fiber field setup
- PSF and transformation sampling

### Post-Processing
- PSF convolution with edge-blanking
- Fiber combination and summation
- Analysis and reporting tools

## Configuration Files

The framework uses YAML configuration files for reproducible simulations:

```yaml
simulation_type: fabry_perot
band: Y
exposure_time: 30.0

source:
  type: fabry_perot
  scaling_factor: 5e9

fibers:
  mode: all
  fibers: all

output:
  directory: "../{band}/"
  filename_template: "{band}_FP_{exposure}s.fits"
```

See `configs/examples/` for template configurations covering all simulation types.

## Migrating from Original Scripts

The framework replaces all original E2E scripts with a unified interface:

| Original Script | New Command |
|----------------|-------------|
| `FF_code/pyechelle_test_ANDES_ff_single_fiber.py` | `andes-sim flat-field --mode single` |
| `FP_code/pyechelle_test_ANDES_fp.py` | `andes-sim fabry-perot --mode all` |
| `PSF/Dkernel.py` | `andes-sim psf-process` |
| `PSF/sumIFU.py` | `andes-sim combine` |

See `ORIGINAL_SCRIPT_MAPPING.md` for complete migration guide.

## Development

### Running Tests

```bash
uv run pytest
uv run pytest --cov=andes_simulator
```

### Code Formatting

```bash
uv run black andes_simulator/
uv run flake8 andes_simulator/
uv run mypy andes_simulator/
```

### Dry Run Testing

Test configurations without running full simulations:

```bash
andes-sim flat-field --band Y --mode single --fiber 1 --dry-run
andes-sim run-config --config my_config.yaml --dry-run
```

## Requirements

- Python ≥ 3.9
- NumPy ≥ 1.20
- SciPy ≥ 1.7  
- Astropy ≥ 5.0
- PyYAML ≥ 6.0
- Click ≥ 8.0

### Optional Dependencies
- pyechelle ≥ 1.0 (for simulations)
- matplotlib ≥ 3.5 (for visualizations)

## License

MIT License - see LICENSE file for details.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make changes with tests
4. Run formatting and tests
5. Submit a pull request

## Support

- GitHub Issues: https://github.com/ivh/ANDES_E2E_scripts/issues
- Documentation: See `ORIGINAL_SCRIPT_MAPPING.md` for migration guide