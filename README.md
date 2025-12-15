# ANDES Simulator

End-to-end simulation framework for the ANDES high-resolution spectrograph at ESO's ELT.

A unified framework for running ANDES E2E spectrograph simulations including flat field calibrations, Fabry-Perot wavelength calibrations, stellar observations, and post-processing.

## Quick Start

```bash
# List available spectral bands
uv run andes-sim list-bands

# Generate flat field calibration (R-band, single fiber)
uv run andes-sim flat-field --band R --mode single --fiber 21 --output-dir ../R/

# Generate Fabry-Perot wavelength calibration
uv run andes-sim fabry-perot --band R --mode single --fiber 21 --output-dir ../R/

# Generate LFC (Laser Frequency Comb) calibration
uv run andes-sim lfc --band R --mode single --fiber 21 --output-dir ../R/
```

## Features

- **Unified Interface**: Single CLI entry point for all simulation types
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

## Command Line Interface

### Calibration Simulations

```bash
# Flat field - various modes
uv run andes-sim flat-field --band R --mode all --output-dir ../R/
uv run andes-sim flat-field --band R --mode single --fiber 21 --output-dir ../R/
uv run andes-sim flat-field --band R --mode even_odd --output-dir ../R/

# Fabry-Perot wavelength calibration
uv run andes-sim fabry-perot --band R --mode single --fiber 21 --flux 100 --output-dir ../R/

# LFC (Laser Frequency Comb) calibration
uv run andes-sim lfc --band R --mode single --fiber 21 --scaling 1e5 --output-dir ../R/

# Stellar spectrum observation
uv run andes-sim spectrum --band R --spectrum SED/star.csv --fiber 21 --output-dir ../R/
```

### Post-Processing

```bash
# Combine individual fiber outputs into single frame
uv run andes-sim combine --band R --input-pattern "R_FP_fiber{fib:02d}_30s.fits" \
    --mode all --output-dir /path/to/R/ --output R_combined.fits

# PSF convolution
uv run andes-sim psf-process --band R --input-pattern "R_FP_fiber{fib:02d}_*.fits" \
    --fwhm 3.2 --output-dir ../R/
```

### Utilities

```bash
# List available bands
uv run andes-sim list-bands

# Dry run (show what would be done)
uv run andes-sim flat-field --band R --mode single --fiber 21 --dry-run

# Use custom HDF model (band inferred from wavelength content)
uv run andes-sim flat-field --hdf HDF/ANDES_full_F18A33_win_jmr_MC_T0019_Rband_p0.hdf --dry-run

# Run from YAML configuration file
uv run andes-sim run-config configs/examples/flat_field_single_fiber.yaml
```

## Python API

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

### LFC (Laser Frequency Comb) Calibrations
- Unresolved emission lines equidistant in velocity
- ~100-150 lines per spectral order
- Configurable flux per line

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

See `andes_simulator/configs/examples/` for template configurations covering all simulation types.

## Spectral Bands

ANDES covers 8 spectral bands:
- **U, B, V** - Blue arms
- **R** - Red arm
- **IZ** - Red/infrared arm
- **Y, J, H** - Near-infrared arms

## Project Structure

```
ANDES_simulator/
├── andes_simulator/           # Modern framework
│   ├── core/                  # Simulation engine
│   ├── sources/               # Source models (FF, FP, stellar)
│   ├── models/                # HDF builders, thermal models
│   ├── postprocess/           # PSF convolution, fiber combination
│   ├── cli/                   # Command-line interface
│   ├── configs/examples/      # Configuration templates
│   └── README.md              # Framework documentation
├── legacy/                    # Archived scripts (deprecated)
│   └── README.md              # Migration guide
├── HDF/                       # Optical models (ZEMAX .hdf files)
├── SED/                       # Spectral data (.csv files)
├── CLAUDE.md                  # AI assistant instructions
└── README.md                  # This file
```

## Data Files

### HDF/ - Optical Models

ZEMAX optical model files (.hdf) containing:
- PSF transformations
- Wavelength-to-pixel mappings
- Fiber geometries
- Diffraction order information

**Files include:**
- NIR arms (Y, J, H): 75-fiber models
- Optical arms (R, IZ): 66-fiber models
- Thermal variant models (T0019, T0028, T0045, T0108)

### SED/ - Spectral Data

Fabry-Perot spectral energy distributions (.csv):
- NIR FP spectrum (finesse 26)
- RIZ FP spectrum (finesse 23)
- UBV FP spectrum (finesse 26)

## Migrating from Legacy Scripts

The framework replaces all original E2E scripts with a unified interface:

| Legacy Script | New Command |
|---------------|-------------|
| `FF_code/pyechelle_test_ANDES_ff_single_fiber.py` | `andes-sim flat-field --mode single` |
| `FP_code/pyechelle_test_ANDES_fp.py` | `andes-sim fabry-perot --mode all` |
| (new) | `andes-sim lfc` |
| `PSF/Dkernel.py` | `andes-sim psf-process` |
| `PSF/sumIFU.py` | `andes-sim combine` |

See `andes_simulator/ORIGINAL_SCRIPT_MAPPING.md` for complete migration guide.

Legacy scripts are archived in `legacy/` directory for reference and validation purposes only.

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
andes-sim simulate configs/examples/my_config.yaml --dry-run
```

## Requirements

- Python ≥ 3.13
- NumPy ≥ 1.20
- SciPy ≥ 1.7
- Astropy ≥ 5.0
- PyYAML ≥ 6.0
- Click ≥ 8.0

### Optional Dependencies
- pyechelle ≥ 0.4.0 (for simulations)
- matplotlib ≥ 3.5 (for visualizations)

## Repository

- GitHub: https://github.com/ivh/ANDES_simulator
- Branch: master

## Support

For questions about:
- **Framework usage**: See andes_simulator/README.md
- **Legacy scripts**: Check legacy/README.md or git history
- **ANDES instrument**: Contact ANDES team
- **Issues**: https://github.com/ivh/ANDES_simulator/issues

## License

MIT License - see LICENSE file for details.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make changes with tests
4. Run formatting and tests
5. Submit a pull request

---

*Last updated: 2025-12-15*
