# ANDES E2E Simulation Suite

End-to-end simulations for the ANDES high-resolution spectrograph at ESO's ELT.

## Quick Start

### Using the Modern Framework (Recommended)

```bash
cd andes_simulator
andes-sim simulate configs/examples/flat_field_all_fibers.yaml
```

See [andes_simulator/README.md](andes_simulator/README.md) for complete documentation.

### Framework Status

⚠️ **Important**: The `andes_simulator` framework has not yet been fully validated against legacy scripts. Before using in production:
- Run comparison tests between legacy and new framework
- Verify outputs match within numerical precision
- Test all required bands and modes

See validation checklist in project cleanup plan.

## What's in This Directory

### andes_simulator/ - Modern Framework

Unified simulation framework that replaces all legacy scripts with:
- Single command-line interface with subcommands
- YAML-based configuration (no hardcoded parameters)
- Batch processing with parallel execution
- Works from any directory (portable paths)
- Complete feature coverage

**Key features:**
- Flat field simulations (all fiber patterns)
- Fabry-Perot wavelength calibration
- Stellar spectrum observations
- HDF model generation from ZEMAX
- Post-processing (PSF convolution, fiber combination)
- Thermal variation testing

### legacy/ - Archived Scripts

Original standalone scripts (15 files) moved here for reference:
- `FF_code/` - Flat field scripts (6 variants)
- `FP_code/` - Fabry-Perot scripts (2 variants)
- `YJH/` - HDF generation for NIR bands
- `PSF/` - Post-processing tools

These scripts are **deprecated** but kept for validation and backward compatibility. See [legacy/README.md](legacy/README.md) for details.

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

## Spectral Bands

ANDES covers 8 spectral bands:
- **U, B, V** - Blue arms
- **R** - Red arm
- **IZ** - Red/infrared arm
- **Y, J, H** - Near-infrared arms

## Common Tasks

### Simulate a flat field (all fibers)

**New framework:**
```bash
cd andes_simulator
andes-sim simulate configs/examples/flat_field_all_fibers.yaml
```

**Legacy:**
```bash
uv run legacy/FF_code/pyechelle_test_ANDES_ff.py
```

### Simulate Fabry-Perot calibration

**New framework:**
```bash
cd andes_simulator
andes-sim simulate configs/examples/fabry_perot_wavelength_cal.yaml
```

**Legacy:**
```bash
uv run legacy/FP_code/pyechelle_test_ANDES_fp.py
```

### Generate HDF model from ZEMAX

**New framework:**
```bash
cd andes_simulator
andes-sim build-hdf configs/examples/hdf_generation.yaml
```

**Legacy:**
```bash
uv run legacy/YJH/MakeHDF_Yband.py
```

### Post-process: Combine fibers with PSF

**New framework:**
```bash
cd andes_simulator
andes-sim postprocess --mode convolve --band R
```

**Legacy:**
```bash
uv run legacy/PSF/Dkernel.py R 5,5,3.2,left
```

## Development

### Dependencies

- Python 3.13+
- PyEchelle (optical simulation engine)
- NumPy, SciPy, Astropy
- YAML configuration support

Install via:
```bash
cd andes_simulator
uv sync
```

### Project Structure

```
src/
├── andes_simulator/           # Modern framework
│   ├── core/                  # Simulation engine
│   ├── sources/               # Source models (FF, FP, stellar)
│   ├── models/                # HDF builders, thermal models
│   ├── postprocess/           # PSF convolution, fiber combination
│   ├── cli/                   # Command-line interface
│   ├── configs/examples/      # Configuration templates
│   └── README.md              # Framework documentation
├── legacy/                    # Archived scripts
│   └── README.md              # Migration guide
├── HDF/                       # Optical models
├── SED/                       # Spectral data
├── CLAUDE.md                  # AI assistant instructions
└── README.md                  # This file
```

## Documentation

- **Framework usage**: [andes_simulator/README.md](andes_simulator/README.md)
- **Migration guide**: [andes_simulator/ORIGINAL_SCRIPT_MAPPING.md](andes_simulator/ORIGINAL_SCRIPT_MAPPING.md)
- **Legacy scripts**: [legacy/README.md](legacy/README.md)
- **AI context**: [CLAUDE.md](CLAUDE.md)

## Repository

Part of the ANDES E2E Scripts repository:
- GitHub: ivh/ANDES_E2E_scripts
- Branch: master

## Support

For questions about:
- **Framework**: See andes_simulator/README.md
- **Legacy scripts**: Check legacy/README.md or git history
- **ANDES instrument**: Contact ANDES team

---

*Last updated: 2025-12-11*
*Reorganized and cleaned up as part of codebase maintenance*
