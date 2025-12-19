# ANDES Simulator

End-to-end simulations for the ANDES high-resolution spectrograph at ESO's ELT, based on
the optical model (ZEMAX files) and [PyEchelle](https://gitlab.com/Stuermer/pyechelle).

This is a mostly vibe-coded (ClaudeCode with sprinkles of Codex and Gemini) command line
interface to allow quick iterations on simulated frames. No guarantees.


## Quick Start

```bash
# List available spectral bands
uv run andes-sim list-bands

# Generate flat field calibration (R-band, single fiber)
uv run andes-sim simulate --band R --source flat --fiber 21

# Generate Fabry-Perot in the IFU for a small wavelength range in H-band
uv run andes-sim simulate --source fp --subslit ifu --wl-min 1600 --wl-max 1603

# Generate LFC (Laser Frequency Comb) in the calibration fibers, double default flux
uv run andes-sim simulate --band R --source lfc --subslit cal --flux 2

# Varying fiber efficiency in SL slit A
uv run andes-sim simulate --source flat --band H --subslit slitA --wl-min 1600 --wl-max 1602 --fib-eff 0.5-0.95
```

## Installation

Clone this repository and `cd` to it.

### Using uv (Recommended)

Install [uv](https://docs.astral.sh/uv/), if you have not already.

```bash
# Install with basic dependencies
uv sync
```

### Using pip

```bash
pip install -e .
```

## Command Line Interface

### List options

```bash
# List available commands
uv run andes-sim --help

# List options for simulation
uv run andes-sim simulate --help
```

### Calibration Simulations

```bash
# Flat field - various subslits
uv run andes-sim simulate --band R --source flat --subslit all --output-dir ../R/
uv run andes-sim simulate --band R --source flat --fiber 21
uv run andes-sim simulate --band R --source flat --subslit even
uv run andes-sim simulate --band R --source flat --subslit slitA
uv run andes-sim simulate --band R --source flat --subslit cal
uv run andes-sim simulate --band Y --source flat --subslit ifu
uv run andes-sim simulate --band H --source flat --subslit ring1

# Fabry-Perot wavelength calibration
uv run andes-sim simulate --band R --source fp --fiber 21 --flux 100 --output-dir ../R/

# LFC (Laser Frequency Comb) calibration
uv run andes-sim simulate --band R --source lfc --fiber 21 --scaling 1e5 --output-dir ../R/

# Stellar spectrum observation
uv run andes-sim simulate --band R --source SED/star.csv --fiber 21 --output-dir ../R/
```

### Post-Processing

Existing simulations can be convolved with kernels and added together.

```bash
# Combine individual fiber outputs into single frame
uv run andes-sim combine --band R --input-pattern "R_FP_fiber{fib:02d}_30s.fits" \
    --subslit all --output-dir /path/to/R/ --output R_combined.fits

# PSF convolution
uv run andes-sim psf-process --band R --input-pattern "R_FP_fiber{fib:02d}_*.fits" \
    --fwhm 3.2 --output-dir ../R/
```

### Other useful commands

```bash
# List available bands
uv run andes-sim list-bands

# Dry run (show what would be done)
uv run andes-sim simulate --band R --source flat --fiber 21 --dry-run

# Use custom HDF model (band inferred from wavelength content)
uv run andes-sim simulate --source flat --hdf HDF/ANDES_full_F18A33_win_jmr_MC_T0019_Rband_p0.hdf --dry-run

# Run from YAML configuration file
uv run andes-sim run-config configs/examples/flat_field_single_fiber.yaml
```

## Data Files

These are not included in the repository. Ask us, if you do not have them already.

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

---

*Last updated: 2025-12-15*
