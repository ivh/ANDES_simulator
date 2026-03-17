# ELT Spectrograph Simulator

End-to-end simulations for ELT spectrographs, based on
optical models (ZEMAX files) and [PyEchelle](https://gitlab.com/Stuermer/pyechelle).

Supported instruments:
- **ANDES** -- high-resolution echelle spectrograph (bands: U, B, V, R, IZ, Y, J, H)
- **MOSAIC** -- multi-object VPH spectrograph (LR: blue, red, J, H; HR: B1, R1, B2, H)

This is a mostly vibe-coded (ClaudeCode with sprinkles of Codex and Gemini) command line
interface to allow quick iterations on simulated frames. No guarantees.


## Quick Start

```bash
# ANDES
uv run andes-sim list-bands
uv run andes-sim simulate --band R --source flat --fiber 21
uv run andes-sim simulate --source fp --subslit ifu --wl-min 1600 --wl-max 1603
uv run andes-sim simulate --band R --source lfc --subslit cal --flux 2
uv run andes-sim simulate --source flat --band H --subslit slitA --wl-min 1600 --wl-max 1602 --fib-eff 0.5-0.95

# MOSAIC
uv run mosaic-sim list-bands
uv run mosaic-sim simulate --band LR-blue --source flat --subslit all
uv run mosaic-sim simulate --band LR-blue --source flat --fiber bundle:1
uv run mosaic-sim simulate --band LR-J --source fp --fiber bundle:1-10 --flux 100
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

Each instrument has its own CLI (`andes-sim`, `mosaic-sim`) with the same subcommands.

```bash
uv run andes-sim --help
uv run mosaic-sim --help
uv run andes-sim simulate --help
```

### ANDES Calibration Simulations

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

# Doppler velocity shift (m/s) - applied to source wavelengths via PyEchelle
uv run andes-sim simulate --band R --source lfc --fiber 21 --velocity-shift 2000
uv run andes-sim simulate --band R --source fp --fiber 21 --velocity-shift data/vel_shifts_R.json

# Constant pixel x-shift - applied via LocalDisturber to the optical model
uv run andes-sim simulate --band R --source lfc --fiber 21 --x-shift 0.5
```

### MOSAIC Simulations

```bash
# VIS LR (980 fibers = 140 bundles x 7, 12288x12288 detector)
uv run mosaic-sim simulate --band LR-blue --source flat --subslit all
uv run mosaic-sim simulate --band LR-red --source flat --fiber bundle:1

# NIR LR (630 fibers = 90 bundles x 7, 4096x4096 detector)
uv run mosaic-sim simulate --band LR-J --source flat --fiber bundle:1-10
uv run mosaic-sim simulate --band LR-H --source fp --fiber bundle:5 --flux 100

# HR modes (VIS HR: 1140 fibers = 60 bundles x 19; NIR HR: 630 fibers)
uv run mosaic-sim simulate --band HR-H --source flat --fiber bundle:1

# Bundle selection: bundle:N or bundle:N-M
uv run mosaic-sim simulate --band LR-blue --source fp --fiber bundle:50-60 --flux 100
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
- ANDES NIR arms (Y, J, H): 75-fiber echelle models
- ANDES Optical arms (R, IZ): 66-fiber echelle models
- ANDES Thermal variant models (T0019, T0028, T0045, T0108)
- MOSAIC VIS LR: 980-fiber VPH models (Blue 390-625nm, Red 595-952nm)
- MOSAIC VIS HR: 1140-fiber VPH models (B1 505-580nm, R1 610-680nm, B2 393-458nm)
- MOSAIC NIR LR: 630-fiber VPH models (J 950-1340nm, H 1430-1800nm)
- MOSAIC NIR HR: 630-fiber VPH model (H 1520-1620nm)

### SED/ - Spectral Data

Fabry-Perot spectral energy distributions (.csv):
- NIR FP spectrum (finesse 26)
- RIZ FP spectrum (finesse 23)
- UBV FP spectrum (finesse 26)

---

*Last updated: 2026-03-17*
