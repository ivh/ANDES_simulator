# ANDES E2E Simulation Instructions

## Framework: andes_simulator

**Status**: Production ready (validated for R-band; other bands expected to work similarly)

### Simulation Commands

```bash
# Flat field calibrations
uv run andes-sim flat-field --band R --mode all --output-dir ../R/
uv run andes-sim flat-field --band R --mode single --fiber 21 --output-dir ../R/
uv run andes-sim flat-field --band R --mode even_odd --output-dir ../R/

# Fabry-Perot wavelength calibration
uv run andes-sim fabry-perot --band R --mode single --fiber 21 --flux 100 --output-dir ../R/

# LFC (Laser Frequency Comb) calibration
uv run andes-sim lfc --band R --mode single --fiber 21 --scaling 1e5 --output-dir ../R/

# Stellar spectrum
uv run andes-sim spectrum --band R --spectrum SED/star.csv --fiber 21 --output-dir ../R/
```

### Post-Processing Commands

```bash
# Combine individual fiber outputs
uv run andes-sim combine --band R --input-pattern "R_FP_fiber{fib:02d}_30s.fits" \
    --mode all --output-dir /absolute/path/to/R/ --output R_combined.fits

# PSF convolution
uv run andes-sim psf-process --band R --input-pattern "R_FP_fiber{fib:02d}_*.fits" \
    --fwhm 3.2 --output-dir ../R/
```

### Key Options

- `--mode`: `all`, `single`, `even_odd`, `first_slit`, `second_slit`
- `--output-dir`: Use absolute paths for post-processing tools
- `--dry-run`: Preview without executing

## Technical Notes

- **PyEchelle**: Uses v0.4.0; must use `max_cpu=1` due to multi-CPU bug
- **Sources**: Each fiber needs individual source object (no shared references)
- **Array shapes**: Config uses (X,Y), FITS/numpy uses (Y,X)
- **LFC**: Lines equidistant in velocity (~33 km/s for R-band), ~100-150 lines/order

## Directory Structure

```
src/
├── andes_simulator/    # Main framework
│   ├── cli/            # Command-line interface
│   ├── core/           # Simulator, config, sources
│   ├── sources/        # FF, FP, LFC, stellar sources
│   └── postprocess/    # combine, psf tools
├── legacy/             # Archived scripts (reference only)
├── HDF/                # ZEMAX optical models (.hdf)
└── SED/                # Spectral data (.csv)
```

## Validation Status

See `VALIDATION_REPORT.md` and `CLEANUP_PLAN.md` for detailed test results.
