# ANDES E2E Simulation Instructions

## Framework: andes_simulator

**Status**: Production ready (validated for R-band; other bands expected to work similarly)

### Simulation Commands

```bash
# Flat field
uv run andes-sim simulate --band R --source flat --subslit all
uv run andes-sim simulate --band R --source flat --fiber 21

# Fabry-Perot
uv run andes-sim simulate --band R --source fp --fiber 21 --flux 100

# LFC
uv run andes-sim simulate --band R --source lfc --subslit cal

# YJH IFU
uv run andes-sim simulate --band Y --source flat --subslit ifu

# Stellar spectrum (CSV path auto-detected as source type)
uv run andes-sim simulate --band R --source SED/star.csv --fiber 21

# Doppler velocity shift (m/s, applied via PyEchelle set_radial_velocities)
uv run andes-sim simulate --band R --source lfc --fiber 21 --velocity-shift 2000
uv run andes-sim simulate --band R --source fp --fiber 21 --velocity-shift data/vel_shifts_R.json

# Pixel x-shift (constant pixel offset via LocalDisturber)
uv run andes-sim simulate --band R --source lfc --fiber 21 --x-shift 0.5
```

### Post-Processing Commands

```bash
# Combine fiber outputs
uv run andes-sim combine --band R --input-pattern "R_FP_fiber{fib:02d}_*.fits" --mode all

# PSF convolution
uv run andes-sim psf-process --band R --input-pattern "R_FP_fiber{fib:02d}_*.fits" --fwhm 3.2
```

### Key Options

- `--source`: Source type (`flat`, `fp`, `lfc`) or path to CSV spectrum file
- `--hdf`: Custom HDF model file (infers band from wavelength content)
- `--subslit`: Fiber selection for simulations
  - All bands: `all`, `single`, `even`, `odd`, `slitA`, `slitB`, `cal`
  - YJH only: `ifu`, `ring0`, `ring1`, `ring2`, `ring3`, `ring4`
- `--mode`: Combination mode for post-processing (`all`, `even_odd`, `slits`, `custom`)
- `--velocity-shift`: Doppler shift in m/s (scalar or JSON file with per-fiber values)
- `--x-shift`: Constant pixel shift (scalar or JSON file), applied via LocalDisturber `d_tx`
- `--output-dir`: Use absolute paths for post-processing tools
- `--dry-run`: Preview without executing

## Technical Notes

- **PyEchelle**: Uses v0.4.0; must use `max_cpu=1` due to multi-CPU bug
- **Sources**: Each fiber needs individual source object (no shared references)
- **Array shapes**: Config uses (X,Y), FITS/numpy uses (Y,X)
- **LFC**: Lines equidistant in velocity (~33 km/s for R-band), ~100-150 lines/order
- **Velocity shift vs x-shift**: `--velocity-shift` uses PyEchelle's `set_radial_velocities` (Doppler, shifts source wavelengths). `--x-shift` uses `LocalDisturber(d_tx=...)` (constant pixel offset in all orders). Due to echelle optics (constant tx_span across orders), both produce nearly uniform pixel shifts (~3% variation across orders for Doppler).

## Directory Structure

```
src/
├── andes_simulator/    # Main package
│   ├── cli/            # Command-line interface
│   ├── core/           # Simulator, config, sources
│   ├── sources/        # Source spectrum generators
│   └── postprocess/    # combine, psf tools
├── HDF/                # ZEMAX optical models (.hdf)
└── SED/                # Spectral data (.csv)
```
