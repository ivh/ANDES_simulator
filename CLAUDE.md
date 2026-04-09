# E2E Simulation Instructions

## Multi-instrument framework

The `andes_simulator` package supports multiple ELT spectrographs via separate CLIs
that share the same PyEchelle-based simulation core:
- `andes-sim` -- ANDES high-resolution echelle (bands: U, B, V, R, IZ, Y, J, H)
- `mosaic-sim` -- MOSAIC multi-object VPH spectrograph (bands: LR-blue, LR-red, LR-J, LR-H, HR-B1, HR-R1, HR-B2, HR-H)

Instrument-specific configs live in `core/andes.py` and `core/mosaic.py`.
The registry in `core/instruments.py` merges them so downstream code is instrument-agnostic.

**Status**: Production ready for ANDES (validated for R-band). MOSAIC basic support.

### ANDES Simulation Commands

```bash
# Flat field
uv run andes-sim simulate --band R --source flat --subslit all
uv run andes-sim simulate --band R --source flat --fiber 21

# Fabry-Perot
uv run andes-sim simulate --band R --source fp --fiber 21 --flux 100

# LFC
uv run andes-sim simulate --band R --source lfc --subslit cal_sl

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

### MOSAIC Simulation Commands

```bash
# Flat field
uv run mosaic-sim simulate --band LR-blue --source flat --subslit all
uv run mosaic-sim simulate --band LR-blue --source flat --fiber bundle:1

# Fabry-Perot
uv run mosaic-sim simulate --band LR-blue --source fp --fiber bundle:1 --flux 100

# Bundle selection (7 fibers/bundle for LR and NIR, 19 fibers/bundle for VIS HR)
uv run mosaic-sim simulate --band LR-red --source flat --fiber bundle:5
uv run mosaic-sim simulate --band LR-J --source flat --fiber bundle:1-10

# NIR HR
uv run mosaic-sim simulate --band HR-H --source flat --fiber bundle:1
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
  - ANDES all bands: `all`, `even`, `odd`, `slitA`, `slitB`, `cal_sl`
  - ANDES YJH only: `cal_ifu`
  - ANDES YJH only: `ifu`, `ring0`, `ring1`, `ring2`, `ring3`, `ring4`
  - MOSAIC: `all`, `even`, `odd`
  - MOSAIC bundles: `bundle:N`, `bundle:N-M` (e.g. `bundle:5`, `bundle:1-10`)
- `--mode`: Combination mode for post-processing (`all`, `even_odd`, `slits`, `custom`)
- `--velocity-shift`: Doppler shift in m/s (scalar or JSON file with per-fiber values)
- `--x-shift`: Constant pixel shift (scalar or JSON file), applied via LocalDisturber `d_tx`
- `--output-dir`: Use absolute paths for post-processing tools
- `--dry-run`: Preview without executing

## Technical Notes

- **PyEchelle**: Uses v0.4.0; must use `max_cpu=1` due to multi-CPU bug
- **Numba cache**: `andes_simulator/__init__.py` provisions a per-process tmpdir as `NUMBA_CACHE_DIR` (with atexit cleanup) to avoid "underlying object has vanished" errors from cache corruption. Single-process runs need no setup. An externally-set `NUMBA_CACHE_DIR` is respected, so parallel wrapper scripts (e.g. `scripts/lfc_allfib_allbands.sh`, `scripts/ifu_star.py`, `scripts/R_starsky.py`) that give each worker subprocess its own cache keep working.
- **CSV sources**: PyEchelle's raytracing passes bare micron floats to `get_counts`; the simulator converts CSV wavelengths to microns automatically. CSV files can include a `# scaling: VALUE` header comment for default flux scaling.
- **Sources**: Each fiber needs individual source object (no shared references)
- **Array shapes**: Config uses (X,Y), FITS/numpy uses (Y,X)
- **LFC**: Lines equidistant in velocity (~33 km/s for R-band), ~100-150 lines/order
- **Velocity shift vs x-shift**: `--velocity-shift` uses PyEchelle's `set_radial_velocities` (Doppler, shifts source wavelengths). `--x-shift` uses `LocalDisturber(d_tx=...)` (constant pixel offset in all orders). Due to echelle optics (constant tx_span across orders), both produce nearly uniform pixel shifts (~3% variation across orders for Doppler).

## Directory Structure

```
src/
├── andes_simulator/    # Main package
│   ├── cli/            # CLIs: andes.py, mosaic.py (entry points), main.py (factory)
│   ├── core/           # andes.py, mosaic.py (instrument configs), instruments.py (registry)
│   ├── sources/        # Source spectrum generators
│   └── postprocess/    # combine, psf tools
├── HDF/                # ZEMAX optical models (.hdf) for all instruments
└── SED/                # Spectral data (.csv)
```
