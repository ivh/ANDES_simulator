# E2E Simulation Instructions

## Multi-instrument framework

The `andes_simulator` package supports multiple ELT spectrographs via separate CLIs
that share the same PyEchelle-based simulation core:
- `andes-sim` -- ANDES high-resolution echelle (bands: U, B, V, R, IZ, Y, J, H;
  plus `Y_iq15`, `J_iq15`, `H_iq15`: single-fiber "master_rotfix_iq15" model variants,
  use `--fiber 1`; H_iq15 has gaps in order coverage)
- `mosaic-sim` -- MOSAIC multi-object VPH spectrograph (bands: B_LR, R_LR, J_LR, H_LR, B1_HR, R1_HR, B2_HR, H_HR)

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

# HCL (ThAr hollow-cathode lamp; NIST line list in SED/linelists/)
uv run andes-sim simulate --band R --source hcl --subslit slitA

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
uv run mosaic-sim simulate --band B_LR --source flat --subslit all
uv run mosaic-sim simulate --band B_LR --source flat --fiber bundle:1

# Fabry-Perot
uv run mosaic-sim simulate --band B_LR --source fp --fiber bundle:1 --flux 100

# Bundle selection (7 fibers/bundle for LR and NIR, 19 fibers/bundle for VIS HR)
uv run mosaic-sim simulate --band R_LR --source flat --fiber bundle:5
uv run mosaic-sim simulate --band J_LR --source flat --fiber bundle:1-10

# NIR HR
uv run mosaic-sim simulate --band H_HR --source flat --fiber bundle:1
```

### EDPS Raw-Frame Commands (ANDES only)

Design record: `src/dpr_summary.md`. One MEF per spectrograph arm (UBV/RIZ/YJH),
one uint16 extension per detector band, ESO classification headers
(DRL spec v1.2 Ch. 4.1). Code in `andes_simulator/raw/`.

```bash
# single raw frame from a DPR.TYPE string (grammar: KIND,A,B per Templates
# Manual v2.0; the calibration fibre C goes into --calfib / ins.calfib)
uv run andes-sim make-raw --arm RIZ --dpr "WAVE,HCL,FP" --calfib OFF --exptime 120 --nexp 2 -o rawdata/
uv run andes-sim make-raw --arm YJH --dpr BIAS --exptime 0 --nexp 10 -o rawdata/
uv run andes-sim make-raw --arm YJH --dpr "FLAT,LAMP" --mode IFU-AO --ifu-scale 16 -o rawdata/
uv run andes-sim make-raw --arm RIZ --dpr "SLITMASK,FP,OFF" --ins-mask M1 --calfib FP -o rawdata/

# only simulate some detectors of the arm (other extensions: detector noise, flagged)
uv run andes-sim make-raw --arm YJH --bands Y --dpr "WAVE,FP,OFF" --calfib FP -o rawdata/

# synthetic night from the canonical plan (versioned in the edps repo)
uv run andes-sim night ~/ANDES/edps/calibration_plan.yaml --sets detector,daily --arms RIZ -o night1/
uv run andes-sim night ~/ANDES/edps/calibration_plan.yaml --include night,science --dry-run

# headers-only night (stub 2x2 extensions; for EDPS classification tests)
uv run andes-sim night ~/ANDES/edps/calibration_plan.yaml --arms RIZ --include night,science --headers-only -o night_hdr/
```

Useful flags: `--jobs N` (parallel slot simulations), `--seed` (reproducible
noise), `--boost` (expectation-cache flux boost, default 10), `--cache-dir`
(default `E2E/simcache/`), `--headers-only`, `--dry-run`.

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

### Raw-frame generation (andes_simulator/raw/)

- **DPR grammar** (`raw/dpr.py`): `<KIND>,<A>,<B>` (SL), `<KIND>,<slit>` (IFU),
  no KIND for science (`OBJECT,SKY`, `OBJECT,WAVE` in TC mode). The calibration
  fibre C is not in DPR.TYPE (Templates Manual v2.0 / ESO-044156): its source
  is the `--calfib` parameter, stamped as `ins.calfib`. BIAS/DARK/LED
  `FLAT,LAMP` are detector-only. Source tokens map to simulator sources
  (LAMP->flat, HCL->hcl, WAVE->fp, SKY/OBJECT/...-> CSVs from SED/ chosen by
  band coverage). Masks M1-M3 = every third fiber.
- **Cache** (`raw/cache.py`, `E2E/simcache/`): stores *boosted expectation*
  images (simulated at boost x flux, divided by boost) so each exposure draws
  fresh Poisson noise. Residual correlated noise is 1/boost of shot variance —
  fine for recipe testing, raise --boost for noise studies. Fiber efficiencies
  are seeded per band (static instrument property, consistent across frames).
- **Detector model** (`raw/detector.py`): PRNU, dark+hot pixels, dead pixels/bad
  columns (all static, seeded per band), Poisson, charge binning, RON, gain,
  bias, uint16 saturation. Parameter values in `core/andes.py` DETECTOR_MODELS
  are placeholders until real detector specs exist. No overscan regions yet
  (geometry undefined in the ADs).
- **Flux levels**: no realistic throughput model — slot images are normalized
  to target peak e- per (KIND, token) in `raw/dpr.py` PEAK_TARGETS_E; header
  EXPTIME is scheduling metadata, not a photon integral.
- **Known gap**: AD2-only patterns (one-aperture wave frames like
  `WAVE,FP,OFF`, the combined flat `FLAT,LAMP,LAMP`) generate fine but are
  unclassified by the edps workflow (reconciliation items 1-2 in
  calibration_plan.yaml, which also carries the proposed OrderDef and
  per-slit flat procedures the cascade needs).

- **PyEchelle**: Uses v0.4.0; must use `max_cpu=1` due to multi-CPU bug
- **Numba cache**: `andes_simulator/__init__.py` provisions a per-process tmpdir as `NUMBA_CACHE_DIR` (with atexit cleanup) to avoid "underlying object has vanished" errors from cache corruption. Single-process runs need no setup. An externally-set `NUMBA_CACHE_DIR` is respected, so parallel wrapper scripts (e.g. `scripts/lfc_allfib_allbands.sh`, `scripts/ifu_star.py`, `scripts/R_starsky.py`) that give each worker subprocess its own cache keep working.
- **CSV sources**: PyEchelle's raytracing passes bare micron floats to `get_counts`; the simulator converts CSV wavelengths to microns automatically. CSV files can include a `# scaling: VALUE` header comment for default flux scaling.
- **Sources**: Each fiber needs individual source object (no shared references)
- **Array shapes**: Config uses (X,Y), FITS/numpy uses (Y,X)
- **LFC**: Lines equidistant in velocity (~33 km/s for R-band), ~100-150 lines/order
- **psf-process caveat**: `--fwhm` is labeled "arcseconds" but the conversion in `postprocess/psf.py` multiplies by the band's `sampling` config value (px per resolution element), so it effectively means FWHM in resolution elements. Also the real sampling varies significantly across the detector, so a single uniform Gaussian kernel is not a realistic PSF model. Revisit before relying on psf-process results.
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
