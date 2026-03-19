# MOSAIC Zemax → PyEchelle HDF Pipeline

## Overview

`make_hdf.py` connects to Zemax OpticStudio via ZOS-API (through ZOSPy), loops
over fibers and orders, samples transformations + Huygens PSFs at ~15 wavelength
points per order, and writes everything to HDF5. PyEchelle then uses these HDF
files to simulate detector images.

Standalone script with PEP 723 headers — `uv run` handles all dependencies.

## Requirements

- Windows with Zemax OpticStudio (license must support ZOS-API)
- Python 3.10 (for pythonnet/ZOSPy compatibility)
- uv

## Automatic fixes in make_hdf.py

Post-generation fixes applied automatically to every HDF:

- **`fix_dispersion_axis()`**: detects if dispersion is along ty instead of tx
  (happens when the Zemax model's spatial axis differs from PyEchelle's
  convention, e.g. NIR model). Swaps tx↔ty and sx↔sy if needed.
- **`unwrap_rotation()`**: fixes rotation angle wrapping (±pi jumps). PyEchelle
  stores rotation near ±pi that can wrap; linear interpolation then passes
  through 0, producing wrong trace positions. Shifts negative values by +2pi.

Other workarounds:

- **`--skip-ccd-check`**: bypasses pyechelle's detector aperture assertion
- **`--wl-range MIN MAX`**: overrides walk-to-detector-edge algorithm
  (needed because `SingleRayNormUnpol` returns error=2 on OpticStudio 17.09)
- **`--config N`**: switches Zemax multi-configuration before tracing
- **`--groove-density N`**: overrides grating groove density (l/mm)
- **`ensure_detector_aperture()`**: auto-adds rectangular aperture to detector
- **`patch_wavelength_range()`**: monkey-patches wavelength range with fixed values
- **Active grating auto-detection**: when multiple surfaces match `--grating-surface`
  name, picks the one with nonzero groove density (MCE sets inactive gratings to 0)
- **Slit axis auto-detection**: checks native field coordinates to determine if the
  slit runs along X or Y, and places fibers accordingly

## All spectrograph modes

### VIS (`MOSAIC_VIS.zmx`, 44 surfaces, 6 multi-configs)

| Config | Name  | Mode | Wavelength (µm) | Grating (l/mm) | Blaze | Status |
|--------|-------|------|-----------------|-----------------|-------|--------|
| 1      | LR-B  | LR   | 0.390–0.625     | 1335            | 19.95°| Done (1.9 GB) |
| 2      | LR-R  | LR   | 0.595–0.952     | 876             | 19.95°| Done (1.9 GB) |
| 3      | HR-B1 | HR   | 0.505–0.580     | 2703            | 47.5° | Done (2.1 GB) |
| 4      | HR-R1 | HR   | 0.610–0.680     | 2278            | 47.5° | Done (2.1 GB) |
| 5      | HR-B2 | HR   | 0.393–0.458     | 3442            | 47.5° | Done (2.1 GB) |
| 6      | HR-R2 | HR   | 0.765–0.890     | 1770            | 47.5° | Done (2.1 GB) |

VIS LR: 140 bundles × 7 fibers = 980, fiber pitch 177 µm, bundle gap 177 µm,
detector 12288×12288 (2×2 mosaic 6K×6K), 15 µm/px. ~4 hours per HDF.

VIS HR: 60 bundles × 19 fibers = 1140, fiber pitch 152 µm, bundle gap 152 µm,
detector 12288×12288, 15 µm/px. ~4 hours per HDF.

### NIR (`Mosaic_2Cam.ZMX`, 71 surfaces, 3 multi-configs)

| Config | Name  | Mode | Wavelength (µm) | VPH (lpm) | CWL (µm) | Status |
|--------|-------|------|-----------------|-----------|-----------|--------|
| 1      | LR-J  | LR   | 0.95–1.34       | 514       | 1.145     | Done (1.2 GB) |
| 2      | LR-H  | LR   | 1.43–1.80       | 451       | 1.615     | Done (1.2 GB) |
| 3      | HR    | HR   | 1.52–1.62       | 1003      | 1.570     | Done (1.2 GB) |

NIR: 90 bundles × 7 fibers = 630, fiber core 150 µm, pitch 203 µm,
bundle gap 203 µm, detector 4096×4096, 15 µm/px. ~20 min per HDF.

**Requires** `MOSAIC130.agf` + `.Bgf` glass catalogs installed in
`Documents/Zemax/GLASSCAT/` (bundled in the NIR design directory).

Two identical spectrographs (NIR-SP1, NIR-SP2), 90 bundles each = 180 total.
The Zemax model covers one spectrograph.

## VIS commands

```bash
# LR Blue (Config 1) — ~4 hours, ~1.9 GB
uv run --python 3.10 make_hdf.py \
    "opticaldesign/MOSAIC-VIS-Optical_Design/MOSAIC_VIS.zmx" \
    MOSAIC_VIS_LR_Blue.hdf \
    --orders -1 0 \
    --grating-surface "VPHG" --blaze 19.95 \
    --name MOSAIC-VIS-LR-Blue \
    --skip-ccd-check --wl-range 0.39 0.625 \
    --layout bundles --nbundles 140 --fibers-per-bundle 7 \
    --fiber-size 177 --bundle-gap 177 \
    --nx 12288 --ny 12288 --pixelsize 15

# LR Red (Config 2) — ~4 hours, ~1.9 GB
uv run --python 3.10 make_hdf.py \
    "opticaldesign/MOSAIC-VIS-Optical_Design/MOSAIC_VIS.zmx" \
    MOSAIC_VIS_LR_Red.hdf \
    --orders -1 0 \
    --grating-surface "VPHG" --blaze 19.95 --config 2 \
    --name MOSAIC-VIS-LR-Red \
    --skip-ccd-check --wl-range 0.595 0.952 \
    --layout bundles --nbundles 140 --fibers-per-bundle 7 \
    --fiber-size 177 --bundle-gap 177 \
    --nx 12288 --ny 12288 --pixelsize 15
```

For quick validation, use `--nbundles 2` (14 fibers, ~6 min).

## NIR commands

```bash
# Install glass catalog first (one-time):
cp opticaldesign/MOSAIC-NIR-Optical_Design/MOSAIC130.{agf,Bgf} \
    ~/Documents/Zemax/GLASSCAT/

# LR-J (Config 1) — ~20 min, ~1.2 GB
uv run --python 3.10 make_hdf.py \
    "opticaldesign/MOSAIC-NIR-Optical_Design/Mosaic_2Cam.ZMX" \
    MOSAIC_NIR_LR_J.hdf \
    --orders -1 0 \
    --grating-surface "VPH grating" --blaze 17.15 --config 1 \
    --name MOSAIC-NIR-LR-J \
    --skip-ccd-check --wl-range 0.95 1.34 \
    --layout bundles --nbundles 90 --fibers-per-bundle 7 \
    --fiber-size 203 --bundle-gap 203 \
    --nx 4096 --ny 4096 --pixelsize 15

# LR-H (Config 2) — ~20 min, ~1.2 GB
uv run --python 3.10 make_hdf.py \
    "opticaldesign/MOSAIC-NIR-Optical_Design/Mosaic_2Cam.ZMX" \
    MOSAIC_NIR_LR_H.hdf \
    --orders -1 0 \
    --grating-surface "VPH grating" --blaze 17.15 --config 2 \
    --name MOSAIC-NIR-LR-H \
    --skip-ccd-check --wl-range 1.43 1.80 \
    --layout bundles --nbundles 90 --fibers-per-bundle 7 \
    --fiber-size 203 --bundle-gap 203 \
    --nx 4096 --ny 4096 --pixelsize 15

# HR (Config 3) — ~20 min, ~1.2 GB
uv run --python 3.10 make_hdf.py \
    "opticaldesign/MOSAIC-NIR-Optical_Design/Mosaic_2Cam.ZMX" \
    MOSAIC_NIR_HR.hdf \
    --orders -1 0 \
    --grating-surface "VPH grating" --blaze 17.15 --config 3 \
    --name MOSAIC-NIR-HR \
    --skip-ccd-check --wl-range 1.52 1.62 \
    --layout bundles --nbundles 90 --fibers-per-bundle 7 \
    --fiber-size 203 --bundle-gap 203 \
    --nx 4096 --ny 4096 --pixelsize 15
```

## Running simulations with the HDF

The ANDES simulator (`C:\Users\Thomas\ANDES_simulator`) also covers MOSAIC:

```bash
# From the ANDES_simulator directory:
uv run --python 3.13 mosaic-sim simulate \
    --source flat --band VIS \
    --hdf /path/to/MOSAIC_VIS_LR_Blue.hdf \
    --output-dir /path/to/output
```

Notes on the simulator:
- Requires Python <3.14 (astropy C extensions don't build on 3.14)
- `--band` must be given when using `--hdf` with a non-default HDF path,
  to avoid band inference that reads missing ANDES HDF files
- MOSAIC VIS band config is in `andes_simulator/core/mosaic.py`
  (wavelength_range, n_fibers, detector_size, etc.)
- Bundle support: `--fiber bundle:N` or `--fiber bundle:N-M`

## Bugs found and fixed

1. **Rotation wrapping**: PyEchelle stores rotation near ±pi that wraps during
   interpolation. Fixed by shifting negative values by +2pi after generation.

2. **Zero groove density (gpmm=0)**: When multiple grating surfaces share the
   same comment name, `find_surface_by_comment` returns the first match, which
   may be inactive (groove density=0) in a different multi-config. Fixed by
   auto-detecting the surface with nonzero groove density.

3. **Swapped dispersion axis**: NIR Zemax model has the slit along X and
   disperses along Y, opposite to VIS convention. Fixed by auto-detecting the
   slit axis from native field coordinates, and swapping tx↔ty + sx↔sy in the
   HDF if dispersion ends up along ty.

## TODO

- [x] VIS LR-B HDF (config 1) — done, 1.9 GB
- [x] VIS LR-R HDF (config 2) — done, 1.9 GB
- [x] NIR LR-J HDF (config 1) — done, 1.2 GB
- [x] NIR LR-H HDF (config 2) — done, 1.2 GB
- [x] NIR HR HDF (config 3) — done, 1.2 GB
- [x] VIS HR-B1 (config 3) — done, 2.1 GB
- [x] VIS HR-R1 (config 4) — done, 2.1 GB
- [x] VIS HR-B2 (config 5) — done, 2.1 GB
- [x] VIS HR-R2 (config 6) — done, 2.1 GB
- [ ] Calibrate NIR fiber pitch (203 µm) by comparing to independent sim
- [ ] Detector choice not finalized (CCD290-99 vs STA5800) — pixel size may change
- [ ] Find or create CRIOGENICI glass catalog (referenced in ZMX but not bundled)
