# Plan: field-dependent PSFs for MOSAIC via coarse Zemax grid + reslit

Handoff for a session on the Zemax (Windows/OpticStudio) machine. Background in
`MakeHDF.md`; tools are `make_hdf.py` and `reslit_hdf.py` in this directory.

## What we found (2026-10-01, on the Mac)

1. **All 9 MOSAIC HDFs have identical PSFs in every fiber** (bitwise, checked
   NIR LR-J, VIS LR-B, VIS HR-B1). Cause is in pyechelle <= 0.4.0:
   `Field.push_to_zos()` adds only the 4 box corners because `DeleteAllFields()`
   keeps field 1, and `get_psf()` runs `huygens_psf(field=1)` -- i.e. always the
   leftover field, presumably the slit centre. Transformations use the corner
   fields and are fine. The ANDES HDFs (other tooling) do vary along the slit
   (R band: |dPSF|/2 up to 0.02 end to end).
   **Fix** (commit 5c1196d, untested on Zemax): `make_hdf.py` monkey-patches
   `push_to_zos` to move field 1 onto the fiber centre. `--test-api` and a
   post-build check print the PSF difference between first and last fiber.

2. **Slit layout no longer needs Zemax.** `reslit_hdf.py` interpolates the
   transformations along the slit (cubic spline of fiber-centre image and linear
   part per micron of field size) and writes an HDF for any positions/fiber sizes.
   PSFs: nearest grid fiber. Leave-one-out against the traced fibers:

   | grid points along slit | NIR LR-J max err | VIS LR-B max err | VIS HR-B1 max err |
   |---|---|---|---|
   | ~20-35 (every 35th fiber) | 0.55 mpx | 1.4 mpx | 1.4 mpx |
   | ~10-18 (every 70th) | 4.3 mpx | 1.4 mpx | 2.4 mpx |
   | 6-8 | 49 mpx | 8.7 mpx | 32 mpx |

   1.4 mpx (VIS) / 0.35 mpx (NIR) is the float32 floor of `translation_x/y`, not
   interpolation error. Coarse 19-point HDF -> full 630-fiber layout reproduces
   the traced NIR LR-J to <= 0.49 mpx.

3. Fiber size is only in the affine `scale_x/y` (no attribute), so it can be changed
   by reslit too. NIR flat with real 150 um cores vs the current 203 um boxes:
   valley/peak between neighbouring fibers 0.14 vs 0.46 -- current NIR HDFs blend
   fibers much more than reality.

4. `make_hdf.py` now takes `--positions FILE` (one slit position in microns per
   line; field size still `--fiber-size`) and writes `slit_position_um` /
   `fiber_size_um` attributes, which `reslit_hdf.py` reads automatically.

## Steps on the Zemax machine

Get the current `make_hdf.py` and `reslit_hdf.py` (git pull, or copy from the Mac
`E2E/src/HDF/`). Commands below assume the same working directory as in
`MakeHDF.md` (with `opticaldesign/` below it).

### 0. Grid position files

Every 35th fiber of the existing layouts, ends included, so grid points coincide
with fibers of the old HDFs (allows direct comparison):

```bash
uv run --python 3.10 --with numpy python -c "
import numpy as np
def bp(nb, fpb, pitch, gap):
    p, y = [], 0.0
    for _ in range(nb):
        for _ in range(fpb): p.append(y); y += pitch
        y += gap
    return np.array(p) - np.mean(p)
for name, lay, step in [('nir', (90, 7, 203, 203), 35), ('vislr', (140, 7, 177, 177), 35),
                        ('vishr', (60, 19, 152, 152), 35), ('nir_dense', (90, 7, 203, 203), 7)]:
    p = bp(*lay); i = np.unique(np.r_[np.arange(0, len(p), step), len(p) - 1])
    np.savetxt(f'grid_{name}.txt', p[i], fmt='%.1f'); print(name, len(i), 'points, old fiber numbers', (i + 1)[:4].tolist(), '...')
"
```

Gives 19 (nir), 29 (vislr), 34 (vishr), 91 (nir_dense) points.

### 1. Check the PSF fix (minutes)

```bash
uv run --python 3.10 make_hdf.py \
    "opticaldesign/MOSAIC-NIR-Optical_Design/Mosaic_2Cam.ZMX" test.hdf \
    --orders -1 0 --grating-surface "VPH grating" --blaze 17.15 --config 1 \
    --name MOSAIC-NIR-LR-J --skip-ccd-check --wl-range 0.95 1.34 \
    --positions grid_nir.txt --fiber-size 203 \
    --nx 4096 --ny 4096 --pixelsize 15 --test-api
```

Expect `PSF fiber 1 vs fiber 19: |diff|/2 = <nonzero>`. If it prints the WARNING,
the field-1 patch did not take (check that `Fields.GetField(1).X/Y` are settable
in OpticStudio 17.09; alternative: delete and re-add fields, or pass the field
number of the centre to `huygens_psf` by patching `InteractiveZEMAX.get_psf`).

### 2. NIR LR-J dense grid: verify + measure PSF variation (~3 min)

Same command without `--test-api`, `--positions grid_nir_dense.txt`, output
`MOSAIC_NIR_LR_J_grid91.hdf`. Then check:

- **Transformations unchanged by the patch**: grid fiber k sits at old fiber
  `7(k-1)+1` (plus the last). Compare `translation_x/y`, `scale_x/y` with the old
  `MOSAIC_NIR_LR_J.hdf` -- expect agreement at the float32 level (~1 mpx). If not,
  the field normalization changed; stop and investigate.
- **PSF variation along the slit**: per wavelength sample, |dPSF|/2 between
  neighbouring grid points and vs the centre fiber. This decides how many grid
  points the PSFs need (transformations only need ~20). reslit uses
  nearest-neighbour PSFs, so the neighbour difference is the reslit PSF error.

### 3. Production grids

All 9 modes, commands as in `MakeHDF.md` but with `--positions grid_<arm>.txt`
(fiber size unchanged: 203 NIR, 177 VIS LR, 152 VIS HR), outputs named
`<old name>_grid.hdf`. Expected time ~1/30 of the full runs (VIS ~8 min, NIR
~1 min each). Use a denser grid if step 2 shows the PSF needs it.

### 4. Bring back to the Mac

Grid HDFs are small (~19-34 PSF sets). On the Mac, expand with e.g.

```bash
uv run python HDF/reslit_hdf.py HDF/MOSAIC_NIR_LR_J_grid.hdf HDF/MOSAIC_NIR_LR_J_core150.hdf \
    --nbundles 90 --fibers-per-bundle 7 --pitch 203 --gap 203 --fiber-size 150
```

## Open items (not for the Zemax session)

- pyechelle `ZEMAX.psfs()` caches PSFs per order **without the fiber key**, and
  our simulator uses one `ZEMAX` object per run: multi-fiber runs (`--subslit all`,
  slitA, ...) use the first fiber's PSF set for all fibers, also for ANDES.
  Single-fiber runs are fine. Needs a patch in the simulator (or upstream).
- Real slit numbers: NIR core 150 um at 203 um pitch (pitch still to calibrate);
  VIS core sizes and bundle gaps unknown -- currently pitch = field size.
- reslit: PSF blending between grid points (centroid-aligned) if nearest-neighbour
  is too coarse; per-fiber wavelength grids (needed for ANDES HDFs).
- `mosaic-sim` takes n_fibers / fibers_per_bundle from `core/mosaic.py`, so a
  reslit HDF with a different fiber count needs a config entry.
