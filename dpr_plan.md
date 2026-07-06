# Plan: simulating EDPS-ready raw frames

Context: the ANDES EDPS workflow (`~/ANDES/edps`) is complete and validated
end-to-end with header-only fake FITS and dummy recipes. The next step is
generating *proper* simulated raw frames for each data type with this
simulator, so the real recipes have realistic pixels to chew on. This note
records the agreed design; see also `~/ANDES/edps/wkf_status.md` (Round 1)
and `~/ANDES/edps/calibration_plan.yaml`.

Revised 2026-07-06 after review: resolved the cache/shot-noise design flaw,
the header-ownership contradiction, and the raw-format question (now follows
DRL spec v1.2 Ch. 4.1); corrected the fiber-list claim; HCL turned out to be
mostly existing machinery.

## Key idea: the DPR grammar is the simulation recipe

The workflow's resolved DPR.TYPE grammar (`~/ANDES/edps/andes/andes_classification.py`)
is `<KIND>,<A>,<C>,<B>`: sub-slit A, calibration fibre C, sub-slit B, e.g.
`WAVE,HCL,FP,FP` = HCL on slit A, FP on the calibration fibre, FP on slit B.
IFU frames use `<KIND>,<slit>[,<calfib>]`, science frames have no KIND
(`OBJECT,FP,SKY`). That maps 1:1 onto this simulator's building blocks: one
per-subslit simulation per slot, then sum.

So: **no script per data type**. One generic raw-frame builder, driven by the
DPR string plus the declarative calibration-plan YAML. Per-type variety lives
in data, not code.

## Raw FITS format (per DRL spec v1.2 Ch. 4.1)

DRL baseline: **one FITS file per spectrograph** (UBV, RIZ, YJH), each
detector arm in a separate image extension; the DRS adapts to either option,
but we follow the baseline. Consequences:

- `SEQ.ARM` stays spectrograph-level (`RIZ`, `YJH`) — matching the existing
  edps test data, so no workflow change at all (this reverses the earlier
  per-detector-arm idea).
- MEF layout: header-only primary HDU carrying all classification keywords;
  one extension per detector band (`EXTNAME` = band name, `BUNIT = 'ADU'`,
  gain/RON keys). Data is uint16 with saturation clipping (16-bit ADCs).
- Filenames contain the spectrograph name (`ANDES_RIZ_0001.fits`) so two
  files taken at the same time cannot collide in DPID (DRL 4.1 requirement).
- A `--bands` filter lets you simulate only some detectors of an arm; the
  other extensions are still present but contain detector-noise-only pixels
  and are flagged `HIERARCH ESO SIM SIMULATED = F`. Recipes under test read
  their extension; the file format stays uniform.
- Prescan/overscan regions: geometry is not defined in any AD yet, so frames
  are exactly detector-sized for now. The detector layer is the single place
  to add overscan later; bias recipes work from the BIAS frames meanwhile.

## Components (implementation: `andes_simulator/raw/`)

1. **`raw/dpr.py`** — grammar parser. DPR string + INS.MODE -> list of
   (fiber list, source token) slots, plus KIND / DPR.CATG / DPR.TECH
   defaults. Slot -> subslit: A->`slitA`, C->`cal_sl`, B->`slitB`; IFU ->
   `ifu` / `cal_ifu`. `BIAS`, `DARK` and LED `FLAT,LAMP` (tech IMAGE) are
   detector-only (no pyechelle). Source tokens: LAMP->flat, FP->fp,
   LFC->lfc, HCL->hcl, OFF/DARK->skip; SKY/OBJECT/FLUX/TELLURIC/RV map to
   CSV spectra chosen per band from `SED/` (first file covering the band
   wavelength range; constant-flux fallback).
   Slit-mask frames (`ins.mask` M1/M2/M3): the mask illuminates every third
   fibre; implemented as an intersection of each slot's fiber list with the
   mask pattern. (Correction to the earlier note: the CLI does *not* expose
   arbitrary fiber lists, but `FiberConfig(mode='custom', fibers=[...])`
   does internally, which is what the builder uses.)

2. **`raw/detector.py`** — detector layer (pyechelle output is pure photon
   counts). Pipeline: expectation e- image x PRNU map -> + dark current
   (incl. hot pixels) -> dead pixels -> Poisson -> charge binning (VIS
   configs) -> read noise -> /gain -> + bias level -> uint16 clip.
   CCD (UBV/RIZ; fast/slow readout with different RON/gain) vs HAWAII4RG
   (YJH; no bias frames, higher dark). Cosmetics (PRNU, hot/dead pixels,
   bad columns) are *static per band* (seeded from the band name) so bad
   pixel maps and flat fields are consistent across frames — essential for
   recipes. Same layer generates BIAS / DARK / LED-flat frames without
   pyechelle (LED = uniform illumination x PRNU, exposure-time series for
   gain/linearity). Parameter values in `core/andes.py` are documented
   placeholders until real detector specs exist.

3. **`raw/cache.py`** — simulation cache, with the shot-noise fix.
   Pyechelle output is a Monte Carlo *realization*, so a naively cached
   image would carry a frozen shot-noise pattern into every frame built
   from it (breaking master-frame stacking statistics and any
   photon-transfer/gain measurement). Instead the cache stores
   **boosted expectation images**: simulate once at `boost` x nominal flux
   (default 10), divide by boost, store as float32 rate image keyed by
   (band, hdf model, source signature, fiber tuple, boost). Each exposure
   then draws fresh Poisson noise from the scaled expectation. Residual
   correlated noise is 1/boost of shot variance — equivalent to the frames
   sharing a common calibration truth at the level of a boost-frame stack.
   Good enough for recipe functional testing; raise `--boost` (or disable
   the cache) for noise-property studies. Total photon budget at boost 10
   is about the same as simulating a 10-frame stack once, and the cached
   slot images are reused across *all* frame types (the FP-on-`cal_sl`
   image appears in slit, LSF, wave and every science frame).
   Fibre efficiencies are seeded per band (not random per run): they are a
   static property of the instrument, and flat-fielding logic breaks if
   flat and wave frames see different fibre throughputs.

4. **Flux normalization**: no realistic end-to-end throughput model exists,
   so per-frame photon levels are set by target peak counts per
   (KIND, source token) — e.g. flats ~40 ke- peak, wave lines ~30 ke-,
   night-sky slots faint — applied as a scale on the cached expectation.
   EXPTIME in the header is scheduling metadata from the plan YAML, not a
   photon integral. Per-band flux hints in the YAML can override later.

5. **`raw/headers.py`** — ESO header stamping: INSTRUME, TELESCOP, MJD-OBS,
   DATE-OBS, EXPTIME, DPR.CATG/TYPE/TECH, SEQ.ARM, INS.MODE,
   DET.BINX/BINY, TPL.START/ID/NEXP/EXPNO, plus the proposed keywords from
   the plan YAML (det.readout, ins.ifu.scale, ins.mask) and simulation
   provenance keys. Ownership resolution: the earlier "header code in
   exactly one place" clashed with keeping `make_test_data.py` in edps for
   fast CI. Resolution: **`calibration_plan.yaml` is the single source of
   truth**, not any code file. E2E's header writer and edps'
   header-only generator are both thin writers checked against the YAML —
   the edps classification tests enforce consistency on their side, the
   night driver consumes the YAML directly on ours. Drift is caught by
   tests, not prevented by shared imports across repos.

6. **`raw/builder.py` + `andes-sim make-raw`**:
   `make-raw --arm RIZ --dpr "WAVE,HCL,FP,FP" --mode SL-UNI --exptime 120
   --nexp 2 -o dir` — parse, fetch/simulate slot expectations via the
   cache, normalize, sum, apply per-exposure detector layer, write one MEF
   per exposure. `--seed` for reproducible noise. Parallelization happens
   at the level of slot simulations (subprocesses; the package already
   provisions per-process `NUMBA_CACHE_DIR`, pyechelle needs `max_cpu=1`).

7. **`raw/night.py` + `andes-sim night <plan.yaml> -o <dir>`**: reads the
   canonical plan (versioned in the edps repo), emits a synthetic
   calibration night (+ optionally night calibrations and science) per
   requested arm. Filters: `--arms`, `--sets` (detector/daily/when-used),
   `--include night,science`, `--bands`. VIS binning/readout configs
   default to 1x1_fast; IFU scale defaults to one value. Replaces
   `make_test_data.py` for real-pixel testing; the header-only generator
   stays in the edps repo for fast CI (see ownership resolution above).

## Sources: what existed and what was missing

- **HCL**: much less work than assumed. `SED/thar.csv` already carries the
  PyReduce ThAr atlas (300-1060 nm, resolved spectrum), and pyechelle 0.4
  ships `ArcLamp`/`pull_catalogue_lines` (NIST line lists via ASDCache,
  ThAr default) plus a generic `LineList` source. `sources/hcl.py` uses
  NIST line lists where available (committed to `SED/linelists/` after
  first fetch, so no network dependency afterwards; needed for J/H where
  the atlas ends), falling back to the atlas spectrum in the VIS. The line
  list used is exported as a catalog file — the DRS needs it as its
  HCL_LINES_TABLE, so simulator and recipe share the same truth.
- **Sky**: `SED/sky_emission_R.csv` (600-780 nm) and `sky_emission_YK.csv`
  (950 nm-) already exist and cover R and YJH; other bands fall back to a
  faint constant continuum until proper skycalc spectra are dropped in.
- **Stars**: `SED/phoenix.csv` (OBJECT, RV), `HR1544.csv` (FLUX),
  `HR1544_transmitted.csv` (TELLURIC), per-band `star_transmitted_*.csv`
  preferred where they cover the band.
- **BIAS / DARK / LED flats**: pure detector-model frames, no pyechelle
  (see detector layer). NIRPS darks in `~/ANDES/E2E/NIRPS_darks/` remain
  an option as realistic NIR noise templates (not wired in yet).

## Open reconciliation items (tracked in calibration_plan.yaml)

The AD2 calibration plan and the DRL spec v1.2 disagree on wave-frame
illumination patterns, order-definition frames, and several keywords
(readout mode, IFU spaxel scale, mask position). The YAML's
`reconciliation` section is the running list; resolve there first, then
update workflow classification and this simulator together. The simulator
follows AD2 (via the YAML) and writes the proposed keywords; nothing here
blocks on the reconciliation.
