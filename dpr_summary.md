# EDPS-ready raw frames: design summary (as built)

The `andes_simulator/raw/` package generates simulated ANDES raw frames
that the EDPS workflow (`~/ANDES/edps`) classifies and organizes like real
instrument data. Two entry points: `andes-sim make-raw` (one template
execution from a DPR.TYPE string) and `andes-sim night` (a full synthetic
night driven by `~/ANDES/edps/calibration_plan.yaml`). Command examples
live in CLAUDE.md; this file records the design decisions and their
reasons. It replaces the original dpr_plan.md (2026-07-06/07).

## The DPR grammar is the simulation recipe

DPR.TYPE follows the Templates Manual E-AND-SW-MAN-06-00-001 v2.0
(Tables 2 and 5), which superseded the DRL spec v1.2 reading — see
reconciliation item 10 in calibration_plan.yaml:

- SL calibrations: `<KIND>,<A>,<B>` (first value = fiber A, second =
  fiber B); KIND one of ORDERDEF, SLITMASK, FLAT, WAVE, STD.
- IFU calibrations: `<KIND>,<slit>`.
- Science: no KIND — `OBJECT,SKY` / `SKY,OBJECT` (TS), `OBJECT,WAVE` /
  `WAVE,OBJECT` (TC, simultaneous reference in the other fiber), plain
  `OBJECT`/`SKY` for IFU.
- The calibration fibre C is NOT in DPR.TYPE (ESO-044156): its source is
  a dedicated keyword, provisionally `ins.calfib` (FP/HCL/LFC/LAMP/OFF,
  name pending the ICD). In make-raw it is the `--calfib` option; the
  builder still simulates the C fibers, so pixels stay physically complete.
- `BIAS`, `DARK` and the LED flat `FLAT,LAMP` (DPR.TECH=IMAGE) are
  detector-only frames that never touch pyechelle.

Each DPR string + calfib resolves to per-band slots (fiber list + source);
one cached pyechelle simulation per slot, summed, then the detector layer.
No script per data type: per-type variety lives in the plan YAML, not code.

## Raw FITS format (DRL spec v1.2 Ch. 4.1)

One MEF per spectrograph (UBV/RIZ/YJH): header-only primary HDU with the
ESO classification keywords, one uint16 image extension per detector band
(EXTNAME = band, gain/RON keys, saturation-clipped). `SEQ.ARM` is
spectrograph-level. Filenames are `ANDES_<ARM>_<DATE-OBS>.fits` with
colons replaced (`ANDES_RIZ_2026-07-07T12_00_00.000.fits`): the arm
disambiguates simultaneous exposures (DPID safety, DRL 4.1), the
timestamp makes names unique per exposure and file age visible;
regenerating with the same `--tpl-start` overwrites the same files.
A `--bands` filter simulates only some detectors; the other extensions get
detector-noise-only pixels flagged `HIERARCH ESO SIM SIMULATED = F`.
`--headers-only` replaces all extensions by 2x2 stubs and skips pyechelle
and the detector model entirely (a full night takes seconds) — enough for
EDPS classification/organization tests. No prescan/overscan regions yet
(geometry undefined in the ADs); the detector layer is the single place to
add them later.

## Night driver

`andes-sim night` is the batch layer over the same builder that serves
`make-raw`: it expands the canonical calibration plan into a complete
night, so the frame lists exist only in the YAML. Three phases:

1. *Planning* (`plan_night`): walk `procedures` (and, with
   `--include night,science`, the night calibrations and observations),
   filter by `--arms`/`--sets`, skip `reference` entries, choose one VIS
   binning/readout config (`--vis-config`) and one IFU scale. Every
   exposure entry becomes a planned run with DPR type, count, exptime
   (LED `varied` expands to the 1-60 s linearity series) and the
   per-exposure keywords (`ins.calfib`, `ins.mask`) from the YAML.
   `--dry-run` prints this table.
2. *Time and grouping* (`run_night`): a synthetic clock (start `--date`,
   default 10:00 UT) advances by exptime + 60 s readout per frame and
   feeds DATE-OBS/MJD-OBS and the filenames. Runs sharing (arm, template,
   mode) share one TPL.START group with continuous TPL.EXPNO: EDPS groups
   by tpl.start, and AD2 spreads e.g. wavelength calibration over six CPs
   of one template whose HCL and FP frames must land in one group
   (reconciliation item 9).
3. *Building*: each run is one `RawFrameBuilder.build()` call — the
   make-raw code path with the cache, detector layer and MEF writing;
   `--bands`, `--jobs`, `--headers-only`, `--seed`, `--boost` pass
   through.

## Simulation cache and the shot-noise fix

Pyechelle output is a Monte Carlo realization, so caching raw simulations
would freeze one shot-noise pattern into every frame built from them —
breaking master-frame stacking statistics and photon-transfer gain
measurement. The cache (`E2E/simcache/`) therefore stores *boosted
expectation* images: each unique (band, model, source, fibers) slot is
simulated once at boost x nominal flux (default 10) and divided by boost;
every exposure then draws fresh Poisson noise from the scaled expectation.
Verified: two exposures sharing a cached slot show variance ratio 0.993 vs
ideal independence; a 10-frame bias stack averages down exactly sqrt(10);
photon transfer on an LED pair recovers the true gain (1.997 vs 2.0).

The frozen MC noise left in the expectation is set by the *simulated*
photon count relative to the shipped peak level (the peak-target
normalization rescales the image, not its statistics). Measured on the
R-band daily set: FP line cores carry a ~2.5 percent static pattern
(identical in all frames, so it cancels in frame-to-frame drift tests but
floors absolute-wavecal accuracy at roughly the m/s level); HCL is ~0.8
percent after its scaling fix (2026-07-07: hcl lost the /20 divisor that
made sense for the equal-intensity LFC but starved the ThAr lines, whose
NIST intensities span ~4 decades below the brightest). The resolved flux
scaling is part of the cache key, so retuning invalidates entries.

Fibre efficiencies and detector cosmetics (PRNU, hot/dead pixels, bad
columns) are seeded per band: they are static instrument properties, and
flat-fielding logic breaks if they vary between frames.

## Flux levels

There is no realistic end-to-end throughput model: slot images are
normalized to target peak counts per (KIND, source token) —
`PEAK_TARGETS_E` in `raw/dpr.py` (flats ~40 ke-, lines ~30 ke-, night-sky
slots faint, twilight sky flats bright). Header EXPTIME is scheduling
metadata from the plan YAML, not a photon integral.

## Header ownership

`calibration_plan.yaml` (versioned in the edps repo) is the single source
of truth linking AD2 procedures, templates, DPR keywords and recipes. Both
header writers — this repo's `raw/headers.py` (real pixels) and the edps
repo's `tests/make_test_data.py` (header-only CI shortcut) — are thin and
checked against it: the edps classification tests enforce consistency on
their side, the night driver consumes the YAML directly on ours.

## Sources

- HCL: NIST ThAr line list in `SED/linelists/ThAr_nist.csv` (20322 lines,
  300-2500 nm, fetched once via pyechelle/ASDCache; covers J/H beyond the
  PyReduce atlas `SED/thar.csv`, which serves as VIS fallback). Exported as
  a line catalog so DRS and simulator share the same wavelength truth.
- Sky: `SED/sky_emission_R.csv` (R) and `sky_emission_YK.csv` (YJH);
  other bands fall back to a faint continuum.
- Stars: `SED/phoenix.csv` (OBJECT, RV), `HR1544.csv` (FLUX),
  `HR1544_transmitted.csv` / `star_transmitted_R.csv` (TELLURIC), chosen
  per band by wavelength coverage.
- The TC-mode `WAVE` science token defaults to the FP etalon (the lamp is
  a template parameter per the manual).
- Slit-mask frames: masks M1-M3 = every third fiber, intersected with the
  slot fiber lists (`ins.mask` keyword, `mask1-3` subslit aliases).

## Validation status (2026-07-07)

- Real-pixel frames (YJH BIAS/DARK/WAVE with HCL+FP) classify to exactly
  their tags with the real edps rules; noise physics verified (see above).
- A full RIZ headers-only night (86 frames + static tables) organizes in
  EDPS into 40/40 complete jobs across the whole cascade (bias through
  science); the rv_std chain organizes complete when targeted (it is a
  leaf, absent from the default science-driven dump).
- EDPS gotcha: its file bookkeeping keys on path — regenerating different
  content under identical paths serves stale records. The timestamped
  filenames avoid this unless `--tpl-start` is pinned; use
  `edps -w andes.andes_wkf -r` after workflow edits.

## Known gaps (tracked in calibration_plan.yaml reconciliation)

- AD2-only patterns stay unclassified by design: the combined flat
  `FLAT,LAMP,LAMP` and the one-aperture wave frames `WAVE,FP,OFF` /
  `WAVE,OFF,FP` / `WAVE,HCL,OFF` / `WAVE,OFF,HCL` (items 1-2). The AD2
  cross-calibration pairs map onto the manual's `WAVE,HCL,FP` /
  `WAVE,FP,HCL` and do classify.
- AD2 lacks order-definition and per-slit flat procedures entirely; the
  plan YAML carries them as `status: proposed` entries (manual templates)
  because the DRS cascade cannot organize without ORDER_TABLE and per-slit
  flat products.
- How the C fibre gets order-defined and flat-fielded is open in all
  documents now that C left the DPR grammar (item 2).
- Detector parameters in `core/andes.py` DETECTOR_MODELS are placeholders;
  psf-process caveats in CLAUDE.md still apply.
