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
spectrograph-level; filenames carry the spectrograph name (DPID safety).
A `--bands` filter simulates only some detectors; the other extensions get
detector-noise-only pixels flagged `HIERARCH ESO SIM SIMULATED = F`.
`--headers-only` replaces all extensions by 2x2 stubs and skips pyechelle
and the detector model entirely (a full night takes seconds) — enough for
EDPS classification/organization tests. No prescan/overscan regions yet
(geometry undefined in the ADs); the detector layer is the single place to
add them later.

## Simulation cache and the shot-noise fix

Pyechelle output is a Monte Carlo realization, so caching raw simulations
would freeze one shot-noise pattern into every frame built from them —
breaking master-frame stacking statistics and photon-transfer gain
measurement. The cache (`E2E/simcache/`) therefore stores *boosted
expectation* images: each unique (band, model, source, fibers) slot is
simulated once at boost x nominal flux (default 10) and divided by boost;
every exposure then draws fresh Poisson noise from the scaled expectation.
Residual correlated noise is 1/boost of the shot variance — fine for
recipe testing; raise `--boost` or disable the cache for noise studies.
Verified: two exposures sharing a cached slot show variance ratio 0.993 vs
ideal independence; a 10-frame bias stack averages down exactly sqrt(10);
photon transfer on an LED pair recovers the true gain (1.997 vs 2.0).

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
- EDPS gotcha: its file bookkeeping keys on path — always generate into a
  fresh directory and use `edps -w andes.andes_wkf -r` after workflow edits.

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
