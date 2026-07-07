"""DPR.TYPE grammar parsing for raw-frame generation.

Grammar per Templates Manual E-AND-SW-MAN-06-00-001 v2.0 (see
~/ANDES/edps/calibration_plan.yaml, reconciliation item 10):
'<KIND>,<A>,<B>' for SL echelle calibrations, '<KIND>,<slit>' for IFU,
no KIND for science ('OBJECT,SKY', 'OBJECT,WAVE' in TC mode, plain
'OBJECT'/'SKY' for IFU). The calibration fibre C is NOT part of DPR.TYPE;
its source is passed separately (calfib) and stamped into the dedicated
ins.calfib keyword. BIAS, DARK and the LED flat 'FLAT,LAMP'
(DPR.TECH=IMAGE) are detector-only frames that never touch pyechelle.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from ..core.instruments import get_instrument_config, get_band_wavelength_range

KINDS = {'BIAS', 'DARK', 'FLAT', 'ORDERDEF', 'WAVE', 'SLITMASK', 'STD'}

SOURCE_TOKENS = {'LAMP', 'FP', 'LFC', 'HCL', 'SKY', 'OBJECT', 'FLUX',
                 'TELLURIC', 'RV', 'WAVE', 'OFF', 'DARK'}

CALFIB_TOKENS = {'FP', 'HCL', 'LFC', 'LAMP', 'OFF'}

DARK_TOKENS = {'OFF', 'DARK'}

# Target peak signal (e-) for the brightest pixels of a slot image. There is
# no realistic end-to-end throughput model, so frame levels are set here;
# EXPTIME in headers is scheduling metadata, not a photon integral.
PEAK_TARGETS_E = {
    'LAMP': 40000.0,
    'FP': 30000.0,
    'LFC': 30000.0,
    'HCL': 30000.0,
    'OBJECT': 25000.0,
    'FLUX': 30000.0,
    'TELLURIC': 30000.0,
    'RV': 25000.0,
    'WAVE': 30000.0,
    'SKY': 1500.0,
}
PEAK_TARGET_SKYFLAT_E = 30000.0  # twilight sky flats are bright

# CSV spectra per token; the first file covering the band is used.
CSV_CANDIDATES = {
    'SKY': ['SED/sky_emission_R.csv', 'SED/sky_emission_YK.csv'],
    'OBJECT': ['SED/phoenix.csv'],
    'RV': ['SED/phoenix.csv'],
    'FLUX': ['SED/HR1544.csv', 'SED/phoenix.csv'],
    'TELLURIC': ['SED/HR1544_transmitted.csv', 'SED/star_transmitted_R.csv'],
}
MIN_BAND_COVERAGE = 0.6

_csv_range_cache: Dict[Path, Tuple[float, float]] = {}


@dataclass
class Slot:
    """One illuminated slot of a raw frame (dark slots are omitted)."""
    name: str          # A, B, C, IFU, CALFIB
    subslit: str       # slitA, slitB, cal_sl, ifu, cal_ifu
    token: str         # source token from the DPR string (or calfib value)
    fibers: List[int]  # resolved 1-based fiber list (after mask intersection)


@dataclass
class DprSpec:
    dpr_type: str
    kind: Optional[str]        # None for science frames
    catg: str                  # CALIB, SCIENCE or TECHNICAL
    tech: str                  # IMAGE, ECHELLE,FIBER or ECHELLE,IFU
    mode: Optional[str]        # SL-UNI, IFU-AO or None (detector-only)
    detector_only: bool
    led: bool                  # detector-only LED flat (uniform illumination)
    calfib: Optional[str] = None   # calibration fibre source (ins.calfib)
    slots: List[Slot] = field(default_factory=list)
    ins_mask: Optional[str] = None


SL_SLOTS = [('A', 'slitA'), ('B', 'slitB')]


def subslit_fibers(band: str, subslit: str) -> List[int]:
    """Resolve a subslit name to its 1-based fiber list for a band."""
    cfg = get_instrument_config(band)
    if subslit in ('slitA', 'slitB', 'cal_sl'):
        slits = cfg.get('fiber_config', {}).get('slits', {})
        key = 'cal_fibers' if subslit == 'cal_sl' else subslit
        if key not in slits:
            raise ValueError(f"Subslit '{subslit}' not available for {band}-band")
        return list(slits[key])
    if subslit in ('ifu', 'cal_ifu'):
        slits = cfg.get('ifu_config', {}).get('slits', {})
        if not slits:
            raise ValueError(f"IFU subslits not available for {band}-band")
        if subslit == 'cal_ifu':
            return list(slits['cal_fibers'])
        fibers = []
        for key in ('ring0', 'ring1', 'ring2', 'ring3', 'ring4'):
            fibers.extend(slits.get(key, []))
        return fibers
    raise ValueError(f"Unknown subslit '{subslit}'")


def mask_fibers(band: str, mask: str) -> List[int]:
    """Fibre-mask pattern: M1/M2/M3 illuminate every third fibre."""
    if mask not in ('M1', 'M2', 'M3'):
        raise ValueError(f"Unknown mask position '{mask}' (expected M1/M2/M3)")
    offset = int(mask[1]) - 1
    n_fibers = get_instrument_config(band)['n_fibers']
    return [f for f in range(1, n_fibers + 1) if (f - 1) % 3 == offset]


def parse_dpr(dpr_type: str, band: Optional[str] = None,
              mode: Optional[str] = None, ins_mask: Optional[str] = None,
              calfib: Optional[str] = None,
              catg: Optional[str] = None, tech: Optional[str] = None) -> DprSpec:
    """Parse a DPR.TYPE string (+ calfib keyword) into a frame specification.

    band is required for frames with slots (fiber lists are band-specific);
    detector-only frames (BIAS, DARK, LED flat) parse without it.
    """
    tokens = [t.strip().upper() for t in dpr_type.split(',') if t.strip()]
    if not tokens:
        raise ValueError("Empty DPR.TYPE")

    kind = tokens[0] if tokens[0] in KINDS else None
    slot_tokens = tokens[1:] if kind else tokens

    calfib = calfib.upper() if calfib else None
    if calfib is not None and calfib not in CALFIB_TOKENS:
        raise ValueError(f"Unknown calfib source '{calfib}' "
                         f"(expected one of {sorted(CALFIB_TOKENS)})")

    if kind in ('BIAS', 'DARK'):
        if slot_tokens:
            raise ValueError(f"{kind} takes no slot tokens: '{dpr_type}'")
        return DprSpec(dpr_type=dpr_type, kind=kind, catg=catg or 'CALIB',
                       tech=tech or 'IMAGE', mode=None,
                       detector_only=True, led=False)

    if kind == 'FLAT' and slot_tokens == ['LAMP'] and mode != 'IFU-AO':
        return DprSpec(dpr_type=dpr_type, kind=kind, catg=catg or 'TECHNICAL',
                       tech=tech or 'IMAGE', mode=None,
                       detector_only=True, led=True)

    unknown = [t for t in slot_tokens if t not in SOURCE_TOKENS]
    if unknown:
        raise ValueError(f"Unknown source token(s) {unknown} in '{dpr_type}'")

    if len(slot_tokens) == 2:
        slot_defs = SL_SLOTS
        calfib_subslit = 'cal_sl'
        mode = mode or 'SL-UNI'
        default_tech = 'ECHELLE,FIBER'
    elif len(slot_tokens) == 1 and mode == 'IFU-AO':
        slot_defs = [('IFU', 'ifu')]
        calfib_subslit = 'cal_ifu'
        default_tech = 'ECHELLE,IFU'
    else:
        raise ValueError(
            f"Cannot interpret DPR.TYPE '{dpr_type}' with mode {mode}: "
            f"expected 2 slot tokens (SL) or 1 with --mode IFU-AO")

    if band is None:
        raise ValueError(f"Band required to resolve fibers for '{dpr_type}'")

    default_catg = 'SCIENCE' if kind is None else 'CALIB'

    mask_set = set(mask_fibers(band, ins_mask)) if ins_mask else None

    def make_slot(name, subslit, token):
        fibers = subslit_fibers(band, subslit)
        if mask_set is not None:
            fibers = [f for f in fibers if f in mask_set]
        if fibers:
            return Slot(name=name, subslit=subslit, token=token, fibers=fibers)
        return None

    slots = []
    for (name, subslit), token in zip(slot_defs, slot_tokens):
        if token in DARK_TOKENS:
            continue
        slot = make_slot(name, subslit, token)
        if slot:
            slots.append(slot)
    if calfib is not None and calfib not in DARK_TOKENS:
        slot = make_slot('C', calfib_subslit, calfib)
        if slot:
            slots.append(slot)

    return DprSpec(dpr_type=dpr_type, kind=kind, catg=catg or default_catg,
                   tech=tech or default_tech, mode=mode, calfib=calfib,
                   detector_only=False, led=False, slots=slots,
                   ins_mask=ins_mask)


def peak_target_e(kind: Optional[str], token: str) -> float:
    """Target peak signal in e- for a slot image."""
    if kind == 'FLAT' and token == 'SKY':
        return PEAK_TARGET_SKYFLAT_E  # twilight sky flat
    return PEAK_TARGETS_E[token]


def _csv_wavelength_range(path: Path) -> Tuple[float, float]:
    """First/last data wavelength (nm) of a two-column CSV."""
    if path in _csv_range_cache:
        return _csv_range_cache[path]
    first = last = None
    with path.open() as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            try:
                wl = float(stripped.split(',')[0])
            except ValueError:
                continue
            if first is None:
                first = wl
            last = wl
    if first is None:
        raise ValueError(f"No data rows in {path}")
    _csv_range_cache[path] = (first, last)
    return first, last


def source_spec_for_token(token: str, band: str,
                          project_root: Path) -> Dict[str, str]:
    """Map a DPR source token to SourceConfig kwargs.

    CSV tokens pick the first candidate file covering enough of the band;
    SKY falls back to a faint constant continuum where no sky spectrum
    exists (the level is set by peak-target normalization anyway). WAVE
    (TC-mode simultaneous reference) defaults to the FP etalon — the lamp
    is a template parameter per the Templates Manual.
    """
    if token == 'LAMP':
        return {'type': 'constant'}
    if token in ('FP', 'WAVE'):
        return {'type': 'fabry_perot'}
    if token == 'LFC':
        return {'type': 'lfc'}
    if token == 'HCL':
        return {'type': 'hcl'}
    if token in CSV_CANDIDATES:
        band_lo, band_hi = get_band_wavelength_range(band, project_root)
        for rel in CSV_CANDIDATES[token]:
            path = project_root / rel
            if not path.exists():
                continue
            lo, hi = _csv_wavelength_range(path)
            overlap = min(hi, band_hi) - max(lo, band_lo)
            if overlap / (band_hi - band_lo) >= MIN_BAND_COVERAGE:
                return {'type': 'csv', 'filepath': rel}
        if token == 'SKY':
            return {'type': 'constant'}
        raise ValueError(
            f"No spectrum in SED/ covers {band}-band for token {token}")
    raise ValueError(f"No source mapping for token '{token}'")
