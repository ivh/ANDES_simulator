"""Hollow-cathode lamp (HCL, ThAr) source for wavelength calibration.

Uses the NIST ThAr line list shipped in SED/linelists/ThAr_nist.csv
(fetched once via pyechelle's NIST/ASDCache machinery; regenerate with
HCLSource.fetch_linelist if needed). Falls back to the resolved PyReduce
ThAr atlas (SED/thar.csv, 300-1060 nm) if the line list is unavailable.
"""

from pathlib import Path
from typing import Any, List, Optional, Tuple

import numpy as np

from pyechelle.sources import CSVSource

from ..core.instruments import get_band_wavelength_range, INSTRUMENTS

LINELIST_REL_PATH = Path('SED/linelists/ThAr_nist.csv')
ATLAS_REL_PATH = Path('SED/thar.csv')


class HCLLineListUnavailable(RuntimeError):
    pass


class HCLSource:
    """ThAr hollow-cathode lamp source with discrete emission lines.

    flux_scale is the photon flux (ph/s) of the brightest line within the
    used wavelength range; all other lines scale by their NIST relative
    intensity.
    """

    def __init__(self,
                 band: str,
                 flux_scale: float = 1e5,
                 project_root: Optional[Path] = None):
        self.band = band
        self.flux_scale = flux_scale
        if project_root is None:
            self.project_root = Path(__file__).parent.parent.parent
        else:
            self.project_root = Path(project_root)

        if band not in INSTRUMENTS:
            raise ValueError(f"Unknown band '{band}'. Available: {list(INSTRUMENTS.keys())}")

    @staticmethod
    def fetch_linelist(project_root: Path, wl_range_um=(0.3, 2.5),
                       elements=('Th', 'Ar')) -> Path:
        """(Re)fetch the NIST line list; needs network on first run."""
        from pyechelle.sources import pull_catalogue_lines

        rows = []
        for elem in elements:
            wl, inten = pull_catalogue_lines(*wl_range_um, elem)
            wl_nm = np.asarray(wl.value if hasattr(wl, 'value') else wl) * 1000.0
            inten = np.asarray(inten, dtype=float)
            ok = np.isfinite(wl_nm) & np.isfinite(inten) & (inten > 0)
            rows.extend(zip(wl_nm[ok], inten[ok], [elem] * int(ok.sum())))
        if not rows:
            raise HCLLineListUnavailable("NIST fetch returned no lines")
        rows.sort()

        out = project_root / LINELIST_REL_PATH
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open('w') as f:
            f.write('# ThAr line list from NIST ASD (via pyechelle/ASDCache), vacuum wavelengths\n')
            f.write('# columns: wavelength_nm, relative_intensity, element\n')
            for w, i, e in rows:
                f.write(f'{w:.6f},{i:.6g},{e}\n')
        return out

    def _load_line_list(self) -> Tuple[np.ndarray, np.ndarray]:
        path = self.project_root / LINELIST_REL_PATH
        if not path.exists():
            raise HCLLineListUnavailable(
                f"{path} missing; run HCLSource.fetch_linelist() (needs network)")
        wavelengths, intensities = [], []
        with path.open() as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith('#'):
                    continue
                parts = stripped.split(',')
                wavelengths.append(float(parts[0]))
                intensities.append(float(parts[1]))
        return np.array(wavelengths), np.array(intensities)

    def get_lines(self, wl_min: Optional[float] = None,
                  wl_max: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray]:
        """Line wavelengths (nm) and photon fluxes (ph/s) within limits."""
        band_lo, band_hi = get_band_wavelength_range(self.band, self.project_root)
        lo = max(band_lo, wl_min) if wl_min is not None else band_lo
        hi = min(band_hi, wl_max) if wl_max is not None else band_hi

        wl, inten = self._load_line_list()
        mask = (wl >= lo) & (wl <= hi)
        wl, inten = wl[mask], inten[mask]
        if wl.size == 0:
            raise ValueError(f"No ThAr lines within {lo}-{hi} nm")
        fluxes = inten / inten.max() * self.flux_scale
        return wl, fluxes

    def _create_hcl_source(self, wl_min: Optional[float] = None,
                           wl_max: Optional[float] = None) -> CSVSource:
        try:
            wavelengths, fluxes = self.get_lines(wl_min, wl_max)
        except HCLLineListUnavailable:
            return self._create_atlas_source(wl_min, wl_max)

        hcl_dir = self.project_root / 'SED' / '.hcl_temp'
        hcl_dir.mkdir(exist_ok=True)
        clip_tag = ""
        if wl_min is not None or wl_max is not None:
            clip_tag = f"_clip_{self._fmt(wl_min)}-{self._fmt(wl_max)}"
        path = hcl_dir / f'HCL_{self.band}_{self.flux_scale:.0e}{clip_tag}.csv'
        np.savetxt(path, np.column_stack([wavelengths, fluxes]),
                   delimiter=',', fmt='%.10e')

        return CSVSource(file_path=str(path), wavelength_units="nm",
                         flux_units="ph/s", list_like=True)

    def _create_atlas_source(self, wl_min: Optional[float],
                             wl_max: Optional[float]) -> CSVSource:
        """Resolved ThAr atlas spectrum (VIS only, 300-1060 nm)."""
        atlas = self.project_root / ATLAS_REL_PATH
        if not atlas.exists():
            raise HCLLineListUnavailable(
                f"Neither {LINELIST_REL_PATH} nor {ATLAS_REL_PATH} available")
        band_lo, band_hi = get_band_wavelength_range(self.band, self.project_root)
        atlas_lo, atlas_hi = 300.0, 1060.0
        if band_lo < atlas_lo - 1 or band_hi > atlas_hi + 1:
            raise HCLLineListUnavailable(
                f"ThAr atlas covers {atlas_lo}-{atlas_hi} nm, not {self.band}-band "
                f"({band_lo:.0f}-{band_hi:.0f} nm); NIST line list required")

        lo = max(band_lo, wl_min) if wl_min is not None else band_lo
        hi = min(band_hi, wl_max) if wl_max is not None else band_hi

        hcl_dir = self.project_root / 'SED' / '.hcl_temp'
        hcl_dir.mkdir(exist_ok=True)
        path = hcl_dir / (f'HCL_atlas_{self.band}_{self.flux_scale:.0e}'
                          f'_{self._fmt(lo)}-{self._fmt(hi)}.csv')
        if not path.exists():
            # atlas is peak-normalized; scale to flux_scale and convert to
            # microns (pyechelle passes bare micron floats to get_counts)
            with atlas.open() as src, path.open('w') as out:
                for line in src:
                    stripped = line.strip()
                    if not stripped or stripped.startswith('#'):
                        continue
                    parts = stripped.split(',')
                    wl = float(parts[0])
                    if wl < lo or wl > hi:
                        continue
                    out.write(f"{wl * 1e-3},{float(parts[1]) * self.flux_scale}\n")
        return CSVSource(file_path=str(path), wavelength_units="um",
                         flux_units="ph/s/AA")

    @staticmethod
    def _fmt(value: Optional[float]) -> str:
        if value is None:
            return "x"
        return f"{value:.3f}".rstrip('0').rstrip('.').replace('.', 'p')

    def get_line_info(self) -> dict:
        wavelengths, fluxes = self.get_lines()
        return {
            'band': self.band,
            'wavelength_range_nm': get_band_wavelength_range(self.band, self.project_root),
            'n_lines': len(wavelengths),
            'flux_scale': self.flux_scale,
        }
