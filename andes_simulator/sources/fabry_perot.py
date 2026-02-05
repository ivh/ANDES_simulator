"""
Synthetic Fabry-Perot source for ANDES wavelength calibrations.

Computes the Airy transmission profile from physical parameters
(finesse, gap thickness, refractive index) on a logarithmic wavelength grid.
"""

from typing import List, Dict, Any, Optional
from pathlib import Path

import numpy as np
from pyechelle.sources import CSVSource, ConstantPhotonFlux

from ..core.instruments import (
    get_band_wavelength_range, INSTRUMENTS,
    BAND_ORDER_ESTIMATES, DEFAULT_FINESSE,
)


class FabryPerotSource:
    """
    Synthetic Fabry-Perot source using Airy transmission function.

    T(lambda) = 1 / (1 + (2F/pi)^2 * sin^2(2*pi*n*d / lambda))
    """

    def __init__(self,
                 band: str,
                 finesse: Optional[float] = None,
                 gap_mm: Optional[float] = None,
                 n_refr: float = 1.0,
                 lines_per_order: int = 100,
                 scaling_factor: float = 1e8,
                 project_root: Optional[Path] = None):
        self.band = band
        self.finesse = finesse if finesse is not None else DEFAULT_FINESSE[band]
        self.n_refr = n_refr
        self.lines_per_order = lines_per_order
        self.scaling_factor = scaling_factor

        if project_root is None:
            self.project_root = Path(__file__).parent.parent.parent
        else:
            self.project_root = project_root

        if band not in INSTRUMENTS:
            raise ValueError(f"Unknown band '{band}'. Available: {list(INSTRUMENTS.keys())}")

        self.dark_source = ConstantPhotonFlux(0.0)

        wl_min_nm, wl_max_nm = get_band_wavelength_range(band, self.project_root)
        self.wl_min = wl_min_nm
        self.wl_max = wl_max_nm
        wl_center = (wl_min_nm + wl_max_nm) / 2.0
        m_center = BAND_ORDER_ESTIMATES[band]

        if gap_mm is not None:
            self.gap_mm = gap_mm
        else:
            # auto-compute gap so ~lines_per_order peaks fit in the central order
            # FSR = lambda / m  =>  FSR = lambda^2 / (2*n*d)
            # => d = lines_per_order * m_center * lambda_center / (2*n)
            # lambda_center in nm, convert to mm: 1 nm = 1e-6 mm
            self.gap_mm = (lines_per_order * m_center * wl_center * 1e-6) / (2.0 * n_refr)

        self.m_center = m_center

    def _generate_spectrum(self,
                           wl_min: Optional[float] = None,
                           wl_max: Optional[float] = None
                           ) -> tuple:
        """
        Generate Airy transmission spectrum on a logarithmic wavelength grid.

        Returns (wavelengths_nm, transmission) arrays.
        """
        lo = wl_min if wl_min is not None else self.wl_min
        hi = wl_max if wl_max is not None else self.wl_max

        samples_per_fsr = max(50, int(5 * self.finesse))
        # wavelength ratio for log-uniform grid
        # FSR in relative terms is ~1/(lines_per_order * m), so step size
        # must be 1/(samples_per_fsr * lines_per_order * m) to resolve peaks
        r = 1.0 + 1.0 / (samples_per_fsr * self.lines_per_order * self.m_center)
        n_samples = int(np.log(hi / lo) / np.log(r)) + 2
        wavelengths = lo * np.power(r, np.arange(n_samples))
        wavelengths = wavelengths[wavelengths <= hi]

        # Airy function: T = 1 / (1 + (2F/pi)^2 * sin^2(delta/2))
        # where delta = 4*pi*n*d / lambda
        # gap in mm, wavelength in nm => convert both to same unit (nm)
        gap_nm = self.gap_mm * 1e6
        coeff = (2.0 * self.finesse / np.pi) ** 2
        phase = 2.0 * np.pi * self.n_refr * gap_nm / wavelengths
        transmission = 1.0 / (1.0 + coeff * np.sin(phase) ** 2)

        return wavelengths, transmission * self.scaling_factor

    def _create_fp_source(self,
                          wl_min: Optional[float] = None,
                          wl_max: Optional[float] = None) -> CSVSource:
        """Create a CSVSource from the synthetic Airy spectrum."""
        wavelengths, flux = self._generate_spectrum(wl_min, wl_max)
        if wavelengths.size == 0:
            raise ValueError(
                f"FP spectrum has no samples within {wl_min}-{wl_max} nm")

        fp_dir = self.project_root / 'SED' / '.fp_temp'
        fp_dir.mkdir(exist_ok=True)

        clip_tag = ""
        if wl_min is not None or wl_max is not None:
            clip_tag = f"_clip_{self._fmt_wl(wl_min)}-{self._fmt_wl(wl_max)}"

        fp_filename = (f'FP_{self.band}_F{self.finesse:.0f}'
                       f'_gap{self.gap_mm:.4f}{clip_tag}.csv')
        fp_path = fp_dir / fp_filename

        np.savetxt(fp_path, np.column_stack([wavelengths, flux]),
                   delimiter=',', fmt='%.10e')

        return CSVSource(
            file_path=str(fp_path),
            wavelength_units="nm",
            flux_units="ph/s",
            list_like=False
        )

    def _fmt_wl(self, value: Optional[float]) -> str:
        if value is None:
            return "x"
        text = f"{value:.3f}".rstrip('0').rstrip('.')
        return text.replace('.', 'p')

    def get_all_fibers_sources(self,
                               n_fibers: int,
                               skip_fibers: Optional[List[int]] = None) -> List[Any]:
        sources = []
        for fiber_idx in range(n_fibers):
            fiber_num = fiber_idx + 1
            if skip_fibers and fiber_num in skip_fibers:
                sources.append(ConstantPhotonFlux(0.0))
            else:
                sources.append(self._create_fp_source())
        return sources

    def get_single_fiber_sources(self,
                                 fiber_num: int,
                                 n_fibers: int,
                                 velocity_shift: Optional[float] = None) -> List[Any]:
        if not (1 <= fiber_num <= n_fibers):
            raise ValueError(f"Fiber number {fiber_num} out of range [1, {n_fibers}]")

        sources = [ConstantPhotonFlux(0.0) for _ in range(n_fibers)]
        sources[fiber_num - 1] = self._create_fp_source()
        return sources

    def get_batch_single_fiber_configs(self,
                                       fiber_nums: List[int],
                                       n_fibers: int,
                                       velocity_shifts: Optional[List[float]] = None) -> Dict[int, Dict]:
        configs = {}
        if velocity_shifts is None:
            velocity_shifts = [None] * len(fiber_nums)
        elif len(velocity_shifts) != len(fiber_nums):
            raise ValueError("Number of velocity shifts must match number of fibers")

        for fiber_num, v_shift in zip(fiber_nums, velocity_shifts):
            sources = self.get_single_fiber_sources(fiber_num, n_fibers, v_shift)
            configs[fiber_num] = {
                'sources': sources,
                'velocity_shift': v_shift,
                'fiber_num': fiber_num,
                'metadata': {
                    'band': self.band,
                    'finesse': self.finesse,
                    'gap_mm': self.gap_mm,
                    'scaling_factor': self.scaling_factor,
                }
            }
        return configs

    def get_spectrum_info(self) -> dict:
        """Summary of the synthetic FP parameters."""
        gap_nm = self.gap_mm * 1e6
        wl_center = (self.wl_min + self.wl_max) / 2.0
        fsr_nm = wl_center ** 2 / (2.0 * self.n_refr * gap_nm)
        fwhm_nm = fsr_nm / self.finesse
        wavelengths, _ = self._generate_spectrum()
        return {
            'band': self.band,
            'finesse': self.finesse,
            'gap_mm': self.gap_mm,
            'n_refr': self.n_refr,
            'fsr_nm_at_center': fsr_nm,
            'fwhm_nm_at_center': fwhm_nm,
            'n_samples': len(wavelengths),
            'wavelength_range_nm': (self.wl_min, self.wl_max),
        }
