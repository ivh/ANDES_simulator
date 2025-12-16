"""
Laser Frequency Comb (LFC) source for ANDES wavelength calibrations.

LFC produces unresolved emission lines that are equidistant in frequency,
which translates to logarithmic spacing in wavelength (constant delta-v).
"""

from typing import List, Any, Optional, Tuple
from pathlib import Path
import tempfile
import numpy as np

from pyechelle.sources import CSVSource, ConstantPhotonFlux


# Wavelength ranges for each band (nm) - approximate usable ranges
BAND_WAVELENGTH_RANGES = {
    'U': (310, 390),
    'B': (390, 490),
    'V': (490, 590),
    'R': (620, 950),
    'IZ': (820, 1000),
    'Y': (1000, 1100),
    'J': (1150, 1350),
    'H': (1450, 1800),
}

# Approximate echelle order numbers at band center (from ZEMAX models)
BAND_ORDER_ESTIMATES = {
    'U': 130,
    'B': 100,
    'V': 85,
    'R': 90,   # orders 81-98
    'IZ': 55,
    'Y': 118,  # orders 109-127
    'J': 99,   # orders 90-108
    'H': 76,   # orders 68-83
}


class LFCSource:
    """
    Laser Frequency Comb source for wavelength calibration.

    Generates emission lines equidistant in velocity (constant delta-lambda/lambda),
    which is the characteristic of LFC sources. The line spacing is chosen to give
    approximately the target number of lines per spectral order.
    """

    def __init__(self,
                 band: str,
                 flux_per_line: float = 1e5,
                 lines_per_order: int = 100,
                 project_root: Optional[Path] = None):
        """
        Initialize LFC source.

        Parameters
        ----------
        band : str
            Spectral band (U, B, V, R, IZ, Y, J, H)
        flux_per_line : float
            Photon flux per line (ph/s)
        lines_per_order : int
            Target number of lines per spectral order
        project_root : Path, optional
            Project root directory
        """
        self.band = band
        self.flux_per_line = flux_per_line
        self.lines_per_order = lines_per_order

        if project_root is None:
            self.project_root = Path(__file__).parent.parent.parent
        else:
            self.project_root = project_root

        if band not in BAND_WAVELENGTH_RANGES:
            raise ValueError(f"Unknown band '{band}'. Available: {list(BAND_WAVELENGTH_RANGES.keys())}")

        self._temp_files: List[Any] = []
        self._cached_source = None

    def _calculate_velocity_spacing(self) -> float:
        """
        Calculate velocity spacing to achieve target lines per order.

        For an echelle spectrograph, each order spans approximately
        delta_lambda/lambda = 1/m where m is the order number.

        Returns
        -------
        float
            Velocity spacing in m/s
        """
        c = 299792458.0  # m/s
        m = BAND_ORDER_ESTIMATES[self.band]
        delta_v_over_c = 1.0 / (self.lines_per_order * m)
        return delta_v_over_c * c

    def _generate_lfc_wavelengths(
        self,
        wl_min: Optional[float] = None,
        wl_max: Optional[float] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate LFC wavelengths with constant velocity spacing.

        Returns
        -------
        wavelengths : ndarray
            Wavelengths in nm
        fluxes : ndarray
            Flux values (all equal to flux_per_line)
        """
        band_wl_min, band_wl_max = BAND_WAVELENGTH_RANGES[self.band]
        delta_v = self._calculate_velocity_spacing()
        c = 299792458.0

        # ratio between consecutive wavelengths: lambda_n+1/lambda_n = 1 + delta_v/c
        ratio = 1.0 + delta_v / c

        wavelengths = []
        wl = band_wl_min
        while wl <= band_wl_max:
            wavelengths.append(wl)
            wl *= ratio

        wavelengths = np.array(wavelengths)
        if wl_min is not None or wl_max is not None:
            mask = np.ones_like(wavelengths, dtype=bool)
            if wl_min is not None:
                mask &= (wavelengths >= wl_min)
            if wl_max is not None:
                mask &= (wavelengths <= wl_max)
            wavelengths = wavelengths[mask]

        fluxes = np.full_like(wavelengths, self.flux_per_line)

        return wavelengths, fluxes

    def _create_lfc_source(
        self,
        wl_min: Optional[float] = None,
        wl_max: Optional[float] = None
    ) -> CSVSource:
        """
        Create an LFC source with discrete emission lines.

        Returns
        -------
        CSVSource
            LFC spectrum source with list_like=True for discrete lines
        """
        wavelengths, fluxes = self._generate_lfc_wavelengths(wl_min, wl_max)
        if wavelengths.size == 0:
            raise ValueError(
                f"LFC spectrum has no lines within {wl_min}-{wl_max} nm"
            )

        # write to temp file
        lfc_dir = self.project_root / 'SED' / '.lfc_temp'
        lfc_dir.mkdir(exist_ok=True)

        clip_tag = ""
        if wl_min is not None or wl_max is not None:
            clip_tag = f"_clip_{self._format_wavelength_tag(wl_min)}-{self._format_wavelength_tag(wl_max)}"

        lfc_filename = f'LFC_{self.band}_{self.lines_per_order}lpo_{self.flux_per_line:.0e}{clip_tag}.csv'
        lfc_path = lfc_dir / lfc_filename

        np.savetxt(lfc_path, np.column_stack([wavelengths, fluxes]),
                   delimiter=',', fmt='%.10e')

        # list_like=True tells PyEchelle these are discrete emission lines
        lfc_source = CSVSource(
            file_path=str(lfc_path),
            wavelength_units="nm",
            flux_units="ph/s",
            list_like=True
        )

        return lfc_source

    def _format_wavelength_tag(self, value: Optional[float]) -> str:
        if value is None:
            return "x"
        text = f"{value:.3f}"
        text = text.rstrip('0').rstrip('.')
        return text.replace('.', 'p')

    def get_all_fibers_sources(self,
                               n_fibers: int,
                               skip_fibers: Optional[List[int]] = None) -> List[Any]:
        """
        Get LFC sources for all fibers.

        Parameters
        ----------
        n_fibers : int
            Total number of fibers
        skip_fibers : List[int], optional
            List of 1-based fiber numbers to keep dark

        Returns
        -------
        List
            Source list with all fibers having LFC illumination
        """
        sources = []
        for fiber_idx in range(n_fibers):
            fiber_num = fiber_idx + 1
            if skip_fibers and fiber_num in skip_fibers:
                sources.append(ConstantPhotonFlux(0.0))
            else:
                sources.append(self._create_lfc_source())

        return sources

    def get_single_fiber_sources(self,
                                 fiber_num: int,
                                 n_fibers: int) -> List[Any]:
        """
        Get LFC source for a single fiber.

        Parameters
        ----------
        fiber_num : int
            1-based fiber number to illuminate
        n_fibers : int
            Total number of fibers

        Returns
        -------
        List
            Source list with only specified fiber having LFC illumination
        """
        if not (1 <= fiber_num <= n_fibers):
            raise ValueError(f"Fiber number {fiber_num} out of range [1, {n_fibers}]")

        sources = [ConstantPhotonFlux(0.0) for _ in range(n_fibers)]
        sources[fiber_num - 1] = self._create_lfc_source()

        return sources

    def get_line_info(self) -> dict:
        """
        Get information about the generated LFC lines.

        Returns
        -------
        dict
            Dictionary with line statistics
        """
        wavelengths, _ = self._generate_lfc_wavelengths()
        delta_v = self._calculate_velocity_spacing()

        return {
            'band': self.band,
            'wavelength_range_nm': BAND_WAVELENGTH_RANGES[self.band],
            'n_lines': len(wavelengths),
            'velocity_spacing_m_s': delta_v,
            'velocity_spacing_km_s': delta_v / 1000,
            'target_lines_per_order': self.lines_per_order,
            'estimated_order_number': BAND_ORDER_ESTIMATES[self.band],
            'flux_per_line': self.flux_per_line,
        }
