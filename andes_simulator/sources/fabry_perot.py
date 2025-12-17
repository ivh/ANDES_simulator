"""
Fabry-Perot source configurations for ANDES wavelength calibrations.

Provides Fabry-Perot etalon spectrum sources for wavelength calibration
with support for velocity shifts and different finesse values.
"""

from typing import List, Dict, Any, Optional
from pathlib import Path
import random

from pyechelle.sources import CSVSource, ConstantPhotonFlux


class FabryPerotSource:
    """
    Manages Fabry-Perot source configurations for wavelength calibrations.
    
    Handles loading of FP spectra with different finesse values for
    different spectral bands and supports velocity shift simulations.
    """
    
    def __init__(self,
                 band: str,
                 scaling_factor: float = 1e8,
                 project_root: Optional[Path] = None):
        """
        Initialize Fabry-Perot source.
        
        Parameters
        ----------
        band : str
            Spectral band (Y, J, H use finesse 26; others use finesse 23)
        scaling_factor : float
            Flux scaling factor for FP spectrum
        project_root : Path, optional
            Project root directory for finding SED files
        """
        self.band = band
        self.scaling_factor = scaling_factor
        
        if project_root is None:
            self.project_root = Path(__file__).parent.parent.parent
        else:
            self.project_root = project_root
        
        # Define FP spectrum files by band group
        self.fp_files = {
            'YJH': 'FP_simulation_YJH_finesse_26.csv',
            'RIZ': 'FP_simulation_RIZ_finesse_23.csv', 
            'UBV': 'FP_simulation_UBV_finesse_23.csv'
        }
        
        self.dark_source = ConstantPhotonFlux(0.0)
        self._base_fp_source = None
    
    def _get_fp_file(self) -> Path:
        """Get the appropriate FP spectrum file for this band."""
        if self.band in ['Y', 'J', 'H']:
            fp_file = self.fp_files['YJH']
        elif self.band in ['R', 'IZ']:
            fp_file = self.fp_files['RIZ']
        elif self.band in ['U', 'B', 'V']:
            fp_file = self.fp_files['UBV']
        else:
            raise ValueError(f"Unknown band '{self.band}' for FP source")
        
        fp_path = self.project_root / 'SED' / fp_file
        if not fp_path.exists():
            raise FileNotFoundError(f"FP spectrum file not found: {fp_path}")
        
        return fp_path
    
    def _create_fp_source(self, velocity_shift: Optional[float] = None) -> CSVSource:
        """
        Create a single FP source with optional velocity shift.

        Parameters
        ----------
        velocity_shift : float, optional
            Velocity shift in m/s for Doppler simulation

        Returns
        -------
        CSVSource
            FP spectrum source
        """
        import pandas as pd

        fp_path = self._get_fp_file()

        # Read CSV and apply scaling factor BEFORE creating CSVSource
        # This preserves units properly (modifying after creation breaks them)
        df = pd.read_csv(fp_path, header=None, names=['wavelength', 'flux'])
        df['flux'] = df['flux'] * self.scaling_factor

        # Create temporary file with scaled data in project SED directory
        # Use a persistent location instead of system temp to avoid cleanup issues
        scaled_fp_dir = self.project_root / 'SED' / '.scaled_fp'
        scaled_fp_dir.mkdir(exist_ok=True)

        # Create unique filename based on band and scaling factor
        scaled_filename = f'FP_{self.band}_scaled_{self.scaling_factor:.0e}.csv'
        scaled_fp_path = scaled_fp_dir / scaled_filename

        # Write scaled data only if not already exists
        if not scaled_fp_path.exists():
            df.to_csv(scaled_fp_path, index=False, header=False)

        # Create CSVSource with scaled data
        # Use ph/s (integrated photons) not ph/s/AA (flux density)
        # list_like=False (default) - FP is a densely sampled continuous spectrum
        fp_source = CSVSource(
            file_path=str(scaled_fp_path),
            wavelength_units="nm",
            flux_units="ph/s",
            list_like=False
        )
        
        # Note: Velocity shifts are typically handled at the spectrograph level
        # through LocalDisturber rather than modifying the source spectrum
        if velocity_shift is not None:
            # This would require more sophisticated wavelength shifting
            # For now, we'll handle this in the simulator setup
            pass
        
        return fp_source
    
    def get_all_fibers_sources(self,
                              n_fibers: int,
                              skip_fibers: Optional[List[int]] = None) -> List[Any]:
        """
        Get FP sources for all fibers.

        Parameters
        ----------
        n_fibers : int
            Total number of fibers
        skip_fibers : List[int], optional
            List of 1-based fiber numbers to keep dark

        Returns
        -------
        List
            Source list with all fibers having FP illumination
        """
        # Create individual FP source objects for each fiber
        # PyEchelle requires separate object instances to avoid state sharing
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
        """
        Get FP source for single fiber with optional velocity shift.

        Parameters
        ----------
        fiber_num : int
            1-based fiber number to illuminate
        n_fibers : int
            Total number of fibers
        velocity_shift : float, optional
            Velocity shift in m/s

        Returns
        -------
        List
            Source list with only specified fiber having FP illumination
        """
        if not (1 <= fiber_num <= n_fibers):
            raise ValueError(f"Fiber number {fiber_num} out of range [1, {n_fibers}]")

        # Create individual dark source objects for each fiber
        # PyEchelle requires separate object instances to avoid state sharing
        sources = [ConstantPhotonFlux(0.0) for _ in range(n_fibers)]
        fp_source = self._create_fp_source(velocity_shift)
        sources[fiber_num - 1] = fp_source

        return sources
    
    def get_random_velocity_shift(self, sigma_ms: float = 100.0) -> float:
        """
        Generate random velocity shift for Doppler testing.
        
        Parameters
        ----------
        sigma_ms : float
            Standard deviation of velocity shifts in m/s
            
        Returns
        -------
        float
            Random velocity shift in m/s
        """
        return random.gauss(0.0, sigma_ms)
    
    def get_batch_single_fiber_configs(self, 
                                     fiber_nums: List[int],
                                     n_fibers: int,
                                     velocity_shifts: Optional[List[float]] = None) -> Dict[int, Dict]:
        """
        Get configuration for batch single-fiber FP simulations.
        
        Parameters
        ----------
        fiber_nums : List[int]
            List of 1-based fiber numbers to simulate
        n_fibers : int
            Total number of fibers
        velocity_shifts : List[float], optional
            Velocity shifts for each fiber (m/s). If None, no shifts applied.
            
        Returns
        -------
        Dict
            Dictionary mapping fiber number to source configuration and metadata
        """
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
                    'fp_file': str(self._get_fp_file()),
                    'scaling_factor': self.scaling_factor
                }
            }
        
        return configs
    
    @property
    def base_fp_source(self) -> CSVSource:
        """Get the base FP source (cached)."""
        if self._base_fp_source is None:
            self._base_fp_source = self._create_fp_source()
        return self._base_fp_source