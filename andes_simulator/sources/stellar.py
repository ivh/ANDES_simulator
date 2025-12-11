"""
Stellar source configurations for ANDES science simulations.

Provides stellar spectrum sources for science observations using
custom CSVSource spectral energy distributions.
"""

from typing import List, Dict, Any, Optional, Union
from pathlib import Path
import astropy.units as u

from pyechelle.sources import CSVSource, ConstantPhotonFlux


class StellarSource:
    """
    Manages stellar source configurations for science simulations.
    
    Handles loading of stellar spectra from CSVSource files with proper
    unit handling and flux scaling.
    """
    
    def __init__(self,
                 spectrum_file: Union[str, Path],
                 scaling_factor: float = 5e3,
                 wavelength_unit: str = "nm",
                 flux_unit: str = "ph/s",
                 project_root: Optional[Path] = None):
        """
        Initialize stellar source.
        
        Parameters
        ----------
        spectrum_file : str or Path
            Path to CSVSource spectrum file (relative to project root or absolute)
        scaling_factor : float
            Flux scaling factor
        wavelength_unit : str
            Wavelength units in spectrum file
        flux_unit : str
            Flux units in spectrum file
        project_root : Path, optional
            Project root directory
        """
        self.spectrum_file = Path(spectrum_file)
        self.scaling_factor = scaling_factor
        self.wavelength_unit = wavelength_unit
        self.flux_unit = flux_unit
        
        if project_root is None:
            self.project_root = Path(__file__).parent.parent.parent
        else:
            self.project_root = project_root
        
        self.dark_source = ConstantPhotonFlux(0.0)
        self._base_stellar_source = None
        
        # Resolve full path to spectrum file
        if not self.spectrum_file.is_absolute():
            self.spectrum_file = self.project_root / self.spectrum_file
        
        if not self.spectrum_file.exists():
            raise FileNotFoundError(f"Spectrum file not found: {self.spectrum_file}")
    
    def _create_stellar_source(self) -> CSVSource:
        """
        Create a stellar spectrum source from CSVSource file.
        
        Returns
        -------
        CSVSource
            Stellar spectrum source
        """
        # Determine flux units for pyechelle
        flux_units = u.Unit(self.flux_unit) if self.flux_unit != "ph/s" else u.Unit("ph/s")
        
        # Create CSVSource source
        stellar_source = CSVSource(
            str(self.spectrum_file),
            wavelength_units=self.wavelength_unit,
            flux_units=flux_units
        )
        
        # Apply scaling factor
        if hasattr(stellar_source, 'data') and 'flux' in stellar_source.data.columns:
            # Use to_numpy() to get raw values without units for scaling
            raw_flux = stellar_source.data["flux"].to_numpy()
            scaled_flux = raw_flux * self.scaling_factor
            stellar_source.data["flux"] = scaled_flux
        
        return stellar_source
    
    def get_single_fiber_sources(self, fiber_num: int, n_fibers: int) -> List[Any]:
        """
        Get stellar source for single fiber.
        
        Parameters
        ----------
        fiber_num : int
            1-based fiber number to illuminate
        n_fibers : int
            Total number of fibers
            
        Returns
        -------
        List
            Source list with only specified fiber having stellar illumination
        """
        if not (1 <= fiber_num <= n_fibers):
            raise ValueError(f"Fiber number {fiber_num} out of range [1, {n_fibers}]")
        
        sources = [self.dark_source] * n_fibers
        stellar_source = self._create_stellar_source()
        sources[fiber_num - 1] = stellar_source
        
        return sources
    
    def get_all_fibers_sources(self, 
                              n_fibers: int,
                              skip_fibers: Optional[List[int]] = None) -> List[Any]:
        """
        Get stellar sources for all fibers (same spectrum).
        
        Parameters
        ----------
        n_fibers : int
            Total number of fibers
        skip_fibers : List[int], optional
            List of 1-based fiber numbers to keep dark
            
        Returns
        -------
        List
            Source list with all fibers having stellar illumination
        """
        stellar_source = self._create_stellar_source()
        sources = [stellar_source] * n_fibers
        
        if skip_fibers:
            for fiber_num in skip_fibers:
                if 1 <= fiber_num <= n_fibers:
                    sources[fiber_num - 1] = self.dark_source
        
        return sources
    
    def get_multi_fiber_sources(self, 
                               fiber_configs: Dict[int, Dict],
                               n_fibers: int) -> List[Any]:
        """
        Get sources for multiple fibers with different spectra.
        
        Parameters
        ----------
        fiber_configs : Dict[int, Dict]
            Dictionary mapping fiber number to spectrum configuration.
            Each config should have 'spectrum_file' and optionally 'scaling_factor'.
        n_fibers : int
            Total number of fibers
            
        Returns
        -------
        List
            Source list with specified fibers having their respective spectra
        """
        sources = [self.dark_source] * n_fibers
        
        for fiber_num, config in fiber_configs.items():
            if not (1 <= fiber_num <= n_fibers):
                raise ValueError(f"Fiber number {fiber_num} out of range [1, {n_fibers}]")
            
            # Create source for this fiber
            spectrum_file = config['spectrum_file']
            scaling = config.get('scaling_factor', self.scaling_factor)
            
            # Create temporary stellar source with this configuration
            temp_source = StellarSource(
                spectrum_file=spectrum_file,
                scaling_factor=scaling,
                wavelength_unit=self.wavelength_unit,
                flux_unit=self.flux_unit,
                project_root=self.project_root
            )
            
            sources[fiber_num - 1] = temp_source._create_stellar_source()
        
        return sources
    
    def get_custom_sources(self, 
                          illuminated_fibers: List[int], 
                          n_fibers: int) -> List[Any]:
        """
        Get stellar sources for custom fiber selection.
        
        Parameters
        ----------
        illuminated_fibers : List[int]
            List of 1-based fiber numbers to illuminate
        n_fibers : int
            Total number of fibers
            
        Returns
        -------
        List
            Source list with specified fibers illuminated
        """
        sources = [self.dark_source] * n_fibers
        stellar_source = self._create_stellar_source()
        
        for fiber_num in illuminated_fibers:
            if not (1 <= fiber_num <= n_fibers):
                raise ValueError(f"Fiber number {fiber_num} out of range [1, {n_fibers}]")
            sources[fiber_num - 1] = stellar_source
        
        return sources
    
    @property
    def base_stellar_source(self) -> CSVSource:
        """Get the base stellar source (cached)."""
        if self._base_stellar_source is None:
            self._base_stellar_source = self._create_stellar_source()
        return self._base_stellar_source
    
    @classmethod
    def from_template(cls, 
                     template_name: str,
                     project_root: Optional[Path] = None,
                     **kwargs) -> 'StellarSource':
        """
        Create stellar source from predefined template.
        
        Parameters
        ----------
        template_name : str
            Name of template spectrum ('sun', 'k_star', 'm_dwarf', etc.)
        project_root : Path, optional
            Project root directory
        **kwargs
            Additional parameters for source creation
            
        Returns
        -------
        StellarSource
            Stellar source configured with template spectrum
        """
        if project_root is None:
            project_root = Path(__file__).parent.parent.parent
        
        # Define template spectra paths
        templates = {
            'sun': 'SED/sun_spectrum.csv',
            'k_star': 'SED/k_star_spectrum.csv', 
            'm_dwarf': 'SED/m_dwarf_spectrum.csv',
            'test': 'SED/test_spectrum.csv'
        }
        
        if template_name not in templates:
            raise ValueError(f"Unknown template '{template_name}'. Available: {list(templates.keys())}")
        
        spectrum_file = templates[template_name]
        
        return cls(
            spectrum_file=spectrum_file,
            project_root=project_root,
            **kwargs
        )