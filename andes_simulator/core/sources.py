"""
Source factory for ANDES simulations.

Centralizes creation of PyEchelle source objects to eliminate duplication
across simulator methods.
"""

import logging
import tempfile
from pathlib import Path
from typing import List, Any, Optional

import pandas as pd
from pyechelle.sources import ConstantPhotonFlux, CSVSource

from .config import SourceConfig
from .instruments import get_sed_path

# Physical constants
SPEED_OF_LIGHT = 299792458.0  # m/s
DEFAULT_BANDWIDTH_AA = 100.0  # Angstroms - typical order bandwidth


class SourceFactory:
    """
    Factory for creating PyEchelle source objects.
    
    Handles all source types (constant, CSV, Fabry-Perot) and manages
    temporary files for scaled spectra.
    """
    
    def __init__(self, project_root: Path):
        """
        Initialize the factory.
        
        Parameters
        ----------
        project_root : Path
            Project root directory for resolving relative paths
        """
        self.project_root = project_root
        self.logger = logging.getLogger(__name__)
        self._temp_files: List[Any] = []  # Track temp files for cleanup
    
    def create_source(self, config: SourceConfig, band: str = None) -> Any:
        """
        Create a single source object based on configuration.
        
        Parameters
        ----------
        config : SourceConfig
            Source configuration
        band : str, optional
            Spectral band (required for Fabry-Perot sources)
            
        Returns
        -------
        Source object (ConstantPhotonFlux or CSVSource)
        """
        if config.type == "constant":
            flux_value = self._convert_flux_units(config.flux, config.flux_unit)
            return ConstantPhotonFlux(flux_value)
        
        elif config.type == "csv":
            return self._create_csv_source(config)
        
        elif config.type == "fabry_perot":
            if band is None:
                raise ValueError("Band required for Fabry-Perot sources")
            return self._create_fabry_perot_source(config, band)
        
        else:
            raise ValueError(f"Unknown source type: {config.type}")
    
    def create_dark_source(self) -> ConstantPhotonFlux:
        """Create a dark (zero flux) source."""
        return ConstantPhotonFlux(0.0)
    
    def create_fiber_sources(
        self,
        config: SourceConfig,
        n_fibers: int,
        illuminated_fibers: List[int],
        band: str = None
    ) -> List[Any]:
        """
        Create source list for all fibers.
        
        Parameters
        ----------
        config : SourceConfig
            Source configuration
        n_fibers : int
            Total number of fibers
        illuminated_fibers : List[int]
            1-based fiber numbers to illuminate
        band : str, optional
            Spectral band (required for Fabry-Perot)
            
        Returns
        -------
        List of source objects, one per fiber
        """
        # Start with all dark sources
        sources = [self.create_dark_source() for _ in range(n_fibers)]
        
        # Illuminate specified fibers
        for fiber_idx in illuminated_fibers:
            if 1 <= fiber_idx <= n_fibers:
                # Create separate source instance for each fiber
                # (PyEchelle may modify sources during simulation)
                sources[fiber_idx - 1] = self.create_source(config, band)
        
        return sources
    
    def create_even_odd_sources(
        self,
        config: SourceConfig,
        n_fibers: int,
        band: str = None
    ) -> dict:
        """
        Create source configurations for even/odd fiber illumination.
        
        Parameters
        ----------
        config : SourceConfig
            Source configuration
        n_fibers : int
            Total number of fibers
        band : str, optional
            Spectral band
            
        Returns
        -------
        Dict with 'even' and 'odd' source lists
        """
        even_sources = [self.create_dark_source() for _ in range(n_fibers)]
        odd_sources = [self.create_dark_source() for _ in range(n_fibers)]
        
        for i in range(n_fibers):
            fiber_num = i + 1
            if fiber_num % 2 == 0:
                even_sources[i] = self.create_source(config, band)
            else:
                odd_sources[i] = self.create_source(config, band)
        
        return {'even': even_sources, 'odd': odd_sources}
    
    def _create_csv_source(self, config: SourceConfig) -> CSVSource:
        """Create a CSV-based source."""
        csv_path = Path(config.filepath)
        if not csv_path.is_absolute():
            csv_path = self.project_root / csv_path
        
        if not csv_path.exists():
            raise FileNotFoundError(f"CSV source file not found: {csv_path}")
        
        source = CSVSource(
            filepath=str(csv_path),
            wavelength_unit=config.wavelength_unit,
            flux_in_photons=(config.flux_unit == "ph/s")
        )
        
        # Apply scaling if specified
        if hasattr(source, 'flux_data') and config.scaling_factor != 1.0:
            source.flux_data *= config.scaling_factor
        
        return source
    
    def _create_fabry_perot_source(self, config: SourceConfig, band: str) -> CSVSource:
        """Create a Fabry-Perot source with scaled spectrum."""
        sed_path = get_sed_path(band, self.project_root)
        
        if not sed_path.exists():
            raise FileNotFoundError(f"Fabry-Perot SED file not found: {sed_path}")
        
        # Load and scale the FP spectrum
        df = pd.read_csv(sed_path, header=None, names=['wavelength', 'flux'])
        df['flux'] *= config.scaling_factor
        
        # Write to temporary file
        temp_file = tempfile.NamedTemporaryFile(
            mode='w', suffix='.csv', delete=False
        )
        df.to_csv(temp_file.name, index=False, header=False)
        temp_file.close()
        self._temp_files.append(temp_file)
        
        return CSVSource(
            file_path=temp_file.name,
            wavelength_units="nm",
            flux_units="ph/s",
            list_like=False
        )
    
    def _convert_flux_units(self, flux: float, flux_unit: str) -> float:
        """
        Convert flux to ph/s/AA as required by ConstantPhotonFlux.
        
        Parameters
        ----------
        flux : float
            Flux value in the specified units
        flux_unit : str
            Unit of flux: "ph/s", "ph/s/AA", or "ph/s/nm"
            
        Returns
        -------
        float
            Flux in ph/s/AA
        """
        if flux_unit == "ph/s/AA":
            return flux
        
        elif flux_unit == "ph/s/nm":
            # 1 nm = 10 AA, so ph/s/nm / 10 = ph/s/AA
            return flux / 10.0
        
        elif flux_unit == "ph/s":
            self.logger.warning(
                f"Converting total flux ({flux} ph/s) to flux density. "
                f"Using assumed bandwidth of {DEFAULT_BANDWIDTH_AA} AA. "
                f"For accurate results, use 'ph/s/AA' or 'ph/s/nm' units."
            )
            return flux / DEFAULT_BANDWIDTH_AA
        
        else:
            raise ValueError(
                f"Unknown flux_unit '{flux_unit}'. "
                f"Supported: 'ph/s', 'ph/s/AA', 'ph/s/nm'"
            )
    
    def cleanup(self) -> None:
        """Clean up temporary files created during source creation."""
        import os
        for temp_file in self._temp_files:
            try:
                os.unlink(temp_file.name)
            except OSError:
                pass
        self._temp_files.clear()
    
    def __del__(self):
        """Ensure cleanup on deletion."""
        self.cleanup()
