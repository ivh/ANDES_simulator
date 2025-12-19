"""
Main ANDES simulation engine.

Coordinates all aspects of ANDES spectrograph simulations including
source setup, instrument configuration, and output generation.
"""

import logging
import numpy as np
from pathlib import Path
from typing import List, Dict, Any, Optional, Union

# Import pyechelle components
from pyechelle.simulator import Simulator
from pyechelle.telescope import Telescope
from pyechelle.spectrograph import ZEMAX, LocalDisturber
from pyechelle.CCD import CCD
from pyechelle.hdfbuilder import HDFBuilder
from pyechelle.spectrograph import InteractiveZEMAX

from .config import SimulationConfig
from .instruments import get_hdf_model_path, BAND_WAVELENGTH_RANGES
from .sources import SourceFactory, SPEED_OF_LIGHT

# Monkey-patch pyechelle to include fiber number in order logging
# Fragile: depends on pyechelle internals. If broken, just remove this block.
import pyechelle.raytracing as _raytracing
_original_prepare_raytracing = _raytracing.prepare_raytracing
_current_fiber = None  # set by our wrapper before calling simulator.run()

def _patched_prepare_raytracing(o, fiber, *args, **kwargs):
    import builtins, io, sys
    # Capture original output
    capture = io.StringIO()
    _original_print = builtins.print
    builtins.print = lambda *a, **k: print(*a, **{**k, 'file': capture})
    try:
        result = _original_prepare_raytracing(o, fiber, *args, **kwargs)
    finally:
        builtins.print = _original_print
    # Reprint with fiber prefix
    line = capture.getvalue().strip()
    if line:
        print(f"Fiber {fiber:2d} {line}")
    return result

_raytracing.prepare_raytracing = _patched_prepare_raytracing


class AndesSimulator:
    """
    Main ANDES simulation engine.
    
    Handles all types of ANDES simulations including flat field calibrations,
    Fabry-Perot wavelength calibrations, stellar observations, and HDF model generation.
    """
    
    def __init__(self, config: SimulationConfig):
        """
        Initialize the simulator with a configuration.
        
        Parameters
        ----------
        config : SimulationConfig
            Complete simulation configuration
        """
        self.config = config
        self.instrument_config = config.instrument_config
        self.project_root = Path(__file__).parent.parent.parent
        
        # Initialize logging
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        
        # Initialize source factory
        self.source_factory = SourceFactory(self.project_root)
        
        # Will be set during simulation setup
        self.simulator = None
        self.sources = None
    
    def cleanup(self) -> None:
        """Clean up resources (temporary files, etc.)."""
        self.source_factory.cleanup()

    def _parse_fib_eff(self, fib_eff: str) -> tuple:
        """
        Parse fiber efficiency specification.

        Parameters
        ----------
        fib_eff : str
            Efficiency value as "0.9" or range as "0.7-0.9"

        Returns
        -------
        tuple
            (min_eff, max_eff) - if equal, use constant; if different, use random
        """
        if '-' in fib_eff:
            parts = fib_eff.split('-')
            if len(parts) != 2:
                raise ValueError(f"Invalid fib_eff format: {fib_eff}")
            return float(parts[0]), float(parts[1])
        else:
            val = float(fib_eff)
            return val, val

    def _apply_fiber_efficiency(self, spec, fib_eff: str) -> None:
        """
        Apply fiber efficiency to spectrograph by setting PyEchelle's internal cache.

        Parameters
        ----------
        spec : ZEMAX or LocalDisturber
            The spectrograph object
        fib_eff : str
            Efficiency specification ("0.9" or "0.7-0.9")
        """
        from pyechelle.efficiency import TabulatedEfficiency, SystemEfficiency

        eff_min, eff_max = self._parse_fib_eff(fib_eff)

        # Get the ZEMAX object (may be wrapped in LocalDisturber)
        if hasattr(spec, 'spectrograph'):
            zemax = spec.spectrograph
        else:
            zemax = spec

        n_fibers = self.instrument_config['n_fibers']
        ccd_index = 1

        # Get wavelength range for this band (convert nm to microns)
        wl_min_nm, wl_max_nm = BAND_WAVELENGTH_RANGES[self.config.band]
        wl_min_um = wl_min_nm / 1000.0
        wl_max_um = wl_max_nm / 1000.0
        wavelengths = np.array([wl_min_um, (wl_min_um + wl_max_um) / 2, wl_max_um])

        self.logger.info(f"Applying fiber efficiency: {fib_eff}")

        # Ensure the efficiency cache dict exists
        if not hasattr(zemax, '_efficiency') or zemax._efficiency is None:
            zemax._efficiency = {}
        if ccd_index not in zemax._efficiency:
            zemax._efficiency[ccd_index] = {}

        for fiber_num in range(1, n_fibers + 1):
            # Determine efficiency value for this fiber
            if eff_min == eff_max:
                eff_value = eff_min
            else:
                eff_value = np.random.uniform(eff_min, eff_max)

            efficiency = np.array([eff_value, eff_value, eff_value])

            # Create efficiency object and set in cache
            fiber_eff_obj = TabulatedEfficiency("FiberEff", wavelengths, efficiency)
            zemax._efficiency[ccd_index][fiber_num] = SystemEfficiency(
                [fiber_eff_obj], "System"
            )
    
    def __del__(self):
        """Ensure cleanup on deletion."""
        self.cleanup()
        
    def setup_simulator(self) -> None:
        """Set up the pyechelle simulator with instrument configuration."""
        # Get HDF model path
        hdf_model = self.config.hdf_model or 'default'
        if hdf_model != 'default' and (hdf_model.endswith('.hdf') or '/' in hdf_model):
            # Already a full path
            hdf_path = Path(hdf_model)
            if not hdf_path.is_absolute():
                hdf_path = self.project_root / hdf_path
        else:
            hdf_path = get_hdf_model_path(
                self.config.band,
                hdf_model,
                self.project_root
            )
        
        if not hdf_path.exists():
            raise FileNotFoundError(f"HDF model not found: {hdf_path}")
        
        self.logger.info(f"Using HDF model: {hdf_path}")
        
        # Create spectrograph with optional velocity shift
        if self.config.velocity_shift is not None:
            tx = self.config.velocity_shift / SPEED_OF_LIGHT
            spec = LocalDisturber(ZEMAX(str(hdf_path)), d_tx=tx)
            self.logger.info(f"Applied velocity shift: {self.config.velocity_shift} m/s")
        else:
            spec = ZEMAX(str(hdf_path))

        # Apply fiber efficiency if specified
        if self.config.fib_eff:
            self._apply_fiber_efficiency(spec, self.config.fib_eff)

        # Create simulator
        self.simulator = Simulator(spec)
        self.simulator.set_ccd(1)
        
        # Set fiber configuration
        fiber_range = range(1, self.instrument_config['n_fibers'] + 1)
        self.simulator.set_fibers(fiber_range)
        
        # Set telescope
        tel_config = self.instrument_config['telescope']
        telescope = Telescope(
            tel_config['primary_diameter'],
            tel_config['central_obstruction']
        )
        self.simulator.set_telescope(telescope)
        
        # Set computation parameters
        self.simulator.set_cuda(self.config.use_cuda)
        self.simulator.max_cpu = self.config.max_cpu
        self.simulator.set_exposure_time(self.config.exposure_time)
        
        self.logger.info(f"Simulator configured for {self.config.band}-band")
    
    def setup_sources(self) -> List[Any]:
        """
        Set up source configuration based on simulation type and fiber configuration.
        
        Returns
        -------
        List
            List of source objects for each fiber
        """
        n_fibers = self.instrument_config['n_fibers']
        illuminated_fibers = self.config.get_fiber_list()
        self.sources = self.source_factory.create_fiber_sources(
            self.config.source,
            n_fibers,
            illuminated_fibers,
            self.config.band,
            wl_min=self.config.wl_min,
            wl_max=self.config.wl_max
        )

        return self.sources
    
    def run_simulation(self, output_path: Optional[Path] = None) -> Any:
        """
        Run the complete simulation.
        
        Parameters
        ----------
        output_path : Path, optional
            Override output path. If None, uses config default.
            
        Returns
        -------
        Any
            Simulation result (usually image data)
        """
        if self.config.simulation_type == "hdf_generation":
            return self.run_hdf_generation(output_path)
        else:
            return self.run_standard_simulation(output_path)
    
    def run_hdf_generation(self, output_path: Optional[Path] = None) -> None:
        """
        Generate HDF model files from ZEMAX.
        
        Parameters
        ----------
        output_path : Path, optional
            Override output path
        """
        self.logger.info(f"Starting HDF generation for {self.config.band}-band")
        
        # Get ZEMAX file path (this would need to be configured)
        band = self.config.band
        if band not in self.instrument_config.get('zemax_files', {}):
            raise ValueError(f"No ZEMAX file configured for band {band}")
        
        zemax_file = self.instrument_config['zemax_files'][band]
        
        # Create InteractiveZEMAX instance
        zmx = InteractiveZEMAX(
            name=f'ANDES_{band}',
            zemax_filepath=zemax_file
        )
        
        # Set grating specifications
        zmx.set_grating(surface='ECHELLE', blaze=76)
        
        # Add CCD information
        detector_size = self.instrument_config['detector_size']
        pixel_size = self.instrument_config['pixel_size']
        zmx.add_ccd(1, CCD(detector_size[0], detector_size[1], pixelsize=pixel_size))
        
        # Set up fibers
        n_fibers = self.instrument_config['n_fibers']
        fiber_size = self.instrument_config['fiber_config']['fiber_size']
        
        # Fiber field positions (vertical distribution)
        y_field = (np.arange(n_fibers) - n_fibers/2 + 0.5) * fiber_size / 1000
        x_field = np.zeros(n_fibers)  # Vertical slit
        
        # Add fibers and set diffraction orders
        diffraction_orders = self.instrument_config.get('diffraction_orders')
        for i in range(n_fibers):
            zmx.add_field(
                x_field[i], y_field[i], 
                fiber_size, fiber_size, 
                shape='circular', 
                name='Science fiber'
            )
            if diffraction_orders:
                zmx.set_orders(1, i+1, diffraction_orders)
        
        # Set PSF settings
        zmx.psf_settings(
            image_delta=3,
            image_sampling="128x128",
            pupil_sampling="64x64"
        )
        
        # Generate HDF file
        if output_path is None:
            output_path = self.config.get_output_path()
        
        self.logger.info(f"Generating HDF file: {output_path}")
        hdf = HDFBuilder(zmx, str(output_path))
        hdf.save_to_hdf(n_transformation_per_order=15, n_psfs_per_order=15)
        hdf.close()
        
        self.logger.info("HDF generation completed")
    
    def run_standard_simulation(self, output_path: Optional[Path] = None) -> Any:
        """
        Run standard pyechelle simulation.
        
        Parameters
        ----------
        output_path : Path, optional
            Override output path
            
        Returns
        -------
        Any
            Simulation result
        """
        # Set up simulator and sources
        self.setup_simulator()
        
        # Standard simulation
        sources = self.setup_sources()
        self.simulator.set_sources(sources)
        
        # Set output path
        if output_path is None:
            output_path = self.config.get_output_path()
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        self.simulator.set_output(str(output_path), overwrite=self.config.output.overwrite)
        
        self.logger.info(f"Output: {output_path}")
        
        # Run simulation
        result = self.simulator.run()
        
        self.logger.info("Simulation completed")
        return result
    
    def run_single_fiber_batch(self, output_dir: Optional[Path] = None) -> Dict[int, Any]:
        """
        Run single fiber simulations for all fibers in the configuration.
        
        Parameters
        ----------
        output_dir : Path, optional
            Override output directory
            
        Returns
        -------
        Dict
            Results indexed by fiber number
        """
        if self.config.fibers.mode != "single":
            raise ValueError("Single fiber batch requires fibers.mode = 'single'")
        
        self.setup_simulator()
        results = {}
        n_fibers = self.instrument_config['n_fibers']
        
        # Get list of fibers to simulate
        if isinstance(self.config.fibers.fibers, list):
            fiber_list = self.config.fibers.fibers
        else:
            fiber_list = range(1, n_fibers + 1)
        
        for fiber_num in fiber_list:
            if fiber_num in self.config.fibers.skip_fibers:
                continue
                
            # Use source factory to create sources for this single fiber
            sources = self.source_factory.create_fiber_sources(
                self.config.source,
                n_fibers,
                [fiber_num],
                self.config.band,
                wl_min=self.config.wl_min,
                wl_max=self.config.wl_max
            )
            self.simulator.set_sources(sources)
            
            # Set output path
            if output_dir is None:
                output_path = self.config.get_output_path(fiber_num=fiber_num)
            else:
                filename = f"{self.config.band}_fiber{fiber_num:02d}_{int(self.config.exposure_time)}s.fits"
                output_path = output_dir / filename
            
            output_path.parent.mkdir(parents=True, exist_ok=True)
            self.simulator.set_output(str(output_path), overwrite=self.config.output.overwrite)
            
            # Run simulation
            result = self.simulator.run()
            results[fiber_num] = result
            
            self.logger.info(f"Fiber {fiber_num} completed: {output_path}")
        
        return results
    
    @classmethod
    def from_config_file(cls, config_path: Union[str, Path]) -> 'AndesSimulator':
        """
        Create simulator from a YAML configuration file.
        
        Parameters
        ----------
        config_path : str or Path
            Path to YAML configuration file
            
        Returns
        -------
        AndesSimulator
            Configured simulator instance
        """
        config = SimulationConfig.from_yaml(config_path)
        return cls(config)
    
    @classmethod
    def quick_flat_field(cls, band: str, fiber_mode: str = "all", **kwargs) -> 'AndesSimulator':
        """
        Quick setup for flat field simulations.
        
        Parameters
        ----------
        band : str
            Spectral band
        fiber_mode : str
            Fiber illumination mode
        **kwargs
            Additional configuration parameters
            
        Returns
        -------
        AndesSimulator
            Configured simulator
        """
        from .config import SimulationConfig, SourceConfig, FiberConfig
        
        config = SimulationConfig(
            simulation_type="flat_field",
            band=band,
            source=SourceConfig(type="constant", flux=0.001),
            fibers=FiberConfig(mode=fiber_mode),
            **kwargs
        )
        return cls(config)
    
    @classmethod
    def quick_fabry_perot(cls, band: str, **kwargs) -> 'AndesSimulator':
        """
        Quick setup for Fabry-Perot simulations.
        
        Parameters
        ----------
        band : str
            Spectral band
        **kwargs
            Additional configuration parameters
            
        Returns
        -------
        AndesSimulator
            Configured simulator
        """
        from .config import SimulationConfig, SourceConfig, FiberConfig
        
        config = SimulationConfig(
            simulation_type="fabry_perot",
            band=band,
            source=SourceConfig(type="fabry_perot", scaling_factor=5e9),
            fibers=FiberConfig(mode="all"),
            **kwargs
        )
        return cls(config)