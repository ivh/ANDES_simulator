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
from pyechelle.sources import ConstantPhotonFlux, CSVSource
from pyechelle.telescope import Telescope
from pyechelle.spectrograph import ZEMAX, LocalDisturber
from pyechelle.CCD import CCD
from pyechelle.hdfbuilder import HDFBuilder
from pyechelle.spectrograph import InteractiveZEMAX

from .config import SimulationConfig
from .instruments import get_instrument_config, get_hdf_model_path, get_sed_path
from .sources import SourceFactory, SPEED_OF_LIGHT


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
    
    def __del__(self):
        """Ensure cleanup on deletion."""
        self.cleanup()
        
    def setup_simulator(self) -> None:
        """Set up the pyechelle simulator with instrument configuration."""
        # Get HDF model path
        hdf_model = self.config.hdf_model or 'default'
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
        
        # Handle even/odd mode specially
        if self.config.fibers.mode == "even_odd":
            return self.source_factory.create_even_odd_sources(
                self.config.source,
                n_fibers,
                self.config.band
            )
        
        # Standard illumination pattern
        self.sources = self.source_factory.create_fiber_sources(
            self.config.source,
            n_fibers,
            illuminated_fibers,
            self.config.band
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
        
        # Handle even/odd mode specially
        if self.config.fibers.mode == "even_odd":
            return self.run_even_odd_simulation(output_path)
        
        # Standard simulation
        sources = self.setup_sources()
        self.simulator.set_sources(sources)
        
        # Set output path
        if output_path is None:
            output_path = self.config.get_output_path()
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        self.simulator.set_output(str(output_path), overwrite=self.config.output.overwrite)
        
        self.logger.info(f"Running {self.config.simulation_type} simulation")
        self.logger.info(f"Output: {output_path}")
        
        # Run simulation
        result = self.simulator.run()
        
        self.logger.info("Simulation completed")
        return result
    
    def run_even_odd_simulation(self, output_path: Optional[Path] = None) -> Dict[str, Any]:
        """
        Run even/odd fiber simulation (creates two output files).

        Parameters
        ----------
        output_path : Path, optional
            Base output path

        Returns
        -------
        Dict
            Results for even and odd fiber configurations
        """
        sources_dict = self.setup_sources()
        results = {}

        for mode in ['even', 'odd']:
            self.logger.info(f"Running {mode} fiber simulation")

            # Create fresh simulator for each run to avoid state contamination
            self.setup_simulator()

            # Set sources for this mode
            self.simulator.set_sources(sources_dict[mode])

            # Set output path
            if output_path is None:
                out_path = self.config.get_output_path(suffix=mode)
            else:
                out_path = output_path.parent / f"{output_path.stem}_{mode}{output_path.suffix}"

            out_path.parent.mkdir(parents=True, exist_ok=True)
            self.simulator.set_output(str(out_path), overwrite=self.config.output.overwrite)

            # Run simulation
            result = self.simulator.run()
            results[mode] = result

            self.logger.info(f"{mode.capitalize()} simulation completed: {out_path}")

        return results
    
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
                
            self.logger.info(f"Simulating fiber {fiber_num}")
            
            # Use source factory to create sources for this single fiber
            sources = self.source_factory.create_fiber_sources(
                self.config.source,
                n_fibers,
                [fiber_num],
                self.config.band
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