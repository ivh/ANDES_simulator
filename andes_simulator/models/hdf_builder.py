"""
HDF model generation for ANDES spectrograph.

Provides automated generation of HDF instrument models from ZEMAX
optical designs with proper fiber field configurations.
"""

import logging
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional

from pyechelle.CCD import CCD
from pyechelle.hdfbuilder import HDFBuilder
from pyechelle.spectrograph import InteractiveZEMAX

from ..core.instruments import get_instrument_config


class AndesHDFBuilder:
    """
    Manages HDF model generation for ANDES spectrograph bands.
    
    Automates the process of creating HDF instrument models from ZEMAX
    optical designs with proper fiber field configurations and diffraction
    order specifications.
    """
    
    def __init__(self, 
                 band: str,
                 project_root: Optional[Path] = None,
                 zemax_data_dir: Optional[Path] = None):
        """
        Initialize HDF builder for a specific band.
        
        Parameters
        ----------
        band : str
            Spectral band (Y, J, H, R, IZ, U, B, V)
        project_root : Path, optional
            Project root directory
        zemax_data_dir : Path, optional
            Directory containing ZEMAX files
        """
        self.band = band
        self.instrument_config = get_instrument_config(band)
        
        if project_root is None:
            self.project_root = Path(__file__).parent.parent.parent
        else:
            self.project_root = project_root
        
        if zemax_data_dir is None:
            self.zemax_data_dir = self.project_root / "ZEMAX"
        else:
            self.zemax_data_dir = zemax_data_dir
        
        self.logger = logging.getLogger(__name__)
        
        # HDF builder instance
        self.hdf_builder = None
        self.zmx = None
    
    def setup_zemax_interface(self, zemax_file: Optional[str] = None) -> None:
        """
        Set up the InteractiveZEMAX interface.
        
        Parameters
        ----------
        zemax_file : str, optional
            ZEMAX file name. If None, uses default from instrument config.
        """
        if zemax_file is None:
            if 'zemax_files' not in self.instrument_config:
                raise ValueError(f"No ZEMAX file configured for band {self.band}")
            zemax_file = self.instrument_config['zemax_files'][self.band]
        
        zemax_path = self.zemax_data_dir / zemax_file
        if not zemax_path.exists():
            self.logger.warning(f"ZEMAX file not found: {zemax_path}")
            self.logger.info("Proceeding with filename only (assuming ZEMAX will find it)")
            zemax_path = zemax_file
        
        self.logger.info(f"Initializing ZEMAX interface: {zemax_path}")
        
        # Create InteractiveZEMAX instance
        self.zmx = InteractiveZEMAX(
            name=f'ANDES_{self.band}',
            zemax_filepath=str(zemax_path)
        )
        
        self.logger.info("ZEMAX interface initialized")
    
    def configure_grating(self, blaze_angle: float = 76.0) -> None:
        """
        Configure echelle grating specifications.
        
        Parameters
        ----------
        blaze_angle : float
            Grating blaze angle in degrees
        """
        if self.zmx is None:
            raise RuntimeError("ZEMAX interface not initialized")
        
        self.logger.info(f"Setting grating blaze angle: {blaze_angle}°")
        self.zmx.set_grating(surface='ECHELLE', blaze=blaze_angle)
    
    def configure_detector(self) -> None:
        """Configure CCD detector specifications."""
        if self.zmx is None:
            raise RuntimeError("ZEMAX interface not initialized")
        
        detector_size = self.instrument_config['detector_size']
        pixel_size = self.instrument_config['pixel_size']
        
        self.logger.info(f"Configuring detector: {detector_size[0]}x{detector_size[1]} pixels, {pixel_size}μm pixel size")
        
        ccd = CCD(detector_size[0], detector_size[1], pixelsize=pixel_size)
        self.zmx.add_ccd(1, ccd)
    
    def configure_fibers(self) -> None:
        """Configure fiber field positions and diffraction orders."""
        if self.zmx is None:
            raise RuntimeError("ZEMAX interface not initialized")
        
        n_fibers = self.instrument_config['n_fibers']
        fiber_size = self.instrument_config['fiber_config']['fiber_size']
        
        self.logger.info(f"Configuring {n_fibers} fibers with {fiber_size}μm size")
        
        # Calculate fiber field positions (vertical distribution)
        y_field = (np.arange(n_fibers) - n_fibers/2 + 0.5) * fiber_size / 1000
        x_field = np.zeros(n_fibers)  # Vertical slit configuration
        
        # Add fibers to ZEMAX
        for i in range(n_fibers):
            self.zmx.add_field(
                x_field[i], y_field[i],
                fiber_size, fiber_size,
                shape='circular',
                name='Science fiber'
            )
            
            # Set diffraction orders if available
            if 'diffraction_orders' in self.instrument_config:
                orders = self.instrument_config['diffraction_orders']
                self.zmx.set_orders(1, i+1, orders)
                self.logger.debug(f"Fiber {i+1}: orders {orders[0]}-{orders[-1]}")
        
        self.logger.info("Fiber configuration completed")
    
    def configure_psf_settings(self, 
                              image_delta: int = 3,
                              image_sampling: str = "128x128",
                              pupil_sampling: str = "64x64") -> None:
        """
        Configure PSF calculation settings.
        
        Parameters
        ----------
        image_delta : int
            Image delta parameter for PSF calculation
        image_sampling : str
            Image plane sampling
        pupil_sampling : str
            Pupil plane sampling
        """
        if self.zmx is None:
            raise RuntimeError("ZEMAX interface not initialized")
        
        self.logger.info(f"Configuring PSF settings: {image_sampling} image, {pupil_sampling} pupil")
        
        self.zmx.psf_settings(
            image_delta=image_delta,
            image_sampling=image_sampling,
            pupil_sampling=pupil_sampling
        )
    
    def build_hdf(self, 
                  output_path: Optional[Path] = None,
                  n_transformation_per_order: int = 15,
                  n_psfs_per_order: int = 15) -> None:
        """
        Build the HDF model file.
        
        Parameters
        ----------
        output_path : Path, optional
            Output path for HDF file. If None, uses default naming.
        n_transformation_per_order : int
            Number of transformation samples per diffraction order
        n_psfs_per_order : int
            Number of PSF samples per diffraction order
        """
        if self.zmx is None:
            raise RuntimeError("ZEMAX interface not initialized")
        
        # Determine output path
        if output_path is None:
            if self.band in ['Y', 'J', 'H']:
                filename = f'ANDES_75fibre_{self.band}.hdf'
            else:
                filename = f'ANDES_123_{self.band}3.hdf'
            output_path = self.project_root / 'HDF' / filename
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        self.logger.info(f"Building HDF model: {output_path}")
        self.logger.info(f"Samples per order: {n_transformation_per_order} transformations, {n_psfs_per_order} PSFs")
        
        # Create HDF builder
        self.hdf_builder = HDFBuilder(self.zmx, str(output_path))
        
        # Build the model (this is the time-consuming step)
        self.logger.warning("HDF generation starting - this will take a long time...")
        self.hdf_builder.save_to_hdf(
            n_transformation_per_order=n_transformation_per_order,
            n_psfs_per_order=n_psfs_per_order
        )
        
        # Close the builder
        self.hdf_builder.close()
        self.logger.info("HDF model generation completed")
    
    def build_complete_model(self, 
                            output_path: Optional[Path] = None,
                            zemax_file: Optional[str] = None,
                            **build_kwargs) -> Path:
        """
        Complete HDF model generation process.
        
        Parameters
        ----------
        output_path : Path, optional
            Output path for HDF file
        zemax_file : str, optional
            ZEMAX file to use
        **build_kwargs
            Additional arguments for HDF building
            
        Returns
        -------
        Path
            Path to generated HDF file
        """
        self.logger.info(f"Starting complete HDF model generation for {self.band}-band")
        
        try:
            # Set up ZEMAX interface
            self.setup_zemax_interface(zemax_file)
            
            # Configure all components
            self.configure_grating()
            self.configure_detector()
            self.configure_fibers()
            self.configure_psf_settings()
            
            # Build HDF model
            self.build_hdf(output_path, **build_kwargs)
            
            # Determine final output path
            if output_path is None:
                if self.band in ['Y', 'J', 'H']:
                    filename = f'ANDES_75fibre_{self.band}.hdf'
                else:
                    filename = f'ANDES_123_{self.band}3.hdf'
                output_path = self.project_root / 'HDF' / filename
            
            self.logger.info(f"HDF model generation completed: {output_path}")
            return output_path
            
        except Exception as e:
            self.logger.error(f"HDF model generation failed: {e}")
            raise
    
    @classmethod
    def build_all_bands(cls, 
                       bands: Optional[List[str]] = None,
                       project_root: Optional[Path] = None,
                       **build_kwargs) -> Dict[str, Path]:
        """
        Build HDF models for multiple bands.
        
        Parameters
        ----------
        bands : List[str], optional
            List of bands to build. If None, builds all configured bands.
        project_root : Path, optional
            Project root directory
        **build_kwargs
            Additional arguments for HDF building
            
        Returns
        -------
        Dict[str, Path]
            Dictionary mapping band names to generated HDF file paths
        """
        if bands is None:
            bands = ['Y', 'J', 'H', 'R', 'IZ', 'U', 'B', 'V']
        
        results = {}
        logger = logging.getLogger(__name__)
        
        for band in bands:
            logger.info(f"Building HDF model for {band}-band")
            try:
                builder = cls(band, project_root)
                output_path = builder.build_complete_model(**build_kwargs)
                results[band] = output_path
            except Exception as e:
                logger.error(f"Failed to build {band}-band model: {e}")
                results[band] = None
        
        return results