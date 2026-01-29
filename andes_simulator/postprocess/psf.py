"""
PSF convolution and processing for ANDES simulations.

Provides tools for applying PSF convolution kernels to simulation outputs
with various edge-blanking effects and kernel configurations.
"""

import numpy as np
import random
from datetime import datetime
from pathlib import Path
from typing import Optional, Tuple
from astropy.io import fits
from scipy.signal import convolve2d
import concurrent.futures
import functools

from ..core.instruments import get_instrument_config


class PSFProcessor:
    """
    Handles PSF convolution processing for ANDES simulation outputs.

    Provides customizable 2D Gaussian kernels with edge-blanking effects
    and supports parallel processing of multiple fiber outputs.
    """

    # Headers to propagate from input files
    PROPAGATE_HEADERS = ['HDFMODEL', 'SIMTYPE', 'SRCTYPE', 'SRCFLUX', 'SRCSCALE',
                         'EXPTIME', 'VSHIFT', 'FIBEFF', 'WL_MIN', 'WL_MAX']

    def __init__(self,
                 band: str,
                 project_root: Optional[Path] = None):
        """
        Initialize PSF processor.

        Parameters
        ----------
        band : str
            Spectral band for detector size and sampling configuration
        project_root : Path, optional
            Project root directory
        """
        self.band = band
        self.instrument_config = get_instrument_config(band)

        if project_root is None:
            self.project_root = Path(__file__).parent.parent.parent
        else:
            self.project_root = project_root

        # Get band-specific parameters
        self.detector_size = self.instrument_config['detector_size']
        self.sampling = self.instrument_config.get('sampling', 1.0)
        self.skip_fibers = self.instrument_config.get('skip_fibers', [])
        self._source_header = None  # Will store headers from first loaded file
    
    def create_gaussian_kernel(self,
                              dimx: int = 4, 
                              dimy: int = 4,
                              sigma: float = 1.0,
                              edge_blank: str = 'random') -> np.ndarray:
        """
        Create a 2D Gaussian kernel with optional edge blanking.
        
        Parameters
        ----------
        dimx : int
            Kernel width
        dimy : int
            Kernel height
        sigma : float
            Gaussian standard deviation
        edge_blank : str
            Edge to zero out ('none', 'top', 'bottom', 'left', 'right', 'random')
            
        Returns
        -------
        np.ndarray
            2D Gaussian kernel with dimensions (dimy, dimx)
        """
        # Create meshgrid for x and y coordinates, centered for even-sized kernels
        y, x = np.ogrid[:dimy, :dimx]
        y = y - (dimy - 1) / 2
        x = x - (dimx - 1) / 2
        
        # Calculate 2D Gaussian
        kernel = np.exp(-(x**2 + y**2) / (2 * sigma**2))
        
        # Apply edge blanking
        if edge_blank == 'random':
            edge = random.choice(['none', 'top', 'bottom', 'left', 'right'])
        else:
            edge = edge_blank
        
        if edge == 'top':
            kernel[0, :] = 0
        elif edge == 'bottom':
            kernel[-1, :] = 0
        elif edge == 'left':
            kernel[:, 0] = 0
        elif edge == 'right':
            kernel[:, -1] = 0
        # 'none' case: no blanking
        
        # Normalize after zeroing
        if np.sum(kernel) > 0:
            kernel = kernel / np.sum(kernel)
        
        return kernel
    
    def convolve_image(self, 
                      image: np.ndarray,
                      kernel: np.ndarray) -> np.ndarray:
        """
        Convolve image with PSF kernel.
        
        Parameters
        ----------
        image : np.ndarray
            Input image data
        kernel : np.ndarray
            PSF kernel
            
        Returns
        -------
        np.ndarray
            Convolved image
        """
        return convolve2d(image, kernel, mode='same')
    
    def process_single_fiber(self, 
                            fiber_num: int,
                            input_pattern: str,
                            kernel_params: Optional[Tuple] = None) -> Optional[np.ndarray]:
        """
        Process a single fiber file with optional PSF convolution.
        
        Parameters
        ----------
        fiber_num : int
            1-based fiber number
        input_pattern : str
            File pattern for finding fiber files (with {fib} placeholder)
        kernel_params : Tuple, optional
            Kernel parameters (dimx, dimy, sigma, edge_blank)
            
        Returns
        -------
        np.ndarray or None
            Processed image data, or None if fiber skipped/not found
        """
        if fiber_num in self.skip_fibers:
            return np.zeros(self.detector_size)
        
        # Find matching files
        from glob import glob
        pattern = str(self.project_root.parent / input_pattern.format(fib=fiber_num))
        matching_files = glob(pattern)
        
        if not matching_files:
            print(f"Warning: No files found for fiber {fiber_num} with pattern: {pattern}")
            return np.zeros(self.detector_size)
        
        # Load image data
        with fits.open(matching_files[0]) as hdul:
            image_data = hdul[0].data
            # Capture headers from first successfully loaded file
            if self._source_header is None:
                self._source_header = {k: hdul[0].header.get(k)
                                       for k in self.PROPAGATE_HEADERS
                                       if k in hdul[0].header}

        # Apply PSF convolution if requested
        if kernel_params is not None:
            kernel = self.create_gaussian_kernel(*kernel_params)
            image_data = self.convolve_image(image_data, kernel)
        
        return image_data
    
    def process_all_fibers(self,
                          input_pattern: str,
                          kernel_params: Optional[Tuple] = None,
                          max_workers: int = 6) -> np.ndarray:
        """
        Process and sum all fiber outputs with optional PSF convolution.
        
        Parameters
        ----------
        input_pattern : str
            File pattern for finding fiber files (with {fib} placeholder)
        kernel_params : Tuple, optional
            Kernel parameters (dimx, dimy, sigma, edge_blank)
        max_workers : int
            Maximum number of parallel workers
            
        Returns
        -------
        np.ndarray
            Combined image from all processed fibers
        """
        n_fibers = self.instrument_config['n_fibers']
        combined_image = np.zeros(self.detector_size)
        
        # Create partial function for parallel processing
        if kernel_params is not None:
            process_func = functools.partial(
                self._process_fiber_with_kernel,
                input_pattern=input_pattern,
                kernel_params=kernel_params
            )
        else:
            process_func = functools.partial(
                self._process_fiber_without_kernel,
                input_pattern=input_pattern
            )
        
        # Process fibers in parallel
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            fiber_nums = range(1, n_fibers + 1)
            for image_data in executor.map(process_func, fiber_nums):
                if image_data is not None:
                    combined_image += image_data
        
        return combined_image
    
    def _process_fiber_with_kernel(self, 
                                  fiber_num: int,
                                  input_pattern: str,
                                  kernel_params: Tuple) -> Optional[np.ndarray]:
        """Helper method for parallel processing with kernel."""
        return self.process_single_fiber(fiber_num, input_pattern, kernel_params)
    
    def _process_fiber_without_kernel(self, 
                                     fiber_num: int,
                                     input_pattern: str) -> Optional[np.ndarray]:
        """Helper method for parallel processing without kernel."""
        return self.process_single_fiber(fiber_num, input_pattern, None)
    
    def save_processed_image(self, 
                            image_data: np.ndarray,
                            output_path: Path,
                            kernel_params: Optional[Tuple] = None) -> None:
        """
        Save processed image to FITS file.
        
        Parameters
        ----------
        image_data : np.ndarray
            Image data to save
        output_path : Path
            Output file path
        kernel_params : Tuple, optional
            Kernel parameters for filename generation
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Create FITS HDU
        hdu = fits.PrimaryHDU(image_data)
        
        # Add kernel information to header if applicable
        if kernel_params is not None:
            dimx, dimy, sigma, edge_blank = kernel_params
            hdu.header['PSFDIMX'] = dimx
            hdu.header['PSFDIMY'] = dimy  
            hdu.header['PSFSIGMA'] = sigma
            hdu.header['PSFBLANK'] = edge_blank
            hdu.header['PSFCONV'] = True
        else:
            hdu.header['PSFCONV'] = False
        
        hdu.header['DATE-OBS'] = (datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S'), 'Observation date')
        hdu.header['INSTRUME'] = ('ANDES', 'Instrument name')
        hdu.header['BAND'] = self.band
        hdu.header['DETSIZE'] = f"{self.detector_size[0]}x{self.detector_size[1]}"

        # Propagate headers from source files
        if self._source_header:
            for key, value in self._source_header.items():
                if value is not None:
                    hdu.header[key] = value

        # Write to file
        hdu.writeto(str(output_path), overwrite=True)
        print(f"Saved processed image: {output_path}")
    
    def generate_output_filename(self,
                                kernel_params: Optional[Tuple] = None,
                                suffix: str = "sum") -> str:
        """
        Generate appropriate output filename.
        
        Parameters
        ----------
        kernel_params : Tuple, optional
            Kernel parameters for filename
        suffix : str
            Base suffix for filename
            
        Returns
        -------
        str
            Generated filename
        """
        if kernel_params is not None:
            dimx, dimy, sigma, edge_blank = kernel_params
            # Convert sigma from pixels back to FWHM in arcsec
            fwhm = sigma * 2.35 / self.sampling
            filename = f"{self.band}_FP_{suffix}_kern{dimx}x{dimy}s{fwhm:.1f}{edge_blank}.fits"
        else:
            filename = f"{self.band}_FP_{suffix}.fits"
        
        return filename
    
    def process_fabry_perot_files(self,
                                 kernel_params: Optional[Tuple] = None,
                                 output_dir: Optional[Path] = None,
                                 input_pattern: Optional[str] = None) -> Path:
        """
        Process Fabry-Perot simulation files with PSF convolution.
        
        Parameters
        ----------
        kernel_params : Tuple, optional
            Kernel parameters (dimx, dimy, fwhm_arcsec, edge_blank)
            FWHM should be in arcseconds and will be converted to sigma in pixels
        output_dir : Path, optional
            Output directory
        input_pattern : str, optional
            Input file pattern
            
        Returns
        -------
        Path
            Path to output file
        """
        if output_dir is None:
            output_dir = self.project_root.parent / self.band
        
        if input_pattern is None:
            input_pattern = f"{self.band}/{self.band}_FP_fiber{{fib:02d}}_shift*.fits"
        
        # Convert FWHM to sigma in pixels if kernel params provided
        if kernel_params is not None:
            dimx, dimy, fwhm_arcsec, edge_blank = kernel_params
            sigma_pixels = fwhm_arcsec * self.sampling / 2.35
            kernel_params = (dimx, dimy, sigma_pixels, edge_blank)
        
        # Process all fibers
        combined_image = self.process_all_fibers(input_pattern, kernel_params)
        
        # Generate output filename and path
        filename = self.generate_output_filename(kernel_params, "IFUsum")
        output_path = output_dir / filename
        
        # Save result
        self.save_processed_image(combined_image, output_path, kernel_params)
        
        return output_path
    
    @classmethod
    def quick_psf_process(cls,
                         band: str,
                         kernel_size: Tuple[int, int] = (4, 4),
                         fwhm_arcsec: float = 3.2,
                         edge_blank: str = 'left',
                         project_root: Optional[Path] = None) -> Path:
        """
        Quick PSF processing with standard parameters.
        
        Parameters
        ----------
        band : str
            Spectral band
        kernel_size : Tuple[int, int]
            Kernel dimensions (dimx, dimy)
        fwhm_arcsec : float
            FWHM in arcseconds
        edge_blank : str
            Edge blanking mode
        project_root : Path, optional
            Project root directory
            
        Returns
        -------
        Path
            Path to output file
        """
        processor = cls(band, project_root)
        kernel_params = (*kernel_size, fwhm_arcsec, edge_blank)
        return processor.process_fabry_perot_files(kernel_params)
    
    def create_kernel_visualization(self, 
                                   kernel_params: Tuple,
                                   output_path: Optional[Path] = None) -> Optional[Path]:
        """
        Create and save visualization of PSF kernel.
        
        Parameters
        ----------
        kernel_params : Tuple
            Kernel parameters (dimx, dimy, sigma, edge_blank)
        output_path : Path, optional
            Output path for visualization
            
        Returns
        -------
        Path or None
            Path to saved visualization, or None if matplotlib not available
        """
        try:
            import matplotlib.pyplot as plt
            
            kernel = self.create_gaussian_kernel(*kernel_params)
            
            plt.figure(figsize=(6, 6))
            plt.imshow(kernel, interpolation='nearest', cmap='viridis')
            plt.colorbar(label='Kernel Value')
            plt.title(f"PSF Kernel: {kernel_params[0]}x{kernel_params[1]}, σ={kernel_params[2]:.2f}, edge={kernel_params[3]}")
            
            if output_path is None:
                output_path = self.project_root / f"psf_kernel_{self.band}.png"
            
            plt.savefig(str(output_path), dpi=150, bbox_inches='tight')
            plt.close()
            
            return output_path
            
        except ImportError:
            print("Matplotlib not available for kernel visualization")
            return None