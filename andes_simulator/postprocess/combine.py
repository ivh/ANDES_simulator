"""
Fiber combination and summation tools for ANDES simulations.

Provides tools for combining individual fiber outputs into
integrated detector frames with various combination modes.
"""

import numpy as np
from pathlib import Path
from typing import Dict, List, Any, Optional
from astropy.io import fits
import concurrent.futures
import functools
from glob import glob

from ..core.instruments import get_instrument_config


class FiberCombiner:
    """
    Handles combination of individual fiber simulation outputs.

    Provides tools for summing fiber outputs with various selection
    criteria and weighting schemes for integrated detector analysis.
    """

    # Headers to propagate from input files to combined output
    PROPAGATE_HEADERS = ['HDFMODEL', 'SIMTYPE', 'SRCTYPE', 'SRCFLUX', 'SRCSCALE',
                         'EXPTIME', 'VSHIFT', 'FIBEFF', 'WL_MIN', 'WL_MAX']

    def __init__(self,
                 band: str,
                 project_root: Optional[Path] = None,
                 input_dir: Optional[Path] = None):
        """
        Initialize fiber combiner.

        Parameters
        ----------
        band : str
            Spectral band for fiber and detector configuration
        project_root : Path, optional
            Project root directory
        input_dir : Path, optional
            Directory to search for input files (defaults to band output dir)
        """
        self.band = band
        self.instrument_config = get_instrument_config(band)

        if project_root is None:
            self.project_root = Path(__file__).parent.parent.parent
        else:
            self.project_root = project_root

        # Input directory for finding fiber files
        if input_dir is None:
            self.input_dir = self.project_root.parent / band
        else:
            self.input_dir = Path(input_dir)

        # Get band-specific parameters
        self.n_fibers = self.instrument_config['n_fibers']
        # Config stores (X, Y) but numpy/FITS uses (Y, X) - swap for array operations
        config_size = self.instrument_config['detector_size']
        self.detector_size = (config_size[1], config_size[0])
        self.skip_fibers = self.instrument_config.get('skip_fibers', [])
        self._source_header = None  # Will store headers from first loaded file
    
    def load_fiber_data(self, 
                       fiber_num: int,
                       input_pattern: str) -> Optional[np.ndarray]:
        """
        Load data for a single fiber.
        
        Parameters
        ----------
        fiber_num : int
            1-based fiber number
        input_pattern : str
            File pattern with {fib} placeholder
            
        Returns
        -------
        np.ndarray or None
            Fiber image data, or None if not found/skipped
        """
        if fiber_num in self.skip_fibers:
            return np.zeros(self.detector_size)
        
        # Create search pattern
        pattern = input_pattern.format(fib=fiber_num)
        if not Path(pattern).is_absolute():
            pattern = str(self.input_dir / pattern)
        
        matching_files = glob(pattern)
        
        if not matching_files:
            print(f"Warning: No files found for fiber {fiber_num} with pattern: {pattern}")
            return np.zeros(self.detector_size)
        
        # Load the first matching file
        try:
            with fits.open(matching_files[0]) as hdul:
                # Capture headers from first successfully loaded file
                if self._source_header is None:
                    self._source_header = {k: hdul[0].header.get(k)
                                           for k in self.PROPAGATE_HEADERS
                                           if k in hdul[0].header}
                return hdul[0].data.copy()
        except Exception as e:
            print(f"Error loading fiber {fiber_num} data: {e}")
            return np.zeros(self.detector_size)
    
    def combine_all_fibers(self,
                          input_pattern: str,
                          weights: Optional[List[float]] = None,
                          max_workers: int = 6) -> np.ndarray:
        """
        Combine all fiber outputs with optional weighting.
        
        Parameters
        ----------
        input_pattern : str
            File pattern for fiber files
        weights : List[float], optional
            Weights for each fiber (must match n_fibers)
        max_workers : int
            Maximum parallel workers
            
        Returns
        -------
        np.ndarray
            Combined detector image
        """
        if weights is not None and len(weights) != self.n_fibers:
            raise ValueError(f"Weights length ({len(weights)}) must match n_fibers ({self.n_fibers})")
        
        combined_image = np.zeros(self.detector_size)
        
        # Create partial function for parallel loading
        load_func = functools.partial(self.load_fiber_data, input_pattern=input_pattern)
        
        # Load all fiber data in parallel
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            fiber_nums = range(1, self.n_fibers + 1)
            fiber_data_list = list(executor.map(load_func, fiber_nums))
        
        # Combine with weights
        for i, fiber_data in enumerate(fiber_data_list):
            if fiber_data is not None:
                weight = weights[i] if weights is not None else 1.0
                combined_image += weight * fiber_data
        
        return combined_image
    
    def combine_fiber_subset(self,
                           fiber_nums: List[int],
                           input_pattern: str,
                           weights: Optional[List[float]] = None) -> np.ndarray:
        """
        Combine specific subset of fibers.
        
        Parameters
        ----------
        fiber_nums : List[int]
            List of 1-based fiber numbers to combine
        input_pattern : str
            File pattern for fiber files
        weights : List[float], optional
            Weights for each specified fiber
            
        Returns
        -------
        np.ndarray
            Combined detector image
        """
        if weights is not None and len(weights) != len(fiber_nums):
            raise ValueError("Weights length must match fiber_nums length")
        
        combined_image = np.zeros(self.detector_size)
        
        for i, fiber_num in enumerate(fiber_nums):
            fiber_data = self.load_fiber_data(fiber_num, input_pattern)
            if fiber_data is not None:
                weight = weights[i] if weights is not None else 1.0
                combined_image += weight * fiber_data
        
        return combined_image
    
    def combine_even_odd_fibers(self, 
                               input_pattern: str) -> Dict[str, np.ndarray]:
        """
        Combine fibers separately for even and odd fiber numbers.
        
        Parameters
        ----------
        input_pattern : str
            File pattern for fiber files
            
        Returns
        -------
        Dict
            Dictionary with 'even' and 'odd' combined images
        """
        even_fibers = [i for i in range(1, self.n_fibers + 1) if i % 2 == 0]
        odd_fibers = [i for i in range(1, self.n_fibers + 1) if i % 2 == 1]
        
        results = {}
        results['even'] = self.combine_fiber_subset(even_fibers, input_pattern)
        results['odd'] = self.combine_fiber_subset(odd_fibers, input_pattern)
        
        return results
    
    def combine_pseudo_slits(self, 
                            input_pattern: str) -> Dict[str, np.ndarray]:
        """
        Combine fibers by pseudo-slit for optical bands.
        
        Parameters
        ----------
        input_pattern : str
            File pattern for fiber files
            
        Returns
        -------
        Dict
            Dictionary with slit-specific combined images
        """
        if self.band in ['Y', 'J', 'H']:
            raise ValueError("Pseudo-slit combination not applicable to YJH bands")
        
        # Get fiber configuration for optical bands
        fiber_config = self.instrument_config['fiber_config']
        slits = fiber_config.get('slits', {})

        results = {}

        if 'slitA' in slits:
            results['first_slit'] = self.combine_fiber_subset(slits['slitA'], input_pattern)

        if 'slitB' in slits:
            results['second_slit'] = self.combine_fiber_subset(slits['slitB'], input_pattern)

        if 'cal_fibers' in slits:
            results['calibration'] = self.combine_fiber_subset(slits['cal_fibers'], input_pattern)
        
        return results
    
    def create_fiber_map(self, 
                        input_pattern: str,
                        output_dir: Optional[Path] = None) -> Path:
        """
        Create a FITS file with each fiber in a separate extension.
        
        Parameters
        ----------
        input_pattern : str
            File pattern for fiber files
        output_dir : Path, optional
            Output directory
            
        Returns
        -------
        Path
            Path to created fiber map FITS file
        """
        if output_dir is None:
            output_dir = self.project_root.parent / self.band
        
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{self.band}_fiber_map.fits"
        
        # Create primary HDU
        primary = fits.PrimaryHDU()
        primary.header['BAND'] = self.band
        primary.header['NFIBERS'] = self.n_fibers
        primary.header['DETSIZE'] = f"{self.detector_size[0]}x{self.detector_size[1]}"
        
        hdul = fits.HDUList([primary])
        
        # Add each fiber as an extension
        for fiber_num in range(1, self.n_fibers + 1):
            fiber_data = self.load_fiber_data(fiber_num, input_pattern)
            
            if fiber_data is not None:
                hdu = fits.ImageHDU(fiber_data, name=f'FIBER_{fiber_num:02d}')
                hdu.header['FIBER'] = fiber_num
                hdu.header['SKIPPED'] = fiber_num in self.skip_fibers
            else:
                # Create empty extension for missing fibers
                hdu = fits.ImageHDU(np.zeros(self.detector_size), name=f'FIBER_{fiber_num:02d}')
                hdu.header['FIBER'] = fiber_num
                hdu.header['SKIPPED'] = True
                hdu.header['MISSING'] = True
            
            hdul.append(hdu)
        
        # Write to file
        hdul.writeto(str(output_path), overwrite=True)
        print(f"Created fiber map: {output_path}")
        
        return output_path
    
    def save_combined_image(self,
                           image_data: np.ndarray,
                           output_path: Path,
                           combination_info: Optional[Dict] = None) -> None:
        """
        Save combined image with metadata.
        
        Parameters
        ----------
        image_data : np.ndarray
            Combined image data
        output_path : Path
            Output file path
        combination_info : Dict, optional
            Information about combination process
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Create HDU with metadata
        hdu = fits.PrimaryHDU(image_data)
        hdu.header['BAND'] = self.band
        hdu.header['NFIBERS'] = self.n_fibers
        hdu.header['DETSIZE'] = f"{self.detector_size[0]}x{self.detector_size[1]}"
        hdu.header['COMBINED'] = True

        # Propagate headers from source files
        if self._source_header:
            for key, value in self._source_header.items():
                if value is not None:
                    hdu.header[key] = value

        if combination_info:
            for key, value in combination_info.items():
                if isinstance(value, (int, float, str, bool)):
                    hdu.header[key.upper()[:8]] = value  # FITS keyword limit

        if self.skip_fibers:
            hdu.header['SKIPFIB'] = ','.join(map(str, self.skip_fibers))
        
        # Write to file
        hdu.writeto(str(output_path), overwrite=True)
        print(f"Saved combined image: {output_path}")
    
    def analyze_fiber_contributions(self, 
                                   input_pattern: str) -> Dict[str, Any]:
        """
        Analyze relative contributions of each fiber.
        
        Parameters
        ----------
        input_pattern : str
            File pattern for fiber files
            
        Returns
        -------
        Dict
            Analysis results including flux statistics
        """
        fiber_stats = {}
        total_flux = 0.0
        
        for fiber_num in range(1, self.n_fibers + 1):
            fiber_data = self.load_fiber_data(fiber_num, input_pattern)
            
            if fiber_data is not None:
                flux = np.sum(fiber_data)
                peak = np.max(fiber_data)
                
                fiber_stats[fiber_num] = {
                    'total_flux': flux,
                    'peak_value': peak,
                    'mean_value': np.mean(fiber_data),
                    'skipped': fiber_num in self.skip_fibers
                }
                
                if fiber_num not in self.skip_fibers:
                    total_flux += flux
        
        # Calculate relative contributions
        for fiber_num in fiber_stats:
            if not fiber_stats[fiber_num]['skipped']:
                fiber_stats[fiber_num]['flux_fraction'] = (
                    fiber_stats[fiber_num]['total_flux'] / total_flux if total_flux > 0 else 0
                )
        
        analysis = {
            'fiber_stats': fiber_stats,
            'total_flux': total_flux,
            'n_active_fibers': len([f for f in range(1, self.n_fibers + 1) if f not in self.skip_fibers]),
            'n_skipped_fibers': len(self.skip_fibers),
            'band': self.band
        }
        
        return analysis
    
    def create_combination_report(self,
                                 input_pattern: str,
                                 output_dir: Optional[Path] = None) -> Path:
        """
        Create a detailed report of fiber combination analysis.
        
        Parameters
        ----------
        input_pattern : str
            File pattern for fiber files
        output_dir : Path, optional
            Output directory
            
        Returns
        -------
        Path
            Path to created report file
        """
        if output_dir is None:
            output_dir = self.project_root.parent / self.band
        
        output_dir.mkdir(parents=True, exist_ok=True)
        report_path = output_dir / f"{self.band}_combination_report.txt"
        
        # Analyze fiber contributions
        analysis = self.analyze_fiber_contributions(input_pattern)
        
        # Write report
        with open(report_path, 'w') as f:
            f.write(f"ANDES {self.band}-band Fiber Combination Report\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Band: {self.band}\n")
            f.write(f"Total fibers: {self.n_fibers}\n")
            f.write(f"Active fibers: {analysis['n_active_fibers']}\n")
            f.write(f"Skipped fibers: {analysis['n_skipped_fibers']}\n")
            f.write(f"Total flux: {analysis['total_flux']:.2e}\n\n")
            
            if self.skip_fibers:
                f.write(f"Skipped fiber numbers: {self.skip_fibers}\n\n")
            
            f.write("Fiber Statistics:\n")
            f.write("-" * 30 + "\n")
            f.write("Fiber  Total Flux     Peak Value    Flux Fraction\n")
            f.write("-" * 50 + "\n")
            
            for fiber_num in sorted(analysis['fiber_stats'].keys()):
                stats = analysis['fiber_stats'][fiber_num]
                if stats['skipped']:
                    f.write(f"{fiber_num:5d}  SKIPPED\n")
                else:
                    f.write(f"{fiber_num:5d}  {stats['total_flux']:10.2e}  "
                           f"{stats['peak_value']:10.2e}  {stats['flux_fraction']:8.4f}\n")
        
        print(f"Created combination report: {report_path}")
        return report_path
    
    @classmethod
    def quick_combine(cls,
                     band: str,
                     input_pattern: str,
                     output_filename: Optional[str] = None,
                     project_root: Optional[Path] = None) -> Path:
        """
        Quick combination of all fibers with default settings.
        
        Parameters
        ----------
        band : str
            Spectral band
        input_pattern : str
            File pattern for fiber files
        output_filename : str, optional
            Output filename
        project_root : Path, optional
            Project root directory
            
        Returns
        -------
        Path
            Path to combined output file
        """
        combiner = cls(band, project_root)
        
        # Combine all fibers
        combined_image = combiner.combine_all_fibers(input_pattern)
        
        # Generate output path
        if output_filename is None:
            output_filename = f"{band}_combined.fits"
        
        output_dir = combiner.project_root.parent / band
        output_path = output_dir / output_filename
        
        # Save result
        combiner.save_combined_image(combined_image, output_path, 
                                    {'combination_mode': 'all_fibers'})
        
        return output_path