"""
Flat field source configurations for ANDES calibrations.

Provides different illumination patterns for flat field calibrations
including single fiber, even/odd patterns, and slit-specific illuminations.
"""

from typing import List, Dict, Any, Optional
from pyechelle.sources import Constant


class FlatFieldSource:
    """
    Manages flat field source configurations for ANDES calibrations.
    
    Supports various illumination patterns:
    - Single fiber illumination for fiber response mapping
    - Even/odd fiber patterns for cross-talk analysis
    - Slit-specific illumination for pseudo-slit characterization
    - Full illumination for detector response
    """
    
    def __init__(self, flux_level: float = 0.001):
        """
        Initialize flat field source.
        
        Parameters
        ----------
        flux_level : float
            Constant flux level for illuminated fibers (photons/s)
        """
        self.flux_level = flux_level
        self.illuminated_source = Constant(flux_level)
        self.dark_source = Constant(0.0)
    
    def get_single_fiber_sources(self, fiber_num: int, n_fibers: int) -> List[Any]:
        """
        Get source configuration for single fiber illumination.
        
        Parameters
        ----------
        fiber_num : int
            1-based fiber number to illuminate
        n_fibers : int
            Total number of fibers
            
        Returns
        -------
        List
            Source list with only specified fiber illuminated
        """
        if not (1 <= fiber_num <= n_fibers):
            raise ValueError(f"Fiber number {fiber_num} out of range [1, {n_fibers}]")
        
        sources = [self.dark_source] * n_fibers
        sources[fiber_num - 1] = self.illuminated_source  # Convert to 0-based
        return sources
    
    def get_all_fibers_sources(self, n_fibers: int, skip_fibers: Optional[List[int]] = None) -> List[Any]:
        """
        Get source configuration with all fibers illuminated.
        
        Parameters
        ----------
        n_fibers : int
            Total number of fibers
        skip_fibers : List[int], optional
            List of 1-based fiber numbers to keep dark
            
        Returns
        -------
        List
            Source list with all fibers illuminated (except skipped ones)
        """
        sources = [self.illuminated_source] * n_fibers
        
        if skip_fibers:
            for fiber_num in skip_fibers:
                if 1 <= fiber_num <= n_fibers:
                    sources[fiber_num - 1] = self.dark_source
        
        return sources
    
    def get_even_odd_sources(self, n_fibers: int, illuminate_even: bool = True) -> List[Any]:
        """
        Get source configuration for even or odd fiber illumination.
        
        Parameters
        ----------
        n_fibers : int
            Total number of fibers
        illuminate_even : bool
            If True, illuminate even-numbered fibers. If False, illuminate odd-numbered.
            
        Returns
        -------
        List
            Source list with even or odd fibers illuminated
        """
        sources = [self.dark_source] * n_fibers
        
        for i in range(n_fibers):
            fiber_num = i + 1  # Convert to 1-based
            is_even = (fiber_num % 2 == 0)
            
            if (illuminate_even and is_even) or (not illuminate_even and not is_even):
                sources[i] = self.illuminated_source
        
        return sources
    
    def get_first_slit_sources(self, n_fibers: int = 66) -> List[Any]:
        """
        Get source configuration for first pseudo-slit illumination (optical bands).
        
        Parameters
        ----------
        n_fibers : int
            Total number of fibers (should be 66 for optical bands)
            
        Returns
        -------
        List
            Source list with first 31 fibers illuminated
        """
        if n_fibers != 66:
            raise ValueError("First slit mode only applicable to optical bands with 66 fibers")
        
        sources = [self.dark_source] * n_fibers
        
        # Illuminate first 31 fibers (indices 0-30)
        for i in range(31):
            sources[i] = self.illuminated_source
        
        return sources
    
    def get_second_slit_sources(self, n_fibers: int = 66) -> List[Any]:
        """
        Get source configuration for second pseudo-slit illumination (optical bands).
        
        Parameters
        ----------
        n_fibers : int
            Total number of fibers (should be 66 for optical bands)
            
        Returns
        -------
        List
            Source list with last 31 fibers illuminated (skipping calibration fibers)
        """
        if n_fibers != 66:
            raise ValueError("Second slit mode only applicable to optical bands with 66 fibers")
        
        sources = [self.dark_source] * n_fibers
        
        # Illuminate fibers 35-65 (indices 34-64), skipping cal fibers 32-34
        for i in range(34, n_fibers):  # Start from index 34 (fiber 35)
            sources[i] = self.illuminated_source
        
        return sources
    
    def get_calibration_sources(self, n_fibers: int = 66) -> List[Any]:
        """
        Get source configuration for calibration fiber illumination pattern.
        
        This creates a specific pattern with calibration fibers illuminated
        and science fibers in specific positions.
        
        Parameters
        ----------
        n_fibers : int
            Total number of fibers (should be 66 for optical bands)
            
        Returns
        -------
        List
            Source list with calibration pattern
        """
        if n_fibers != 66:
            raise ValueError("Calibration mode only applicable to optical bands with 66 fibers")
        
        sources = [self.dark_source] * n_fibers
        
        # Based on original pyechelle_test_ANDES_ff_calib.py pattern:
        # science_b = [Constant(0)] * 31 (fibers 1-31 dark)
        # dark_b = [Constant(0)] (fiber 32 dark)  
        # sim_cal_b = [science] (fiber 33 illuminated)
        # sim_cal_a = [science] (fiber 34 illuminated)
        # dark_a = [Constant(0)] (fiber 35 dark)
        # science_a = [Constant(0)] * 31 (fibers 36-66 dark)
        
        # Illuminate calibration fibers 33 and 34 (indices 32 and 33)
        sources[32] = self.illuminated_source  # Fiber 33
        sources[33] = self.illuminated_source  # Fiber 34
        
        return sources
    
    def get_custom_sources(self, illuminated_fibers: List[int], n_fibers: int) -> List[Any]:
        """
        Get source configuration for custom fiber illumination pattern.
        
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
        
        for fiber_num in illuminated_fibers:
            if 1 <= fiber_num <= n_fibers:
                sources[fiber_num - 1] = self.illuminated_source
            else:
                raise ValueError(f"Fiber number {fiber_num} out of range [1, {n_fibers}]")
        
        return sources