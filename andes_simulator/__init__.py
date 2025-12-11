"""
ANDES E2E Simulation Framework

A consolidated framework for running ANDES spectrograph simulations
including flat field calibrations, Fabry-Perot wavelength calibrations,
stellar observations, and post-processing.
"""

__version__ = "1.0.0"
__author__ = "ANDES E2E Team"

from .core.simulator import AndesSimulator
from .core.config import SimulationConfig
from .core.instruments import get_instrument_config

__all__ = [
    'AndesSimulator',
    'SimulationConfig', 
    'get_instrument_config'
]