"""
ANDES E2E Simulation Framework

A consolidated framework for running ANDES spectrograph simulations
including flat field calibrations, Fabry-Perot wavelength calibrations,
stellar observations, and post-processing.
"""

__version__ = "1.0.0"
__author__ = "ANDES E2E Team"

__all__ = [
    'AndesSimulator',
    'SimulationConfig',
    'get_instrument_config'
]

def __getattr__(name):
    """Lazy import heavy modules only when accessed."""
    if name == 'AndesSimulator':
        from .core.simulator import AndesSimulator
        return AndesSimulator
    elif name == 'SimulationConfig':
        from .core.config import SimulationConfig
        return SimulationConfig
    elif name == 'get_instrument_config':
        from .core.instruments import get_instrument_config
        return get_instrument_config
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")