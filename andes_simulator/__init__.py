"""
ANDES E2E Simulation Framework

A consolidated framework for running ANDES spectrograph simulations
including flat field calibrations, Fabry-Perot wavelength calibrations,
stellar observations, and post-processing.
"""

# Give this process its own numba cache dir to avoid "underlying object has
# vanished" errors from cache corruption (seen both with parallel workers
# sharing a cache and with sequential multi-fiber runs reusing it).
# Respect an externally-set NUMBA_CACHE_DIR so wrapper scripts that manage
# their own per-subprocess caches keep working.
# Must run before any pyechelle/numba import — this module is reached first
# via the CLI entry points.
import os as _os
if "NUMBA_CACHE_DIR" not in _os.environ:
    import tempfile as _tempfile
    import atexit as _atexit
    import shutil as _shutil
    _numba_cache = _tempfile.mkdtemp(prefix="numba_andes_")
    _os.environ["NUMBA_CACHE_DIR"] = _numba_cache
    _atexit.register(_shutil.rmtree, _numba_cache, ignore_errors=True)

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