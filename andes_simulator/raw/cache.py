"""Simulation cache for raw-frame building.

Pyechelle output is a Monte Carlo realization, so caching raw simulations
would freeze one shot-noise pattern into every frame built from them.
Instead the cache stores *boosted expectation* images: each unique
(band, model, source, fibers) slot is simulated once at boost x nominal
flux and divided by boost; exposures then draw fresh Poisson noise from
the scaled expectation (see dpr_summary.md). Residual correlated noise is
1/boost of the shot variance.

Simulations seed numpy's global RNG from the band name before running so
per-fiber efficiencies (drawn in AndesSimulator._apply_fiber_efficiency)
are a static instrument property, identical across all cached slots of a
band. Flat-fielding logic would break otherwise.
"""

import hashlib
import json
import logging
import tempfile
import zlib
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from astropy.io import fits

logger = logging.getLogger(__name__)

DEFAULT_BOOST = 10.0
REF_EXPOSURE_S = 1.0

_SIM_TYPE_FOR_SOURCE = {
    'constant': 'flat_field',
    'fabry_perot': 'fabry_perot',
    'lfc': 'lfc',
    'hcl': 'hcl',
    'csv': 'spectrum',
}


def _band_rng_seed(band: str) -> int:
    return zlib.crc32(f"ANDES-fibeff-{band}".encode())


class SimCache:
    def __init__(self, cache_dir: Path, project_root: Path,
                 boost: float = DEFAULT_BOOST):
        self.cache_dir = Path(cache_dir)
        self.project_root = Path(project_root)
        self.boost = boost
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def signature(self, band: str, source_kwargs: Dict, fibers: Sequence[int],
                  hdf: Optional[str] = None) -> Dict:
        sig = {
            'band': band,
            'hdf': hdf or 'default',
            'source': dict(sorted(source_kwargs.items())),
            'fibers': list(fibers),
            'boost': self.boost,
            'ref_exposure': REF_EXPOSURE_S,
        }
        filepath = source_kwargs.get('filepath')
        if filepath:
            p = self.project_root / filepath
            if p.exists():
                st = p.stat()
                sig['source_file_stat'] = [int(st.st_size), int(st.st_mtime)]
        return sig

    def path_for(self, sig: Dict) -> Path:
        digest = hashlib.sha1(
            json.dumps(sig, sort_keys=True).encode()).hexdigest()[:16]
        return self.cache_dir / f"{sig['band']}_{digest}.fits"

    def load(self, sig: Dict) -> Optional[np.ndarray]:
        path = self.path_for(sig)
        if not path.exists():
            return None
        with fits.open(path) as hdul:
            return hdul[0].data.astype(np.float32)

    def resolve_fibers(self, band: str, fibers: Sequence[int]) -> Tuple[int, ...]:
        """Drop the band's dead fibers and canonicalize the list."""
        from ..core.instruments import get_instrument_config
        skip = set(get_instrument_config(band).get('skip_fibers', []))
        return tuple(sorted(set(fibers) - skip))

    def get_or_simulate(self, band: str, source_kwargs: Dict,
                        fibers: Sequence[int],
                        hdf: Optional[str] = None) -> np.ndarray:
        """Expectation rate image (e-/s per unbinned pixel, float32)."""
        fibers = self.resolve_fibers(band, fibers)
        if not fibers:
            raise ValueError(f"No live fibers left for {band}-band slot")
        sig = self.signature(band, source_kwargs, fibers, hdf)
        cached = self.load(sig)
        if cached is not None:
            logger.info("cache hit: %s", self.path_for(sig).name)
            return cached
        logger.info("cache miss: simulating %s %s on %d fibers (boost %g)",
                    band, source_kwargs.get('type'), len(fibers), self.boost)
        rate = _simulate_slot(self.project_root, band, source_kwargs,
                              fibers, self.boost, hdf)
        self.store(sig, rate)
        return rate

    def store(self, sig: Dict, rate: np.ndarray) -> Path:
        path = self.path_for(sig)
        hdu = fits.PrimaryHDU(rate.astype(np.float32))
        hdu.header['KEYJSON'] = json.dumps(sig, sort_keys=True)
        hdu.header['BUNIT'] = 'electron/s'
        tmp = path.with_suffix('.tmp.fits')
        hdu.writeto(tmp, overwrite=True)
        tmp.rename(path)
        return path


def _simulate_slot(project_root: Path, band: str, source_kwargs: Dict,
                   fibers: Tuple[int, ...], boost: float,
                   hdf: Optional[str]) -> np.ndarray:
    """Run one boosted pyechelle simulation, return the rate image."""
    from ..core.config import SimulationConfig, SourceConfig, FiberConfig, OutputConfig
    from ..core.simulator import AndesSimulator
    from ..cli.utils import resolve_source_scaling

    source_type = source_kwargs['type']
    scaling, use_file_scaling = resolve_source_scaling(
        source_type, band, flux=boost, user_scaling=None)
    kwargs = dict(source_kwargs)
    kwargs['scaling_factor'] = scaling
    if source_type == 'constant':
        kwargs['flux'] = scaling
    if use_file_scaling is not None:
        kwargs['use_file_scaling'] = use_file_scaling

    with tempfile.TemporaryDirectory(prefix='andes_raw_') as tmpdir:
        config = SimulationConfig(
            simulation_type=_SIM_TYPE_FOR_SOURCE[source_type],
            band=band,
            exposure_time=REF_EXPOSURE_S,
            hdf_model=hdf,
            source=SourceConfig(**kwargs),
            fibers=FiberConfig(mode='custom', fibers=list(fibers)),
            output=OutputConfig(directory=tmpdir, filename='slot.fits'),
        )
        # deterministic per-band fiber efficiencies (static instrument property)
        np.random.seed(_band_rng_seed(band))
        simulator = AndesSimulator(config)
        try:
            simulator.run_simulation()
        finally:
            simulator.cleanup()
        with fits.open(Path(tmpdir) / 'slot.fits') as hdul:
            data = hdul[0].data.astype(np.float32)

    return data / (boost * REF_EXPOSURE_S)


def _ensure_entry(args) -> str:
    """Worker for parallel cache filling; returns the cache file name."""
    cache_dir, project_root, boost, band, source_kwargs, fibers, hdf = args
    cache = SimCache(Path(cache_dir), Path(project_root), boost)
    cache.get_or_simulate(band, source_kwargs, fibers, hdf)
    return cache.path_for(
        cache.signature(band, source_kwargs,
                        cache.resolve_fibers(band, fibers), hdf)).name


def ensure_many(cache: SimCache, entries: List[Dict], jobs: int = 1) -> None:
    """Fill cache for many slots, optionally in parallel subprocesses.

    Each worker process gets its own NUMBA_CACHE_DIR via the package
    __init__; pyechelle itself must stay at max_cpu=1.
    """
    todo = {}
    for e in entries:
        fibers = cache.resolve_fibers(e['band'], e['fibers'])
        if not fibers:
            continue
        sig = cache.signature(e['band'], e['source_kwargs'], fibers,
                              e.get('hdf'))
        path = cache.path_for(sig)
        if not path.exists():
            todo[path] = (str(cache.cache_dir), str(cache.project_root),
                          cache.boost, e['band'], e['source_kwargs'],
                          fibers, e.get('hdf'))
    if not todo:
        return
    if jobs <= 1 or len(todo) == 1:
        for args in todo.values():
            _ensure_entry(args)
        return
    import concurrent.futures as cf
    import multiprocessing as mp
    ctx = mp.get_context('spawn')
    with cf.ProcessPoolExecutor(max_workers=jobs, mp_context=ctx) as pool:
        for name in pool.map(_ensure_entry, todo.values()):
            logger.info("cached %s", name)
