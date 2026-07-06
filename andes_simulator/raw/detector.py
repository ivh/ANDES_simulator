"""Detector layer: photon-expectation images to realistic raw ADU frames.

Pipeline: expectation e- x PRNU -> + dark current (incl. hot pixels) ->
dead pixels -> Poisson -> charge binning -> read noise -> gain -> bias ->
uint16 with saturation clipping.

Cosmetics (PRNU, hot/dead pixels, bad columns) are static per band, seeded
from the band name, so bad-pixel maps and flat fields are consistent across
frames. Maps are generated per row-block from spawned RNG streams instead of
being held in memory (the VIS detectors are 9216x9232).
"""

import zlib
from typing import Optional, Tuple

import numpy as np

from ..core.andes import SPECTROGRAPHS, DETECTOR_MODELS
from ..core.instruments import get_instrument_config

BLOCK_ROWS = 1024
SATURATION_ADU = 65535
LED_LEVEL_E_S = 800.0


def spectrograph_for_band(band: str) -> Tuple[str, dict]:
    """The spectrograph (arm) a detector band belongs to."""
    for arm, cfg in SPECTROGRAPHS.items():
        if band in cfg['bands']:
            return arm, cfg
    raise ValueError(f"Band '{band}' not part of any ANDES spectrograph arm")


def _band_seed(band: str) -> int:
    return zlib.crc32(f"ANDES-detector-{band}".encode())


class DetectorModel:
    """Per-band detector model with static cosmetics."""

    def __init__(self, band: str, readout: Optional[str] = None):
        self.band = band
        self.arm, arm_cfg = spectrograph_for_band(band)
        self.detector_type = arm_cfg['detector_type']
        self.has_bias = arm_cfg['has_bias']
        self.params = DETECTOR_MODELS[self.detector_type]

        self.readout = readout or self.params['default_readout']
        if self.readout not in self.params['readout_modes']:
            raise ValueError(
                f"Unknown readout mode '{self.readout}' for {self.detector_type} "
                f"(available: {list(self.params['readout_modes'])})")
        mode = self.params['readout_modes'][self.readout]
        self.ron_e = mode['ron_e']
        self.gain_e_adu = mode['gain_e_adu']

        size_xy = get_instrument_config(band)['detector_size']
        self.shape = (size_xy[1], size_xy[0])  # numpy (Y, X)

        self._seed = _band_seed(band)
        rng = np.random.default_rng(self._seed)
        n_bad = self.params['n_bad_columns']
        self.bad_columns = sorted(
            rng.choice(self.shape[1], size=n_bad, replace=False).tolist()) if n_bad else []

    def _block_maps(self, y0: int, n_rows: int):
        """Static PRNU / dark / dead maps for rows [y0, y0+n_rows)."""
        block_idx = y0 // BLOCK_ROWS
        rng = np.random.default_rng(
            np.random.SeedSequence(entropy=self._seed, spawn_key=(block_idx,)))
        shape = (n_rows, self.shape[1])
        p = self.params

        prnu = rng.normal(1.0, p['prnu_rms'], shape).astype(np.float32)

        dark = np.full(shape, p['dark_e_s'], dtype=np.float32)
        hot = rng.random(shape) < p['hot_pixel_frac']
        if hot.any():
            lo, hi = p['hot_dark_e_s']
            dark[hot] = np.exp(rng.uniform(np.log(lo), np.log(hi), int(hot.sum())))

        dead = rng.random(shape) < p['dead_pixel_frac']
        for col in self.bad_columns:
            dead[:, col] = True

        return prnu, dark, dead

    def apply(self, expectation_e: Optional[np.ndarray], exptime: float,
              rng: np.random.Generator, binx: int = 1, biny: int = 1) -> np.ndarray:
        """Turn a photon-expectation image (e-, unbinned) into a raw frame.

        expectation_e None means no light (bias/dark frames). rng drives the
        per-exposure noise; cosmetics stay fixed.
        """
        ny, nx = self.shape
        if expectation_e is not None and expectation_e.shape != self.shape:
            raise ValueError(
                f"Expectation shape {expectation_e.shape} != detector {self.shape}")
        if ny % biny or nx % binx:
            raise ValueError(f"Detector {self.shape} not divisible by binning "
                             f"{binx}x{biny}")

        out = np.empty((ny // biny, nx // binx), dtype=np.uint16)
        ron_adu = self.ron_e / self.gain_e_adu

        for y0 in range(0, ny, BLOCK_ROWS):
            n_rows = min(BLOCK_ROWS, ny - y0)
            prnu, dark, dead = self._block_maps(y0, n_rows)

            lam = dark * exptime
            if expectation_e is not None:
                lam = lam + expectation_e[y0:y0 + n_rows].astype(np.float32) * prnu
            lam[dead] = 0.0
            np.clip(lam, 0.0, None, out=lam)

            electrons = rng.poisson(lam.astype(np.float64))
            if binx > 1 or biny > 1:
                electrons = electrons.reshape(
                    n_rows // biny, biny, nx // binx, binx).sum(axis=(1, 3))

            adu = electrons / self.gain_e_adu + self.params['bias_adu']
            adu += rng.normal(0.0, ron_adu, electrons.shape)
            np.rint(adu, out=adu)
            np.clip(adu, 0, SATURATION_ADU, out=adu)
            out[y0 // biny:(y0 + n_rows) // biny] = adu.astype(np.uint16)

        return out

    def led_expectation(self, exptime: float,
                        level_e_s: float = LED_LEVEL_E_S) -> np.ndarray:
        """Uniform LED illumination (PRNU applied later in apply())."""
        return np.full(self.shape, level_e_s * exptime, dtype=np.float32)
