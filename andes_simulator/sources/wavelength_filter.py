"""Wavelength-filtered wrapper for sources lacking spectral files."""

from __future__ import annotations

from typing import Optional
import numpy as np


class WavelengthFilteredSource:
    """Wrap a PyEchelle source to enforce wavelength limits."""

    def __init__(self, source, wl_min: Optional[float] = None, wl_max: Optional[float] = None):
        self._source = source
        self.wl_min = wl_min
        self.wl_max = wl_max

        # Propagate common attributes used elsewhere
        self.list_like = getattr(source, 'list_like', False)
        self.name = getattr(source, 'name', 'filtered_source')

    def get_counts(self, wl, integration_time):
        counts = self._source.get_counts(wl, integration_time)
        mask = self._build_mask(wl)
        return self._apply_mask(counts, mask)

    def get_counts_rv(self, wl, integration_time, rv=0):
        if hasattr(self._source, 'get_counts_rv'):
            counts = self._source.get_counts_rv(wl, integration_time, rv)
        else:
            counts = self._source.get_counts(wl, integration_time)
        mask = self._build_mask(wl)
        return self._apply_mask(counts, mask)

    # ------------------------------------------------------------------
    def _build_mask(self, wl):
        if isinstance(wl, (list, tuple)):
            return [self._build_mask(item) for item in wl]
        if isinstance(wl, np.ndarray) and wl.dtype == object:
            return [self._build_mask(item) for item in wl.tolist()]

        wl_nm = self._to_nanometers(wl)
        mask = np.ones_like(wl_nm, dtype=bool)
        if self.wl_min is not None:
            mask &= (wl_nm >= self.wl_min)
        if self.wl_max is not None:
            mask &= (wl_nm <= self.wl_max)
        return mask

    def _to_nanometers(self, wl):
        if hasattr(wl, 'to'):
            return wl.to('nm').value
        wl_arr = np.asarray(wl)
        if wl_arr.size == 0:
            return wl_arr
        sample = wl_arr.flat[0]
        if sample < 10:  # assume microns
            return wl_arr * 1000.0
        return wl_arr

    def _apply_mask(self, counts, mask):
        if isinstance(mask, list):
            if isinstance(counts, np.ndarray) and counts.dtype == object:
                masked = [self._apply_mask(c, m) for c, m in zip(counts.tolist(), mask)]
                return np.array(masked, dtype=object)
            if isinstance(counts, list):
                return [self._apply_mask(c, m) for c, m in zip(counts, mask)]
            if isinstance(counts, tuple):
                return tuple(self._apply_mask(c, m) for c, m in zip(counts, mask))

        mask_arr = np.asarray(mask)
        if mask_arr.ndim == 0:
            mask_arr = np.ones_like(np.asarray(counts), dtype=bool)

        if hasattr(counts, 'value'):
            data = np.asarray(counts.value)
            masked = np.where(mask_arr, data, 0)
            counts.value[:] = masked
            return counts

        counts_arr = np.asarray(counts)
        masked = np.where(mask_arr, counts_arr, 0)
        return masked
