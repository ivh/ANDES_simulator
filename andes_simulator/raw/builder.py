"""Raw frame builder: DPR string -> EDPS-ready MEF raw frames.

One template execution = one builder call: parse the DPR grammar, fetch or
simulate the slot expectation images through the cache, normalize to target
peak levels, sum, and per exposure apply the detector layer and write one
MEF per spectrograph arm (DRL v1.2 Ch. 4.1 layout).
"""

import logging
import re
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from ..core.andes import SPECTROGRAPHS
from .cache import SimCache, ensure_many
from .detector import DetectorModel
from .dpr import parse_dpr, peak_target_e, source_spec_for_token, DprSpec
from .headers import build_raw_hdul

logger = logging.getLogger(__name__)

READOUT_OVERHEAD_S = 60.0
NORM_PERCENTILE = 99.9


class RawFrameBuilder:
    def __init__(self, project_root: Path, cache: SimCache, output_dir: Path,
                 seed: Optional[int] = None, jobs: int = 1,
                 headers_only: bool = False):
        self.project_root = Path(project_root)
        self.cache = cache
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.rng_master = np.random.default_rng(seed)
        self.jobs = jobs
        # headers-only: stub 2x2 extensions, no pyechelle, no detector model;
        # enough for classification/organization testing of the workflow
        self.headers_only = headers_only
        self._counters: Dict[str, int] = {}
        self._detectors: Dict[str, DetectorModel] = {}

    def _detector(self, band: str, readout: Optional[str]) -> DetectorModel:
        key = f"{band}:{readout}"
        if key not in self._detectors:
            self._detectors[key] = DetectorModel(band, readout=readout)
        return self._detectors[key]

    def _next_filename(self, arm: str) -> Path:
        if arm not in self._counters:
            pattern = re.compile(rf"ANDES_{arm}_(\d+)\.fits$")
            existing = [int(m.group(1)) for f in self.output_dir.glob(f"ANDES_{arm}_*.fits")
                        if (m := pattern.match(f.name))]
            self._counters[arm] = max(existing, default=0)
        self._counters[arm] += 1
        return self.output_dir / f"ANDES_{arm}_{self._counters[arm]:04d}.fits"

    def slot_entries(self, arm: str, dpr_type: str, mode: Optional[str],
                     bands: Optional[Sequence[str]] = None,
                     ins_mask: Optional[str] = None,
                     calfib: Optional[str] = None) -> List[Dict]:
        """Cache entries needed for a frame (for parallel pre-filling)."""
        arm_cfg = SPECTROGRAPHS[arm]
        sim_bands = [b for b in arm_cfg['bands'] if bands is None or b in bands]
        entries = []
        for band in sim_bands:
            spec = parse_dpr(dpr_type, band=band, mode=mode, ins_mask=ins_mask,
                             calfib=calfib)
            if spec.detector_only:
                continue
            for slot in spec.slots:
                entries.append({
                    'band': band,
                    'source_kwargs': source_spec_for_token(
                        slot.token, band, self.project_root),
                    'fibers': slot.fibers,
                })
        return entries

    def _expectation_for_band(self, band: str, spec: DprSpec,
                              detector: DetectorModel,
                              exptime: float) -> Optional[np.ndarray]:
        if spec.detector_only:
            return detector.led_expectation(exptime) if spec.led else None

        total = None
        for slot in spec.slots:
            source_kwargs = source_spec_for_token(slot.token, band,
                                                  self.project_root)
            rate = self.cache.get_or_simulate(band, source_kwargs, slot.fibers)
            lit = rate[rate > 0]
            if lit.size == 0:
                logger.warning("slot %s (%s) produced no photons on %s",
                               slot.name, slot.token, band)
                continue
            peak = float(np.percentile(lit, NORM_PERCENTILE))
            scale = peak_target_e(spec.kind, slot.token) / peak
            contrib = rate * np.float32(scale)
            total = contrib if total is None else total + contrib
        return total

    def build(self, arm: str, dpr_type: str,
              mode: Optional[str] = None,
              exptime: float = 10.0,
              nexp: int = 1,
              bands: Optional[Sequence[str]] = None,
              tpl_start: Optional[datetime] = None,
              obs_start: Optional[datetime] = None,
              tpl_id: Optional[str] = None,
              tpl_nexp: Optional[int] = None,
              tpl_expno_start: int = 1,
              catg: Optional[str] = None,
              tech: Optional[str] = None,
              ins_mask: Optional[str] = None,
              calfib: Optional[str] = None,
              ifu_scale: Optional[int] = None,
              binx: int = 1, biny: int = 1,
              readout: Optional[str] = None,
              extra_keywords: Optional[Dict] = None) -> List[Path]:
        """Generate nexp raw frames for one template execution on one arm."""
        arm_cfg = SPECTROGRAPHS[arm]
        arm_bands = arm_cfg['bands']
        sim_bands = [b for b in arm_bands if bands is None or b in bands]
        if not sim_bands:
            raise ValueError(f"No bands of arm {arm} selected ({bands})")
        if arm_cfg['detector_type'] != 'CCD':
            binx = biny = 1

        # reference spec (grammar/keywords are band-independent)
        ref_spec = parse_dpr(dpr_type, band=sim_bands[0], mode=mode,
                             ins_mask=ins_mask, calfib=calfib,
                             catg=catg, tech=tech)

        expectations: Dict[str, Optional[np.ndarray]] = {}
        detectors: Dict[str, DetectorModel] = {}
        if not self.headers_only:
            entries = self.slot_entries(arm, dpr_type, mode, sim_bands,
                                        ins_mask, calfib)
            if entries:
                ensure_many(self.cache, entries, jobs=self.jobs)

            for band in arm_bands:
                detectors[band] = self._detector(band, readout)
                if band in sim_bands:
                    spec = parse_dpr(dpr_type, band=band, mode=mode,
                                     ins_mask=ins_mask, calfib=calfib,
                                     catg=catg, tech=tech)
                    expectations[band] = self._expectation_for_band(
                        band, spec, detectors[band], exptime)
                else:
                    expectations[band] = None

        tpl_start = tpl_start or datetime.now(timezone.utc)
        tpl_start_str = tpl_start.strftime('%Y-%m-%dT%H:%M:%S')

        keywords = {
            'dpr.catg': ref_spec.catg,
            'dpr.type': ref_spec.dpr_type,
            'dpr.tech': ref_spec.tech,
            'ins.mode': ref_spec.mode,
            'det.binx': binx,
            'det.biny': biny,
            'tpl.start': tpl_start_str,
            'tpl.id': tpl_id,
            'tpl.nexp': tpl_nexp if tpl_nexp is not None else nexp,
        }
        if arm_cfg['detector_type'] == 'CCD':
            keywords['det.readout'] = (readout or
                                       self._detector(arm_bands[0], readout).readout)
        if not ref_spec.detector_only:
            keywords['ins.calfib'] = ref_spec.calfib or 'OFF'
        if ins_mask:
            keywords['ins.mask'] = ins_mask
        if ifu_scale is not None and ref_spec.mode == 'IFU-AO':
            keywords['ins.ifu.scale'] = ifu_scale
        if extra_keywords:
            keywords.update(extra_keywords)

        paths = []
        obs_time = obs_start or tpl_start
        for i in range(nexp):
            images, ext_meta = {}, {}
            for band in arm_bands:
                if self.headers_only:
                    images[band] = np.zeros((2, 2), dtype=np.uint16)
                    ext_meta[band] = {'simulated': False}
                    continue
                det = detectors[band]
                rng = np.random.default_rng(self.rng_master.integers(2**63))
                images[band] = det.apply(expectations[band], exptime, rng,
                                         binx=binx, biny=biny)
                ext_meta[band] = {
                    'gain': det.gain_e_adu,
                    'ron': det.ron_e,
                    'simulated': band in sim_bands,
                }

            keywords['tpl.expno'] = tpl_expno_start + i
            hdul = build_raw_hdul(arm, images, ext_meta, obs_time, exptime,
                                  keywords)
            path = self._next_filename(arm)
            hdul.writeto(path, overwrite=True)
            paths.append(path)
            logger.info("wrote %s (%s, %s)", path.name, ref_spec.dpr_type,
                        ref_spec.catg)
            obs_time = obs_time + timedelta(
                seconds=exptime + READOUT_OVERHEAD_S)

        return paths
