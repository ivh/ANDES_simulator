"""Night driver: calibration_plan.yaml -> a synthetic night of raw frames.

Reads the canonical plan (versioned in the edps repo), plans the exposure
sequence and drives RawFrameBuilder. Template executions sharing the same
template name on the same arm share one TPL.START group (the workflow
groups by tpl.start and needs e.g. all wave exposures in one group;
reconciliation item 9 in the plan YAML).
"""

import logging
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import yaml

from ..core.andes import SPECTROGRAPHS
from .builder import RawFrameBuilder, READOUT_OVERHEAD_S

logger = logging.getLogger(__name__)

LED_EXPTIME_SERIES = [1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 60.0]
DEFAULT_EXPTIME_S = 120.0
DEFAULT_SCIENCE_EXPTIME_S = 300.0
DEFAULT_N = 2


@dataclass
class PlannedExposures:
    """One homogeneous run of exposures within a template execution."""
    arm: str
    dpr_type: str
    n: int
    exptime: float
    mode: Optional[str] = None
    catg: Optional[str] = None
    tech: Optional[str] = None
    tpl_id: Optional[str] = None
    group: Optional[str] = None      # frames with equal group share tpl.start
    ins_mask: Optional[str] = None
    ifu_scale: Optional[int] = None
    binx: int = 1
    biny: int = 1
    readout: Optional[str] = None
    extra_keywords: Dict = field(default_factory=dict)


def load_plan(path: Path) -> Dict:
    with open(path) as f:
        return yaml.safe_load(f)


def _exptimes_for(exp: Dict, n: int) -> List[float]:
    raw = exp.get('exptime_s', DEFAULT_EXPTIME_S)
    if isinstance(raw, (int, float)):
        return [float(raw)] * n
    if raw == 'varied':
        return [LED_EXPTIME_SERIES[i % len(LED_EXPTIME_SERIES)] for i in range(n)]
    return [DEFAULT_EXPTIME_S] * n


def _n_for(exp: Dict) -> int:
    n = exp.get('n', DEFAULT_N)
    return n if isinstance(n, int) else DEFAULT_N


def _vis_config(plan: Dict, name: str) -> Dict:
    for cfg in plan.get('setups', {}).get('vis_configs', []):
        if cfg.get('name') == name:
            return cfg
    raise ValueError(f"VIS config '{name}' not in plan setups.vis_configs")


def plan_night(plan: Dict,
               arms: Optional[Sequence[str]] = None,
               sets: Optional[Sequence[str]] = None,
               include: Sequence[str] = (),
               vis_config: str = '1x1_fast',
               ifu_scale: int = 16) -> List[PlannedExposures]:
    """Expand the plan into a concrete exposure sequence."""
    known_arms = [a for a in plan.get('setups', {}).get('arms', {})
                  if a in SPECTROGRAPHS]
    requested = [a for a in (arms or known_arms) if a in known_arms]
    if not requested:
        raise ValueError(f"No usable arms among {arms} (known: {known_arms})")

    vis = _vis_config(plan, vis_config)
    planned: List[PlannedExposures] = []

    def vis_kwargs(arm: str) -> Dict:
        if SPECTROGRAPHS[arm]['detector_type'] != 'CCD':
            return {}
        return {'binx': int(vis.get('det.binx', 1)),
                'biny': int(vis.get('det.biny', 1)),
                'readout': vis.get('det.readout')}

    for proc in plan.get('procedures', []):
        if 'reference' in proc:
            continue
        if sets and proc.get('set') not in sets:
            continue
        applies = proc.get('applies_to', {})
        proc_arms = applies.get('arms', known_arms)
        modes = applies.get('modes', [])
        mode = modes[0] if modes else None
        per_scale = applies.get('per') == 'ifu_scale'
        dpr = proc.get('dpr', {})
        tpl_id = proc.get('template')

        for arm in requested:
            if arm not in proc_arms:
                continue
            for exp in proc.get('exposures', []):
                n = _n_for(exp)
                exptimes = _exptimes_for(exp, n)
                kw = exp.get('keywords', {})
                for exptime, count in _runs(exptimes):
                    planned.append(PlannedExposures(
                        arm=arm, dpr_type=exp['type'], n=count,
                        exptime=exptime, mode=mode,
                        catg=dpr.get('catg'),
                        tech=exp.get('tech') or dpr.get('tech'),
                        tpl_id=tpl_id, group=f"{arm}:{tpl_id}:{mode}",
                        ins_mask=kw.get('ins.mask'),
                        ifu_scale=ifu_scale if per_scale else None,
                        **vis_kwargs(arm)))

    if 'night' in include:
        for entry in plan.get('night_calibrations', []):
            if entry.get('status') == 'upgrade' or not entry.get('exposures'):
                continue
            dpr = entry.get('dpr', {})
            tpl_id = entry.get('template')
            for arm in requested:
                for exp in entry.get('exposures', []):
                    tech = exp.get('tech') or dpr.get('tech')
                    mode = 'IFU-AO' if tech and 'IFU' in tech else 'SL-UNI'
                    if mode == 'IFU-AO' and arm != 'YJH':
                        continue
                    n = _n_for(exp)
                    for exptime, count in _runs(_exptimes_for(exp, n)):
                        planned.append(PlannedExposures(
                            arm=arm, dpr_type=exp['type'], n=count,
                            exptime=exptime, mode=mode,
                            catg=dpr.get('catg'), tech=tech,
                            tpl_id=tpl_id,
                            group=f"{arm}:{tpl_id or entry.get('name')}:{mode}",
                            **vis_kwargs(arm)))

    if 'science' in include:
        for entry in plan.get('observations', []):
            dpr = entry.get('dpr', {})
            tech = dpr.get('tech', '')
            mode = 'IFU-AO' if 'IFU' in tech else 'SL-UNI'
            tpl_id = (entry.get('templates') or [entry.get('name')])[0]
            for arm in requested:
                if mode == 'IFU-AO' and arm != 'YJH':
                    continue
                for exp in entry.get('exposures', []):
                    n = exp.get('n', 1)
                    n = n if isinstance(n, int) else 1
                    planned.append(PlannedExposures(
                        arm=arm, dpr_type=exp['type'], n=n,
                        exptime=DEFAULT_SCIENCE_EXPTIME_S, mode=mode,
                        catg=dpr.get('catg'), tech=exp.get('tech') or dpr.get('tech'),
                        tpl_id=tpl_id,
                        group=f"{arm}:{entry.get('name')}:{mode}",
                        **vis_kwargs(arm)))

    return planned


def _runs(exptimes: List[float]):
    """Collapse an exposure-time list into (exptime, count) runs."""
    runs = []
    for t in exptimes:
        if runs and runs[-1][0] == t:
            runs[-1][1] += 1
        else:
            runs.append([t, 1])
    return [(t, c) for t, c in runs]


def describe(planned: List[PlannedExposures]) -> str:
    lines = []
    total = 0
    for p in planned:
        total += p.n
        lines.append(
            f"  {p.arm:4s} {p.n:3d} x {p.exptime:7.1f}s  {p.dpr_type:24s} "
            f"{p.mode or '-':7s} {p.catg or 'auto':9s} tpl={p.tpl_id or '-'}")
    lines.append(f"  total: {total} frames in {len(planned)} exposure runs")
    return "\n".join(lines)


def run_night(plan_path: Path, output_dir: Path, builder: RawFrameBuilder,
              planned: List[PlannedExposures],
              start: Optional[datetime] = None,
              bands: Optional[Sequence[str]] = None) -> List[Path]:
    clock = start or datetime.now(timezone.utc).replace(
        hour=10, minute=0, second=0, microsecond=0)

    group_starts: Dict[str, datetime] = {}
    group_expno: Dict[str, int] = {}
    group_nexp: Dict[str, int] = {}
    for p in planned:
        group_nexp[p.group] = group_nexp.get(p.group, 0) + p.n

    written: List[Path] = []
    for p in planned:
        tpl_start = group_starts.setdefault(p.group, clock)
        expno = group_expno.get(p.group, 1)
        paths = builder.build(
            arm=p.arm, dpr_type=p.dpr_type, mode=p.mode,
            exptime=p.exptime, nexp=p.n, bands=bands,
            tpl_start=tpl_start, obs_start=clock,
            tpl_id=p.tpl_id, tpl_nexp=group_nexp[p.group],
            tpl_expno_start=expno,
            catg=p.catg, tech=p.tech, ins_mask=p.ins_mask,
            ifu_scale=p.ifu_scale, binx=p.binx, biny=p.biny,
            readout=p.readout, extra_keywords=p.extra_keywords)
        written.extend(paths)
        group_expno[p.group] = expno + p.n
        clock = clock + timedelta(
            seconds=p.n * (p.exptime + READOUT_OVERHEAD_S))
    return written
