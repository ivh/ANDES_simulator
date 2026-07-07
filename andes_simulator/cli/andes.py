"""ANDES CLI entry point.

Adds the ANDES-specific raw-frame commands (make-raw, night) on top of the
shared instrument CLI. See src/dpr_summary.md for the raw-frame design.
"""

from datetime import datetime, timezone
from pathlib import Path

import click

from .main import create_cli
from ..core.andes import ANDES_BANDS, SUBSLIT_CHOICES, SPECTROGRAPHS

cli = create_cli('ANDES', ANDES_BANDS, SUBSLIT_CHOICES)


def _project_root(ctx) -> Path:
    return ctx.obj.get('project_root') or Path(__file__).parent.parent.parent


def _parse_time(value):
    if value is None:
        return None
    dt = datetime.fromisoformat(value)
    return dt if dt.tzinfo else dt.replace(tzinfo=timezone.utc)


def _parse_csv(value):
    return [v.strip() for v in value.split(',') if v.strip()] if value else None


def raw_options(f):
    f = click.option('--seed', type=int, default=None,
                     help='Seed for per-exposure noise (reproducible frames)')(f)
    f = click.option('--boost', type=float, default=10.0, show_default=True,
                     help='Expectation-cache flux boost (see dpr_summary.md)')(f)
    f = click.option('--cache-dir', type=click.Path(path_type=Path), default=None,
                     help='Simulation cache directory (default: E2E/simcache)')(f)
    f = click.option('--jobs', type=int, default=1, show_default=True,
                     help='Parallel subprocesses for slot simulations')(f)
    f = click.option('--bands', type=str, default=None,
                     help='Comma list of detector bands to simulate; other '
                          'extensions get detector noise only')(f)
    f = click.option('-o', '--output-dir', type=click.Path(path_type=Path),
                     default=Path('.'), show_default=True)(f)
    f = click.option('--headers-only', is_flag=True,
                     help='Stub 2x2 extensions, no pyechelle/detector model; '
                          'for workflow classification tests')(f)
    f = click.option('--dry-run', is_flag=True)(f)
    return f


def _make_builder(ctx, output_dir, cache_dir, boost, seed, jobs,
                  headers_only=False):
    from ..raw.builder import RawFrameBuilder
    from ..raw.cache import SimCache

    root = _project_root(ctx)
    cache_dir = cache_dir or root.parent / 'simcache'
    cache = SimCache(cache_dir, root, boost=boost)
    return RawFrameBuilder(root, cache, output_dir, seed=seed, jobs=jobs,
                           headers_only=headers_only)


# --- make-raw ---

@cli.command('make-raw')
@click.option('--arm', type=click.Choice(list(SPECTROGRAPHS)), required=True,
              help='Spectrograph arm (one MEF per arm, DRL 4.1)')
@click.option('--dpr', 'dpr_type', required=True,
              help='DPR.TYPE string, e.g. "WAVE,HCL,FP" or BIAS')
@click.option('--mode', type=click.Choice(['SL-UNI', 'IFU-AO']), default=None,
              help='INS.MODE (default: SL-UNI for 2-slot types)')
@click.option('--calfib', type=click.Choice(['FP', 'HCL', 'LFC', 'LAMP', 'OFF']),
              default=None,
              help='Calibration fibre source (ins.calfib keyword; default OFF)')
@click.option('--exptime', type=float, default=10.0, show_default=True)
@click.option('--nexp', type=int, default=1, show_default=True)
@click.option('--tpl-start', type=str, default=None, help='ISO time (default now)')
@click.option('--tpl-id', type=str, default=None, help='Template name')
@click.option('--catg', type=str, default=None, help='DPR.CATG override')
@click.option('--tech', type=str, default=None, help='DPR.TECH override')
@click.option('--ins-mask', type=click.Choice(['M1', 'M2', 'M3']), default=None)
@click.option('--ifu-scale', type=int, default=None)
@click.option('--binx', type=int, default=1, show_default=True)
@click.option('--biny', type=int, default=1, show_default=True)
@click.option('--readout', type=str, default=None,
              help='Readout mode (CCD: fast/slow)')
@raw_options
@click.pass_context
def make_raw(ctx, arm, dpr_type, mode, calfib, exptime, nexp, tpl_start, tpl_id,
             catg, tech, ins_mask, ifu_scale, binx, biny, readout,
             seed, boost, cache_dir, jobs, bands, output_dir, headers_only,
             dry_run):
    """Generate EDPS-ready raw frame(s) from a DPR.TYPE string."""
    from ..raw.dpr import parse_dpr

    bands_list = _parse_csv(bands)
    arm_bands = SPECTROGRAPHS[arm]['bands']
    if bands_list:
        bad = [b for b in bands_list if b not in arm_bands]
        if bad:
            raise click.UsageError(f"Bands {bad} not in arm {arm} ({arm_bands})")

    if dry_run:
        click.echo(f"Dry run - would generate {nexp} frame(s):")
        click.echo(f"  Arm: {arm} (bands: {bands_list or arm_bands})")
        for band in (bands_list or arm_bands):
            spec = parse_dpr(dpr_type, band=band, mode=mode,
                             ins_mask=ins_mask, calfib=calfib,
                             catg=catg, tech=tech)
            if spec.detector_only:
                what = 'LED flat' if spec.led else 'detector-only'
                click.echo(f"  {band}: {what} ({spec.catg}, {spec.tech})")
            else:
                for slot in spec.slots:
                    click.echo(f"  {band} slot {slot.name}: {slot.token} on "
                               f"{slot.subslit} ({len(slot.fibers)} fibers)")
        return

    builder = _make_builder(ctx, output_dir, cache_dir, boost, seed, jobs,
                            headers_only)
    paths = builder.build(
        arm=arm, dpr_type=dpr_type, mode=mode, exptime=exptime, nexp=nexp,
        bands=bands_list, tpl_start=_parse_time(tpl_start), tpl_id=tpl_id,
        catg=catg, tech=tech, ins_mask=ins_mask, calfib=calfib,
        ifu_scale=ifu_scale, binx=binx, biny=biny, readout=readout)
    for p in paths:
        click.echo(f"wrote {p}")


# --- night ---

@cli.command('night')
@click.argument('plan', type=click.Path(exists=True, path_type=Path))
@click.option('--arms', type=str, default=None,
              help='Comma list of arms (UBV,RIZ,YJH; default: all in plan)')
@click.option('--sets', type=str, default=None,
              help='Comma list of calibration sets (detector,daily,when-used)')
@click.option('--include', type=str, default='',
              help='Extra frame groups: night,science')
@click.option('--vis-config', default='1x1_fast', show_default=True,
              help='VIS binning/readout config name from the plan')
@click.option('--ifu-scale', type=int, default=16, show_default=True)
@click.option('--date', type=str, default=None,
              help='UTC start time ISO (default: today 10:00)')
@raw_options
@click.pass_context
def night(ctx, plan, arms, sets, include, vis_config, ifu_scale, date,
          seed, boost, cache_dir, jobs, bands, output_dir, headers_only,
          dry_run):
    """Generate a synthetic (calibration) night from a calibration plan YAML."""
    from ..raw.night import load_plan, plan_night, describe, run_night

    plan_data = load_plan(plan)
    planned = plan_night(
        plan_data, arms=_parse_csv(arms), sets=_parse_csv(sets),
        include=_parse_csv(include) or (), vis_config=vis_config,
        ifu_scale=ifu_scale)

    click.echo(f"Planned exposure sequence from {plan}:")
    click.echo(describe(planned))
    if dry_run:
        return

    builder = _make_builder(ctx, output_dir, cache_dir, boost, seed, jobs,
                            headers_only)
    written = run_night(plan, output_dir, builder, planned,
                        start=_parse_time(date), bands=_parse_csv(bands))
    click.echo(f"wrote {len(written)} frames to {output_dir}")


def main():
    cli()
