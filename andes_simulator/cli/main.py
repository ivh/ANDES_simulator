"""
CLI factory for spectrograph simulation frameworks.

Provides create_cli() which builds a Click CLI parameterized by
instrument-specific band lists and subslit choices.
"""

import sys
from pathlib import Path
from typing import Optional, List

try:
    import click
except ImportError:
    print("Error: click package required for CLI. Install with: pip install click")
    sys.exit(1)

from .utils import (
    run_simulation_command,
    build_config_from_options,
    format_dry_run_output,
    setup_logging
)


def create_cli(instrument_name: str, bands: List[str], subslit_choices: List[str]):
    """Build a Click CLI for a specific instrument.

    Parameters
    ----------
    instrument_name : str
        Instrument name for help text (e.g. 'ANDES', 'MOSAIC')
    bands : list of str
        Valid band names for this instrument
    subslit_choices : list of str
        Valid subslit/fiber mode names
    """

    def validate_fiber_spec(ctx, param, value):
        if value is None:
            return 'all'
        try:
            return int(value)
        except ValueError:
            pass
        if value.startswith('bundle:'):
            return value
        if value not in subslit_choices:
            raise click.BadParameter(
                f"Must be fiber number, bundle:N, bundle:N-M, "
                f"or one of: {', '.join(subslit_choices)}")
        return value

    def common_options(f):
        f = click.option('--dry-run', is_flag=True, help='Show what would be done without running')(f)
        f = click.option('--output-dir', type=click.Path(path_type=Path), help='Output directory')(f)
        f = click.option('--exposure', default=1.0, type=float, help='Exposure time in seconds')(f)
        f = click.option('--hdf', type=click.Path(exists=True, path_type=Path),
                         help='HDF model file (infers band if --band not given)')(f)
        f = click.option('--wl-min', type=float, help='Minimum wavelength in nm')(f)
        f = click.option('--wl-max', type=float, help='Maximum wavelength in nm')(f)
        f = click.option('--fib-eff', type=str, default='0.9-0.95', show_default=True,
                         help='Fiber efficiency: single value (0.9) or range (0.7-0.9)')(f)
        return f

    def resolve_band_and_hdf(band, hdf, project_root, wl_min=None, wl_max=None):
        from ..core.instruments import infer_band_from_hdf, infer_band_from_wavelengths, get_hdf_model_path

        if hdf:
            inferred_band = infer_band_from_hdf(hdf, restrict_to=bands)
            if band and band != inferred_band:
                raise click.UsageError(
                    f"--band {band} conflicts with HDF file (contains {inferred_band} data)")
            return inferred_band, str(hdf)

        if not band:
            if wl_min is not None or wl_max is not None:
                try:
                    band = infer_band_from_wavelengths(wl_min, wl_max, restrict_to=bands)
                    click.echo(f"Inferred band: {band} (from wavelength limits)")
                except ValueError as e:
                    raise click.UsageError(str(e))
            else:
                raise click.UsageError("Either --band, --hdf, or wavelength limits (--wl-min/--wl-max) required")

        if project_root is None:
            project_root = Path(__file__).parent.parent.parent
        default_hdf = get_hdf_model_path(band, 'default', project_root)
        return band, str(default_hdf) if default_hdf.exists() else None

    def subslit_options(f):
        f = click.option('--subslit', '--fiber', 'fiber_spec', default='all',
                         callback=validate_fiber_spec, is_eager=True,
                         help='Fiber number (e.g. 21) or mode (all, slitA, ifu, ...)')(f)
        return f

    def flux_options(f):
        f = click.option('--scaling', type=float, default=None,
                         help='Base scaling factor (default: per-band)')(f)
        f = click.option('--flux', default=1.0, type=float,
                         help='Flux multiplier (multiplied with scaling)')(f)
        return f

    # --- CLI group ---

    @click.group()
    @click.option('--verbose', '-v', is_flag=True, help='Enable verbose logging')
    @click.option('--project-root', type=click.Path(exists=True, path_type=Path),
                  help='Project root directory (auto-detected if not specified)')
    @click.pass_context
    def cli(ctx, verbose, project_root):
        """E2E Simulation Framework."""
        setup_logging(verbose)
        ctx.ensure_object(dict)
        ctx.obj['verbose'] = verbose
        ctx.obj['project_root'] = project_root

    cli.help = f"{instrument_name} E2E Simulation Framework"

    # --- simulate ---

    @cli.command()
    @click.option('--band', type=click.Choice(bands),
                  help='Spectral band (inferred from --hdf or wavelengths if not given)')
    @click.option('--source', 'source_spec', required=True, type=str,
                  help='Source type: flat, fp, lfc, or path to CSV file')
    @subslit_options
    @flux_options
    @common_options
    @click.option('--output-name', type=str, help='Output filename (overrides default)')
    @click.option('--velocity-shift', type=str, default=None,
                  help='Velocity shift: value in m/s or path to JSON offsets file')
    @click.option('--x-shift', type=str, default=None,
                  help='Pixel x-shift: value in pixels or path to JSON offsets file')
    @click.option('--finesse', type=float, default=None,
                  help='FP finesse (default: band-dependent)')
    @click.option('--fp-gap', type=float, default=None,
                  help='FP gap thickness in mm (default: auto-computed for ~100 lines/order)')
    @click.pass_context
    def simulate(ctx, band, source_spec, fiber_spec, flux, scaling, exposure,
                 output_dir, output_name, hdf, wl_min, wl_max, fib_eff, velocity_shift,
                 x_shift, finesse, fp_gap, dry_run):
        """Run detector simulation with specified source.

        Source types:
          flat    - Constant flux (flat field calibration)
          fp      - Fabry-Perot calibration spectrum
          lfc     - Laser Frequency Comb calibration
          *.csv   - Custom spectrum from CSV file
        """
        source_spec_lower = source_spec.lower()
        if source_spec_lower == 'flat':
            source_type = 'constant'
            simulation_type = 'flat_field'
            spectrum_path = None
        elif source_spec_lower == 'fp':
            source_type = 'fabry_perot'
            simulation_type = 'fabry_perot'
            spectrum_path = None
        elif source_spec_lower == 'lfc':
            source_type = 'lfc'
            simulation_type = 'lfc'
            spectrum_path = None
        elif source_spec.endswith('.csv') or Path(source_spec).exists():
            source_type = 'csv'
            simulation_type = 'spectrum'
            spectrum_path = Path(source_spec)
            if not spectrum_path.exists():
                raise click.BadParameter(f"Spectrum file not found: {source_spec}")
        else:
            raise click.BadParameter(
                f"Unknown source '{source_spec}'. Use: flat, fp, lfc, or path to CSV")

        band, hdf_path = resolve_band_and_hdf(band, hdf, ctx.obj['project_root'], wl_min, wl_max)

        if isinstance(fiber_spec, int):
            fiber_mode = 'single'
            fiber = fiber_spec
            custom_fibers = None
        elif isinstance(fiber_spec, str) and fiber_spec.startswith('bundle:'):
            from ..core.instruments import get_instrument_config
            inst_cfg = get_instrument_config(band)
            bundle_size = inst_cfg.get('bundle_size')
            n_bundles = inst_cfg.get('n_bundles')
            if not bundle_size:
                raise click.UsageError(
                    f"bundle: syntax not supported for {band}-band (no bundle_size in config)")
            bundle_str = fiber_spec.split(':', 1)[1]
            if '-' in bundle_str:
                lo, hi = bundle_str.split('-', 1)
                bundles = range(int(lo), int(hi) + 1)
            else:
                bundles = [int(bundle_str)]
            custom_fibers = []
            for b in bundles:
                if b < 1 or (n_bundles and b > n_bundles):
                    raise click.UsageError(f"Bundle {b} out of range (1-{n_bundles})")
                start = (b - 1) * bundle_size + 1
                custom_fibers.extend(range(start, start + bundle_size))
            fiber_mode = 'custom'
            fiber = None
            click.echo(f"Bundle {bundle_str} -> fibers {custom_fibers}")
        else:
            fiber_mode = fiber_spec
            fiber = None
            custom_fibers = None

        velocity_shift_file = None
        velocity_shift_value = None
        if velocity_shift is not None:
            try:
                velocity_shift_value = float(velocity_shift)
            except ValueError:
                import json
                vshift_path = Path(velocity_shift)
                if not vshift_path.exists():
                    raise click.BadParameter(
                        f"Not a number and file not found: {velocity_shift}",
                        param_hint="--velocity-shift")
                with open(vshift_path) as f:
                    offsets = json.load(f)
                velocity_shift_file = str(vshift_path)
                if fiber is not None:
                    key = str(fiber)
                    if key not in offsets:
                        raise click.UsageError(f"Fiber {fiber} not found in {vshift_path}")
                    velocity_shift_value = float(offsets[key])
                else:
                    # Multi-fiber mode: build per-fiber RV list, defer to after config
                    velocity_shift_value = offsets  # pass raw dict, resolved below

        x_shift_file = None
        x_shift_value = None
        if x_shift is not None:
            try:
                x_shift_value = float(x_shift)
            except ValueError:
                import json
                xshift_path = Path(x_shift)
                if not xshift_path.exists():
                    raise click.BadParameter(
                        f"Not a number and file not found: {x_shift}",
                        param_hint="--x-shift")
                if fiber is None:
                    raise click.UsageError(
                        "JSON x-shift file requires --fiber N (single fiber mode)")
                with open(xshift_path) as f:
                    offsets = json.load(f)
                key = str(fiber)
                if key not in offsets:
                    raise click.UsageError(f"Fiber {fiber} not found in {xshift_path}")
                x_shift_value = float(offsets[key])
                x_shift_file = str(xshift_path)

        # For multi-fiber JSON velocity shifts, pass None initially and resolve after
        vshift_for_config = None if isinstance(velocity_shift_value, dict) else velocity_shift_value

        sim_config = build_config_from_options(
            simulation_type=simulation_type,
            band=band,
            exposure=exposure,
            source_type=source_type,
            fiber_mode=fiber_mode,
            output_dir=output_dir,
            output_name=output_name,
            fiber=custom_fibers if custom_fibers else fiber,
            flux=flux,
            scaling=scaling,
            velocity_shift=vshift_for_config,
            x_shift=x_shift_value,
            hdf=hdf_path,
            wl_min=wl_min,
            wl_max=wl_max,
            fib_eff=fib_eff,
            spectrum_path=spectrum_path,
            finesse=finesse,
            fp_gap=fp_gap,
        )

        # Resolve per-fiber velocity shifts from JSON dict
        if isinstance(velocity_shift_value, dict):
            n_fibers = sim_config.n_fibers
            rv_list = [0.0] * n_fibers
            for k, v in velocity_shift_value.items():
                idx = int(k) - 1  # 1-based fiber number to 0-based index
                if 0 <= idx < n_fibers:
                    rv_list[idx] = float(v)
            sim_config.velocity_shift = rv_list

        run_simulation_command(
            sim_config,
            dry_run,
            lambda: format_dry_run_output(sim_config, velocity_shift_file=velocity_shift_file,
                                          velocity_shift_fiber=fiber,
                                          x_shift_file=x_shift_file,
                                          x_shift_fiber=fiber),
            f"Simulation completed ({source_spec})"
        )

    # --- generate-offsets ---

    @cli.command()
    @click.option('--band', required=True, type=click.Choice(bands),
                  help='Spectral band')
    @subslit_options
    @click.option('--rms', required=True, type=float,
                  help='RMS of velocity offsets in m/s')
    @click.option('--seed', default=None, type=int,
                  help='Random seed for reproducibility')
    @click.option('-o', '--output', required=True, type=click.Path(path_type=Path),
                  help='Output JSON file path')
    def generate_offsets(band, fiber_spec, rms, seed, output):
        """Generate random per-fiber velocity offsets and save to JSON."""
        import json
        import numpy as np
        from ..core.config import SimulationConfig, SourceConfig, FiberConfig, OutputConfig

        if isinstance(fiber_spec, int):
            fiber_mode = 'single'
            fibers = [fiber_spec]
        else:
            fiber_mode = fiber_spec
            fibers = "all"

        config = SimulationConfig(
            simulation_type='flat_field',
            band=band,
            source=SourceConfig(),
            fibers=FiberConfig(mode=fiber_mode, fibers=fibers),
            output=OutputConfig(),
        )
        fiber_list = config.get_fiber_list()

        rng = np.random.default_rng(seed)
        offsets = {fib: round(float(rng.normal(0, rms)), 2) for fib in fiber_list}

        output.parent.mkdir(parents=True, exist_ok=True)
        with open(output, 'w') as f:
            json.dump(offsets, f, indent=2)

        click.echo(f"Generated {len(offsets)} offsets (rms={rms} m/s, seed={seed}) -> {output}")
        for fib, v in sorted(offsets.items()):
            click.echo(f"  fiber {fib:2d}: {v:+8.2f} m/s")

    # --- generate-hdf ---

    @cli.command()
    @click.option('--band', required=True, type=click.Choice(bands),
                  help='Spectral band')
    @click.option('--zemax-file', type=str, help='ZEMAX file name (uses default if not specified)')
    @click.option('--output', type=click.Path(path_type=Path), help='Output HDF file path')
    @click.option('--n-transform', default=15, type=int, help='Transformations per order')
    @click.option('--n-psf', default=15, type=int, help='PSFs per order')
    @click.option('--dry-run', is_flag=True, help='Show what would be done without running')
    @click.pass_context
    def generate_hdf(ctx, band, zemax_file, output, n_transform, n_psf, dry_run):
        """Generate HDF instrument model files from ZEMAX."""
        from ..models.hdf_builder import AndesHDFBuilder

        if dry_run:
            click.echo("Dry run - would execute:")
            click.echo(f"  Band: {band}")
            click.echo(f"  ZEMAX file: {zemax_file or 'default'}")
            click.echo(f"  Output: {output or 'default path'}")
            click.echo(f"  Samples: {n_transform} transforms, {n_psf} PSFs per order")
            click.echo("  WARNING: This would take a very long time!")
            return

        click.echo(f"Generating HDF model for {band}-band...")
        click.echo("WARNING: This will take a very long time!")

        if not click.confirm("Do you want to continue?"):
            return

        builder = AndesHDFBuilder(band, ctx.obj['project_root'])

        try:
            output_path = builder.build_complete_model(
                output_path=output,
                zemax_file=zemax_file,
                n_transformation_per_order=n_transform,
                n_psfs_per_order=n_psf
            )
            click.echo(f"HDF model generated: {output_path}")
        except Exception as e:
            click.echo(f"Error generating HDF model: {e}", err=True)
            sys.exit(1)

    # --- psf-process ---

    @cli.command()
    @click.option('--band', required=True, type=click.Choice(bands),
                  help='Spectral band')
    @click.option('--input-pattern', required=True, type=str,
                  help='Input file pattern')
    @click.option('--kernel-size', default='4,4', type=str, help='Kernel size as "width,height"')
    @click.option('--fwhm', default=3.2, type=float, help='FWHM in arcseconds')
    @click.option('--edge-blank', default='left',
                  type=click.Choice(['none', 'top', 'bottom', 'left', 'right', 'random']),
                  help='Edge blanking mode')
    @click.option('--output-dir', type=click.Path(path_type=Path), help='Output directory')
    @click.option('--visualize', is_flag=True, help='Create kernel visualization')
    @click.option('--dry-run', is_flag=True, help='Show what would be done without running')
    @click.pass_context
    def psf_process(ctx, band, input_pattern, kernel_size, fwhm, edge_blank, output_dir, visualize, dry_run):
        """Apply PSF convolution to simulation outputs."""
        from ..postprocess.psf import PSFProcessor

        try:
            dimx, dimy = map(int, kernel_size.split(','))
        except ValueError:
            raise click.BadParameter("kernel-size must be 'width,height' format")

        if dry_run:
            click.echo("Dry run - would execute:")
            click.echo(f"  Band: {band}")
            click.echo(f"  Input pattern: {input_pattern}")
            click.echo(f"  Kernel: {dimx}x{dimy}, FWHM={fwhm}\", edge={edge_blank}")
            click.echo(f"  Output dir: {output_dir or 'default'}")
            return

        processor = PSFProcessor(band, ctx.obj['project_root'])
        kernel_params = (dimx, dimy, fwhm, edge_blank)

        try:
            output_path = processor.process_fabry_perot_files(
                kernel_params=kernel_params,
                output_dir=output_dir,
                input_pattern=input_pattern
            )
            click.echo(f"PSF processing completed: {output_path}")

            if visualize:
                viz_path = processor.create_kernel_visualization(kernel_params)
                if viz_path:
                    click.echo(f"Kernel visualization saved: {viz_path}")

        except Exception as e:
            click.echo(f"Error in PSF processing: {e}", err=True)
            sys.exit(1)

    # --- combine ---

    @cli.command()
    @click.option('--band', required=True, type=click.Choice(bands),
                  help='Spectral band')
    @click.option('--input-pattern', required=True, type=str,
                  help='Input file pattern')
    @click.option('--mode', default='all',
                  type=click.Choice(['all', 'all_but_dark', 'even_odd', 'slits', 'custom']),
                  help='Combination mode')
    @click.option('--fibers', type=str, help='Fiber list for custom mode (e.g., "1,5,10-15")')
    @click.option('--output', type=str, help='Output filename')
    @click.option('--output-dir', type=click.Path(path_type=Path), help='Output directory')
    @click.option('--report', is_flag=True, help='Generate combination report')
    @click.option('--dry-run', is_flag=True, help='Show what would be done without running')
    @click.pass_context
    def combine(ctx, band, input_pattern, mode, fibers, output, output_dir, report, dry_run):
        """Combine individual fiber outputs."""
        from ..postprocess.combine import FiberCombiner

        if dry_run:
            click.echo("Dry run - would execute:")
            click.echo(f"  Band: {band}")
            click.echo(f"  Input pattern: {input_pattern}")
            click.echo(f"  Mode: {mode}")
            if fibers:
                click.echo(f"  Fibers: {fibers}")
            return

        combiner = FiberCombiner(band, ctx.obj['project_root'], input_dir=output_dir)

        try:
            if mode in ('all', 'all_but_dark'):
                skip_dark = (mode == 'all_but_dark')
                combined_image = combiner.combine_all_fibers(input_pattern, skip_dark=skip_dark)
                output_filename = output or f"{band}_combined_{mode}.fits"

            elif mode == 'even_odd':
                results = combiner.combine_even_odd_fibers(input_pattern)
                for eo_mode, image_data in results.items():
                    eo_filename = output or f"{band}_combined_{eo_mode}.fits"
                    if output_dir:
                        eo_path = output_dir / eo_filename
                    else:
                        eo_path = combiner.project_root.parent / band / eo_filename
                    combiner.save_combined_image(image_data, eo_path, {'mode': eo_mode})
                    click.echo(f"Saved {eo_mode} combination: {eo_path}")

                if report:
                    report_path = combiner.create_combination_report(input_pattern, output_dir)
                    click.echo(f"Report created: {report_path}")
                return

            elif mode == 'slits':
                results = combiner.combine_pseudo_slits(input_pattern)
                for slit_name, image_data in results.items():
                    slit_filename = output or f"{band}_combined_{slit_name}.fits"
                    if output_dir:
                        slit_path = output_dir / slit_filename
                    else:
                        slit_path = combiner.project_root.parent / band / slit_filename
                    combiner.save_combined_image(image_data, slit_path, {'mode': slit_name})
                    click.echo(f"Saved {slit_name} combination: {slit_path}")

                if report:
                    report_path = combiner.create_combination_report(input_pattern, output_dir)
                    click.echo(f"Report created: {report_path}")
                return

            elif mode == 'custom':
                if not fibers:
                    raise click.BadParameter("--fibers required for custom mode")

                fiber_nums = []
                for part in fibers.split(','):
                    if '-' in part:
                        start, end = map(int, part.split('-'))
                        fiber_nums.extend(range(start, end + 1))
                    else:
                        fiber_nums.append(int(part))

                combined_image = combiner.combine_fiber_subset(fiber_nums, input_pattern)
                output_filename = output or f"{band}_combined_custom.fits"

            if mode in ['all', 'all_but_dark', 'custom']:
                if output_dir:
                    output_path = output_dir / output_filename
                else:
                    output_path = combiner.project_root.parent / band / output_filename

                combiner.save_combined_image(combined_image, output_path, {'mode': mode})
                click.echo(f"Combined image saved: {output_path}")

            if report:
                report_path = combiner.create_combination_report(input_pattern, output_dir)
                click.echo(f"Report created: {report_path}")

        except Exception as e:
            click.echo(f"Error in fiber combination: {e}", err=True)
            sys.exit(1)

    # --- run-config ---

    @cli.command()
    @click.option('--config', required=True, type=click.Path(exists=True, path_type=Path),
                  help='YAML configuration file')
    @click.option('--dry-run', is_flag=True, help='Show what would be done without running')
    @click.pass_context
    def run_config(ctx, config, dry_run):
        """Run simulation from YAML configuration file."""
        from ..core.config import SimulationConfig

        sim_config = SimulationConfig.from_yaml(config)

        run_simulation_command(
            sim_config,
            dry_run,
            lambda: format_dry_run_output(sim_config, [f"Config: {config}"]),
            f"Simulation completed: {sim_config.simulation_type}"
        )

    # --- create-templates ---

    @cli.command()
    @click.option('--output-dir', default='configs/examples', type=click.Path(path_type=Path),
                  help='Output directory for template files')
    @click.pass_context
    def create_templates(ctx, output_dir):
        """Create template configuration files."""
        from ..core.config import create_template_configs

        if not output_dir.is_absolute():
            if ctx.obj['project_root']:
                output_dir = ctx.obj['project_root'] / output_dir
            else:
                output_dir = Path.cwd() / output_dir

        try:
            create_template_configs(output_dir)
            click.echo(f"Template configurations created in: {output_dir}")
        except Exception as e:
            click.echo(f"Error creating templates: {e}", err=True)
            sys.exit(1)

    # --- list-bands ---

    @cli.command()
    @click.pass_context
    def list_bands(ctx):
        """List available spectral bands."""
        click.echo(f"Available {instrument_name} bands:")
        for band in bands:
            click.echo(f"  {band}")

    return cli
