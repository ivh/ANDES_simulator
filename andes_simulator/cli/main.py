"""
Main command-line interface for ANDES simulation framework.

Provides unified CLI with subcommands for all simulation types
including flat field, Fabry-Perot, stellar observations, and post-processing.
"""

import sys
from pathlib import Path
from typing import Optional
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

# Lightweight band list for CLI validation (avoids importing instruments module)
ANDES_BANDS = ['U', 'B', 'V', 'R', 'IZ', 'Y', 'J', 'H']

SUBSLIT_CHOICES = ['all', 'even', 'odd', 'slitA', 'slitB', 'cal',
                   'ifu', 'ring0', 'ring1', 'ring2', 'ring3', 'ring4']


def validate_fiber_spec(ctx, param, value):
    """Validate --subslit/--fiber value: either int or valid mode string."""
    if value is None:
        return 'all'
    # Try to parse as integer (fiber number)
    try:
        return int(value)
    except ValueError:
        pass
    # Validate as mode string
    if value not in SUBSLIT_CHOICES:
        raise click.BadParameter(
            f"Must be fiber number or one of: {', '.join(SUBSLIT_CHOICES)}")
    return value


def common_options(f):
    """Common options for simulation commands."""
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


def resolve_band_and_hdf(
    band: Optional[str],
    hdf: Optional[Path],
    project_root: Optional[Path],
    wl_min: Optional[float] = None,
    wl_max: Optional[float] = None
) -> tuple:
    """Resolve band and HDF model path, inferring band from HDF or wavelengths if needed."""
    from ..core.instruments import infer_band_from_hdf, infer_band_from_wavelengths, get_hdf_model_path

    if hdf:
        inferred_band = infer_band_from_hdf(hdf)
        if band and band != inferred_band:
            raise click.UsageError(
                f"--band {band} conflicts with HDF file (contains {inferred_band} data)")
        return inferred_band, str(hdf)

    if not band:
        # Try to infer from wavelength limits
        if wl_min is not None or wl_max is not None:
            try:
                band = infer_band_from_wavelengths(wl_min, wl_max)
                click.echo(f"Inferred band: {band} (from wavelength limits)")
            except ValueError as e:
                raise click.UsageError(str(e))
        else:
            raise click.UsageError("Either --band, --hdf, or wavelength limits (--wl-min/--wl-max) required")

    # Use default HDF for band
    if project_root is None:
        project_root = Path(__file__).parent.parent.parent
    default_hdf = get_hdf_model_path(band, 'default', project_root)
    return band, str(default_hdf) if default_hdf.exists() else None


def subslit_options(f):
    """Subslit/fiber selection options."""
    f = click.option('--subslit', '--fiber', 'fiber_spec', default='all',
                     callback=validate_fiber_spec, is_eager=True,
                     help='Fiber number (e.g. 21) or mode (all, slitA, ifu, ...)')(f)
    return f


def flux_options(default_scaling=1e5):
    """Flux/scaling options with configurable default scaling."""
    def decorator(f):
        f = click.option('--scaling', type=float, default=default_scaling,
                         help='Base scaling factor')(f)
        f = click.option('--flux', default=1.0, type=float,
                         help='Flux multiplier (multiplied with scaling)')(f)
        return f
    return decorator


@click.group()
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose logging')
@click.option('--project-root', type=click.Path(exists=True, path_type=Path), 
              help='Project root directory (auto-detected if not specified)')
@click.pass_context
def cli(ctx, verbose, project_root):
    """ANDES E2E Simulation Framework - Unified simulation and analysis tools."""
    setup_logging(verbose)
    
    # Store context
    ctx.ensure_object(dict)
    ctx.obj['verbose'] = verbose
    ctx.obj['project_root'] = project_root


@cli.command()
@click.option('--band', type=click.Choice(ANDES_BANDS),
              help='Spectral band (inferred from --hdf or wavelengths if not given)')
@click.option('--source', 'source_spec', required=True, type=str,
              help='Source type: flat, fp, lfc, or path to CSV file')
@subslit_options
@flux_options(default_scaling=1e5)
@common_options
@click.option('--output-name', type=str, help='Output filename (overrides default)')
@click.option('--velocity-shift', type=float, help='Velocity shift in m/s')
@click.option('--finesse', type=float, default=None,
              help='FP finesse (default: band-dependent)')
@click.option('--fp-gap', type=float, default=None,
              help='FP gap thickness in mm (default: auto-computed for ~100 lines/order)')
@click.pass_context
def simulate(ctx, band, source_spec, fiber_spec, flux, scaling, exposure,
             output_dir, output_name, hdf, wl_min, wl_max, fib_eff, velocity_shift,
             finesse, fp_gap, dry_run):
    """Run detector simulation with specified source.

    Source types:
      flat    - Constant flux (flat field calibration)
      fp      - Fabry-Perot calibration spectrum
      lfc     - Laser Frequency Comb calibration
      *.csv   - Custom spectrum from CSV file

    Examples:
      andes-sim simulate --band R --source flat
      andes-sim simulate --band R --source fp --fiber 21 --velocity-shift 100
      andes-sim simulate --band R --source lfc --subslit cal
      andes-sim simulate --band R --source SED/star.csv --fiber 21
    """
    # Resolve source type
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

    # Parse fiber_spec: int means single fiber, string means mode
    if isinstance(fiber_spec, int):
        fiber_mode = 'single'
        fiber = fiber_spec
    else:
        fiber_mode = fiber_spec
        fiber = None

    sim_config = build_config_from_options(
        simulation_type=simulation_type,
        band=band,
        exposure=exposure,
        source_type=source_type,
        fiber_mode=fiber_mode,
        output_dir=output_dir,
        output_name=output_name,
        fiber=fiber,
        flux=flux,
        scaling=scaling,
        velocity_shift=velocity_shift,
        hdf=hdf_path,
        wl_min=wl_min,
        wl_max=wl_max,
        fib_eff=fib_eff,
        spectrum_path=spectrum_path,
        finesse=finesse,
        fp_gap=fp_gap,
    )

    run_simulation_command(
        sim_config,
        dry_run,
        lambda: format_dry_run_output(sim_config),
        f"Simulation completed ({source_spec})"
    )


@cli.command()
@click.option('--band', required=True, type=click.Choice(ANDES_BANDS),
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
    
    # Create HDF builder
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


@cli.command()
@click.option('--band', required=True, type=click.Choice(ANDES_BANDS),
              help='Spectral band')
@click.option('--input-pattern', required=True, type=str,
              help='Input file pattern (e.g., "{band}_FP_fiber{fib:02d}_shift*.fits")')
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

    # Parse kernel size
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
    
    # Create PSF processor
    processor = PSFProcessor(band, ctx.obj['project_root'])
    
    # Set up kernel parameters
    kernel_params = (dimx, dimy, fwhm, edge_blank)
    
    try:
        # Process files
        output_path = processor.process_fabry_perot_files(
            kernel_params=kernel_params,
            output_dir=output_dir,
            input_pattern=input_pattern
        )
        click.echo(f"PSF processing completed: {output_path}")
        
        # Create visualization if requested
        if visualize:
            viz_path = processor.create_kernel_visualization(kernel_params)
            if viz_path:
                click.echo(f"Kernel visualization saved: {viz_path}")
        
    except Exception as e:
        click.echo(f"Error in PSF processing: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option('--band', required=True, type=click.Choice(ANDES_BANDS),
              help='Spectral band')
@click.option('--input-pattern', required=True, type=str,
              help='Input file pattern (e.g., "{band}_FP_fiber{fib:02d}_shift*.fits")')
@click.option('--mode', default='all',
              type=click.Choice(['all', 'even_odd', 'slits', 'custom']),
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
    
    # Create fiber combiner (use output_dir as input_dir for finding files)
    combiner = FiberCombiner(band, ctx.obj['project_root'], input_dir=output_dir)
    
    try:
        if mode == 'all':
            combined_image = combiner.combine_all_fibers(input_pattern)
            output_filename = output or f"{band}_combined_all.fits"
            
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
            if band in ['Y', 'J', 'H']:
                click.echo("Slit mode not applicable to YJH bands", err=True)
                return
            
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
            
            # Parse fiber list (simple implementation)
            fiber_nums = []
            for part in fibers.split(','):
                if '-' in part:
                    start, end = map(int, part.split('-'))
                    fiber_nums.extend(range(start, end + 1))
                else:
                    fiber_nums.append(int(part))
            
            combined_image = combiner.combine_fiber_subset(fiber_nums, input_pattern)
            output_filename = output or f"{band}_combined_custom.fits"
        
        # Save result (for all, custom modes)
        if mode in ['all', 'custom']:
            if output_dir:
                output_path = output_dir / output_filename
            else:
                output_path = combiner.project_root.parent / band / output_filename
            
            combiner.save_combined_image(combined_image, output_path, {'mode': mode})
            click.echo(f"Combined image saved: {output_path}")
        
        # Generate report if requested
        if report:
            report_path = combiner.create_combination_report(input_pattern, output_dir)
            click.echo(f"Report created: {report_path}")
            
    except Exception as e:
        click.echo(f"Error in fiber combination: {e}", err=True)
        sys.exit(1)


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


@cli.command()
@click.pass_context
def list_bands(ctx):
    """List available spectral bands."""
    bands = ANDES_BANDS
    click.echo("Available spectral bands:")
    for band in bands:
        click.echo(f"  {band}")


def main():
    """Main entry point for the CLI."""
    cli()


if __name__ == '__main__':
    main()