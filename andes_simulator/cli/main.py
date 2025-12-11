"""
Main command-line interface for ANDES simulation framework.

Provides unified CLI with subcommands for all simulation types
including flat field, Fabry-Perot, stellar observations, and post-processing.
"""

import sys
import logging
from pathlib import Path
from typing import Optional, List

try:
    import click
except ImportError:
    print("Error: click package required for CLI. Install with: pip install click")
    sys.exit(1)

# Heavy imports moved to lazy loading within commands to speed up --help

# Lightweight band list for CLI validation (avoids importing instruments module)
ANDES_BANDS = ['U', 'B', 'V', 'R', 'IZ', 'Y', 'J', 'H']


@click.group()
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose logging')
@click.option('--project-root', type=click.Path(exists=True, path_type=Path), 
              help='Project root directory (auto-detected if not specified)')
@click.pass_context
def cli(ctx, verbose, project_root):
    """ANDES E2E Simulation Framework - Unified simulation and analysis tools."""
    # Set up logging
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # Store context
    ctx.ensure_object(dict)
    ctx.obj['verbose'] = verbose
    ctx.obj['project_root'] = project_root


@cli.command()
@click.option('--band', required=True, type=click.Choice(ANDES_BANDS), 
              help='Spectral band')
@click.option('--mode', default='all', 
              type=click.Choice(['all', 'single', 'even_odd', 'first_slit', 'second_slit', 'calib']),
              help='Fiber illumination mode')
@click.option('--fiber', type=int, help='Specific fiber number (for single mode)')
@click.option('--flux', default=0.001, type=float, help='Flux level for flat field')
@click.option('--exposure', default=30.0, type=float, help='Exposure time in seconds')
@click.option('--output-dir', type=click.Path(path_type=Path), help='Output directory')
@click.option('--config', type=click.Path(exists=True, path_type=Path), 
              help='YAML configuration file')
@click.option('--dry-run', is_flag=True, help='Show what would be done without running')
@click.pass_context
def flat_field(ctx, band, mode, fiber, flux, exposure, output_dir, config, dry_run):
    """Generate flat field calibration frames."""
    from ..core.simulator import AndesSimulator
    from ..core.config import SimulationConfig, SourceConfig, FiberConfig, OutputConfig

    if config:
        # Load from configuration file
        sim_config = SimulationConfig.from_yaml(config)
    else:
        # Create configuration from command line options
        if mode == 'single' and fiber is None:
            raise click.BadParameter("--fiber required for single mode")
        
        fibers = [fiber] if mode == 'single' else "all"
        
        sim_config = SimulationConfig(
            simulation_type="flat_field",
            band=band,
            exposure_time=exposure,
            source=SourceConfig(type="constant", flux=flux, flux_unit="ph/s/AA"),
            fibers=FiberConfig(mode=mode, fibers=fibers),
            output=OutputConfig(directory=str(output_dir) if output_dir else "../{band}/")
        )
    
    if dry_run:
        click.echo("Dry run - would execute:")
        click.echo(f"  Band: {sim_config.band}")
        click.echo(f"  Mode: {sim_config.fibers.mode}")
        click.echo(f"  Exposure: {sim_config.exposure_time}s")
        click.echo(f"  Flux: {sim_config.source.flux}")
        return
    
    # Create and run simulator
    simulator = AndesSimulator(sim_config)
    
    if mode == 'single' and isinstance(sim_config.fibers.fibers, list) and len(sim_config.fibers.fibers) == 1:
        # Single fiber mode
        result = simulator.run_simulation()
        click.echo(f"Flat field simulation completed for fiber {sim_config.fibers.fibers[0]}")
    elif mode == 'even_odd':
        # Even/odd mode creates two files
        results = simulator.run_simulation()
        click.echo("Even/odd flat field simulations completed")
    else:
        # Standard modes
        result = simulator.run_simulation()
        click.echo(f"Flat field simulation completed: {mode} mode")


@cli.command()
@click.option('--band', required=True, type=click.Choice(ANDES_BANDS),
              help='Spectral band')
@click.option('--mode', default='all',
              type=click.Choice(['all', 'single']),
              help='Fiber illumination mode')
@click.option('--fiber', type=int, help='Specific fiber number (for single mode)')
@click.option('--velocity-shift', type=float, help='Velocity shift in m/s')
@click.option('--scaling', default=5e9, type=float, help='FP flux scaling factor')
@click.option('--exposure', default=30.0, type=float, help='Exposure time in seconds')
@click.option('--output-dir', type=click.Path(path_type=Path), help='Output directory')
@click.option('--config', type=click.Path(exists=True, path_type=Path),
              help='YAML configuration file')
@click.option('--dry-run', is_flag=True, help='Show what would be done without running')
@click.pass_context
def fabry_perot(ctx, band, mode, fiber, velocity_shift, scaling, exposure, output_dir, config, dry_run):
    """Generate Fabry-Perot wavelength calibration frames."""
    from ..core.simulator import AndesSimulator
    from ..core.config import SimulationConfig, SourceConfig, FiberConfig, OutputConfig

    if config:
        sim_config = SimulationConfig.from_yaml(config)
    else:
        if mode == 'single' and fiber is None:
            raise click.BadParameter("--fiber required for single mode")
        
        fibers = [fiber] if mode == 'single' else "all"
        
        sim_config = SimulationConfig(
            simulation_type="fabry_perot",
            band=band,
            exposure_time=exposure,
            velocity_shift=velocity_shift,
            source=SourceConfig(type="fabry_perot", scaling_factor=scaling),
            fibers=FiberConfig(mode=mode, fibers=fibers),
            output=OutputConfig(directory=str(output_dir) if output_dir else "../{band}/")
        )
    
    if dry_run:
        click.echo("Dry run - would execute:")
        click.echo(f"  Band: {sim_config.band}")
        click.echo(f"  Mode: {sim_config.fibers.mode}")
        click.echo(f"  Velocity shift: {sim_config.velocity_shift} m/s" if sim_config.velocity_shift else "  No velocity shift")
        click.echo(f"  Scaling: {sim_config.source.scaling_factor}")
        return
    
    # Create and run simulator
    simulator = AndesSimulator(sim_config)
    result = simulator.run_simulation()
    click.echo("Fabry-Perot simulation completed")


@cli.command()
@click.option('--band', required=True, type=click.Choice(ANDES_BANDS),
              help='Spectral band')
@click.option('--spectrum', required=True, type=click.Path(exists=True, path_type=Path),
              help='CSV spectrum file')
@click.option('--fiber', required=True, type=int, help='Fiber number to illuminate')
@click.option('--scaling', default=5e3, type=float, help='Spectrum flux scaling factor')
@click.option('--exposure', default=30.0, type=float, help='Exposure time in seconds')
@click.option('--output-dir', type=click.Path(path_type=Path), help='Output directory')
@click.option('--config', type=click.Path(exists=True, path_type=Path),
              help='YAML configuration file')
@click.option('--dry-run', is_flag=True, help='Show what would be done without running')
@click.pass_context
def spectrum(ctx, band, spectrum, fiber, scaling, exposure, output_dir, config, dry_run):
    """Generate stellar spectrum observations."""
    from ..core.simulator import AndesSimulator
    from ..core.config import SimulationConfig, SourceConfig, FiberConfig, OutputConfig

    if config:
        sim_config = SimulationConfig.from_yaml(config)
    else:
        from ..core.config import SourceConfig, FiberConfig, OutputConfig
        
        sim_config = SimulationConfig(
            simulation_type="spectrum",
            band=band,
            exposure_time=exposure,
            source=SourceConfig(
                type="csv",
                filepath=str(spectrum),
                scaling_factor=scaling
            ),
            fibers=FiberConfig(mode="single", fibers=[fiber]),
            output=OutputConfig(directory=str(output_dir) if output_dir else "../{band}/")
        )
    
    if dry_run:
        click.echo("Dry run - would execute:")
        click.echo(f"  Band: {sim_config.band}")
        click.echo(f"  Spectrum: {sim_config.source.filepath}")
        click.echo(f"  Fiber: {fiber}")
        click.echo(f"  Scaling: {sim_config.source.scaling_factor}")
        return
    
    # Create and run simulator
    simulator = AndesSimulator(sim_config)
    result = simulator.run_simulation()
    click.echo(f"Spectrum simulation completed for fiber {fiber}")


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
    
    # Create fiber combiner
    combiner = FiberCombiner(band, ctx.obj['project_root'])
    
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
    from ..core.simulator import AndesSimulator
    from ..core.config import SimulationConfig

    # Load configuration
    sim_config = SimulationConfig.from_yaml(config)
    
    if dry_run:
        click.echo("Dry run - would execute:")
        click.echo(f"  Config: {config}")
        click.echo(f"  Type: {sim_config.simulation_type}")
        click.echo(f"  Band: {sim_config.band}")
        click.echo(f"  Fibers: {sim_config.fibers.mode}")
        return
    
    # Create and run simulator
    simulator = AndesSimulator(sim_config)
    
    try:
        result = simulator.run_simulation()
        click.echo(f"Simulation completed: {sim_config.simulation_type}")
    except Exception as e:
        click.echo(f"Simulation failed: {e}", err=True)
        sys.exit(1)


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