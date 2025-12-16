"""
CLI utilities for ANDES simulation commands.

Provides shared functionality to reduce code duplication across CLI commands.
"""

import sys
import logging
from pathlib import Path
from typing import Optional, Callable, Any

try:
    import click
except ImportError:
    print("Error: click package required for CLI. Install with: pip install click")
    sys.exit(1)


def run_simulation_command(
    config,
    dry_run: bool,
    dry_run_callback: Callable[[], None],
    success_message: str = "Simulation completed"
) -> Any:
    """
    Execute a simulation with standard dry-run and error handling.
    
    Parameters
    ----------
    config : SimulationConfig
        The simulation configuration
    dry_run : bool
        If True, only show what would be done
    dry_run_callback : callable
        Function to call for dry-run output (should use click.echo)
    success_message : str
        Message to display on successful completion
        
    Returns
    -------
    Any
        Simulation result, or None if dry_run
    """
    if dry_run:
        dry_run_callback()
        return None
    
    # Lazy import to speed up --help
    from ..core.simulator import AndesSimulator
    
    try:
        simulator = AndesSimulator(config)
        result = simulator.run_simulation()
        click.echo(success_message)
        return result
    except Exception as e:
        click.echo(f"Simulation failed: {e}", err=True)
        sys.exit(1)


def build_config_from_options(
    simulation_type: str,
    band: str,
    exposure: float,
    source_type: str,
    fiber_mode: str,
    output_dir: Optional[Path] = None,
    output_name: Optional[str] = None,
    fiber: Optional[int] = None,
    flux: float = 1.0,
    scaling: Optional[float] = None,
    velocity_shift: Optional[float] = None,
    spectrum_path: Optional[Path] = None,
    flux_unit: str = "ph/s/AA",
    hdf: Optional[str] = None,
    wl_min: Optional[float] = None,
    wl_max: Optional[float] = None,
    fib_eff: Optional[str] = None
):
    """
    Build a SimulationConfig from CLI options.
    
    This centralizes the config-building logic shared across commands.
    
    Parameters
    ----------
    simulation_type : str
        Type of simulation (flat_field, fabry_perot, spectrum)
    band : str
        Spectral band
    exposure : float
        Exposure time in seconds
    source_type : str
        Source type (constant, fabry_perot, csv)
    fiber_mode : str
        Fiber illumination mode
    output_dir : Path, optional
        Output directory
    fiber : int, optional
        Specific fiber number for single mode
    flux : float
        Flux level
    scaling : float, optional
        Scaling factor (overrides flux-based calculation if provided)
    velocity_shift : float, optional
        Velocity shift in m/s
    spectrum_path : Path, optional
        Path to spectrum file for CSV sources
    flux_unit : str
        Unit for flux values
        
    Returns
    -------
    SimulationConfig
        Complete simulation configuration
    """
    from ..core.config import SimulationConfig, SourceConfig, FiberConfig, OutputConfig
    
    # Determine fibers based on mode
    if fiber_mode == 'single':
        if fiber is None:
            raise click.BadParameter("--fiber required for single mode")
        fibers = [fiber]
    else:
        fibers = "all"
    
    # Build source config
    source_kwargs = {'type': source_type, 'flux': flux, 'flux_unit': flux_unit}
    if source_type in ('fabry_perot', 'constant'):
        # For FP and constant sources, total scaling = flux × scaling
        effective_scaling = flux * (scaling if scaling is not None else 1.0)
        source_kwargs['scaling_factor'] = effective_scaling
        if source_type == 'constant':
            # For constant sources, the flux IS the scaling factor
            source_kwargs['flux'] = effective_scaling
    elif scaling is not None:
        source_kwargs['scaling_factor'] = scaling
    if spectrum_path is not None:
        source_kwargs['filepath'] = str(spectrum_path)
    
    # Build output config (default to current working directory)
    output_directory = str(output_dir) if output_dir else str(Path.cwd())

    return SimulationConfig(
        simulation_type=simulation_type,
        band=band,
        exposure_time=exposure,
        velocity_shift=velocity_shift,
        hdf_model=hdf,
        wl_min=wl_min,
        wl_max=wl_max,
        fib_eff=fib_eff,
        source=SourceConfig(**source_kwargs),
        fibers=FiberConfig(mode=fiber_mode, fibers=fibers),
        output=OutputConfig(directory=output_directory, filename=output_name)
    )


def format_dry_run_output(config, extra_lines: Optional[list] = None) -> None:
    """
    Print standardized dry-run output.
    
    Parameters
    ----------
    config : SimulationConfig
        The simulation configuration
    extra_lines : list, optional
        Additional lines to print
    """
    click.echo("Dry run - would execute:")
    click.echo(f"  Type: {config.simulation_type}")
    click.echo(f"  Band: {config.band}")
    if config.hdf_model:
        click.echo(f"  HDF model: {config.hdf_model}")
    click.echo(f"  Subslit: {config.fibers.mode}")
    click.echo(f"  Exposure: {config.exposure_time}s")
    
    if config.source.type == "constant":
        click.echo(f"  Flux: {config.source.flux}")
    elif config.source.type == "fabry_perot":
        click.echo(f"  Scaling: {config.source.scaling_factor:.2e}")
    elif config.source.type == "lfc":
        click.echo(f"  LFC flux per line: {config.source.scaling_factor:.2e} ph/s")
    elif config.source.type == "csv":
        click.echo(f"  Spectrum: {config.source.filepath}")
        click.echo(f"  Scaling: {config.source.scaling_factor}")
    
    if config.velocity_shift:
        click.echo(f"  Velocity shift: {config.velocity_shift} m/s")
    if config.wl_min or config.wl_max:
        wl_range = f"{config.wl_min or '...'}-{config.wl_max or '...'} nm"
        click.echo(f"  Wavelength range: {wl_range}")
    if config.fib_eff:
        click.echo(f"  Fiber efficiency: {config.fib_eff}")
    if config.output.filename:
        click.echo(f"  Output: {config.output.filename}")

    if extra_lines:
        for line in extra_lines:
            click.echo(f"  {line}")


def setup_logging(verbose: bool) -> None:
    """Configure logging based on verbosity."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
