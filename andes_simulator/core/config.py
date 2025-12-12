"""
Configuration management for ANDES simulations.

Handles loading and validation of simulation configurations from YAML files
and command-line parameters.
"""

import yaml
import logging
from pathlib import Path
from typing import Dict, List, Any, Union, Optional
from dataclasses import dataclass, field

from .instruments import get_instrument_config, get_all_bands


@dataclass
class OutputConfig:
    """Configuration for simulation outputs."""
    directory: str = "../{band}/"
    filename_template: str = "{band}_{type}_{exposure}s.fits"
    overwrite: bool = True


@dataclass
class SourceConfig:
    """Configuration for source spectrum."""
    type: str = "constant"  # constant, csv, fabry_perot
    flux: float = 1.0
    scaling_factor: float = 1.0
    filepath: Optional[str] = None  # For CSV sources
    wavelength_unit: str = "nm"
    flux_unit: str = "ph/s"


@dataclass
class FiberConfig:
    """Configuration for fiber illumination."""
    mode: str = "all"  # all, single, even_odd, first_slit, second_slit, custom
    fibers: Union[str, List[int]] = "all"  # "all", fiber numbers, or ranges
    skip_fibers: List[int] = field(default_factory=list)


@dataclass
class PSFConfig:
    """Configuration for PSF convolution."""
    enabled: bool = False
    kernel_size: tuple = (4, 4)  # (dimx, dimy)
    sigma: float = 1.0
    edge_blank: str = "random"  # none, top, bottom, left, right, random


@dataclass
class SimulationConfig:
    """Main simulation configuration."""
    # Simulation type and basic parameters
    simulation_type: str = "flat_field"  # flat_field, fabry_perot, spectrum, hdf_generation
    band: str = "Y"
    exposure_time: float = 30.0
    
    # Instrument configuration
    hdf_model: Optional[str] = None  # Use default if None
    use_cuda: bool = False
    max_cpu: int = 1  # Use performance cores only (M4 has 10 perf + 4 efficiency cores)
    
    # Source configuration
    source: SourceConfig = field(default_factory=SourceConfig)
    
    # Fiber configuration
    fibers: FiberConfig = field(default_factory=FiberConfig)
    
    # Special features
    velocity_shift: Optional[float] = None  # m/s for Doppler shifts
    thermal_model: Optional[str] = None     # For thermal variations
    
    # Post-processing
    psf: PSFConfig = field(default_factory=PSFConfig)
    combine_fibers: bool = False
    
    # Output configuration
    output: OutputConfig = field(default_factory=OutputConfig)
    
    # Validation and derived properties
    _instrument_config: Optional[Dict] = field(default=None, init=False)
    
    def __post_init__(self):
        """Validate configuration after initialization."""
        self.validate()
        self._instrument_config = get_instrument_config(self.band)
    
    def validate(self):
        """Validate the configuration parameters."""
        # Validate band
        if self.band not in get_all_bands():
            raise ValueError(f"Invalid band '{self.band}'. Available: {get_all_bands()}")
        
        # Validate simulation type
        valid_sim_types = ["flat_field", "fabry_perot", "spectrum", "hdf_generation"]
        if self.simulation_type not in valid_sim_types:
            raise ValueError(f"Invalid simulation_type '{self.simulation_type}'. Available: {valid_sim_types}")
        
        # Validate fiber mode
        valid_fiber_modes = ["all", "single", "even_odd", "first_slit", "second_slit", "custom"]
        if self.fibers.mode not in valid_fiber_modes:
            raise ValueError(f"Invalid fiber mode '{self.fibers.mode}'. Available: {valid_fiber_modes}")
        
        # Validate source type
        valid_source_types = ["constant", "csv", "fabry_perot"]
        if self.source.type not in valid_source_types:
            raise ValueError(f"Invalid source type '{self.source.type}'. Available: {valid_source_types}")
        
        # Check for required parameters based on simulation type
        if self.simulation_type == "spectrum" and self.source.type == "csv":
            if self.source.filepath is None:
                raise ValueError("CSV filepath required for spectrum simulations")
        
        if self.fibers.mode == "single" and isinstance(self.fibers.fibers, str):
            raise ValueError("Specific fiber number required for single fiber mode")
    
    @property
    def instrument_config(self) -> Dict:
        """Get the instrument configuration for this band."""
        if self._instrument_config is None:
            self._instrument_config = get_instrument_config(self.band)
        return self._instrument_config
    
    @property
    def n_fibers(self) -> int:
        """Get the number of fibers for this band."""
        return self.instrument_config['n_fibers']
    
    @property
    def detector_size(self) -> tuple:
        """Get the detector size for this band."""
        return self.instrument_config['detector_size']
    
    def get_fiber_list(self) -> List[int]:
        """
        Get the list of fiber numbers to illuminate based on configuration.
        
        Returns
        -------
        List[int]
            List of 1-based fiber numbers to illuminate
        """
        n_fibers = self.n_fibers
        
        if self.fibers.mode == "all":
            if isinstance(self.fibers.fibers, str) and self.fibers.fibers == "all":
                fibers = list(range(1, n_fibers + 1))
            elif isinstance(self.fibers.fibers, list):
                fibers = self.fibers.fibers
            else:
                raise ValueError("Invalid fiber specification for 'all' mode")
                
        elif self.fibers.mode == "single":
            if isinstance(self.fibers.fibers, int):
                fibers = [self.fibers.fibers]
            elif isinstance(self.fibers.fibers, list) and len(self.fibers.fibers) == 1:
                fibers = self.fibers.fibers
            else:
                raise ValueError("Single fiber mode requires exactly one fiber number")
                
        elif self.fibers.mode == "even_odd":
            # This mode requires special handling in the simulator
            # Return all fibers, but the simulator will handle even/odd logic
            fibers = list(range(1, n_fibers + 1))
            
        elif self.fibers.mode == "first_slit":
            if self.band in ['Y', 'J', 'H']:
                raise ValueError("Slit modes not applicable to YJH bands")
            fibers = list(range(1, 32))  # First 31 fibers
            
        elif self.fibers.mode == "second_slit":
            if self.band in ['Y', 'J', 'H']:
                raise ValueError("Slit modes not applicable to YJH bands")
            fibers = list(range(35, n_fibers + 1))  # Skip cal fibers 32-34
            
        elif self.fibers.mode == "custom":
            if isinstance(self.fibers.fibers, list):
                fibers = self.fibers.fibers
            else:
                raise ValueError("Custom mode requires explicit fiber list")
        else:
            raise ValueError(f"Unknown fiber mode: {self.fibers.mode}")
        
        # Remove skipped fibers
        if self.fibers.skip_fibers:
            fibers = [f for f in fibers if f not in self.fibers.skip_fibers]
        
        return fibers
    
    def get_output_path(self, fiber_num: Optional[int] = None, suffix: str = "") -> Path:
        """
        Generate output file path based on configuration.
        
        Parameters
        ----------
        fiber_num : int, optional
            Fiber number for single-fiber outputs
        suffix : str, optional
            Additional suffix for filename
            
        Returns
        -------
        Path
            Complete output file path
        """
        # Format directory path
        directory = self.output.directory.format(
            band=self.band,
            type=self.simulation_type
        )
        
        # Format filename
        filename_parts = {
            'band': self.band,
            'type': self.simulation_type,
            'exposure': int(self.exposure_time),
            'fiber': f"{fiber_num:02d}" if fiber_num is not None else "",
            'mode': self.fibers.mode
        }
        
        if suffix:
            filename_parts['suffix'] = suffix
        
        # Create filename based on simulation type
        if self.simulation_type == "flat_field":
            if self.fibers.mode == "single" and fiber_num is not None:
                filename = f"{self.band}_FF_fiber{fiber_num:02d}_{int(self.exposure_time)}s.fits"
            elif self.fibers.mode == "even_odd":
                filename = f"{self.band}_FF_{suffix}_{int(self.exposure_time)}s.fits"
            else:
                filename = f"{self.band}_FF_{self.fibers.mode}_{int(self.exposure_time)}s.fits"
        elif self.simulation_type == "fabry_perot":
            if fiber_num is not None:
                shift_str = f"_shift{suffix}" if suffix else ""
                filename = f"{self.band}_FP_fiber{fiber_num:02d}{shift_str}.fits"
            else:
                filename = f"{self.band}_FP_{int(self.exposure_time)}s.fits"
        elif self.simulation_type == "spectrum":
            if fiber_num is not None:
                filename = f"{self.band}_spectrum_fiber{fiber_num:02d}.fits"
            else:
                filename = f"{self.band}_spectrum_{int(self.exposure_time)}s.fits"
        elif self.simulation_type == "hdf_generation":
            filename = f"ANDES_75fibre_{self.band}.hdf" if self.band in ['Y', 'J', 'H'] else f"ANDES_123_{self.band}3.hdf"
        else:
            # Use template
            filename = self.output.filename_template.format(**filename_parts)
        
        return Path(directory) / filename
    
    @classmethod
    def from_yaml(cls, yaml_path: Union[str, Path]) -> 'SimulationConfig':
        """
        Load configuration from a YAML file.
        
        Parameters
        ---------- 
        yaml_path : str or Path
            Path to YAML configuration file
            
        Returns
        -------
        SimulationConfig
            Loaded configuration object
        """
        yaml_path = Path(yaml_path)
        if not yaml_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {yaml_path}")
        
        with open(yaml_path, 'r') as f:
            data = yaml.safe_load(f)
        
        return cls.from_dict(data)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'SimulationConfig':
        """
        Create configuration from a dictionary.
        
        Parameters
        ----------
        data : dict
            Configuration dictionary
            
        Returns
        -------
        SimulationConfig
            Configuration object
        """
        # Handle nested configurations
        if 'source' in data:
            data['source'] = SourceConfig(**data['source'])
        if 'fibers' in data:
            data['fibers'] = FiberConfig(**data['fibers'])
        if 'psf' in data:
            data['psf'] = PSFConfig(**data['psf'])
        if 'output' in data:
            data['output'] = OutputConfig(**data['output'])
        
        return cls(**data)
    
    def to_yaml(self, yaml_path: Union[str, Path]) -> None:
        """
        Save configuration to a YAML file.
        
        Parameters
        ----------
        yaml_path : str or Path
            Output path for YAML file
        """
        yaml_path = Path(yaml_path)
        yaml_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Convert to dictionary, handling dataclasses
        data = {}
        for field_name, field_obj in self.__dataclass_fields__.items():
            value = getattr(self, field_name)
            if hasattr(value, '__dataclass_fields__'):
                # Convert nested dataclass to dict
                data[field_name] = {k: getattr(value, k) for k in value.__dataclass_fields__.keys()}
            elif not field_name.startswith('_'):
                data[field_name] = value
        
        with open(yaml_path, 'w') as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)


def create_template_configs(output_dir: Path) -> None:
    """
    Create template configuration files for different simulation types.
    
    Parameters
    ----------
    output_dir : Path
        Directory to save template files
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Flat field single fiber template
    ff_single = SimulationConfig(
        simulation_type="flat_field",
        band="Y",
        exposure_time=1.0,
        source=SourceConfig(type="constant", flux=0.001),
        fibers=FiberConfig(mode="single", fibers=[1]),
        output=OutputConfig(
            directory="../{band}/",
            filename_template="{band}_FF_fiber{fiber:02d}_{exposure}s.fits"
        )
    )
    ff_single.to_yaml(output_dir / "flat_field_single_fiber.yaml")
    
    # Flat field even/odd template
    ff_even_odd = SimulationConfig(
        simulation_type="flat_field",
        band="R",
        exposure_time=30.0,
        source=SourceConfig(type="constant", flux=0.001),
        fibers=FiberConfig(mode="even_odd"),
        output=OutputConfig(directory="../{band}/")
    )
    ff_even_odd.to_yaml(output_dir / "flat_field_even_odd.yaml")
    
    # Fabry-Perot template
    fp_config = SimulationConfig(
        simulation_type="fabry_perot",
        band="Y",
        exposure_time=30.0,
        source=SourceConfig(type="fabry_perot", scaling_factor=5e9),
        fibers=FiberConfig(mode="all"),
    )
    fp_config.to_yaml(output_dir / "fabry_perot_all_fibers.yaml")
    
    # Spectrum simulation template
    spectrum_config = SimulationConfig(
        simulation_type="spectrum",
        band="J",
        exposure_time=30.0,
        source=SourceConfig(
            type="csv",
            filepath="SED/star_spectrum.csv",
            scaling_factor=5e3
        ),
        fibers=FiberConfig(mode="single", fibers=[33])
    )
    spectrum_config.to_yaml(output_dir / "spectrum_simulation.yaml")
    
    # HDF generation template
    hdf_config = SimulationConfig(
        simulation_type="hdf_generation",
        band="Y",
        source=SourceConfig(type="zemax"),  # Special type for HDF generation
        output=OutputConfig(directory="HDF/")
    )
    hdf_config.to_yaml(output_dir / "hdf_generation.yaml")
    
    logging.info(f"Created template configurations in {output_dir}")
