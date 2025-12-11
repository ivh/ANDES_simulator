"""
Batch processing utilities for ANDES simulations.

Provides tools for running multiple simulations in parallel
and managing large-scale simulation campaigns.
"""

import logging
import yaml
from pathlib import Path
from typing import List, Dict, Any, Optional, Union
import concurrent.futures
from dataclasses import asdict

from ..core.simulator import AndesSimulator
from ..core.config import SimulationConfig


class BatchRunner:
    """
    Manages batch execution of ANDES simulations.
    
    Supports parallel execution of multiple simulation configurations
    with progress tracking and error handling.
    """
    
    def __init__(self, 
                 project_root: Optional[Path] = None,
                 max_workers: int = 4):
        """
        Initialize batch runner.
        
        Parameters
        ----------
        project_root : Path, optional
            Project root directory
        max_workers : int
            Maximum number of parallel simulation workers
        """
        if project_root is None:
            self.project_root = Path(__file__).parent.parent.parent
        else:
            self.project_root = project_root
        
        self.max_workers = max_workers
        self.logger = logging.getLogger(__name__)
        
        # Track execution state
        self.results = {}
        self.errors = {}
    
    def run_single_simulation(self, 
                             config: Union[SimulationConfig, Dict, Path],
                             simulation_id: Optional[str] = None) -> Dict[str, Any]:
        """
        Run a single simulation from configuration.
        
        Parameters
        ----------
        config : SimulationConfig, Dict, or Path
            Simulation configuration
        simulation_id : str, optional
            Unique identifier for this simulation
            
        Returns
        -------
        Dict
            Simulation result with metadata
        """
        # Convert config to SimulationConfig if needed
        if isinstance(config, Path):
            sim_config = SimulationConfig.from_yaml(config)
            config_source = str(config)
        elif isinstance(config, dict):
            sim_config = SimulationConfig.from_dict(config)
            config_source = "dict"
        else:
            sim_config = config
            config_source = "object"
        
        if simulation_id is None:
            simulation_id = f"{sim_config.band}_{sim_config.simulation_type}_{id(config)}"
        
        self.logger.info(f"Starting simulation: {simulation_id}")
        
        try:
            # Create and run simulator
            simulator = AndesSimulator(sim_config)
            result = simulator.run_simulation()
            
            # Package results
            result_data = {
                'id': simulation_id,
                'status': 'success',
                'config': asdict(sim_config),
                'config_source': config_source,
                'result': result,
                'output_path': str(sim_config.get_output_path()) if hasattr(sim_config, 'get_output_path') else None
            }
            
            self.logger.info(f"Completed simulation: {simulation_id}")
            return result_data
            
        except Exception as e:
            self.logger.error(f"Failed simulation {simulation_id}: {e}")
            return {
                'id': simulation_id,
                'status': 'error',
                'error': str(e),
                'config_source': config_source
            }
    
    def run_batch(self, 
                  configs: List[Union[SimulationConfig, Dict, Path]],
                  simulation_ids: Optional[List[str]] = None,
                  parallel: bool = True) -> Dict[str, Any]:
        """
        Run multiple simulations in batch.
        
        Parameters
        ----------
        configs : List
            List of simulation configurations
        simulation_ids : List[str], optional
            Unique identifiers for each simulation
        parallel : bool
            Whether to run simulations in parallel
            
        Returns
        -------
        Dict
            Batch execution results
        """
        if simulation_ids is None:
            simulation_ids = [None] * len(configs)
        elif len(simulation_ids) != len(configs):
            raise ValueError("Number of simulation IDs must match number of configs")
        
        self.logger.info(f"Starting batch execution: {len(configs)} simulations")
        
        results = []
        
        if parallel and self.max_workers > 1:
            # Parallel execution
            with concurrent.futures.ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                future_to_id = {
                    executor.submit(self.run_single_simulation, config, sim_id): sim_id
                    for config, sim_id in zip(configs, simulation_ids)
                }
                
                for future in concurrent.futures.as_completed(future_to_id):
                    result = future.result()
                    results.append(result)
                    
                    if result['status'] == 'success':
                        self.logger.info(f"✓ {result['id']}")
                    else:
                        self.logger.error(f"✗ {result['id']}: {result.get('error', 'Unknown error')}")
        else:
            # Sequential execution
            for config, sim_id in zip(configs, simulation_ids):
                result = self.run_single_simulation(config, sim_id)
                results.append(result)
                
                if result['status'] == 'success':
                    self.logger.info(f"✓ {result['id']}")
                else:
                    self.logger.error(f"✗ {result['id']}: {result.get('error', 'Unknown error')}")
        
        # Summarize results
        successful = [r for r in results if r['status'] == 'success']
        failed = [r for r in results if r['status'] == 'error']
        
        batch_result = {
            'total': len(configs),
            'successful': len(successful),
            'failed': len(failed),
            'success_rate': len(successful) / len(configs) if configs else 0,
            'results': results,
            'successful_ids': [r['id'] for r in successful],
            'failed_ids': [r['id'] for r in failed]
        }
        
        self.logger.info(f"Batch completed: {len(successful)}/{len(configs)} successful")
        
        return batch_result
    
    def run_from_config_directory(self, 
                                 config_dir: Path,
                                 pattern: str = "*.yaml",
                                 parallel: bool = True) -> Dict[str, Any]:
        """
        Run all configurations from a directory.
        
        Parameters
        ----------
        config_dir : Path
            Directory containing YAML configuration files
        pattern : str
            File pattern to match
        parallel : bool
            Whether to run in parallel
            
        Returns
        -------
        Dict
            Batch execution results
        """
        config_files = list(config_dir.glob(pattern))
        
        if not config_files:
            raise ValueError(f"No configuration files found in {config_dir} matching {pattern}")
        
        self.logger.info(f"Found {len(config_files)} configuration files")
        
        # Use filename (without extension) as simulation ID
        simulation_ids = [f.stem for f in config_files]
        
        return self.run_batch(config_files, simulation_ids, parallel)
    
    def create_fiber_sweep_configs(self, 
                                  base_config: Union[SimulationConfig, Dict],
                                  band: str,
                                  fiber_range: Optional[List[int]] = None) -> List[SimulationConfig]:
        """
        Create configurations for sweeping through individual fibers.
        
        Parameters
        ----------
        base_config : SimulationConfig or Dict
            Base configuration to modify
        band : str
            Spectral band
        fiber_range : List[int], optional
            List of fiber numbers to simulate. If None, uses all fibers.
            
        Returns
        -------
        List[SimulationConfig]
            List of single-fiber configurations
        """
        if isinstance(base_config, dict):
            base_config = SimulationConfig.from_dict(base_config)
        
        from ..core.instruments import get_instrument_config
        
        if fiber_range is None:
            n_fibers = get_instrument_config(band)['n_fibers']
            fiber_range = list(range(1, n_fibers + 1))
        
        configs = []
        
        for fiber_num in fiber_range:
            # Create modified config for this fiber
            config_dict = asdict(base_config)
            config_dict['band'] = band
            config_dict['fibers']['mode'] = 'single'
            config_dict['fibers']['fibers'] = [fiber_num]
            
            # Modify output path to include fiber number
            if 'output' in config_dict:
                output_config = config_dict['output']
                if 'filename_template' in output_config:
                    template = output_config['filename_template']
                    if '{fiber}' not in template:
                        # Add fiber number to filename
                        base, ext = template.rsplit('.', 1)
                        template = f"{base}_fiber{{fiber:02d}}.{ext}"
                        output_config['filename_template'] = template
            
            fiber_config = SimulationConfig.from_dict(config_dict)
            configs.append(fiber_config)
        
        return configs
    
    def create_thermal_sweep_configs(self,
                                   base_config: Union[SimulationConfig, Dict],
                                   band: str,
                                   t_numbers: List[Union[str, int]]) -> List[SimulationConfig]:
        """
        Create configurations for thermal model sweep.
        
        Parameters
        ----------
        base_config : SimulationConfig or Dict
            Base configuration
        band : str
            Spectral band (R or IZ)
        t_numbers : List[str or int]
            List of thermal state T-numbers
            
        Returns
        -------
        List[SimulationConfig]
            List of thermal variant configurations
        """
        if isinstance(base_config, dict):
            base_config = SimulationConfig.from_dict(base_config)
        
        from ..models.thermal import ThermalModelManager
        
        thermal_manager = ThermalModelManager(band)
        configs = []
        
        for t_num in t_numbers:
            try:
                model_name = thermal_manager.select_thermal_model_by_t_number(t_num)
                
                # Create modified config
                config_dict = asdict(base_config)
                config_dict['band'] = band
                config_dict['hdf_model'] = model_name
                config_dict['thermal_state'] = t_num
                
                # Modify output path to include thermal state
                if 'output' in config_dict:
                    output_config = config_dict['output']
                    if 'filename_template' in output_config:
                        template = output_config['filename_template']
                        if '{thermal}' not in template:
                            base, ext = template.rsplit('.', 1)
                            t_str = f"T{t_num:04d}" if isinstance(t_num, int) else str(t_num)
                            template = f"{base}_{t_str}.{ext}"
                            output_config['filename_template'] = template
                
                thermal_config = SimulationConfig.from_dict(config_dict)
                configs.append(thermal_config)
                
            except ValueError as e:
                self.logger.warning(f"Skipping T-number {t_num}: {e}")
        
        return configs
    
    def save_batch_report(self, 
                         batch_result: Dict[str, Any],
                         output_path: Path) -> None:
        """
        Save batch execution report to file.
        
        Parameters
        ----------
        batch_result : Dict
            Result from run_batch
        output_path : Path
            Output file path
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Create summary report
        report = {
            'batch_summary': {
                'total_simulations': batch_result['total'],
                'successful': batch_result['successful'],
                'failed': batch_result['failed'],
                'success_rate': batch_result['success_rate']
            },
            'successful_simulations': batch_result['successful_ids'],
            'failed_simulations': batch_result['failed_ids'],
            'detailed_results': batch_result['results']
        }
        
        # Save as YAML
        with open(output_path, 'w') as f:
            yaml.dump(report, f, default_flow_style=False, sort_keys=False)
        
        self.logger.info(f"Batch report saved: {output_path}")
    
    @classmethod
    def quick_fiber_sweep(cls,
                         band: str,
                         simulation_type: str = "fabry_perot",
                         fiber_range: Optional[List[int]] = None,
                         **config_kwargs) -> Dict[str, Any]:
        """
        Quick setup and execution of fiber sweep.
        
        Parameters
        ----------
        band : str
            Spectral band
        simulation_type : str
            Type of simulation
        fiber_range : List[int], optional
            Fibers to simulate
        **config_kwargs
            Additional configuration parameters
            
        Returns
        -------
        Dict
            Batch execution results
        """
        from ..core.config import SimulationConfig, SourceConfig, FiberConfig
        
        # Create base configuration
        if simulation_type == "fabry_perot":
            source_config = SourceConfig(type="fabry_perot", scaling_factor=5e9)
        elif simulation_type == "flat_field":
            source_config = SourceConfig(type="constant", flux=0.001)
        else:
            raise ValueError(f"Unsupported simulation type for quick sweep: {simulation_type}")
        
        base_config = SimulationConfig(
            simulation_type=simulation_type,
            band=band,
            source=source_config,
            fibers=FiberConfig(mode="single"),
            **config_kwargs
        )
        
        # Create batch runner and execute
        runner = cls()
        configs = runner.create_fiber_sweep_configs(base_config, band, fiber_range)
        
        return runner.run_batch(configs)