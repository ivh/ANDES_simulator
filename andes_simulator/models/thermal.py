"""
Thermal model management for ANDES simulations.

Handles selection and management of different thermal states
and optical configurations for studying instrument variations.
"""

import logging
from pathlib import Path
from typing import List, Dict, Any, Optional, Union

from ..core.instruments import get_instrument_config


class ThermalModelManager:
    """
    Manages thermal model variations for ANDES simulations.
    
    Handles selection of different HDF models representing various
    thermal states and mechanical configurations of the instrument.
    """
    
    def __init__(self, 
                 band: str,
                 project_root: Optional[Path] = None):
        """
        Initialize thermal model manager.
        
        Parameters
        ----------
        band : str
            Spectral band (typically R or IZ for thermal studies)
        project_root : Path, optional
            Project root directory
        """
        self.band = band
        self.instrument_config = get_instrument_config(band)
        
        if project_root is None:
            self.project_root = Path(__file__).parent.parent.parent
        else:
            self.project_root = project_root
        
        self.logger = logging.getLogger(__name__)
        
        # Check if this band has thermal variants
        if 'thermal_variants' not in self.instrument_config.get('hdf_models', {}):
            self.logger.warning(f"No thermal variants configured for {band}-band")
            self.thermal_models = []
        else:
            self.thermal_models = self.instrument_config['hdf_models']['thermal_variants']
    
    def get_available_thermal_models(self) -> List[str]:
        """
        Get list of available thermal model names.
        
        Returns
        -------
        List[str]
            List of thermal model identifiers
        """
        return self.thermal_models.copy()
    
    def get_thermal_model_info(self) -> Dict[str, Dict[str, Any]]:
        """
        Get detailed information about available thermal models.
        
        Returns
        -------
        Dict
            Dictionary with thermal model metadata
        """
        info = {}
        
        for model_name in self.thermal_models:
            # Extract information from model name
            t_number = self._extract_t_number(model_name)
            model_type = self._classify_model_type(model_name)
            
            info[model_name] = {
                't_number': t_number,
                'model_type': model_type,
                'hdf_path': self.get_thermal_model_path(model_name),
                'description': self._get_model_description(model_name)
            }
        
        return info
    
    def _extract_t_number(self, model_name: str) -> Optional[str]:
        """Extract T-number (thermal state identifier) from model name."""
        import re
        match = re.search(r'T(\d{4})', model_name)
        return f"T{match.group(1)}" if match else None
    
    def _classify_model_type(self, model_name: str) -> str:
        """Classify the type of thermal model."""
        if 'MC' in model_name:
            return 'monte_carlo'
        elif 'VM' in model_name:
            return 'vignetting_map'
        else:
            return 'standard'
    
    def _get_model_description(self, model_name: str) -> str:
        """Get human-readable description of thermal model."""
        descriptions = {
            'ANDES_full_F18A33_win_jmr_MC_T0019_Rband_p0': 
                'R-band thermal state T0019 with Monte Carlo analysis',
            'Andes_full_F18A33_win_jmr_MC_T0108_Rband_P0_cfg1':
                'R-band thermal state T0108 with Monte Carlo analysis, configuration 1',
            'Andes_F18A33_VM246aa_win_jmr9_MC_T0045_IZband_P0_cf1':
                'IZ-band thermal state T0045 with vignetting map and Monte Carlo',
            'Andes_full_F18A33_win_jmr_MC_T0028_IZband_P0':
                'IZ-band thermal state T0028 with Monte Carlo analysis'
        }
        
        return descriptions.get(model_name, f"Thermal model: {model_name}")
    
    def get_thermal_model_path(self, model_name: str) -> Path:
        """
        Get full path to a thermal model HDF file.
        
        Parameters
        ----------
        model_name : str
            Thermal model name
            
        Returns
        -------
        Path
            Full path to HDF file
        """
        if model_name not in self.thermal_models:
            raise ValueError(f"Unknown thermal model '{model_name}' for {self.band}-band")
        
        hdf_path = self.project_root / 'HDF' / f'{model_name}.hdf'
        return hdf_path
    
    def select_thermal_model_by_t_number(self, t_number: Union[str, int]) -> str:
        """
        Select thermal model by T-number.
        
        Parameters
        ----------
        t_number : str or int
            T-number (e.g., 'T0019', 'T0108', or just 19, 108)
            
        Returns
        -------
        str
            Full thermal model name
        """
        # Normalize T-number format
        if isinstance(t_number, int):
            t_str = f"T{t_number:04d}"
        elif isinstance(t_number, str):
            if not t_number.startswith('T'):
                t_str = f"T{int(t_number):04d}"
            else:
                t_str = t_number
        else:
            raise ValueError(f"Invalid T-number format: {t_number}")
        
        # Find matching model
        for model_name in self.thermal_models:
            if t_str in model_name:
                self.logger.info(f"Selected thermal model {model_name} for {t_str}")
                return model_name
        
        raise ValueError(f"No thermal model found for {t_str} in {self.band}-band")
    
    def get_model_comparison_set(self) -> Dict[str, List[str]]:
        """
        Get sets of thermal models suitable for comparison studies.
        
        Returns
        -------
        Dict
            Dictionary grouping models by comparison type
        """
        r_band_models = [m for m in self.thermal_models if 'Rband' in m]
        iz_band_models = [m for m in self.thermal_models if 'IZband' in m]
        
        comparison_sets = {}
        
        if len(r_band_models) > 1:
            comparison_sets['r_band_thermal'] = r_band_models
        
        if len(iz_band_models) > 1:
            comparison_sets['iz_band_thermal'] = iz_band_models
        
        # Group by analysis type
        mc_models = [m for m in self.thermal_models if 'MC' in m]
        if len(mc_models) > 1:
            comparison_sets['monte_carlo'] = mc_models
        
        return comparison_sets
    
    def validate_thermal_model(self, model_name: str) -> bool:
        """
        Validate that a thermal model exists and is accessible.
        
        Parameters
        ----------
        model_name : str
            Thermal model name
            
        Returns
        -------
        bool
            True if model is valid and accessible
        """
        if model_name not in self.thermal_models:
            self.logger.error(f"Thermal model '{model_name}' not configured for {self.band}-band")
            return False
        
        hdf_path = self.get_thermal_model_path(model_name)
        if not hdf_path.exists():
            self.logger.error(f"Thermal model file not found: {hdf_path}")
            return False
        
        return True
    
    def create_thermal_study_config(self, 
                                   t_numbers: List[Union[str, int]],
                                   base_config: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Create configuration list for thermal variation study.
        
        Parameters
        ----------
        t_numbers : List[str or int]
            List of T-numbers to study
        base_config : Dict
            Base simulation configuration
            
        Returns
        -------
        List[Dict]
            List of configurations for each thermal state
        """
        configs = []
        
        for t_num in t_numbers:
            try:
                model_name = self.select_thermal_model_by_t_number(t_num)
                
                # Create config for this thermal state
                config = base_config.copy()
                config['hdf_model'] = model_name
                config['thermal_state'] = t_num
                
                # Update output path to include thermal state
                if 'output' in config:
                    output_config = config['output'].copy()
                    if 'filename_template' in output_config:
                        template = output_config['filename_template']
                        # Insert thermal state into filename
                        if '{thermal}' not in template:
                            # Add thermal state before file extension
                            base, ext = template.rsplit('.', 1)
                            template = f"{base}_T{t_num:04d}.{ext}" if isinstance(t_num, int) else f"{base}_{t_num}.{ext}"
                            output_config['filename_template'] = template
                    config['output'] = output_config
                
                configs.append(config)
                
            except ValueError as e:
                self.logger.warning(f"Skipping T-number {t_num}: {e}")
        
        return configs
    
    @classmethod
    def get_all_thermal_models(cls, project_root: Optional[Path] = None) -> Dict[str, 'ThermalModelManager']:
        """
        Get thermal model managers for all bands with thermal variants.
        
        Parameters
        ----------
        project_root : Path, optional
            Project root directory
            
        Returns
        -------
        Dict
            Dictionary mapping band names to thermal model managers
        """
        thermal_bands = ['R', 'IZ']  # Bands with thermal model variants
        managers = {}
        
        for band in thermal_bands:
            try:
                manager = cls(band, project_root)
                if manager.thermal_models:  # Only include if thermal models exist
                    managers[band] = manager
            except Exception as e:
                logging.getLogger(__name__).warning(f"Could not create thermal manager for {band}: {e}")
        
        return managers