# Original Script to New Framework Mapping

This document maps the original ANDES E2E scripts to the new consolidated framework.

## Script Equivalents

### HDF Generation Scripts

| Original Script | New Framework Command | Configuration File |
|----------------|----------------------|-------------------|
| `YJH/MakeHDF_Yband.py` | `andes-sim generate-hdf --band Y` | `configs/examples/hdf_generation.yaml` |
| `YJH/MakeHDF_Jband.py` | `andes-sim generate-hdf --band J` | Modified for J-band |
| `YJH/MakeHDF_Hband.py` | `andes-sim generate-hdf --band H` | Modified for H-band |

### Flat Field Scripts

| Original Script | New Framework Command | Configuration File |
|----------------|----------------------|-------------------|
| `FF_code/pyechelle_test_ANDES_ff_single_fiber.py` | `andes-sim flat-field --band Y --mode single --fiber N` | `configs/examples/flat_field_single_fiber.yaml` |
| `FF_code/pyechelle_test_ANDES_ff_even_odd.py` | `andes-sim flat-field --band R --mode even_odd` | `configs/examples/flat_field_even_odd.yaml` |
| `FF_code/pyechelle_test_ANDES_ff_calib.py` | `andes-sim flat-field --band IZ --mode calib` | Custom config needed |
| `FF_code/pyechelle_test_ANDES_ff_first_slit.py` | `andes-sim flat-field --band IZ --mode first_slit` | `configs/examples/flat_field_slits.yaml` |
| `FF_code/pyechelle_test_ANDES_ff_second_slit.py` | `andes-sim flat-field --band IZ --mode second_slit` | Modified slit config |
| `FF_code/pyechelle_test_ANDES_ff.py` | `andes-sim flat-field --band IZ --mode all` | Standard all-fiber config |

### Fabry-Perot Scripts

| Original Script | New Framework Command | Configuration File |
|----------------|----------------------|-------------------|
| `FP_code/pyechelle_test_ANDES_fp.py` | `andes-sim fabry-perot --band Y --mode all` | `configs/examples/fabry_perot_all_fibers.yaml` |
| `FP_code/pyechelle_test_ANDES_fp_single_fiber.py` | `andes-sim fabry-perot --band J --mode single --fiber 33 --velocity-shift 100` | `configs/examples/fabry_perot_single_fiber.yaml` |

### Spectrum Simulation Scripts

| Original Script | New Framework Command | Configuration File |
|----------------|----------------------|-------------------|
| `pyechelle_test_ANDES.py` | `andes-sim spectrum --band J --spectrum SED/star.csv --fiber 33` | `configs/examples/spectrum_simulation.yaml` |

### Thermal Variation Scripts

| Original Script | New Framework Command | Configuration File |
|----------------|----------------------|-------------------|
| `Rband_varyZemax.py` | `andes-sim run-config --config thermal_config.yaml` | `configs/examples/thermal_variation.yaml` |

### Post-Processing Scripts

| Original Script | New Framework Command | Configuration File |
|----------------|----------------------|-------------------|
| `PSF/Dkernel.py` | `andes-sim psf-process --band Y --kernel-size 4,4 --fwhm 3.2 --edge-blank left` | CLI parameters |
| `PSF/sumIFU.py` | `andes-sim combine --band Y --mode all --input-pattern "Y_FP_fiber{fib:02d}_*.fits"` | CLI parameters |

## Migration Examples

### Running Single Fiber Flat Field (Original)
```bash
cd YJH
python ../FF_code/pyechelle_test_ANDES_ff_single_fiber.py
# (with hardcoded parameters in script)
```

### Running Single Fiber Flat Field (New Framework)
```bash
# Command line approach
andes-sim flat-field --band Y --mode single --fiber 1 --flux 0.001 --exposure 1

# Configuration file approach  
andes-sim run-config --config configs/examples/flat_field_single_fiber.yaml

# Batch approach for all fibers
andes-sim run-config --config configs/batch/all_fibers_ff.yaml
```

### Running Fabry-Perot with Velocity Shift (Original)
```bash
python FP_code/pyechelle_test_ANDES_fp_single_fiber.py J 33 100
```

### Running Fabry-Perot with Velocity Shift (New Framework)
```bash
andes-sim fabry-perot --band J --mode single --fiber 33 --velocity-shift 100
```

### Running PSF Processing (Original)
```bash
python PSF/Dkernel.py Y 5,5,3.2,left
```

### Running PSF Processing (New Framework)
```bash
andes-sim psf-process --band Y --kernel-size 5,5 --fwhm 3.2 --edge-blank left --input-pattern "Y_FP_fiber{fib:02d}_shift*.fits"
```

## Advanced Features Not in Original Scripts

### Batch Processing
```bash
# Run all fibers for flat field
andes-sim flat-field --band Y --mode single --config batch_single_fibers.yaml

# Run thermal variations
andes-sim run-config --config thermal_sweep_config.yaml
```

### Configuration Management
```bash
# Create template configurations
andes-sim create-templates

# Validate configuration
andes-sim run-config --config my_config.yaml --dry-run
```

### Post-Processing Integration
```bash
# Generate, then process in one workflow
andes-sim fabry-perot --band Y --mode all
andes-sim psf-process --band Y --kernel-size 4,4 --fwhm 3.2
andes-sim combine --band Y --mode all
```

## Python API Usage

### Original Style
```python
# Multiple separate scripts with hardcoded parameters
```

### New Framework Style
```python
from andes_simulator import AndesSimulator, SimulationConfig

# Quick creation
sim = AndesSimulator.quick_flat_field('Y', 'single', fiber=1)
result = sim.run_simulation()

# Full configuration
config = SimulationConfig.from_yaml('my_config.yaml')
sim = AndesSimulator(config)
result = sim.run_simulation()

# Batch processing
from andes_simulator.scripts.batch_runner import BatchRunner
runner = BatchRunner()
results = runner.quick_fiber_sweep('Y', 'fabry_perot')
```

## Benefits of New Framework

1. **Unified Interface**: Single entry point for all simulation types
2. **Configuration Management**: YAML-based configuration with validation
3. **Batch Processing**: Parallel execution of multiple simulations
4. **Portable Paths**: Works from any directory without hardcoded paths
5. **Extensible**: Easy to add new simulation types and features
6. **Reproducible**: Configuration files ensure reproducibility
7. **Post-Processing Integration**: Built-in PSF and fiber combination tools
8. **Error Handling**: Better error reporting and recovery
9. **Documentation**: Self-documenting configuration and CLI help
10. **Testing**: Dry-run capability for validation without execution