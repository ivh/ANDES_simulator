# Legacy ANDES E2E Scripts

## Status: Archived (Reference Only)

These scripts represent the original ANDES E2E simulation implementation. They have been **archived** and replaced by the modern **andes_simulator** framework.

## Why Archived?

The legacy scripts had several limitations:
- **Code duplication**: ~60-70% of code was duplicated across 15 files
- **Hardcoded parameters**: All settings embedded in scripts, not configurable
- **Fragile path management**: Required running from specific directories
- **No batch processing**: Manual loops required for multiple simulations
- **Inconsistent interfaces**: Different parameter handling across scripts

## Migration to andes_simulator

The new framework provides:
- Single unified command-line interface
- YAML-based configuration
- Batch processing with parallel execution
- Portable (works from any directory)
- Complete feature parity with legacy scripts

For detailed migration instructions, see:
- `../andes_simulator/README.md` - Framework documentation
- `../andes_simulator/ORIGINAL_SCRIPT_MAPPING.md` - Script-by-script mapping

## Directory Structure

```
legacy/
├── FF_code/                  # Flat field scripts (6 files)
├── FP_code/                  # Fabry-Perot scripts (2 files)
├── YJH/                      # HDF generation scripts (3 files)
├── PSF/                      # Post-processing scripts (2 files)
├── pyechelle_test_ANDES.py  # Stellar spectrum simulation
└── Rband_varyZemax.py        # Thermal variation testing
```

## Using Legacy Scripts (If Needed)

If you need to run legacy scripts for comparison or testing:

1. **Set working directory**: Scripts expect to run from specific locations
2. **Check hardcoded paths**: May need editing for your environment
3. **Install dependencies**: Same as andes_simulator (PyEchelle, numpy, etc.)

### Example: Running a flat field simulation

```bash
cd /path/to/ANDES/E2E/src/legacy/FF_code
python pyechelle_test_ANDES_ff.py
```

### Example: Fabry-Perot single fiber

```bash
cd /path/to/ANDES/E2E/src/legacy/FP_code
python pyechelle_test_ANDES_fp_single_fiber.py R 10 5.0
```

## Quick Migration Guide

| Legacy Script | New Framework Command |
|--------------|----------------------|
| `FF_code/pyechelle_test_ANDES_ff.py` | `andes-sim simulate configs/examples/flat_field_all_fibers.yaml` |
| `FF_code/pyechelle_test_ANDES_ff_single_fiber.py` | `andes-sim simulate configs/examples/flat_field_single_fiber.yaml` |
| `FP_code/pyechelle_test_ANDES_fp.py` | `andes-sim simulate configs/examples/fabry_perot_all_fibers.yaml` |
| `PSF/Dkernel.py` | `andes-sim postprocess --mode convolve` |
| `PSF/sumIFU.py` | `andes-sim postprocess --mode sum` |

For complete mappings and all scripts, see `../andes_simulator/ORIGINAL_SCRIPT_MAPPING.md`

## Framework Validation Status

⚠️ **Important**: The andes_simulator framework has not yet been fully validated against these legacy scripts. Before relying exclusively on the new framework:

1. Run comparison tests (see validation plan in cleanup documentation)
2. Verify outputs match within numerical precision
3. Test all bands and modes you need
4. Check performance is acceptable

Once validated, these legacy scripts can be safely deprecated.

## Support

For questions about:
- **Legacy scripts**: Contact original authors or check git history
- **New framework**: See `../andes_simulator/README.md`
- **Migration issues**: Consult `../andes_simulator/ORIGINAL_SCRIPT_MAPPING.md`

---

*Last updated: 2025-12-11*
*Archived as part of E2E codebase cleanup*
