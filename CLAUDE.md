# ANDES E2E Simulation Instructions

## Primary Tool: andes_simulator

The recommended way to run E2E simulations is using the modern **andes_simulator** framework:

```bash
cd andes_simulator
andes-sim simulate configs/examples/fabry_perot_single_fiber.yaml
```

**Framework Status**: ⚠️ Not yet validated against legacy scripts.

**Documentation**: See `andes_simulator/README.md` for complete usage guide.

## Legacy Scripts (Archived)

The original standalone scripts have been moved to `legacy/` directory. They are kept for reference and validation purposes only.

To run legacy scripts (if needed):
```bash
uv run legacy/FP_code/pyechelle_test_ANDES_fp_single_fiber.py J 33 0
```

See `legacy/README.md` for migration guide and script mappings.

## Directory Structure

```
src/
├── andes_simulator/    # Modern framework (preferred)
├── legacy/             # Archived original scripts
├── HDF/                # ZEMAX optical models
└── SED/                # Spectral energy distributions
```
