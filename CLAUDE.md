# ANDES E2E Simulation Instructions

## Primary Tool: andes_simulator

The recommended way to run E2E simulations is using the modern **andes_simulator** framework.

**Framework Status**: Under active development.

### Available Commands

```bash
# Flat field calibrations
uv run andes-sim flat-field --band R --mode all --output-dir ../R/

# Fabry-Perot wavelength calibration
uv run andes-sim fabry-perot --band R --mode single --fiber 21 --output-dir ../R/

# Laser Frequency Comb calibration
uv run andes-sim lfc --band R --mode single --fiber 21 --output-dir ../R/

# Stellar spectrum
uv run andes-sim spectrum --band R --spectrum SED/star.csv --fiber 21 --output-dir ../R/

# List all bands
uv run andes-sim list-bands
```

### LFC (Laser Frequency Comb)

The `lfc` command generates wavelength calibration frames with unresolved emission lines equidistant in velocity (constant delta-lambda/lambda), mimicking real laser frequency combs:

```bash
uv run andes-sim lfc --band R --mode single --fiber 21 --scaling 1e5 --output-dir ../R/
```

Options:
- `--scaling`: Flux per line in ph/s (default 1e5)
- `--flux`: Brightness multiplier (default 1.0)
- Lines are spaced ~30-50 km/s apart, giving ~100-150 lines per spectral order

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
