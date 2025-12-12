# Fiber Efficiency in PyEchelle HDF Files

## The Warning Message

When running simulations, you may see:
```
root - WARNING - No spectrograph efficiency data found for fiber 5.
```

This means PyEchelle cannot find fiber efficiency data in the HDF model file.

## What PyEchelle Expects

PyEchelle looks for fiber efficiency data stored in the HDF file at:
- **Path**: `CCD_{ccd_index}/fiber_{fiber_num}` (as an attribute)
- **Attribute name**: `"efficiency"`
- **Format**: A tuple of `(wavelength_array, efficiency_array)`
  - Wavelengths in microns
  - Efficiency values between 0 and 1
  - Both arrays must have the same length

### Code Reference

From PyEchelle's ZEMAX class (`pyechelle/spectrograph.py`):

```python
def get_efficiency(self, fiber: int, ccd_index: int) -> SystemEfficiency:
    try:
        # Try to load grating efficiency from HDF attributes
        ge = GratingEfficiency(...)
    except KeyError:
        ge = ConstantEfficiency("Spectrograph", eff=1.0)

    try:
        # Try to load fiber efficiency from HDF attributes
        self._efficiency[ccd_index][fiber] = SystemEfficiency(
            [
                ge,
                TabulatedEfficiency(
                    "System",
                    *self.h5f[f"CCD_{ccd_index}/fiber_{fiber}"].attrs["efficiency"],
                ),
            ],
            "System",
        )
    except KeyError:
        # This is where the warning comes from
        logging.warning(f"No spectrograph efficiency data found for fiber {fiber}.")
        self._efficiency[ccd_index][fiber] = SystemEfficiency([ge], "System")
```

## Expected HDF Structure

```
HDF File
└── CCD_1/
    ├── Spectrograph (group)
    │   └── attrs: blaze, gpmm
    ├── fiber_1/
    │   ├── attrs: field_shape, efficiency
    │   ├── order109/ (transformation data)
    │   ├── psf_order_109/ (PSF data)
    │   └── ...
    ├── fiber_2/
    │   └── attrs: efficiency
    └── ...
```

## What Happens Without Efficiency Data

If the `efficiency` attribute is missing, PyEchelle:
1. Issues a warning for each fiber
2. Falls back to grating efficiency only
3. Continues simulation (doesn't crash)
4. May produce less realistic results

## How to Add Efficiency Data

### Option 1: Manually Patch Existing HDF Files

Use h5py to add efficiency data to existing HDF files:

```python
import h5py
import numpy as np

# Open HDF file in read-write mode
with h5py.File('src/HDF/ANDES_123_R3.hdf', 'r+') as h5f:
    for fiber_num in range(1, 67):  # 66 fibers for optical bands
        # Define wavelength range for this band
        # R-band example: 0.6-0.8 microns
        wavelengths = np.linspace(0.6, 0.8, 50)

        # Define efficiency values
        # Option A: Constant efficiency
        efficiency = np.ones(50) * 0.85  # 85% constant

        # Option B: Wavelength-dependent efficiency
        # efficiency = 0.5 + 0.4 * np.exp(-((wavelengths - 0.7)**2) / 0.01)

        # Save as attribute (tuple of two arrays)
        fiber_group = h5f[f'CCD_1/fiber_{fiber_num}']
        fiber_group.attrs['efficiency'] = (wavelengths, efficiency)

print("Efficiency data added to all fibers")
```

### Option 2: Add During HDF Generation from ZEMAX

Modify the HDF builder code to include efficiency when creating HDF files from ZEMAX models. This requires modifying PyEchelle's `HDFBuilder` class or post-processing the generated HDF.

## Wavelength Ranges by Band

When adding efficiency data, use appropriate wavelength ranges:

| Band | Wavelength Range (microns) | Wavelength Range (nm) |
|------|---------------------------|----------------------|
| U    | 0.30 - 0.40              | 300 - 400            |
| B    | 0.40 - 0.50              | 400 - 500            |
| V    | 0.50 - 0.60              | 500 - 600            |
| R    | 0.60 - 0.80              | 600 - 800            |
| IZ   | 0.80 - 1.00              | 800 - 1000           |
| Y    | 0.95 - 1.15              | 950 - 1150           |
| J    | 1.15 - 1.35              | 1150 - 1350          |
| H    | 1.50 - 1.80              | 1500 - 1800          |

## Sources of Efficiency Data

Efficiency data can come from:

1. **ZEMAX ray-tracing**: Transmission calculations from optical model
2. **Fiber manufacturer specifications**: Typical fiber transmission curves
3. **Measured calibration data**: From actual instrument calibrations
4. **Constant values**: Simple approximation (e.g., 0.85 for all wavelengths)
5. **Model curves**: Gaussian or polynomial fits to expected performance

## Example: Inspect Existing Efficiency Data

The Y-band file `ANDES_Y01_wFiberEff.hdf` contains efficiency data. Inspect it:

```python
import h5py
import matplotlib.pyplot as plt

with h5py.File('src/HDF/ANDES_Y01_wFiberEff.hdf', 'r') as h5f:
    # Check if fiber 1 has efficiency data
    fiber_group = h5f['CCD_1/fiber_1']

    if 'efficiency' in fiber_group.attrs:
        wavelengths, efficiency = fiber_group.attrs['efficiency']

        print(f"Number of wavelength points: {len(wavelengths)}")
        print(f"Wavelength range: {wavelengths[0]:.3f} - {wavelengths[-1]:.3f} microns")
        print(f"Efficiency range: {efficiency.min():.3f} - {efficiency.max():.3f}")

        # Plot efficiency curve
        plt.plot(wavelengths * 1000, efficiency)  # Convert to nm
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Efficiency')
        plt.title('Fiber 1 Efficiency Curve')
        plt.grid(True)
        plt.show()
    else:
        print("No efficiency data found")
```

## Quick Fix for Testing

To eliminate warnings quickly without detailed efficiency curves:

```python
import h5py
import numpy as np

# Add constant 85% efficiency to all fibers in R-band
with h5py.File('src/HDF/ANDES_123_R3.hdf', 'r+') as h5f:
    wavelengths = np.array([0.6, 0.7, 0.8])  # R-band in microns
    efficiency = np.array([0.85, 0.85, 0.85])  # Constant 85%

    for fiber_num in range(1, 67):  # 66 fibers
        h5f[f'CCD_1/fiber_{fiber_num}'].attrs['efficiency'] = (wavelengths, efficiency)
```

## Available HDF Files

In `src/HDF/`:
- `ANDES_Y01_wFiberEff.hdf` - **WITH** fiber efficiency (Y-band)
- `ANDES_Y01.hdf` - without fiber efficiency
- `ANDES_75fibre_Y.hdf`, `ANDES_75fibre_J.hdf`, `ANDES_75fibre_H.hdf` - IR bands
- `ANDES_123_R3.hdf`, `ANDES_123_IZ3.hdf` - optical bands

## Configuration

To use HDF files with fiber efficiency in your simulations:

```yaml
# In your config YAML
hdf_model: "with_fiber_eff"  # Instead of "default"
```

Currently only Y-band has a `with_fiber_eff` variant configured in `instruments.py`.

## Notes

- Efficiency is optional but provides more realistic simulations
- Each fiber can have different efficiency (important for multi-fiber instruments)
- If you don't have measured data, constant efficiency is better than none
- The simulator will work without efficiency data, just with warnings
