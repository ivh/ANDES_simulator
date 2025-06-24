from pyechelle.simulator import Simulator
from pyechelle.sources import CSVSource as CSV
from pyechelle.sources import ConstantPhotonFlux as Constant
from pyechelle.telescope import Telescope
from pyechelle.spectrograph import ZEMAX
import os, sys
import astropy.units as u
from pathlib import Path

t_exp = 30  # sec
spec_scale= 5E3

def main(arm, idx, specname, shift=None):
    script_dir = Path(__file__).parent
    project_root = script_dir
    
    if arm in {"Y", "J", "H"}:
        n_fibers = 75
        hdfname = f'ANDES_75fibre_{arm}'
    elif arm in {"U", "B", "V"}:
        n_fibers = 66
        hdfname = f'ANDES_123_{arm}3'
    elif arm in {"R", "IZ"}:
        n_fibers = 66
        hdfname = f'ANDES_123_{arm}3'
    else:
        print(f'Error: Unknown arm {arm}')
        sys.exit(1)

    if not (0 <= idx < n_fibers):
        print(f'Error: fiber_idx must be between 1 and {n_fibers}')
        sys.exit(1)

    hdf_path = project_root / 'HDF' / hdfname
    spec = ZEMAX(str(hdf_path))
    sim = Simulator(spec)
    sim.set_ccd(1)

    set_of_fibers = range(1, n_fibers + 1)
    sim.set_fibers(set_of_fibers)
    sim.set_telescope(Telescope(39.3, 4.09))

    # flux in photons is set to True to avoid a bug in pyechelle,
    # the FP spectra actually doesn't have physical units
    # t_exp = coeff_FP * t_exp_nominal was set to cope with this
    inspec = CSV(specname, wavelength_units="nm", flux_units=u.Unit("ph/s"))
    # Use to_numpy() to get the raw values without units
    raw_flux = inspec.data["flux"].to_numpy()
    # Multiply by scalar
    scaled_flux = raw_flux * spec_scale
    # Assign back to the dataframe (pandas will handle the units)
    inspec.data["flux"] = scaled_flux
    
    dark = Constant(0.00)
    # Reset all fibers to dark
    fibers = [dark] * n_fibers
    # Set the current fiber to flat
    fibers[idx] = inspec

    list_of_sources = fibers

    sim.set_sources(list_of_sources)
    sim.set_exposure_time(t_exp)
    sim.set_output(os.pardir + '/' + arm + '/' +
                   f"{arm}_FP_fiber{idx+1:02d}.fits", overwrite=True)
    sim.set_cuda(False)
    sim.max_cpu = 1
    sim.run()

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print('Usage: python pyechelle_test_ANDES_fp_single_fiber.py <arm> <fiber_idx> <specname>')
        sys.exit(1)
    arm = sys.argv[1]
    try:
        idx = int(sys.argv[2]) - 1
    except ValueError:
        print('Error: fiber_idx must be an integer')
        sys.exit(1)
    
    specname = (sys.argv[3])
    
    main(arm, idx, specname)
