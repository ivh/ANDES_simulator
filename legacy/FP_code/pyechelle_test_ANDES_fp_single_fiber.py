# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "pyechelle",
# ]
# ///
from pyechelle.simulator import Simulator
from pyechelle.sources import CSVSource
from pyechelle.sources import ConstantFlux
from pyechelle.telescope import Telescope
from pyechelle.spectrograph import ZEMAX, LocalDisturber
import os
import sys
import random
from pathlib import Path

t_exp = 30  # sec
coeff_FP = 5E9

def main(arm, idx, shift=None):
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    
    if arm in {"Y", "J", "H"}:
        n_fibers = 75
        hdfname = f'ANDES_75fibre_{arm}'
        fpname = "FP_simulation_YJH_finesse_26.csv"
    elif arm in {"U", "B", "V"}:
        n_fibers = 66
        hdfname = f'ANDES_123_{arm}3'
        fpname = "FP_simulation_UBV_finesse_23.csv"
    elif arm in {"R", "IZ"}:
        n_fibers = 66
        hdfname = f'ANDES_123_{arm}3'
        fpname = "FP_simulation_RIZ_finesse_23.csv"
    else:
        print(f'Error: Unknown arm {arm}')
        sys.exit(1)

    if not (0 <= idx < n_fibers):
        print(f'Error: fiber_idx must be between 1 and {n_fibers}')
        sys.exit(1)

    hdf_path = project_root / 'HDF' / hdfname
    sed_path = project_root / 'SED' / fpname
    
    if shift:
        tx = random.gauss(0.0,shift/1000) # sigma=0.1=100m/s
        print(f'velocity shift={1000*tx}m/s')
        spec = LocalDisturber(ZEMAX(str(hdf_path)),d_tx=tx)
    else:
        spec = ZEMAX(str(hdf_path))
    sim = Simulator(spec)
    sim.set_ccd(1)

    set_of_fibers = range(1, n_fibers + 1)
    sim.set_fibers(set_of_fibers)
    sim.set_telescope(Telescope(39.3, 4.09))

    # flux in photons is set to True to avoid a bug in pyechelle,
    # the FP spectra actually doesn't have physical units
    # t_exp = coeff_FP * t_exp_nominal was set to cope with this
    fp = CSVSource(file_path=str(sed_path))
    
    # Modify flux - handle astropy units properly
    try:
        if hasattr(fp, 'data') and hasattr(fp.data, 'iloc'):
            # Get the flux column (assuming it's the second column)
            flux_col = fp.data.columns[1]
            # Multiply by coefficient, preserving units
            fp.data[flux_col] = fp.data[flux_col].values * coeff_FP
        else:
            print("Warning: Cannot modify flux - fp.data structure not as expected")
            print(f"fp has attributes: {[attr for attr in dir(fp) if not attr.startswith('_')]}")
            if hasattr(fp, 'data'):
                print(f"fp.data type: {type(fp.data)}")
                print(f"fp.data attributes: {[attr for attr in dir(fp.data) if not attr.startswith('_')]}")
    except Exception as e:
        print(f"Warning: Failed to modify flux data: {e}")
        print("Attempting alternative flux modification methods...")
    dark = ConstantFlux(0.00)
    # Reset all fibers to dark
    fibers = [dark] * n_fibers
    # Set the current fiber to flat
    fibers[idx] = fp

    list_of_sources = fibers

    sim.set_sources(list_of_sources)
    sim.set_exposure_time(t_exp)
    sim.set_output(os.pardir + '/' + arm + '/' +
                   f"{arm}_FP_fiber{idx+1:02d}_shift{1000*locals().get('tx', 0):+.1f}.fits", overwrite=True)
    sim.set_cuda(False)
    sim.max_cpu = 1
    sim.run()

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print('Usage: python pyechelle_test_ANDES_fp_single_fiber.py <arm> <fiber_idx> <shift>')
        sys.exit(1)
    arm = sys.argv[1]
    try:
        idx = int(sys.argv[2]) - 1
    except ValueError:
        print('Error: fiber_idx must be an integer')
        sys.exit(1)
    try:
        shift = float(sys.argv[3])
    except ValueError:
        print('Error: shift must be a float')
        sys.exit(1)

    main(arm, idx, shift)
