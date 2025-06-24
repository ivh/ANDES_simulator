from pyechelle.simulator import Simulator
from pyechelle.sources import CSV
from pyechelle.telescope import Telescope
from pyechelle.spectrograph import ZEMAX
import os
from pathlib import Path

arm = 'Y'
prefix = 'ANDES_75fibre'
t_exp_nominal = 30  # sec
coeff_FP = 100000000
t_exp = coeff_FP * t_exp_nominal

script_dir = Path(__file__).parent
project_root = script_dir.parent
hdf_path = project_root / 'HDF' / f'{prefix}_{arm}'
sim = Simulator(ZEMAX(str(hdf_path)))

sim.set_ccd(1)

if arm in {"Y", "J", "H"}:
    n_fibers = 75

if arm in {"R", "IZ", "U", "B", "V"}:
    n_fibers = 66

set_of_fibers = range(1, n_fibers + 1)

sim.set_fibers(set_of_fibers)

sim.set_telescope(Telescope(39.3, 4.09))

# flux in photons is set to True to avoid a bug in pyechelle,
# the FP spectra actually doesn't have physical units
# t_exp = coeff_FP * t_exp_nominal was set to cope with this
sed_path = project_root / 'SED' / 'FP_simulation_YJH_finesse_26.csv'
fp = CSV(filepath=str(sed_path),
         wavelength_unit="nm", flux_in_photons=True)

fibers = [fp] * n_fibers

list_of_sources = fibers

sim.set_sources(list_of_sources)

sim.set_exposure_time(t_exp)
sim.set_output(os.pardir + '/' + arm + '/' +
               prefix + '_FP_' + str(t_exp_nominal) + 's.fits', overwrite=True)
sim.set_cuda(False)
sim.max_cpu=4
image = sim.run()
