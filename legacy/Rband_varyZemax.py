# /// script
# requires-python = "==3.9.6"
# dependencies = [
#     "pyechelle",
# ]
# ///

from pyechelle.simulator import Simulator
from pyechelle.sources import CSV
from pyechelle.sources import Constant
from pyechelle.telescope import Telescope
from pyechelle.spectrograph import ZEMAX
import os
import sys
from pathlib import Path

t_exp = 30  # sec
coeff_FP = 5E9
n_fibers = 66
arm='R'
cpus=1

#check argv
if len(sys.argv) != 3:
    print('Usage: python Rband_varyZemax.py <tnum> <fibidx>')
    sys.exit(1)
tnum=int(sys.argv[1])
fibidx=int(sys.argv[2])

tnum=f'T{tnum:04d}'
zemaxfiles = ['ANDES_full_F18A33_win_jmr_MC_T0019_Rband_p0.hdf',
               'Andes_full_F18A33_win_jmr_MC_T0108_Rband_P0_cfg1.hdf',
               'Andes_F18A33_VM246aa_win_jmr9_MC_T0045_IZband_P0_cf1.hdf',
               'Andes_full_F18A33_win_jmr_MC_T0028_IZband_P0.hdf']

for zemaxfile in zemaxfiles:
    if tnum in zemaxfile:
        break

script_dir = Path(__file__).parent
project_root = script_dir
hdf_path = project_root / 'HDF' / zemaxfile
sim = Simulator(ZEMAX(str(hdf_path)))
sim.set_ccd(1)

set_of_fibers = range(1, n_fibers + 1)
sim.set_fibers(set_of_fibers)
sim.set_telescope(Telescope(39.3, 4.09))

# flux in photons is set to True to avoid a bug in pyechelle,
# the FP spectra actually doesn't have physical units
# t_exp = coeff_FP * t_exp_nominal was set to cope with this
sed_path = project_root / 'SED' / 'FP_simulation_RIZ_finesse_23.csv'
fp = CSV(filepath=str(sed_path),
        wavelength_unit="nm", flux_in_photons=True)
fp.flux_data *= coeff_FP
dark = Constant(0.00)
# Reset all fibers to dark
fibers = [dark] * n_fibers
fibers[fibidx - 1] = fp

sim.set_sources(fibers)
sim.set_exposure_time(t_exp)
sim.set_output(os.pardir + '/' + arm + '/' +
                f"{arm}_{tnum}_fib{fibidx:02d}.fits", overwrite=True)
sim.set_cuda(False)
sim.max_cpu = cpus
sim.run()
