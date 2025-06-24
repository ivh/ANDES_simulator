from pyechelle.simulator import Simulator
from pyechelle.sources import Constant
from pyechelle.telescope import Telescope
from pyechelle.spectrograph import ZEMAX
import os

arm = 'IZ'
version = '04_P0'
t_exp = 30  # sec

sim = Simulator(ZEMAX("ANDES_" + arm + version + "_wFiberEff"))

if arm in {"R", "IZ", "U", "B", "V"}:
    n_fibers = 66
else:
    print("The code does not work with YJH. Check the number of fibers for each pseudo-slit.")
    exit()

sim.set_ccd(1)

set_of_fibers = range(1, n_fibers + 1)

sim.set_fibers(set_of_fibers)

sim.set_telescope(Telescope(39.3, 4.09))

science = Constant(0.001)  # 0.001

science_b = [Constant(0)] * 31

dark_b = [Constant(0)]

sim_cal_b = [Constant(0)]
sim_cal_a = [Constant(0)]

dark_a = [Constant(0)]

sceince_a = [science] * 31

list_of_sources = science_b + dark_b + sim_cal_b + sim_cal_a + dark_a + sceince_a

sim.set_sources(list_of_sources)

sim.set_exposure_time(t_exp)
sim.set_output(os.path.dirname(os.getcwd()) + '/' + arm + '/' + arm + version + '/' + arm + version + '_DRS_FF_first_slit_' + str(t_exp) + 's.fits', overwrite=True)

sim.set_cuda(True)
image = sim.run()
