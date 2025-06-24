from pyechelle.simulator import Simulator
from pyechelle.sources import Constant
from pyechelle.telescope import Telescope
from pyechelle.spectrograph import ZEMAX
import os

arm = 'Y'
version = '01'
t_exp = 1  # sec

sim = Simulator(ZEMAX("ANDES_" + arm + version + "_wFiberEff"))

sim.set_ccd(1)

if arm in {"Y", "J", "H"}:
    n_fibers = 75

if arm in {"R", "IZ", "U", "B", "V"}:
    n_fibers = 66

set_of_fibers = range(1, n_fibers + 1)

sim.set_fibers(set_of_fibers)

sim.set_telescope(Telescope(39.3, 4.09))

flat = Constant(0.001)  # 0.001
dark = Constant(0)

for idx in range(n_fibers):
    # Reset all fibers to dark
    fibers = [dark] * n_fibers
    # Set the current fiber to flat
    fibers[idx] = flat

    list_of_sources = fibers

    sim.set_sources(list_of_sources)

    sim.set_exposure_time(t_exp)
    sim.set_output(os.path.dirname(os.getcwd()) + '/' + arm + '/' + arm + version + '/' +
                   arm + version + '_DRS_FF_fiber_' + str(idx) + '_' + str(t_exp) + 's.fits', overwrite=True)
    sim.set_cuda(True)
    image = sim.run()
