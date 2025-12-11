from pyechelle.simulator import Simulator
from pyechelle.sources import Constant
from pyechelle.telescope import Telescope
from pyechelle.spectrograph import ZEMAX
import os

arm = 'R'
version = '04_P0'
t_exp = 30  # sec

odd_even = 1  # 0 = even | 1 = odd fibes are obscured

if odd_even:
    odd_even_str = 'even'  # odd_even = 1 -> True when idx fiber is odd, and odd fibers are obscured
else:
    odd_even_str = 'odd'  # odd_even = 0

sim = Simulator(ZEMAX("ANDES_" + arm + version + "_wFiberEff"))

if arm in {"Y", "J", "H"}:
    n_fibers = 75

if arm in {"R", "IZ", "U", "B", "V"}:
    n_fibers = 66

sim.set_ccd(1)
set_of_fibers = range(1, n_fibers + 1)

sim.set_fibers(set_of_fibers)

sim.set_telescope(Telescope(39.3, 4.09))

science = Constant(0.001)  # 0.001

fibers = [science] * n_fibers

for idx, fiber in enumerate(fibers):
    if idx % 2 == odd_even:
        fibers[idx] = Constant(0)

list_of_sources = fibers

sim.set_sources(list_of_sources)

sim.set_exposure_time(t_exp)
sim.set_output(os.path.dirname(os.getcwd()) + '/' + arm + '/' + arm + version + '/' + arm + version + '_DRS_FF_' + odd_even_str + '_' + str(t_exp) + 's.fits', overwrite=True)
sim.set_cuda(True)
image = sim.run()
