import xarray as xr
import numpy as np
from mafredo.rao import Rao
from mafredo.hyddb1 import Hyddb1

def test_load_dhyd():
    hyd = Hyddb1()
    hyd.load_from(r'files/barge_100_30_4.dhyd')
    rao = hyd._force[0]

    assert rao.n_frequencies == 28
    assert rao.n_wave_directions == 9



def test_Rao_read_nc():

    test = Rao()
    test.wave_force_from_capytaine(r"files/capytaine.nc", "Roll")

    test.add_symmetry_xz()
    test.regrid_omega(np.linspace(0,4,100))
    test.regrid_direction(np.linspace(0, 360, 5))


    print(test.get_value(omega = 0.11, wave_direction=30))

