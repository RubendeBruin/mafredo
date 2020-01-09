import xarray as xr
import numpy as np
from mafredo.rao import Rao

def test_Rao_read_nc():

    test = Rao()
    test.wave_force_from_capytaine(r"files/capytaine.nc", "Roll")

    test.add_symmetry_xz()
    test.regrid_omega(np.linspace(0,4,100))
    test.regrid_direction(np.linspace(0, 360, 5))


    print(test.get_value(omega = 0.11, wave_direction=30))

