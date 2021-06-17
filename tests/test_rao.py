import xarray as xr
import numpy as np

from mafredo.helpers import MotionMode
from mafredo.rao import Rao
from mafredo.hyddb1 import Hyddb1

def test_load_dhyd():
    hyd = Hyddb1.create_from(r'files/barge_100_30_4.dhyd')
    rao = hyd._force[0]

    assert rao.n_frequencies == 28
    assert rao.n_wave_directions == 9

def test_Rao_read_nc():

    test = Rao()
    test.wave_force_from_capytaine(r"files/capytaine.nc", MotionMode.ROLL)

    test.add_symmetry_xz()
    test.regrid_omega(np.linspace(0,4,100))
    test.regrid_direction(np.linspace(0, 360, 5))

    print(test.get_value(omega = 0.11, wave_direction=30))

def test_rao_set_data():

    test = Rao()

    # create some dummy data
    headings = np.arange(0,181,step=45)
    n_headings = len(headings)

    omegas = np.linspace(0,2,num=35)
    n_omegas = len(omegas)

    amplitudes = np.random.random((n_headings, n_omegas))
    phases = np.zeros((n_headings, n_omegas))

    test.set_data(headings, omegas, amplitudes, phases)
