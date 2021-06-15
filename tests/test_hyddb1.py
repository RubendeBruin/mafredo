import numpy as np
import pytest
from numpy.testing import assert_almost_equal

from mafredo.hyddb1 import Hyddb1



def gimme():
    hyd = Hyddb1()
    hyd.load_from(r'files/barge_100_30_4.dhyd')
    return hyd

def test_load_dhyd():
    hyd = Hyddb1()
    hyd.load_from(r'files/barge_100_30_4.dhyd')
    hyd.plot()

    print(hyd._mass)

def test_load_nc():
    hyd = Hyddb1()
    hyd.load_from_capytaine(r"files/capytaine.nc")

    omega = 0.01

    mass = hyd.amass(omega=omega)
    damping = hyd.damping(omega=omega)
    force = hyd.force(omega=omega, wave_direction=90)

def test_n_frequencies():
    hyd = gimme()
    assert hyd.n_frequencies == 28

def test_n_frequencies_error():
    hyd = Hyddb1()
    with pytest.raises(ValueError):
        assert hyd.n_wave_directions == 9

def test_n_wave_directions():
    hyd = gimme()
    assert hyd.n_wave_directions == 9

def test_wave_directions():
    hyd = gimme()

    from numpy.testing import assert_almost_equal

    assert_almost_equal(hyd.wave_directions , [  0.  , 22.5 , 45. ,  67.5,  90.,  112.5, 135. , 157.5 ,180. ])

def test_n_wave_directions_error():
    hyd = Hyddb1()
    with pytest.raises(ValueError):
        assert hyd.n_wave_directions == 9

def test_write_hyd_file():
    hyd = gimme()
    hyd.to_hyd_file(r'c:\data\test.hyd', hydrostatics=None)

def test_add_addedmass():
    hyd = Hyddb1()
    hyd.set_amass(omega = 1, m6x6 = np.ones((6,6)))

    # by default omega=0 is set to zeros
    # check if interpolation to omega = 0.2 yields 0.2
    actual02 = hyd.amass(0.2)
    expected = 0.2 * np.ones((6,6))

    assert_almost_equal(actual02, expected)

def test_add_damping():
    hyd = Hyddb1()
    hyd.set_damping(omega = 1, m6x6 = np.ones((6,6)))

    actual03 = hyd.damping(0.3)
    expected = 0.3 * np.ones((6,6))

    assert_almost_equal(actual03, expected)

def test_read_hyd():
    hyd = Hyddb1()
    hyd.load_from_hyd(r'files/barge.hyd')

    hyd.plot()