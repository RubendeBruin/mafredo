import xarray as xr
import numpy as np
from mafredo.rao import Rao
from mafredo.hyddb1 import Hyddb1

def test_extrapolate():

    test = Rao()
    test.wave_force_from_capytaine(r"files/capytaine.nc", "Roll")

    test.add_symmetry_xz()
    test.regrid_omega(np.linspace(0,4,100))
    test.regrid_direction(np.linspace(0, 360, 5))

    value = test.get_value(omega = 0, wave_direction=359)

    assert not np.isnan(value)


def test_extrapolate_on_addedmass():
    hyd = Hyddb1.create_from_capytaine(r"files/capytaine.nc")
    # hyd.add_frequencies([0])

    assert not np.any(np.isnan(hyd.amass(0)))


def test_extrapolate_on_damping():
    hyd = Hyddb1.create_from_capytaine(r"files/capytaine.nc")
    # hyd.add_frequency(0)
    assert not np.any(np.isnan(hyd.damping(0)))

def test_extrapolate_on_force():
    hyd = Hyddb1.create_from_capytaine(r"files/capytaine.nc")
    hyd.add_frequency(0)
    assert not np.any(np.isnan(hyd.force(0,90)))

def test_extrapolate_on_force_dir():
    hyd = Hyddb1.create_from_capytaine(r"files/capytaine.nc")
    hyd.add_direction(17)
    assert not np.any(np.isnan(hyd.force(1,17)))  # automatically added


