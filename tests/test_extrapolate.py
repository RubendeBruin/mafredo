import xarray as xr
import numpy as np
from mafredo.rao import Rao
from mafredo.helpers import MotionMode
from mafredo.hyddb1 import Hyddb1
from numpy.testing import assert_allclose

def test_extrapolate():

    test = Rao.create_from_capytaine_wave_force(r"files/capytaine.nc", MotionMode.ROLL)

    test.expand_symmetry_xz()
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

def test_extrapolate_below_heading_range():
    rao = Rao()
    rao._data = xr.Dataset({
            'amplitude': (['wave_direction', 'omega'], [[10.,10],[0,0]]),
            'phase': (['wave_direction', 'omega'], np.zeros((2,2),dtype=float)),
                    },
            coords={'wave_direction': [10,350.],
                    'omega': [0,4.],
                    }
        )

    rao._data

    rao.add_direction(0)

    rao._data.sel(omega=4)

    val = rao.get_value(0, 0)  # expect 5.0

    assert_allclose(val, 5.0)





