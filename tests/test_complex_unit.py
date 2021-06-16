from mafredo.rao import _complex_unit_add,_complex_unit_delete,_complex_unit_to_phase
import xarray as xr
import numpy as np
from numpy.testing import assert_allclose

def test_raw():
    a = xr.Dataset({
            'amplitude': (['wave_direction', 'omega'], np.zeros((2,2))),
            'phase': (['wave_direction', 'omega'], np.zeros((2,2))),
                    },
            coords={'wave_direction': [0,180],
                    'omega': [0,4],
                    }
        )

    a['phase'].values = [[np.pi, 2*np.pi],
                         [0.5*np.pi,  -0.5*np.pi]]

    _complex_unit_add(a)

    print(a['complex_unit'].values)

    expected = [[-1, 1],
                [1j, -1j]]

    assert_allclose(a['complex_unit'].values, expected)

def test_there_and_back():
    a = xr.Dataset({
        'amplitude': (['wave_direction', 'omega'], np.zeros((2, 2))),
        'phase': (['wave_direction', 'omega'], np.zeros((2, 2))),
    },
        coords={'wave_direction': [0, 180],
                'omega': [0, 4],
                }
    )

    expected = [[np.pi, 2 * np.pi],
                         [0.5 * np.pi, -0.5 * np.pi]]

    a['phase'].values = expected
    _complex_unit_add(a)
    _complex_unit_to_phase(a)

    # to avoid the 2*pi phase jumps, simply use the cos and the sine

    assert_allclose(np.cos(a['phase'].values),
                    np.cos(expected))
    assert_allclose(np.sin(a['phase'].values),
                    np.sin(expected))

def test_cleanup():
    a = xr.Dataset({
        'amplitude': (['wave_direction', 'omega'], np.zeros((2, 2))),
        'phase': (['wave_direction', 'omega'], np.zeros((2, 2))),
    },
        coords={'wave_direction': [0, 180],
                'omega': [0, 4],
                }
    )
    _complex_unit_add(a)
    a = _complex_unit_delete(a)

    print(a.variables)

    assert('complex_unit' not in a.variables)





