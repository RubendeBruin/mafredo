import xarray as xr
from mafredo.helpers import dof_names_to_numbers
from numpy.testing import assert_almost_equal

def test_order_dofs():
    ds = xr.open_dataarray(r'files/barge_100_30_4.dhyd', group="mass", engine='netcdf4')
    fixed = dof_names_to_numbers(ds)

    expected = [[1189.62457,0.00000,0.00000,-0.00000,61686.01166,0.00000],
[0.00000,4523.11256,0.00000,-28139.12490,0.00000,0.00000],
[-0.00000,-0.00000,65963.95710,0.00000,-0.00000,-0.00000],
[-0.00000,-23770.56834,-0.00000,1272439.41489,-0.00000,-0.00000],
[54743.93852,0.00000,-0.00000,0.00000,29171788.90566,0.00000],
[0.00000,0.00000,0.00000,-0.00000,-0.00000,3122606.88051]]

    A0 = fixed.isel(omega=0)

    assert_almost_equal(A0, expected, decimal=5)

