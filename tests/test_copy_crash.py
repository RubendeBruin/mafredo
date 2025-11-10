import pytest
import xarray as xr
from mafredo import Hyddb1, Rao, MotionMode


def test_copy_crash(data_path):
    filename = data_path / "barge_100_30_4.dhyd"
    hyd = Hyddb1.create_from(filename)

    hyd2 = hyd.copy()  # This should not crash

    hyd2.regrid_omega([0.1, 0.2, 0.3])  # Regridding should work without issues

    _hyd3 = Hyddb1.create_from(filename)


def test_copy_crash_minimal(
    data_path,
):
    filename = data_path / "barge_100_30_4.dhyd"

    with xr.open_dataarray(filename, group="damping", engine="h5netcdf") as ds:
        damping = ds

    # make a deep-copy and interpolate the copy
    new_damping = damping.copy(deep=True)
    new_damping = new_damping.interp(omega=[0.1, 0.2, 0.3], method="linear")

    # opening "damping" again - ok
    with xr.open_dataarray(filename, group="damping", engine="h5netcdf") as ds:
        damping = ds
    #
    # # opening "mass" again - crashes
    # with xr.open_dataarray(filename, group="mass", engine="netcdf4") as ds:
    #     _mass = ds
