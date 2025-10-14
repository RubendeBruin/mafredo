import numpy as np
import pytest
import os
import tempfile
from numpy.testing import assert_almost_equal

from mafredo.hyddb1 import Hyddb1


def gimme(data_path):
    hyd = Hyddb1.create_from(data_path / "barge_100_30_4.dhyd")
    return hyd


def test_check_dimensions(data_path):
    hyd = gimme(data_path=data_path)
    hyd._check_dimensions()


def test_load_dhyd(data_path):
    hyd = Hyddb1.create_from(data_path / "barge_100_30_4.dhyd")
    hyd.plot()

    print(hyd._mass)


def test_load_nc(data_path):
    hyd = Hyddb1.create_from_capytaine(data_path / "capytaine.nc")

    omega = 0.01

    _mass = hyd.amass(omega=omega)
    _damping = hyd.damping(omega=omega)
    _force = hyd.force(omega=omega, wave_direction=90)


# def test_save_dhyd_to_no_complex():
#     hyd = Hyddb1.create_from(r'files/barge_100_30_4.dhyd')
#     rao = hyd.force_rao(0)
#     r = rao.to_xarray_nocomplex()
#     print(r)


def test_save_dhyd(data_path):
    hyd = Hyddb1.create_from(data_path / "barge_100_30_4.dhyd")

    # find a temporary folder to work in
    import tempfile
    from pathlib import Path

    tempdir = Path(tempfile.gettempdir())

    file = tempdir / "test_write_dhyd_and_read_it_again.dhyd"

    hyd.save_as(file)

    copy = Hyddb1.create_from(file)

    hyd.assert_allclose_to(copy)

    # assert_allclose(hyd._mass.values, copy._mass.values, rtol=1e-3, atol=1)
    # assert_allclose(hyd._damping.values, copy._damping.values, rtol=1e-3, atol=1)
    #
    # for i in range(6):
    #     F0 = hyd._force[i]
    #     Fc = copy._force[i]
    #
    #     assert_allclose(F0._data.amplitude.values, Fc._data.amplitude.values, rtol=1e-3, atol=1)


def test_n_frequencies(data_path):
    hyd = gimme(data_path=data_path)
    assert hyd.n_frequencies == 28


def test_n_frequencies_error():
    """Check that an error is raised if the RAOs have un-equal frequencies when n_frequencies is requested"""
    hyd = Hyddb1()
    hyd._force[1].add_frequency(0.176)
    with pytest.raises(ValueError):
        assert hyd.n_frequencies == 9


def test_n_wave_directions(data_path):
    hyd = gimme(data_path=data_path)
    assert hyd.n_wave_directions == 9


def test_wave_directions(data_path):
    hyd = gimme(data_path=data_path)

    from numpy.testing import assert_almost_equal

    assert_almost_equal(
        hyd.wave_directions, [0.0, 22.5, 45.0, 67.5, 90.0, 112.5, 135.0, 157.5, 180.0]
    )


def test_n_wave_directions_error():
    """Check that an error is raised if the RAOs have un-equal directions when n_wave_directions is requested"""
    hyd = Hyddb1()
    hyd._force[1].add_direction(17)
    with pytest.raises(ValueError):
        assert hyd.n_wave_directions == 9


def test_write_hyd_file(data_path):
    hyd = gimme(data_path=data_path)

    # Create temporary file
    with tempfile.NamedTemporaryFile(suffix=".hyd", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        # Write to the temporary file
        hyd.to_hyd_file(tmp_path, hydrostatics=None)

        # Check if it exists
        assert os.path.exists(tmp_path), "HYD file was not created"
        assert os.path.getsize(tmp_path) > 0, "HYD file is empty"

    finally:
        # Clean up the file
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def test_add_addedmass():
    hyd = Hyddb1()
    hyd.set_amass(omega=1, m6x6=np.ones((6, 6)))

    # by default omega=0 is set to zeros
    # check if interpolation to omega = 0.2 yields 0.2
    actual02 = hyd.amass(0.2)
    expected = 0.2 * np.ones((6, 6))

    assert_almost_equal(actual02, expected)


def test_add_damping():
    hyd = Hyddb1()
    hyd.set_damping(omega=1, m6x6=np.ones((6, 6)))

    actual03 = hyd.damping(0.3)
    expected = 0.3 * np.ones((6, 6))

    assert_almost_equal(actual03, expected)


def test_read_hyd(data_path):
    hyd = Hyddb1.create_from_hyd(data_path / "barge.hyd")

    hyd.plot(do_show=False)


def test_interpolate_amass(data_path):
    hyd = Hyddb1.create_from(data_path / "barge_100_30_4.dhyd")
    omegas = np.linspace(0.1, 4, 100)

    hyd.amass(omegas)


def test_interpolate_damping(data_path):
    hyd = Hyddb1.create_from(data_path / "barge_100_30_4.dhyd")
    omegas = np.linspace(0.1, 4, 100)

    hyd.damping(omegas)


def test_interpolate_damping_single(data_path):
    hyd = Hyddb1.create_from(data_path / "barge_100_30_4.dhyd")
    _omegas = np.linspace(0.1, 4, 100)

    hyd.damping(0.05)


def test_interpolate_forces(data_path):
    hyd = Hyddb1.create_from(data_path / "barge_100_30_4.dhyd")
    _omegas = np.linspace(0.1, 4, 100)

    hyd.force(0.3, wave_direction=90)
    hyd.force(0.301, wave_direction=90)
