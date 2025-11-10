"""
Tests for to_dict() and from_dict() methods
"""

import numpy as np
import json
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from numpy.testing import assert_array_almost_equal

from mafredo.rao import Rao
from mafredo.hyddb1 import Hyddb1
from mafredo.helpers import MotionMode, Symmetry


def _xr_var(d: dict, name: str):
    """Return np.array of a Dataset variable's data from xarray's to_dict() layout."""
    return np.array(d["data_vars"][name]["data"])


def _xr_coord(d: dict, name: str):
    """Return np.array of a Dataset coord's data from xarray's to_dict() layout."""
    return np.array(d["coords"][name]["data"])


def test_rao_to_dict():
    """Test RAO to_dict method"""
    # Create test data
    directions = np.array([0, 90, 180])
    omegas = np.array([0.5, 1.0, 1.5])
    amplitude = np.random.random((3, 3))  # [omega, direction]
    phase = np.random.random((3, 3)) * np.pi

    # Create RAO
    rao = Rao.create_from_data(directions, omegas, amplitude, phase, MotionMode.HEAVE)

    # Convert to dict
    rao_dict = rao.to_dict()

    # This unit has been updated as the unit test did not correspond to the actual structure of xarray to dict

    # Verify structure
    assert (
        "data_vars" in rao_dict
        and "coords" in rao_dict
        and "dims" in rao_dict
        and "attrs" in rao_dict
        and "mode" in rao_dict
    )

    data_vars = rao_dict["data_vars"]
    coords = rao_dict["coords"]
    dims = rao_dict["dims"]
    attrs = rao_dict["attrs"]
    mode = rao_dict["mode"]

    assert "amplitude" in data_vars
    assert "phase" in data_vars
    assert "omega" in coords
    assert "wave_direction" in coords
    assert "omega" in dims
    assert "wave_direction" in dims

    amplitude_dict = data_vars["amplitude"]
    phase_dict = data_vars["phase"]
    omega_dict = coords["omega"]
    wave_direction_dict = coords["wave_direction"]

    assert "data" in amplitude_dict
    assert "data" in phase_dict
    assert "data" in omega_dict
    assert "data" in wave_direction_dict

    amplitude_data = amplitude_dict["data"]
    phase_data = phase_dict["data"]
    omega_data = omega_dict["data"]
    wave_direction_data = wave_direction_dict["data"]

    # Verify data types (should be lists for JSON compatibility)
    assert isinstance(mode, int)
    assert isinstance(amplitude_data, list)
    assert isinstance(phase_data, list)
    assert isinstance(wave_direction_data, list)
    assert isinstance(omega_data, list)
    assert isinstance(attrs, dict)

    # Verify data values
    assert_array_almost_equal(np.array(amplitude_data), amplitude)
    assert_array_almost_equal(np.array(phase_data), phase)
    assert_array_almost_equal(np.array(wave_direction_data), directions)
    assert_array_almost_equal(np.array(omega_data), omegas)
    assert mode == MotionMode.HEAVE.value
    assert not attrs


def test_rao_from_dict():
    """Test RAO from_dict method"""
    # Create test data
    directions = np.array([0, 45, 90, 135, 180])
    omegas = np.array([0.2, 0.5, 1.0, 1.5, 2.0])
    amplitude = np.random.random((5, 5))  # [omega, direction]
    phase = np.random.random((5, 5)) * np.pi

    # Create original RAO
    original_rao = Rao.create_from_data(
        directions, omegas, amplitude, phase, MotionMode.ROLL
    )

    # Convert to dict and back
    rao_dict = original_rao.to_dict()
    restored_rao = Rao.from_dict(rao_dict)

    # Verify data integrity
    assert_array_almost_equal(
        original_rao.wave_directions, restored_rao.wave_directions
    )
    assert_array_almost_equal(original_rao.omega, restored_rao.omega)
    assert_array_almost_equal(
        original_rao["amplitude"].values, restored_rao["amplitude"].values
    )
    assert_array_almost_equal(
        original_rao["phase"].values, restored_rao["phase"].values
    )
    assert original_rao.mode == restored_rao.mode


def test_rao_dict_json_serialization():
    """Test that RAO dict can be serialized to/from JSON"""
    # Create test data
    directions = np.array([0, 90, 180, 270])
    omegas = np.array([0.1, 0.5, 1.0])
    amplitude = np.random.random((3, 4))
    phase = np.random.random((3, 4)) * np.pi

    # Create RAO
    rao = Rao.create_from_data(directions, omegas, amplitude, phase, MotionMode.PITCH)

    # Convert to dict and serialize to JSON
    rao_dict = rao.to_dict()
    json_str = json.dumps(rao_dict)

    # Deserialize from JSON and create RAO
    loaded_dict = json.loads(json_str)
    restored_rao = Rao.from_dict(loaded_dict)

    # Verify data integrity
    assert_array_almost_equal(rao["amplitude"].values, restored_rao["amplitude"].values)
    assert_array_almost_equal(rao["phase"].values, restored_rao["phase"].values)
    assert rao.mode == restored_rao.mode


def test_rao_dict_with_none_mode():
    """Test RAO dict methods with None mode"""
    directions = np.array([0, 180])
    omegas = np.array([1.0])
    amplitude = np.array([[1.0, 2.0]])
    phase = np.array([[0.0, np.pi / 2]])

    # Create RAO with None mode
    rao = Rao.create_from_data(directions, omegas, amplitude, phase, None)

    # Convert to dict and back
    rao_dict = rao.to_dict()
    restored_rao = Rao.from_dict(rao_dict)

    # if mode is none, to_dict will not assign a key mode to the dict, hence it should not exist in the dict
    assert "mode" not in rao_dict

    # however, the Roa class does creates a mode attribute that is None is mode was not present
    assert rao.mode is None
    assert restored_rao.mode is None


def test_hyddb1_to_dict():
    """Test Hyddb1 to_dict method"""
    # Create a simple Hyddb1 object
    hyd = Hyddb1()

    # Convert to dict
    hyd_dict = hyd.to_dict()

    # Verify structure
    assert "mass" in hyd_dict
    assert "damping" in hyd_dict
    assert "force" in hyd_dict
    assert "symmetry" in hyd_dict
    assert "modes" in hyd_dict

    mass_dict = hyd_dict["mass"]

    # Verify mass structure
    assert "data" in mass_dict
    assert "omega" in mass_dict["coords"]
    assert "radiating_dof" in mass_dict["coords"]
    assert "influenced_dof" in mass_dict["coords"]

    # Verify damping structure
    damping_dict = hyd_dict["damping"]
    assert "data" in damping_dict
    assert "omega" in damping_dict["coords"]

    # Verify force structure (should be list of 6 RAO dicts)
    assert isinstance(hyd_dict["force"], list)
    assert len(hyd_dict["force"]) == 6

    # Verify each force RAO dict has proper structure
    for force_dict in hyd_dict["force"]:
        assert (
            "data_vars" in force_dict
            and "coords" in force_dict
            and "dims" in force_dict
            and "attrs" in force_dict
            and "mode" in force_dict
        )

        data_vars = force_dict["data_vars"]
        coords = force_dict["coords"]
        dims = force_dict["dims"]
        attrs = force_dict["attrs"]
        mode = force_dict["mode"]

        assert "amplitude" in data_vars
        assert "phase" in data_vars
        assert "wave_direction" in coords
        assert "omega" in coords
        assert "wave_direction" in dims
        assert "omega" in dims
        assert not attrs
        assert mode is None

        assert dims["omega"] == 2
        assert dims["wave_direction"] == 2


def test_hyddb1_from_dict():
    """Test Hyddb1 from_dict method"""
    import xarray as xr

    # Create test Hyddb1 with realistic data
    hyd = Hyddb1()

    # Set test frequencies
    omegas = np.array([0.5, 1.0, 1.5])

    # Create test matrices
    test_mass = np.random.random((3, 6, 6)) * 1000
    test_damping = np.random.random((3, 6, 6)) * 100

    hyd._mass = xr.DataArray(
        test_mass,
        coords={
            "omega": omegas,
            "radiating_dof": [0, 1, 2, 3, 4, 5],
            "influenced_dof": [0, 1, 2, 3, 4, 5],
        },
        dims=["omega", "radiating_dof", "influenced_dof"],
    )

    hyd._damping = xr.DataArray(
        test_damping,
        coords={
            "omega": omegas,
            "radiating_dof": [0, 1, 2, 3, 4, 5],
            "influenced_dof": [0, 1, 2, 3, 4, 5],
        },
        dims=["omega", "radiating_dof", "influenced_dof"],
    )

    # Create force RAOs
    directions = np.array([0, 90, 180])
    for i, mode in enumerate(
        [
            MotionMode.SURGE,
            MotionMode.SWAY,
            MotionMode.HEAVE,
            MotionMode.ROLL,
            MotionMode.PITCH,
            MotionMode.YAW,
        ]
    ):
        amplitude = np.random.random((3, 3))
        phase = np.random.random((3, 3)) * np.pi
        hyd._force[i] = Rao.create_from_data(directions, omegas, amplitude, phase, mode)

    hyd._symmetry = Symmetry.XZ

    # Convert to dict and back
    hyd_dict = hyd.to_dict()
    restored_hyd = Hyddb1.from_dict(hyd_dict)

    # Verify data integrity
    assert_array_almost_equal(hyd._mass.values, restored_hyd._mass.values)
    assert_array_almost_equal(hyd._damping.values, restored_hyd._damping.values)
    assert hyd._symmetry == restored_hyd._symmetry
    assert len(hyd._force) == len(restored_hyd._force)

    # Check each force RAO
    for i in range(6):
        assert_array_almost_equal(
            hyd._force[i]["amplitude"].values,
            restored_hyd._force[i]["amplitude"].values,
        )
        assert_array_almost_equal(
            hyd._force[i]["phase"].values, restored_hyd._force[i]["phase"].values
        )
        assert hyd._force[i].mode == restored_hyd._force[i].mode


def test_hyddb1_dict_json_serialization():
    """Test that Hyddb1 dict can be serialized to/from JSON"""
    import xarray as xr

    # Create a simple Hyddb1 object
    hyd = Hyddb1()

    # Set some basic data
    omegas = np.array([0.5, 1.0])
    test_mass = np.random.random((2, 6, 6)) * 100
    test_damping = np.random.random((2, 6, 6)) * 10

    hyd._mass = xr.DataArray(
        test_mass,
        coords={
            "omega": omegas,
            "radiating_dof": [0, 1, 2, 3, 4, 5],
            "influenced_dof": [0, 1, 2, 3, 4, 5],
        },
        dims=["omega", "radiating_dof", "influenced_dof"],
    )

    hyd._damping = xr.DataArray(
        test_damping,
        coords={
            "omega": omegas,
            "radiating_dof": [0, 1, 2, 3, 4, 5],
            "influenced_dof": [0, 1, 2, 3, 4, 5],
        },
        dims=["omega", "radiating_dof", "influenced_dof"],
    )

    # Create simple force RAOs
    directions = np.array([0, 180])
    for i, mode in enumerate(
        [
            MotionMode.SURGE,
            MotionMode.SWAY,
            MotionMode.HEAVE,
            MotionMode.ROLL,
            MotionMode.PITCH,
            MotionMode.YAW,
        ]
    ):
        amplitude = np.random.random((2, 2))
        phase = np.random.random((2, 2)) * np.pi
        hyd._force[i] = Rao.create_from_data(directions, omegas, amplitude, phase, mode)

    # Convert to dict and serialize to JSON
    hyd_dict = hyd.to_dict()
    json_str = json.dumps(hyd_dict)

    # Deserialize from JSON and create Hyddb1
    loaded_dict = json.loads(json_str)
    restored_hyd = Hyddb1.from_dict(loaded_dict)

    # Verify data integrity
    assert_array_almost_equal(hyd._mass.values, restored_hyd._mass.values)
    assert_array_almost_equal(hyd._damping.values, restored_hyd._damping.values)
    assert hyd._symmetry == restored_hyd._symmetry
