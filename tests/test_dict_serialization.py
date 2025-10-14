import numpy as np
import json
from numpy.testing import assert_allclose
from mafredo import Rao, Hyddb1, MotionMode, Symmetry


def test_rao_dict_serialization():
    """Test RAO to_dict and from_dict methods"""
    # Create test data
    directions = np.array([0, 90, 180])
    omegas = np.array([0.5, 1.0, 1.5])
    amplitude = np.random.random((3, 3))  # [omega, direction]
    phase = np.random.random((3, 3)) * np.pi

    # Create original RAO
    original_rao = Rao.create_from_data(
        directions, omegas, amplitude, phase, MotionMode.HEAVE
    )

    # Convert to dict and JSON
    rao_dict = original_rao.to_dict()
    json_str = json.dumps(rao_dict)

    # Restore from JSON
    loaded_dict = json.loads(json_str)
    restored_rao = Rao.from_dict(loaded_dict)

    # Verify data integrity
    assert_allclose(original_rao.wave_directions, restored_rao.wave_directions)
    assert_allclose(original_rao.omega, restored_rao.omega)
    assert_allclose(original_rao["amplitude"].values, restored_rao["amplitude"].values)
    assert_allclose(original_rao["phase"].values, restored_rao["phase"].values)
    assert original_rao.mode == restored_rao.mode


def test_hyddb1_dict_serialization():
    """Test Hyddb1 to_dict and from_dict methods"""
    # Create a simple Hyddb1 object with test data
    hyd = Hyddb1()

    # Set test frequencies and create test data
    omegas = np.array([0.5, 1.0, 1.5])
    test_mass = np.random.random((3, 6, 6)) * 1000
    test_damping = np.random.random((3, 6, 6)) * 100

    import xarray as xr

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

    # Create force RAOs for each DOF
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

    # Set symmetry
    hyd._symmetry = Symmetry.XZ

    # Convert to dict and JSON
    hyd_dict = hyd.to_dict()
    json_str = json.dumps(hyd_dict)

    # Restore from JSON
    loaded_dict = json.loads(json_str)
    restored_hyd = Hyddb1.from_dict(loaded_dict)

    # Verify data integrity
    assert_allclose(hyd._mass.values, restored_hyd._mass.values)
    assert_allclose(hyd._damping.values, restored_hyd._damping.values)
    assert hyd._symmetry == restored_hyd._symmetry
    assert len(hyd._force) == len(restored_hyd._force)

    # Check each force RAO
    for i in range(6):
        assert_allclose(
            hyd._force[i]["amplitude"].values,
            restored_hyd._force[i]["amplitude"].values,
        )
        assert_allclose(
            hyd._force[i]["phase"].values, restored_hyd._force[i]["phase"].values
        )
        assert hyd._force[i].mode == restored_hyd._force[i].mode


if __name__ == "__main__":
    test_rao_dict_serialization()
    test_hyddb1_dict_serialization()
    print("All tests passed!")
