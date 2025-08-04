"""
Tests for to_dict() and from_dict() methods
"""

import numpy as np
import json
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from numpy.testing import assert_allclose, assert_array_almost_equal

from mafredo.rao import Rao
from mafredo.hyddb1 import Hyddb1
from mafredo.helpers import MotionMode, Symmetry


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
    
    # Verify structure
    assert 'amplitude' in rao_dict
    assert 'phase' in rao_dict
    assert 'wave_direction' in rao_dict
    assert 'omega' in rao_dict
    assert 'mode' in rao_dict
    
    # Verify data types (should be lists for JSON compatibility)
    assert isinstance(rao_dict['amplitude'], list)
    assert isinstance(rao_dict['phase'], list)
    assert isinstance(rao_dict['wave_direction'], list)
    assert isinstance(rao_dict['omega'], list)
    assert isinstance(rao_dict['mode'], int)
    
    # Verify data values
    assert_array_almost_equal(np.array(rao_dict['amplitude']), amplitude)
    assert_array_almost_equal(np.array(rao_dict['phase']), phase)
    assert_array_almost_equal(np.array(rao_dict['wave_direction']), directions)
    assert_array_almost_equal(np.array(rao_dict['omega']), omegas)
    assert rao_dict['mode'] == MotionMode.HEAVE.value


def test_rao_from_dict():
    """Test RAO from_dict method"""
    # Create test data
    directions = np.array([0, 45, 90, 135, 180])
    omegas = np.array([0.2, 0.5, 1.0, 1.5, 2.0])
    amplitude = np.random.random((5, 5))  # [omega, direction]
    phase = np.random.random((5, 5)) * np.pi
    
    # Create original RAO
    original_rao = Rao.create_from_data(directions, omegas, amplitude, phase, MotionMode.ROLL)
    
    # Convert to dict and back
    rao_dict = original_rao.to_dict()
    restored_rao = Rao.from_dict(rao_dict)
    
    # Verify data integrity
    assert_array_almost_equal(original_rao.wave_directions, restored_rao.wave_directions)
    assert_array_almost_equal(original_rao.omega, restored_rao.omega)
    assert_array_almost_equal(original_rao['amplitude'].values, restored_rao['amplitude'].values)
    assert_array_almost_equal(original_rao['phase'].values, restored_rao['phase'].values)
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
    assert_array_almost_equal(rao['amplitude'].values, restored_rao['amplitude'].values)
    assert_array_almost_equal(rao['phase'].values, restored_rao['phase'].values)
    assert rao.mode == restored_rao.mode


def test_rao_dict_with_none_mode():
    """Test RAO dict methods with None mode"""
    directions = np.array([0, 180])
    omegas = np.array([1.0])
    amplitude = np.array([[1.0, 2.0]])
    phase = np.array([[0.0, np.pi/2]])
    
    # Create RAO with None mode
    rao = Rao.create_from_data(directions, omegas, amplitude, phase, None)
    
    # Convert to dict and back
    rao_dict = rao.to_dict()
    restored_rao = Rao.from_dict(rao_dict)
    
    assert rao_dict['mode'] is None
    assert restored_rao.mode is None


def test_hyddb1_to_dict():
    """Test Hyddb1 to_dict method"""
    # Create a simple Hyddb1 object
    hyd = Hyddb1()
    
    # Convert to dict
    hyd_dict = hyd.to_dict()
    
    # Verify structure
    assert 'mass' in hyd_dict
    assert 'damping' in hyd_dict
    assert 'force' in hyd_dict
    assert 'symmetry' in hyd_dict
    assert 'modes' in hyd_dict
    
    # Verify mass structure
    assert 'data' in hyd_dict['mass']
    assert 'omega' in hyd_dict['mass']
    assert 'radiating_dof' in hyd_dict['mass']
    assert 'influenced_dof' in hyd_dict['mass']
    
    # Verify damping structure
    assert 'data' in hyd_dict['damping']
    assert 'omega' in hyd_dict['damping']
    
    # Verify force structure (should be list of 6 RAO dicts)
    assert isinstance(hyd_dict['force'], list)
    assert len(hyd_dict['force']) == 6
    
    # Verify each force RAO dict has proper structure
    for force_dict in hyd_dict['force']:
        assert 'amplitude' in force_dict
        assert 'phase' in force_dict
        assert 'wave_direction' in force_dict
        assert 'omega' in force_dict
        assert 'mode' in force_dict


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
    for i, mode in enumerate([MotionMode.SURGE, MotionMode.SWAY, MotionMode.HEAVE,
                              MotionMode.ROLL, MotionMode.PITCH, MotionMode.YAW]):
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
        assert_array_almost_equal(hyd._force[i]['amplitude'].values, 
                                  restored_hyd._force[i]['amplitude'].values)
        assert_array_almost_equal(hyd._force[i]['phase'].values, 
                                  restored_hyd._force[i]['phase'].values)
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
    for i, mode in enumerate([MotionMode.SURGE, MotionMode.SWAY, MotionMode.HEAVE,
                              MotionMode.ROLL, MotionMode.PITCH, MotionMode.YAW]):
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
