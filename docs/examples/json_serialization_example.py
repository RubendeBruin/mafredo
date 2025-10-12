#!/usr/bin/env python3
"""
Example demonstrating JSON/dict serialization of mafredo objects

This example shows how to:
1. Create RAO and Hyddb1 objects
2. Convert them to dictionaries (JSON-serializable format)
3. Save/load to/from JSON files
4. Combine hydrodynamic data with other data in a single JSON file
"""

import json
import numpy as np
from mafredo import Rao, Hyddb1, MotionMode, Symmetry


def create_example_rao():
    """Create an example RAO object"""
    directions = np.array([0, 45, 90, 135, 180])
    omegas = np.linspace(0.1, 2.0, 20)

    # Create some example amplitude and phase data
    amplitude = np.outer(np.exp(-omegas / 2), np.ones(len(directions))) + np.outer(
        np.ones(len(omegas)), np.sin(np.radians(directions)) * 0.3
    )
    phase = np.outer(omegas * 0.5, np.cos(np.radians(directions)))

    return Rao.create_from_data(directions, omegas, amplitude, phase, MotionMode.HEAVE)


def create_example_hyddb1():
    """Create an example Hyddb1 object"""
    import xarray as xr

    hyd = Hyddb1()
    omegas = np.linspace(0.1, 2.0, 10)

    # Create example mass and damping matrices
    mass_data = np.zeros((len(omegas), 6, 6))
    damping_data = np.zeros((len(omegas), 6, 6))

    for i, omega in enumerate(omegas):
        # Simple diagonal matrices with frequency-dependent values
        mass_data[i] = np.diag([1000, 1000, 1000, 50000, 50000, 50000]) * (
            1 + omega * 0.1
        )
        damping_data[i] = np.diag([100, 100, 100, 5000, 5000, 5000]) * omega

    hyd._mass = xr.DataArray(
        mass_data,
        coords={
            "omega": omegas,
            "radiating_dof": [0, 1, 2, 3, 4, 5],
            "influenced_dof": [0, 1, 2, 3, 4, 5],
        },
        dims=["omega", "radiating_dof", "influenced_dof"],
    )

    hyd._damping = xr.DataArray(
        damping_data,
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
        amplitude = np.random.random((len(omegas), len(directions))) * 1000
        phase = np.random.random((len(omegas), len(directions))) * np.pi
        hyd._force[i] = Rao.create_from_data(directions, omegas, amplitude, phase, mode)

    hyd._symmetry = Symmetry.XZ
    return hyd


def example_rao_json():
    """Example: RAO to/from JSON"""
    print("=== RAO JSON Serialization Example ===")

    # Create an example RAO
    rao = create_example_rao()
    print(
        f"Original RAO: {rao.n_frequencies} frequencies, {rao.n_wave_directions} directions"
    )

    # Convert to dictionary
    rao_dict = rao.to_dict()
    print(f"Converted to dict with keys: {list(rao_dict.keys())}")

    # Save to JSON file
    with open("example_rao.json", "w") as f:
        json.dump(rao_dict, f, indent=2)
    print("Saved to example_rao.json")

    # Load from JSON file
    with open("example_rao.json", "r") as f:
        loaded_dict = json.load(f)

    restored_rao = Rao.from_dict(loaded_dict)
    print(
        f"Restored RAO: {restored_rao.n_frequencies} frequencies, {restored_rao.n_wave_directions} directions"
    )
    print(f"Mode: {restored_rao.mode}")

    # Verify they're the same
    np.testing.assert_array_almost_equal(
        rao["amplitude"].values, restored_rao["amplitude"].values
    )
    print("✅ Data integrity verified!")


def example_hyddb1_json():
    """Example: Hyddb1 to/from JSON"""
    print("\n=== Hyddb1 JSON Serialization Example ===")

    # Create an example Hyddb1
    hyd = create_example_hyddb1()
    print(
        f"Original Hyddb1: {hyd.n_frequencies} frequencies, {hyd.n_wave_directions} directions"
    )

    # Convert to dictionary
    hyd_dict = hyd.to_dict()
    print(f"Converted to dict with keys: {list(hyd_dict.keys())}")

    # Save to JSON file
    with open("example_hyddb1.json", "w") as f:
        json.dump(hyd_dict, f, indent=2)
    print("Saved to example_hyddb1.json")

    # Load from JSON file
    with open("example_hyddb1.json", "r") as f:
        loaded_dict = json.load(f)

    restored_hyd = Hyddb1.from_dict(loaded_dict)
    print(
        f"Restored Hyddb1: {restored_hyd.n_frequencies} frequencies, {restored_hyd.n_wave_directions} directions"
    )
    print(f"Symmetry: {restored_hyd._symmetry}")

    # Verify they're the same
    np.testing.assert_array_almost_equal(hyd._mass.values, restored_hyd._mass.values)
    print("✅ Data integrity verified!")


def example_combined_data():
    """Example: Combine hydrodynamic data with other data in a single JSON file"""
    print("\n=== Combined Data Example ===")

    # Create hydrodynamic data
    hyd = create_example_hyddb1()

    # Create a combined data structure that includes hydrodynamic data plus other information
    combined_data = {
        "project_info": {
            "name": "Example Vessel",
            "description": "Demonstration of combined data storage",
            "created_by": "mafredo example",
            "date": "2025-08-04",
        },
        "vessel_properties": {
            "length": 100.0,
            "beam": 30.0,
            "draft": 4.0,
            "displacement": 12000.0,
        },
        "environmental_conditions": {
            "water_depth": 50.0,
            "wave_spectrum": "JONSWAP",
            "significant_wave_height": 2.5,
            "peak_period": 8.0,
        },
        "hydrodynamic_data": hyd.to_dict(),  # Embed the hydrodynamic data
        "analysis_settings": {
            "frequency_range": [0.1, 2.0],
            "directions": [0, 90, 180],
            "symmetry": "XZ",
        },
    }

    # Save combined data to JSON
    with open("combined_vessel_data.json", "w") as f:
        json.dump(combined_data, f, indent=2)
    print("Saved combined data to combined_vessel_data.json")

    # Load and extract hydrodynamic data
    with open("combined_vessel_data.json", "r") as f:
        loaded_combined = json.load(f)

    print(f"Project: {loaded_combined['project_info']['name']}")
    print(f"Vessel length: {loaded_combined['vessel_properties']['length']} m")

    # Extract and restore hydrodynamic data
    restored_hyd = Hyddb1.from_dict(loaded_combined["hydrodynamic_data"])
    print(f"Hydrodynamic data: {restored_hyd.n_frequencies} frequencies")

    print("✅ Combined data storage and retrieval successful!")


if __name__ == "__main__":
    example_rao_json()
    example_hyddb1_json()
    example_combined_data()

    print("\n🎉 All examples completed successfully!")
    print("\nGenerated files:")
    print("- example_rao.json")
    print("- example_hyddb1.json")
    print("- combined_vessel_data.json")
