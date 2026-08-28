from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless backend for CI/pytest
import matplotlib.pyplot as plt
import numpy as np

from mafredo.hyddb1 import Hyddb1


def test_example(data_path: Path) -> None:
    # Load hydrodynamics and extract a RAO
    hyd = Hyddb1.create_from_capytaine(data_path / "capytaine.nc")
    rao = hyd.force_rao(0)

    rao_ampl = rao["amplitude"]
    rao_phase = rao["phase"]

    # Basic plot (do not show; just ensure plotting works headless)
    plt.figure()
    rao_ampl.plot()
    plt.title("input")

    # Add resolution in omega direction
    rao.add_frequency(0)
    rao.add_frequency(4)
    rao.add_frequency([1, 2, 1.7, 0.6])

    # Plot after frequency additions
    plt.figure()
    rao_ampl.plot()
    plt.title("Added more frequencies")

    # Add symmetry
    rao.expand_symmetry_xz()

    plt.figure()
    rao_ampl.plot()
    plt.title("Added symmetry")

    # Regrid directions
    rao.regrid_direction(np.linspace(0, 360, 80))

    plt.figure()
    rao_ampl.plot()
    plt.title("Regridded directions")

    plt.figure()
    rao_phase.plot()
    plt.title("Regridded directions (phase)")
