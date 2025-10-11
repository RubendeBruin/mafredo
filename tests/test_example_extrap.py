from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless backend voor CI/pytest
import matplotlib.pyplot as plt
import numpy as np

from mafredo.hyddb1 import Hyddb1


def test_example(data_path: Path) -> None:
    # get a RAO (data onder tests/files)
    hyd = Hyddb1.create_from_capytaine(data_path / "capytaine.nc")
    rao = hyd.force_rao(0)

    # basisplot (niet tonen, alleen genereren)
    plt.figure()
    rao["amplitude"].plot()
    plt.title("input")

    # add resolution in omega direction
    rao.add_frequency(0)
    rao.add_frequency(4)
    rao.add_frequency([1, 2, 1.7, 0.6])

    # plot na extra frequenties
    plt.figure()
    rao["amplitude"].plot()
    plt.title("Added more frequencies")

    # add symmetry
    rao.expand_symmetry_xz()

    plt.figure()
    rao["amplitude"].plot()
    plt.title("Added symmetry")

    # regrid directions
    plt.figure()
    rao.regrid_direction(np.linspace(0, 360, 80))
    rao["amplitude"].plot()
    plt.title("Regridded directions")

    plt.figure()
    rao["phase"].plot()
    plt.title("Regridded directions (phase)")

    # geen plt.show(); figures worden automatisch opgeruimd door pytest run
