from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")  # headless backend for CI/pytest
import matplotlib.pyplot as plt
import numpy as np

from mafredo.hyddb1 import Hyddb1


def _dataset(obj: Any):
    """Return the underlying xarray.Dataset from Rao-like object."""
    # Prefer a public 'data' attribute; fall back to a private '_data' if needed.
    if hasattr(obj, "data"):
        return obj.data
    if hasattr(obj, "_data"):
        return obj._data  # noqa: SLF001 - test-only access
    raise AttributeError("Object does not expose an xarray Dataset via .data or ._data")


def test_example(data_path: Path) -> None:
    # Load hydrodynamics and extract a RAO
    hyd = Hyddb1.create_from_capytaine(data_path / "capytaine.nc")
    rao = hyd.force_rao(0)

    ds = _dataset(rao)  # xarray.Dataset

    # Basic plot (do not show; just ensure plotting works headless)
    plt.figure()
    ds["amplitude"].plot()
    plt.title("input")

    # Add resolution in omega direction
    rao.add_frequency(0)
    rao.add_frequency(4)
    rao.add_frequency([1, 2, 1.7, 0.6])

    # Plot after frequency additions
    plt.figure()
    _dataset(rao)["amplitude"].plot()
    plt.title("Added more frequencies")

    # Add symmetry
    rao.expand_symmetry_xz()

    plt.figure()
    _dataset(rao)["amplitude"].plot()
    plt.title("Added symmetry")

    # Regrid directions
    rao.regrid_direction(np.linspace(0, 360, 80))

    plt.figure()
    _dataset(rao)["amplitude"].plot()
    plt.title("Regridded directions")

    plt.figure()
    _dataset(rao)["phase"].plot()
    plt.title("Regridded directions (phase)")
