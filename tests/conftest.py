from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: tests must never open windows or block on plt.show()

import pytest

ROOT = Path(__file__).resolve().parent.parent  # repo root
TEST_DATA = ROOT / "tests" / "files"


@pytest.fixture(scope="session")
def data_path() -> Path:
    return TEST_DATA


@pytest.fixture(autouse=True)
def _close_figures():
    """Close all pyplot figures after each test.

    mafredo's own figures are object-oriented and never enter pyplot, but
    xarray's .plot() and mafredo.show() do register figures; this keeps the
    registry empty between tests.
    """
    yield
    import matplotlib.pyplot as plt

    plt.close("all")
