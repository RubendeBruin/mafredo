from mafredo.helpers import *
from numpy.testing import assert_almost_equal
import numpy as np

def test_wavelength_deep():
    omega = 0.1
    L = wavelength(omega)

    expected = 2*np.pi*9.81 / omega**2

    assert_almost_equal(L, expected)


def test_wavelength_shallow():
    omega = 0.1
    waterdepth = 30
    L = wavelength(omega, waterdepth=waterdepth)

    # omega^2 == g * k * tank(k*h)
    # k = 2*pi / wavelength

    k = 2*np.pi / L
    omega2 = 9.81 * k * np.tanh(k*waterdepth)

    assert_almost_equal(omega**2, omega2)