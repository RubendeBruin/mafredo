import numpy as np
import pytest
from numpy.testing import assert_allclose


from mafredo.hyddb1 import Hyddb1



def test_write_hyd_and_read_it_again():
    hyd = Hyddb1.create_from(r'files/barge_100_30_4.dhyd')

    # find a temporary folder to work in
    import tempfile
    from pathlib import Path
    tempdir = Path(tempfile.gettempdir())

    file = tempdir / 'test_write_hyd_and_read_it_again.hyd'

    hyd.to_hyd_file(file, hydrostatics=None)

    copy = Hyddb1.create_from_hyd(file)

    assert_allclose(hyd._mass.values, copy._mass.values, rtol=1e-3, atol = 1 )
    assert_allclose(hyd._damping.values, copy._damping.values, rtol=1e-3, atol = 1 )

    for i in range(6):
        F0 = hyd._force[i]
        Fc = copy._force[i]

        assert_allclose(F0._data.amplitude.values, Fc._data.amplitude.values.transpose(), rtol=1e-3, atol=1)


    #
    # hyd.plot(phase=False)
    # copy.plot(phase=False)