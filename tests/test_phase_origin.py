import os
import tempfile

import numpy as np
import pytest
from numpy.testing import assert_almost_equal

from mafredo.hyddb1 import Hyddb1


def test_read_old(data_path):
    hyd = Hyddb1.create_from(data_path / "barge_100_30_4.dhyd")  # does not have phase origin data

    assert hyd.phase_origin == (0,0)


def test_copy(data_path):
    hyd = Hyddb1.create_from(data_path / "barge_100_30_4.dhyd")  # does not have phase origin data

    hyd._phase_origin = (3,14)

    h2 = hyd.copy()
    assert h2.phase_origin == (3,14)


def test_save_reload(data_path):
    hyd = Hyddb1.create_from(data_path / "barge_100_30_4.dhyd")  # does not have phase origin data

    hyd._phase_origin = (3, 14)

    # use a temporary file (context)
    with tempfile.TemporaryDirectory() as tmpdir:
        filename = tmpdir + '/temp.dhyd'
        hyd.save_as(filename)

        h2 = Hyddb1.create_from(filename)

        assert h2.phase_origin == (3,14)
