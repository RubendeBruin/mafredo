# -*- coding: utf-8 -*-
from importlib.metadata import PackageNotFoundError, version
from mafredo.helpers import FrequencyUnit, Symmetry, MotionMode
from mafredo.hyddb1 import Hyddb1
from mafredo.rao import Rao

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = version(dist_name)
except PackageNotFoundError:
    __version__ = "unknown"
finally:
    del PackageNotFoundError, version


__all__ = ["Hyddb1", "Symmetry", "Rao", "FrequencyUnit", "MotionMode"]
