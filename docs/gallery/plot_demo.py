"""
Read hyd and plot
==================

Minimum example

"""

from pathlib import Path

import matplotlib.pyplot as plt

from mafredo import Hyddb1

file_name = Path(__file__).parent / "barge.hyd"
vessel = Hyddb1.create_from_hyd(filename=file_name)
vessel.plot()

plt.show()
