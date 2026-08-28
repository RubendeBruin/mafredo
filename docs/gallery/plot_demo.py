"""
Read hyd and plot
==================

Minimum example

"""

from pathlib import Path

from mafredo import Hyddb1, FrequencyUnit, show

file_name = Path(__file__).parent / "barge.hyd"
vessel = Hyddb1.create_from_hyd(filename=file_name)
show(vessel.plot(unit=FrequencyUnit.seconds, xlim=20))
