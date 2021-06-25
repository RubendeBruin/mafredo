"""
Read hyd and plot
==================

Minimum example

"""

from mafredo import *

# First create the object
vessel = Hyddb1.create_from_hyd('barge.hyd')

# and finally plot
vessel.plot()