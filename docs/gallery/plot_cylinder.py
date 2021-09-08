"""
Read hyd and plot
==================

Minimum example

"""

from mafredo import *



# First create the object
cylinder = Hyddb1.create_from_capytaine('open_cylinder.nc')
cylinder_refined = Hyddb1.create_from_capytaine('open_cylinder_2.nc')
cylinder_further_refined = Hyddb1.create_from_capytaine('open_cylinder_3.nc')


# TODO: implement
# merge
cylinder.add(cylinder_refined)
cylinder.add(cylinder_further_refined)

# and finally plot
cylinder.plot()