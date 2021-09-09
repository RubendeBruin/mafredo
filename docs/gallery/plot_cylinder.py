"""
Merge datasets along frequency
===============================

Example of how to merge refined hydrodynamic databases.

In this example a model was run multiple times with different frequency grids.

RAO data is merged using xarrays 'merge' method while the added-mass and damping are merged
using xarrays 'concat' method along the omega dimension.

"""

from mafredo import *

try:

    # First create the object
    cylinder = Hyddb1.create_from_capytaine('open_cylinder.nc')
    cylinder_refined = Hyddb1.create_from_capytaine('open_cylinder_2.nc')
    cylinder_further_refined = Hyddb1.create_from_capytaine('open_cylinder_3.nc')

    # merge
    cylinder.add(cylinder_refined)
    cylinder.add(cylinder_further_refined)

    # and finally plot
    cylinder.plot()

except:
    pass # read-the-docs