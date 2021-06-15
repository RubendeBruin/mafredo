"""
Example of extrapolation of values

** This is not the easiest way to do this **

"""



from mafredo import *
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

hyd = Hyddb1.create_from_capytaine(r"capytaine.nc")

rao = hyd.force_rao(0)

omegas = np.linspace(0,10,100)

plt.subplot(211)
rao['amplitude'].plot()

for omega in omegas:
    rao.add_frequency(omega)


rao.add_symmetry_xz()

rao.regrid_direction(np.linspace(0,360,80))

plt.subplot(212)
rao['amplitude'].plot()
plt.show()

rao._data['amplitude'].plot()