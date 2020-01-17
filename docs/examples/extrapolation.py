# Example of extrapolation of values

from mafredo import hyddb1, rao
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

hyd = hyddb1.Hyddb1()
hyd.load_from_capytaine(r"capytaine.nc")

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