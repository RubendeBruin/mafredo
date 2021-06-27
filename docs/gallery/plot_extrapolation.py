"""
Interpolation and symmetry
===========================

Example of extrapolation of values

"""

try:

    from mafredo import *
    import xarray as xr
    import numpy as np
    import matplotlib.pyplot as plt

    # get a RAO
    hyd = Hyddb1.create_from_capytaine(r"capytaine.nc")
    rao = hyd.force_rao(0)

    plt.figure()
    rao['amplitude'].plot()
    plt.title('input')

    # add resolution in omega direction

    omegas = np.linspace(0,4,100)
    for omega in omegas:
        rao.add_frequency(omega)

    # add symmetry
    plt.figure()
    rao['amplitude'].plot()
    plt.title('Re-grid frequencies')

    rao.expand_symmetry_xz()

    plt.figure()
    rao['amplitude'].plot()
    plt.title('Added symmetry')

    # regrid directions
    plt.figure()
    rao.regrid_direction(np.linspace(0,360,80))
    rao['amplitude'].plot()
    plt.title('Regridded directions')

    import matplotlib.pyplot as plt
    # plt.show()
except:
    pass # read-the-docs

