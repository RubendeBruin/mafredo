from mafredo import *
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

def test_example():
    # get a RAO
    hyd = Hyddb1.create_from_capytaine(r"./files/capytaine.nc")
    rao = hyd.force_rao(0)

    plt.figure()
    rao['amplitude'].plot()
    plt.title('input')

    # add resolution in omega direction

    rao.add_frequency(0)
    rao.add_frequency(4)
    rao.add_frequency([1,2,1.7,0.6])

    # add symmetry
    plt.figure()
    rao['amplitude'].plot()
    plt.title('Added more frequencies')


    rao.expand_symmetry_xz()

    plt.figure()
    rao['amplitude'].plot()
    plt.title('Added symmetry')


    # regrid directions
    plt.figure()
    rao.regrid_direction(np.linspace(0,360,80))
    rao['amplitude'].plot()
    plt.title('Regridded directions')

    plt.figure()
    rao['phase'].plot()
    plt.title('Regridded directions')

    # plt.show()

