from mafredo import hyddb1
import numpy as np
import matplotlib.pyplot as plt


filename = r"c:\data\temp.nc"

hyd = hyddb1.Hyddb1()

print('reading data')

hyd.load_from_capytaine(r"capytaine.nc")
hyd.save_as(filename)

h2 = hyddb1.Hyddb1()
h2.load_from(filename)