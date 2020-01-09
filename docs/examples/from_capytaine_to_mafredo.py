from mafredo import hyddb1
import numpy as np
import matplotlib.pyplot as plt


filename = r"c:\data\temp.nc"

hyd = hyddb1.Hyddb1()

print('reading data')

hyd.load_from_capytaine(r"cheetah.nc")
hyd.save_as("cheetah.dhyd")

