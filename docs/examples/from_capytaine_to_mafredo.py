from mafredo import hyddb1
import numpy as np
import matplotlib.pyplot as plt


filename = 'cheetah'

hyd = hyddb1.Hyddb1()

print('reading data')

hyd.create_from_capytaine(r"{}.nc".format(filename))

plt.subplots(2,3)
for i in range(6):
    plt.subplot(2,3,i+1)
    rao = hyd.force_rao(i)
    rao['amplitude'].plot()

plt.show()


hyd.save_as("{}.dhyd".format(filename))

