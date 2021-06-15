"""
Wavelength
===========

Wave-length for a given frequency can be calculated using .. function::helpers.wavelength.

"""

import numpy as np
from mafredo.helpers import wavelength
import matplotlib.pyplot as plt

T = np.linspace(1,20)
omega = 2*np.pi / T

for wd in (5,10,20, 40,100,0):
    wave_lengths = [wavelength(om, waterdepth=wd) for om in omega]
    if wd == 0:
        name = 'deep water'
    else:
        name = f'waterdepth = {wd}m'

    plt.plot(T, wave_lengths, label=name)
plt.legend()
plt.xlabel('Wave period [s]')
plt.ylabel('Wave length [m]')
plt.ylim(0,200)
plt.grid()





