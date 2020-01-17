import numpy as np
import xarray as xr
from scipy.optimize import fsolve

def wavelength(omega, waterdepth = 0):
    """Returns the wave-length for this frequency [rad/s] and waterdepth.

    In deep water the wave-length is 2*pi*g / omega^2

    In shallow water the following equations need to be solved numerically:

            1.  k = 2*pi / wavelength
            2.  c = sqrt((g/k) * tanh(k*waterdepth))
            3.  wavelength = T * c = (2*pi/omega) * c


    Args:
        omega : wave-frequency in rad/s
        waterdepth : waterdepth in [m], use 0 for infinite waterdepth
    """


    if waterdepth==0:
        return 2*np.pi*9.81 / omega**2


    def error(guess_wavelength):
        k = 2*np.pi / guess_wavelength
        c = np.sqrt((9.81/k) * np.tanh(k*waterdepth))
        wavelength = 2*np.pi*c / omega

        return guess_wavelength - wavelength

    x = fsolve(error, 2*np.pi*9.81 / omega**2)
    return x

def expand_omega_dim_const(dataset, new_omega):
    """Expands the omega axis of dataset to cover the range of new_omega. Extrapolation is done by repeating the nearest values (ie: keep constant)

    Returns:
        expanded dataset
    """

    if min(new_omega) < min(dataset['omega'].values):
        new_minimum = dataset.sel(omega=min(dataset['omega'].values))  # get the lowest entry
        new_minimum['omega'] = min(new_omega)  # change coordinate
        dataset = xr.concat([dataset, new_minimum], dim='omega')  # and concat

    if max(new_omega) > max(dataset['omega'].values):
        new_max = dataset.sel(omega=max(dataset['omega'].values))  # get the highest entry
        new_max['omega'] = max(new_omega)  # change the coordinate
        dataset = xr.concat([dataset, new_max], dim='omega')  # and concat


    return dataset


def expand_direction_to_full_range(dataset):
    """Adds entries at value+360 or values-360 if the wave_direction in the dataset do not span the 0...360 interval

    Returns:
        expanded dataset
    """
    headings = dataset.coords['wave_direction'].values

    if max(headings) - min(headings) < 360:

        if max(headings) < 360:
            head_first = headings[0]
            head360 = dataset.sel(wave_direction=head_first)
            head360.coords['wave_direction'].values += 360
            dataset = xr.concat([dataset, head360], dim='wave_direction')

        if min(headings) > 0:
            head_last = headings[-1]
            head360 = dataset.sel(wave_direction=head_last)
            head360.coords['wave_direction'].values -= 360
            dataset = xr.concat([dataset, head360], dim='wave_direction')

    return dataset

if __name__ == "__main__":

    l_inf = wavelength(0.4)

    depths = np.linspace(1,300,100)
    lamb = []
    for d in depths:
        lamb.append(wavelength(0.4,d))

    import matplotlib.pyplot as plt
    plt.plot(depths, lamb)

    plt.plot([0,300],[l_inf, l_inf])
    plt.show()



