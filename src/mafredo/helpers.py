import numpy as np
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



