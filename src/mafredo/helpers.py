import numpy as np
import xarray as xr
from scipy.optimize import fsolve
from enum import Enum

class Symmetry(Enum):
    No = 0
    XZ = 1
    YZ = 2
    XZ_and_YZ = 3
    Circular = 4

class MotionMode(Enum):
    SURGE = 0
    SWAY = 1
    HEAVE = 2
    ROLL = 3
    PITCH = 4
    YAW = 5

def MotionModeToStr(mode: MotionMode) -> str:
    """These are used in the netcdf databases"""

    if mode == MotionMode.SURGE:
        cmode = 'Surge'
    elif mode == MotionMode.SWAY:
        cmode = 'Sway'
    elif mode == MotionMode.HEAVE:
        cmode = 'Heave'
    elif mode == MotionMode.ROLL:
        cmode = 'Roll'
    elif mode == MotionMode.PITCH:
        cmode = 'Pitch'
    elif mode == MotionMode.YAW:
        cmode = 'Yaw'
    else:
        raise ValueError('Unknonwn mode')

    return cmode

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

def dof_names_to_numbers(ds):
    """Converts the names of the DOFS to numbers. This is unfortunately required to make sure
    that the order of the dofs makes sense when retrieving added mass or damping matrices."""

    # check if needed
    if not isinstance(ds.influenced_dof.values[0], int):
        if isinstance(ds.radiating_dof.values[0], int):
            # already done
            return ds

    # re-name and order

    def names_to_ind(dof_names):
        dof_names[dof_names == 'Surge'] = 0
        dof_names[dof_names == 'Sway'] = 1
        dof_names[dof_names == 'Heave'] = 2
        dof_names[dof_names == 'Roll'] = 3
        dof_names[dof_names == 'Pitch'] = 4
        dof_names[dof_names == 'Yaw'] = 5
        return dof_names

    dof_names = names_to_ind(ds.influenced_dof.values)
    ds = ds.assign_coords(influenced_dof=dof_names)

    dof_names = names_to_ind(ds.radiating_dof.values)
    ds = ds.assign_coords(radiating_dof=dof_names)

    ds = ds.sortby('radiating_dof')
    ds = ds.sortby('influenced_dof')

    return ds

# def fix_order_dofs(m):
#     """M can have a single omega, or multiple"""
#
#     modes = ('Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw')
#
#     try:
#         n_omega = m['omega'].shape[0]
#     except:
#         n_omega = 1
#
#     if n_omega == 1:
#         r = np.zeros((6, 6), dtype=float)
#         for i, m1 in enumerate(modes):
#             for j,m2 in enumerate(modes):
#                 r[i,j] = m.sel(
#                     influenced_dof=m1, radiating_dof=m2)
#         return r
#     else:
#
#         r = np.zeros((6,6,n_omega))
#
#         for i, m1 in enumerate(modes):
#             for j,m2 in enumerate(modes):
#                 r[i,j,:] = m.sel(
#                     influenced_dof=m1, radiating_dof=m2)
#
#         return r

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


def f10(number, tol=1e-12):
    """Make a length-10 string representing the given number

    We use {:10g} for numbers (general format) which converts the numbers in the best way:
            '{:10g}'.format(1325123551512511.0)  --> '1.32512e+15'  (11 chars)
            '{:10g}'.format(2.0) --> '         2'
            '{:10g}'.format(2) --> '         2'
            '{:10g}'.format(-432.0) --> '      -432'
            '{:10g}'.format(-1.1641532182693481e-10) --> '-1.16415e-10' (12 chars)

    Numbers with an absolute value below tol are exported as 0

    """

    # scientific notation

    if abs(number) < tol:
        return '         0'

    for formatter in ['{:10g}','{:.6g}','{:.5g}','{:.4g}','{:.3g}','{:.2g}']:
        s = formatter.format(number)

        if len(s) == 10:
            return s
        if len(s) < 10:
            return (10-len(s))*' '+s

    raise ValueError(f'Can not convert number {number} to a string with length 10')

if __name__ == "__main__":
    pass





