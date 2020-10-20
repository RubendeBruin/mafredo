"""
Hydrodynamic Database 1st Order

Units:

 added mass : mT, mT*m
 damping: kN/(m/s) or kN*m/(rad/s)

Getting data:

 .force -> gives force RAOs
 .damping -> gives damping matrix
 .amass -> gives added mass matrix

Matrix definitions:

The Mass matrix is constructed as follows:
 M(i,j) is the force or moment that results in dof_i as result of a unit acceleration of dof_j
  i = influenced dof
  j = exited dof

Loading data:

  .load_from_capytaine --> import from capytaine .nc database (netCDF)
  .




"""
import xarray as xr
import numpy as np
from mafredo.rao import Rao
from mafredo.helpers import expand_omega_dim_const,expand_direction_to_full_range

class Hyddb1(object):

    def __init__(self):

        self._force = []
        self._mass = xr.Dataset()
        self._damping = xr.Dataset()

        self._modes = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']

        self._kg_to_mt = 1/1000
        self._N_to_kN = 1/1000

    def save_as(self, filename):
        self._mass.to_netcdf(filename, mode="w", group="mass")
        self._damping.to_netcdf(filename, mode="a", group="damping")
        for i in range(6):
            self._force[i].to_xarray_nocomplex().to_netcdf(filename, mode="a", group=self._modes[i])

    def load_from(self, filename):
        with xr.open_dataarray(filename, group="mass") as ds:
            self._mass = ds
        with xr.open_dataarray(filename, group="damping") as ds:
            self._damping = ds

        self._force = list()
        for i in range(6):
            with xr.open_dataset(filename, group=self._modes[i]) as ds:
                r = Rao()
                r.from_xarray_nocomplex(ds, self._modes[i])
                self._force.append(r)

    def load_from_capytaine(self, filename):

        from capytaine.io.xarray import merge_complex_values
        dataset = merge_complex_values(xr.open_dataset(filename))
        dataset['wave_direction'] *= 180 / np.pi  # convert rad/s to deg

        self._force.clear()
        for mode in self._modes:
            r = Rao()
            r.wave_force_from_capytaine(filename, mode)
            r.scale(self._N_to_kN)
            self._force.append(r)

        self._damping = dataset['radiation_damping'] * self._N_to_kN
        self._mass = dataset['added_mass'] * self._kg_to_mt

        # apply scaling


    def _order_dofs(self, m):
        """M can have a single omega, or multiple"""

        try:
            n_omega = m['omega'].shape[0]
        except:
            n_omega = 1

        if n_omega == 1:
            r = np.zeros((6, 6), dtype=float)
            for i, m1 in enumerate(self._modes):
                for j,m2 in enumerate(self._modes):
                    r[i,j] = m.sel(
                        influenced_dof=m1, radiating_dof=m2)
            return r
        else:

            r = np.zeros((6,6,n_omega))

            for i, m1 in enumerate(self._modes):
                for j,m2 in enumerate(self._modes):
                    r[i,j,:] = m.sel(
                        influenced_dof=m1, radiating_dof=m2)

            return r


    def amass(self, omega):
        """Returns the added mass matrix for given frequency or frequencies.
        Linear interpolated is applied if needed"""

        m = self._mass.interp(omega=omega)
        r = self._order_dofs(m)

        return r

    def damping(self, omega):
        """Returns the damping matrix for given frequency or frequencies.
                Linear interpolated is applied if needed"""

        m = self._damping.interp(omega=omega)
        r = self._order_dofs(m)

        return r

    def force(self, omega, wave_direction):
        """Returns the force vector for given omega/wave-direction"""

        r = np.zeros(6,dtype=complex)
        for i in range(6):
            r[i] = self._force[i].get_value(omega=omega, wave_direction=wave_direction)

        return r

    def force_rao(self, mode):
        """Return a reference to the internal force rao object

        Args:
            mode : 0...5 for surge...yaw
        """

        return self._force[mode]

    @property
    def frequencies(self):
        return self._mass['omega'].values


    def regrid_omega(self, new_omega):

        exp_damping = expand_omega_dim_const(self._damping, new_omega)
        exp_mass = expand_omega_dim_const(self._mass, new_omega)

        self._damping = exp_damping.interp(omega=new_omega, method= 'linear')
        self._mass = exp_mass.interp(omega=new_omega, method = 'linear')

        for i in range(6):
            self._force[i].regrid_omega(new_omega)

    def regrid_direction(self,new_headings):

        # added mass and damping only depend on omega, not on wave-direction
        for i in range(6):
            self._force[i].regrid_direction(new_headings)

    def add_direction(self, wave_direction):
        # added mass and damping only depend on omega, not on wave-direction
        for i in range(6):
            self._force[i].add_direction(wave_direction)


    def add_frequency(self, omega):
        """Adds a frequency to the database by linear interpolation"""
        self.add_frequencies([omega])

    def add_frequencies(self, omegas):
        """Adds more frequencies to the database by linear interpolation"""

        omegas = np.array([*omegas, *self.frequencies])
        omegas = np.unique(omegas)
        omegas.sort()

        self.regrid_omega(omegas)





if __name__ == "__main__":

    import matplotlib.pyplot as plt

    hyd = Hyddb1()
    hyd.load_from_capytaine(r"files/capytaine.nc")

    disp = 100*30*5
    omega = 0.01

    mass = hyd.amass(omega=omega)
    damping = hyd.damping(omega=omega)
    force = hyd.force(omega=omega, wave_direction=90)


