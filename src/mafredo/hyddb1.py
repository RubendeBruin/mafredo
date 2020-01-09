"""
Hydrodynamic Database 1st Order


Units:

 added mass : mT, mT*m
 damping: kN/(m/s) or kN*m/(rad/s)



Matrix definitions:

The Mass matrix is constructed as follows:
 M(i,j) is the force or moment that results in dof_i as result of a unit acceleration of dof_j
  i = influenced dof
  j = exited dof



"""
import xarray as xr
import numpy as np
from mafredo.rao import Rao

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

    @property
    def frequencies(self):
        return self._mass['omega'].values

    def add_frequencies(self, omegas):
        """Adds more frequencies to the database by linear interpolation"""

        omegas = np.array([*omegas, *self.frequencies])
        omegas = np.unique(omegas)
        omegas.sort()

        self._damping = self._damping.interp(omega = omegas)
        self._mass = self._mass.interp(omega=omegas)
        for i in range(6):
            self._force[i].regrid_omega(omegas)




if __name__ == "__main__":

    import matplotlib.pyplot as plt

    hyd = Hyddb1()
    hyd.load_from_capytaine(r"files/capytaine.nc")

    disp = 100*30*5
    omega = 0.01

    mass = hyd.amass(omega=omega)
    damping = hyd.damping(omega=omega)
    force = hyd.force(omega=omega, wave_direction=90)


