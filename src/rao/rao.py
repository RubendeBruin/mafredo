"""
RAO is a data-object for dealing with RAO-type data being
amplitude and phase
as function of
heading and frequency

This class is created to provides methods for the correct interpolation of this type of data which means
- the amplitude and phase are interpolated separately (interpolation of complex numbers result in incorrect amplitude)
- continuity in heading is considered (eg: interpolation of heading 355 from heading 345 and 0)

The dimensions of the dataset are:
- omega [rad/s]
- wave_direction [deg]
- amplitude [any]
- phase [radians]

An attribute "mode"

The amplitude and phase can be anything. But typically would be one of the following:
- A motion RAO  (result of frequency domain response calculation)
- A force or moment RAO (result of diffraction analysis)
- A response spectrum with relative phase angles (result of a motion RAO combined with a wave-spectrum)


It is suggested to define the wave_direction as the direction of wave propagation relative to the X-axis. So heading 90 is propagation along the y-axis.
"""

import xarray as xr
import numpy as np

__license__ = "mpl2"


class Rao(object):

    def __init__(self):

        self._data = xr.Dataset()
        self._data.coords['mode'] = 'not_set'

    def wave_force_from_capytaine(self, filename, mode):
        """
        Reads hydrodynamic data from a netCFD file created with capytaine and copies the
        data for the requested mode into the object.

        Args:
            filename: .nc file to read from
            mode: Name of the mode to read ('Heave' 'Pitch' ... 'Sway' 'Yaw')

        Returns:
            None

        Examples:
            test = Rao()
            test.wave_force_from_capytaine(r"capytaine.nc", "Heave")

        """

        from capytaine.io.xarray import merge_complex_values
        dataset = merge_complex_values(xr.open_dataset(filename))

        if 'excitation_force' not in dataset:
            dataset['excitation_force'] = dataset['Froude_Krylov_force'] + dataset['diffraction_force']

        da = dataset['excitation_force'].sel(influenced_dof = mode)

        self._data['amplitude'] = abs(da)
        self._data['complex_unit'] = xr.DataArray(da/abs(da), dims=('omega','wave_direction'))

        self._data.coords['mode'] = mode

    def regrid_omega(self,new_omega):
        new_amp = self._data['amplitude'].interp(omega = new_omega, method='linear')
        new_cu  = self._data['complex_unit'].interp(omega = new_omega, method='linear')

        self._data = xr.Dataset()
        self._data['complex_unit'] = new_cu
        self._data['amplitude'] = new_amp

    def regrid_heading(self, new_headings):

        # repeat the zero heading at the zero + 360

        headings = self._data.coords['wave_direction'].values

        repeat = self._data

        if max(headings) - min(headings) < 360:

            if max(headings) < 360:
                head_first = headings[0]
                head360 = self._data.sel(wave_direction=head_first)
                head360.coords['wave_direction'].values += 360
                repeat = xr.concat([repeat, head360], dim='wave_direction')

            if min(headings) > 0:

                head_last = headings[-1]
                head360 = self._data.sel(wave_direction=head_last)
                head360.coords['wave_direction'].values -= 360
                repeat = xr.concat([head360, repeat], dim='wave_direction')

        new_amp = repeat['amplitude'].interp(wave_direction=new_headings, method='linear')
        new_cu = repeat['complex_unit'].interp(wave_direction=new_headings, method='linear')

        self._data = xr.Dataset()
        self._data['complex_unit'] = new_cu
        self._data['amplitude'] = new_amp

    def add_symmetry_xz(self):
        """Appends equivalent headings considering xz symmetry to the dataset.

        That is:
        The RAO for heading = a is identical to the RAO for heading = -a

        except that for sway, roll and yaw a sign change will be applied (phase shift of pi)

        """

        try:
            mode = self._data.coords['mode']
        except:
            raise ValueError('Mode coordinate not present - Mode coordinate should be set to surge/sway/...yaw ')
        mode = mode.item().upper()

        if mode in ['SWAY','ROLL','YAW']:
            opposite = True
        elif mode in ['SURGE','HEAVE','PITCH']:
            opposite = False
        else:
            raise ValueError('Mode coordinate should be set to surge/sway/...yaw but is set to {}'.format(mode))

        headings = self._data.coords['wave_direction'].values

        for heading in headings:

            heading_copy = np.mod(-heading, 360)
            if heading_copy in headings:
                continue

            sym = self._data.sel(wave_direction=heading)
            sym.coords['wave_direction'].values = heading_copy

            if opposite:
                sym['complex_unit'] = -sym['complex_unit']

            self._data = xr.concat([self._data, sym], dim='wave_direction')

        self._data = self._data.sortby('wave_direction')

    def __getitem__(self, key):
        if key == "phase":
            return xr.DataArray(np.angle(self._data['complex_unit']),dims=('omega','wave_direction'))
        else:
            return self._data[key]

    def __str__(self):
        return str(self._data)


if __name__ == "__main__":

    import matplotlib.pyplot as plt

    test = Rao()
    test.wave_force_from_capytaine(r"C:\data\python\rao\docs\examples\capytaine.nc", "Pitch")

    test.add_symmetry_xz()
    #
    # test.regrid_omega(np.linspace(0,4,100))
    # test.regrid_heading(np.linspace(0,360,5))

    test['amplitude'].plot()

    plt.show()

