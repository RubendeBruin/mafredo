
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from mafredo.rao import Rao
from mafredo.helpers import expand_omega_dim_const,expand_direction_to_full_range, dof_names_to_numbers, MotionMode, Symmetry, MotionModeToStr




class Hyddb1(object):
    """
    Hydrodynamic Database 1st Order
    ================================

    This class contains all information for first order floating bodies. That is:

    - Added mass
    - Damping
    - Wave-forces

    Typically all of the above are obtained from BEM package such as capytaine, wamit, diffrac, orcawave, etc.

    Getting data in


    The following methods exist for loading data from supported formats:

    Load from mafredo's own format (netcdf)

    .. method:: create_from

    Load from capytaine

    .. method:: create_from_capytaine

    Load from MARIN .hyd file

    .. method:: create_from_hyd

    From data

    .. method:: create_from_data(everything)

    .. method:: set_data(everything)
    .. method:: set_amass(omega, m6x6)
    .. method:: set_damping(omega, m6x6)


    Plotting

    .. method:: plot


    Getting data out

    It is possible to obtain the data directly from the _mass , _damping and _force xarrays.
    But my be easier to use one of the following instead:

    .. method:: force
    .. method:: damping
    .. method:: amass


    Saving to file

    To .hyd format, hydrostatics to be supplied optionally

    .. method:: to_hyd_file



    """

    def __init__(self):

        self._force = [Rao() for i in range(6)]

        """
        _mass and _damping are (named) DataArrays with with dimensions 'omega','radiating_dof' and 'influenced_dof'
        
        They are initialized to 0 for omega=0
        
        """
        self._mass = xr.DataArray(np.zeros((1,6,6)),
                                             coords = {"omega" : [0.],
                                                       "radiating_dof" : [0,1,2,3,4,5],
                                                       "influenced_dof" : [0,1,2,3,4,5]},
                                             dims = ['omega', 'radiating_dof', 'influenced_dof'],
                                             )

        self._damping = xr.DataArray(np.zeros((1,6,6)),
                                     coords = {"omega" : [0.],
                                               "radiating_dof" : [0,1,2,3,4,5],
                                               "influenced_dof" : [0,1,2,3,4,5]},
                                     dims = ['omega', 'radiating_dof', 'influenced_dof'],
                                     )

        self._modes = ('Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw')

        self._symmetry = Symmetry.No

        self._kg_to_mt = 1/1000
        self._N_to_kN = 1/1000

    @property
    def n_frequencies(self):
        """Returns the number of frequencies of the mass, damping and force.
        Raises a ValueError if mass, force or damping have unequal number of frequencies."""
        n_freq_rao = self._force[0].n_frequencies
        n_freq_mass = len(self._mass.omega)
        n_freq_damping = len(self._damping.omega)

        if n_freq_mass != n_freq_damping:
            raise ValueError('Mass and damping have different number of frequencies')
        if n_freq_mass != n_freq_rao:
            raise ValueError('Mass and force have different number of frequencies')

        return n_freq_rao

    @property
    def n_wave_directions(self):
        """Returns the number of wave directions in the first of the RAOs"""

        dirs = [rao.n_wave_directions for rao in self._force]
        uds = np.unique(dirs)
        if len(uds) != 1:
            raise ValueError(f'Force RAOs have different amounts of directions: {uds}')

        return int(uds[0])

    @property
    def wave_directions(self):
        return self._force[0].wave_directions

    @property
    def symmetry(self):
        return self._symmetry

    @symmetry.setter
    def symmetry(self, value):
        if isinstance(value, Symmetry):
            self._symmetry = value


    def save_as(self, filename):
        """Saves the contents of the database using the netcdf format.

        See Also:
            load_from
        """
        self._mass.to_netcdf(filename, mode="w", group="mass")
        self._damping.to_netcdf(filename, mode="a", group="damping")

        for i, mode in enumerate(MotionMode):
            self._force[i].to_xarray_nocomplex().to_netcdf(filename, mode="a", group=MotionModeToStr(mode))

    @staticmethod
    def create_from(filename):
        """Loads hydrodynamic data from a netcdf4 file, for example one as saved using save_as.

        See Also:
            save_as
        """
        R = Hyddb1()

        with xr.open_dataarray(filename, group="mass", engine='netcdf4') as ds:
            R._mass = dof_names_to_numbers(ds)
        with xr.open_dataarray(filename, group="damping", engine='netcdf4') as ds:
            R._damping = dof_names_to_numbers(ds)

        R._force = list()
        for mode in MotionMode:
            with xr.open_dataset(filename, group=MotionModeToStr(mode), engine='netcdf4') as ds:
                r = Rao()
                r.from_xarray_nocomplex(ds, mode)
                R._force.append(r)

        return R

    @staticmethod
    def create_from_capytaine(filename):
        """Loads hydrodynamic data from a  dataset produced with capytaine.

        - Wave forces,
        - radiation_damping and
        - added_mass are read.

        See Also:
            Rao.wave_force_from_capytaine
        """
        R = Hyddb1()

        from capytaine.io.xarray import merge_complex_values
        dataset = merge_complex_values(xr.open_dataset(filename))
        # wave_direction = dataset['wave_direction'] * 180 / np.pi  # convert rad/s to deg
        # dataset = dataset.assign_coords(wave_direction=wave_direction)

        R._force.clear()

        for mode in MotionMode:
            r = Rao()
            r.wave_force_from_capytaine(filename, mode)
            r.scale(R._N_to_kN)
            R._force.append(r)

        R._damping = dof_names_to_numbers(dataset['radiation_damping'] * R._N_to_kN)
        R._mass = dof_names_to_numbers(dataset['added_mass'] * R._kg_to_mt)

        return R

    @staticmethod
    def create_from_hyd(filename):
        """Load from the MARIN .hyd format

        .Hyd files are databases containing both hydrodynamics and hydrostatics.
        Both are expressed about the same origin.
        We only read the hydrodynamcis.

        Definition: https://mods.marin.nl/download/attachments/13139976/HYDFILE_Description_Nov_2012.pdf?version=1&modificationDate=1453716460000&api=v2

        Assumptions:
            Headings are repeated in same order for each omega

        Args:
            filename: file to read
        """
        not_parsed = []

        omegas = []
        admas = []
        admasses = []
        bdamp = []
        bdamps = []
        wdirs = []

        famp = []
        famps = []  # iOmega, iHeading, iMode
        feps = []
        fepss = []  # iOmega, iHeading, iMode

        hyd_info = dict()

        R = Hyddb1()

        with open(filename,'r') as f:
            for line in f.readlines():
                line = line.strip('\n') # remove newline characters

                # split into sections with length 10
                keyword = line[:10]

                if 'IDENT ' in keyword:
                    not_parsed.append(line)
                    continue # skip ident lines

                values = [line[10*i+10:10*i+20] for i in range(7)]

                if 'PARA ' in keyword:
                    n_freq = int(values[0])
                    n_head = int(values[1])
                    sym = int(values[2])

                    if sym==0:
                        R.symmetry = Symmetry.No
                    elif sym==1:
                        R.symmetry = Symmetry.XZ
                    elif sym==2:
                        R.symmetry = Symmetry.XZ_and_YZ
                    else:
                        raise ValueError(f'Unknonwn symmetry value {sym}')

                elif 'REFS ' in keyword:
                    hyd_info['water_depth'] = float(values[0])
                    hyd_info['body_draft'] = float(values[1])
                    hyd_info['waterline'] = float(values[2])

                elif 'SPRING ' in keyword:
                    hyd_info['disp_m3'] = float(values[0])
                    hyd_info['Awl_m2'] = float(values[1])
                    hyd_info['COFX_m'] = float(values[2])
                    hyd_info['COBX_m'] = float(values[3])
                    hyd_info['KMT_m'] = float(values[4])
                    hyd_info['KML_m'] = float(values[5])

                elif 'OMEGA ' in keyword:
                    omega = float(values[0])
                    omegas.append(omega)
                elif 'ADMAS ' in keyword:
                    admas.append([float(v) for v in values[:6]])
                    if len(admas) == 6:
                        admasses.append(np.array(admas, dtype=float))
                        admas = []
                elif 'BDAMP ' in keyword:
                    bdamp.append([float(v) for v in values[:6]])
                    if len(bdamp) == 6:
                        bdamps.append(np.array(bdamp, dtype=float))
                        bdamp = []

                elif 'WDIR ' in keyword:
                    wdirs.append(float(values[0]))

                elif 'FAMP ' in keyword:
                    famp.append([float(v) for v in values[:6]])

                    if len(famp) == n_head:
                        famps.append(np.array(famp, dtype=float))
                        famp = []

                elif 'FEPS ' in keyword:
                    feps.append([float(v) for v in values[:6]])

                    if len(feps) == n_head:
                        fepss.append(np.array(feps, dtype=float))
                        feps = []

                else:
                    not_parsed.append(line)

        # now we have
        # - wdirs
        # - omegas
        #
        # admasses and bdamps [ iOmega, iDof, iDof ]
        # famps and fepss     [ iOmega, iHeading, iMode ]
        #

        famps = np.array(famps)
        fepss = np.array(fepss)

        famps = np.swapaxes(famps, 0, 2) # iMode, iHeading, iOmega
        fepss = np.swapaxes(fepss, 0, 2)  # iMode, iHeading, iOmega

        # cut headings axis to only the first unique entries
        wdirs = wdirs[:n_head]



        R.set_data(omega=omegas,
                      added_mass=admasses,
                      damping=bdamps,
                      directions = wdirs,
                      force_amps = famps,
                      force_phase_rad = (np.pi/180)*fepss)

        hyd_info['not parsed'] = not_parsed
        R.hyd_reader_info = hyd_info

        return R

    # def _order_dofs(self, m):
    #     """M can have a single omega, or multiple"""
    #
    #     try:
    #         n_omega = m['omega'].shape[0]
    #     except:
    #         n_omega = 1
    #
    #     if n_omega == 1:
    #         r = np.zeros((6, 6), dtype=float)
    #         for i, m1 in enumerate(self._modes):
    #             for j,m2 in enumerate(self._modes):
    #                 r[i,j] = m.sel(
    #                     influenced_dof=m1, radiating_dof=m2)
    #         return r
    #     else:
    #
    #         r = np.zeros((6,6,n_omega))
    #
    #         for i, m1 in enumerate(self._modes):
    #             for j,m2 in enumerate(self._modes):
    #                 r[i,j,:] = m.sel(
    #                     influenced_dof=m1, radiating_dof=m2)
    #
    #         return r


    def amass(self, omega):
        """Returns the added mass matrix for given frequency or frequencies.
        Linear interpolated is applied if needed"""

        if omega not in self._mass.omega:
            self.add_frequency(omega)

        m = self._mass.sel(omega=omega)
        # r = self._order_dofs(m)

        return m.values

    def _insert_6x6(self, xarr, omega, m6x6):
        """Helper for set_mass and set_damping"""
        if omega in xarr.omega:
            xarr.loc[{'omega':omega}] = m6x6
        else:
            dummy = xarr.isel(omega=0)
            dummy.values = m6x6
            dummy.omega.values = float(omega)
            xarr = xr.concat([xarr, dummy], dim='omega')

        return xarr

    def set_amass(self, omega, m6x6):
        """Sets the added-mass matrix for a given omega"""
        self._mass = self._insert_6x6(self._mass, omega, m6x6)

    def set_damping(self, omega, m6x6):
        """Sets the damping-mass matrix for a given omega"""
        self._damping = self._insert_6x6(self._damping, omega, m6x6)

    @staticmethod
    def create_from_data(self, omega,
                 added_mass,
                 damping,
                 directions,
                 force_amps,
                 force_phase_rad
                 ):
        """Creates a new database using the provided data.

        Args:
            omega : common omega vector [rad/s]
            added_mass : added mass components : [iOmega, iRadating_dof, iInfluenced_dof]
            damping    : damping components : [iOmega, iRadating_dof, iInfluenced_dof]
            directions : wave directions for wave-forces [degrees, coming from]
            force_amps : wave forces [iMode (0..5) , iDirection, iOmega]
            force_phase_rad : wave force phase in rad [iMode (0..5) , iDirection, iOmega]
        """

        r = Hyddb1()
        r.set_data(omega,
                 added_mass,
                 damping,
                 directions,
                 force_amps,
                 force_phase_rad)
        return r


    def set_data(self, omega,
                 added_mass,
                 damping,
                 directions,
                 force_amps,
                 force_phase_rad
                 ):
        """Sets all internal data for added mass, damping and wave-forces.

        Args:
            omega : common omega vector [rad/s]
            added_mass : added mass components : [iOmega, iRadating_dof, iInfluenced_dof]
            damping    : damping components : [iOmega, iRadating_dof, iInfluenced_dof]
            directions : wave directions for wave-forces [degrees, coming from]
            force_amps : wave forces [iMode (0..5) , iDirection, iOmega]
            force_phase_rad : wave force phase in rad [iMode (0..5) , iDirection, iOmega]


        See Also: create_from_data
        """

        self._mass = xr.DataArray(added_mass,
                                             coords = {"omega" : omega,
                                                       "radiating_dof" : [0,1,2,3,4,5],
                                                       "influenced_dof" : [0,1,2,3,4,5]},
                                             dims = ['omega', 'radiating_dof', 'influenced_dof'],
                                             )

        self._damping = xr.DataArray(damping,
                                  coords={"omega": omega,
                                          "radiating_dof": [0,1,2,3,4,5],
                                          "influenced_dof": [0,1,2,3,4,5]},
                                  dims=['omega', 'radiating_dof', 'influenced_dof'],
                                  )

        self._force = []

        for iMode in range(6):
            amps = force_amps[iMode]
            phases = force_phase_rad[iMode]

            rao = Rao()
            rao.set_data(directions=directions,
                         omegas=omega,
                         amplitude=amps,
                         phase=phases)
            self._force.append(rao)





    def damping(self, omega):
        """Returns the damping matrix for given frequency or frequencies.
                Linear interpolated is applied if needed"""

        if omega not in self._damping.omega:
            self.add_frequency(omega)

        m = self._damping.interp(omega=omega)

        return m.values

    def force(self, omega, wave_direction):
        """Returns the force vector for given omega/wave-direction"""

        r = np.zeros(6,dtype=complex)
        for i in range(6):
            r[i] = self._force[i].get_value(omega=omega, wave_direction=wave_direction)

        return r

    def force_rao(self, mode : MotionMode):
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


    def to_hyd_file(self, filename, hydrostatics):
        """Export the database to a .hyd file.

        The exported file will be a single-body, first order .hyd database.

        .Hyd files are databases containing both hydrodynamics and hydrostatics.
        Both are expressed about the same origin.

        Definition: https://mods.marin.nl/download/attachments/13139976/HYDFILE_Description_Nov_2012.pdf?version=1&modificationDate=1453716460000&api=v2

        A .hyd file contains hydrostatic information as well. This hydrostatic information needs to be provided in the "hydrostatic" argument.

         - water_depth  : waterdepth (not hydrostatics)
         - body_draft   : draft of body (m)
         - waterline    : Z-coordinate of waterline wrt hydrodynamic origin of body (m)
         - disp_m3      : displacement of body (m3)
         - Awl_m2       : waterline plane area (m2)
         - COFX_m       : X-coordinate of centre of flotation (m)
         - COBX_m       : X-coordinate of centre of buoyancy (m)
         - KMT_m        : transverse metacentric height KMT measured from keel (m)
         - KML_m        : longitudinal metacentric height KML measured from keel (m)

        Args:
            filename : file to write to
            hydrostatics: A dictionary containing the hydrostatic properties. Optionally define None for all values 0.

        Returns:
            None

        """

        # checks : frequencies for added mass, damping and wave-frequencies need to be the same

        n_omega_waves = self._force


        if hydrostatics is None:
            hydrostatics = dict()
            hydrostatics['water_depth'] = 0
            hydrostatics['body_draft'] = 0
            hydrostatics['waterline'] = 0
            hydrostatics['disp_m3'] = 0
            hydrostatics['Awl_m2'] = 0
            hydrostatics['COFX_m'] = 0
            hydrostatics['COBX_m'] = 0
            hydrostatics['KMT_m'] = 0
            hydrostatics['KML_m'] = 0

        """ The file consists of 80 character records; each record is divided into 8 sections of 10 
            characters. The first section is reserved for a (compulsory) keyword. The remaining 
            seven sections are reserved for (free format) data.
            
            
            
        """
        from mafredo.helpers import f10

        def fixed_format(ident, sections):
            secs = [f10(s, tol=1e-6) for s in sections]
            return '{:10s}'.format(ident) + ''.join(secs) + '\n'

        with open(filename, 'wt') as f:

            # Fist 15 ident lines
            f.write('IDENT     *** .hyd file exported by MaFreDo        ***\n')
            for i in range(14):
                f.write('IDENT     \n')

            #
            f.write(fixed_format('REFS',
                                 [hydrostatics['water_depth'],
                                  hydrostatics['body_draft'],
                                  hydrostatics['waterline'],
                                  0]))
            f.write(fixed_format('SPRING',
                                 [hydrostatics['disp_m3'],
                                  hydrostatics['Awl_m2'],
                                  hydrostatics['COFX_m'],
                                  hydrostatics['COBX_m'],
                                  hydrostatics['KMT_m'],
                                  hydrostatics['KML_m']]))

            if self.symmetry == Symmetry.No or self.symmetry == Symmetry.Circular:
                sym = 0
            elif self.symmetry == Symmetry.XZ:
                sym = 1
            elif self.symmetry == Symmetry.XZ_and_YZ:
                sym = 2
            else:
                raise ValueError('Unsupported symmetry type')

            f.write(fixed_format('PARA',
                                 [self.n_frequencies,
                                  self.n_wave_directions,
                                  sym]))                      # Symmetry by default set to None

            for i_omega in range(self.n_frequencies):
                omega = self._mass.omega.values[i_omega]
                f.write(fixed_format('OMEGA',[omega]))

                ams = self.amass(omega)

                for row in ams:
                    f.write(fixed_format('ADMAS', row))

                # get added mass matrix for this omega
                damp = self.damping(omega)
                for row in damp:
                    f.write(fixed_format('BDAMP', row))

                # wave-forces
                for dir in self.wave_directions:
                    f.write(fixed_format('WDIR',[dir]))
                    force = self.force(omega,dir)

                    # force is complex
                    amp = np.abs(force)
                    phase = np.rad2deg(np.angle(force))

                    f.write(fixed_format('FAMP', amp))
                    f.write(fixed_format('FEPS', phase))

            f.write(fixed_format('PARA2',
                                 [0, 0]))
            f.write('END\n')

    def plot(self, adm=True, damp=True, amp=True, phase=True, do_show=True):
        """Produces a plot of the contents of the database

        Args:
            adm: plot added mass
            damp: plot damping
            amp: plot force amplitudes
            phase: plot force phases
            do_show : do plt.show()

        Returns:
            figure handles


        """

        import matplotlib.pyplot as plt

        figs = []

        # --- RAO amplitudes

        if amp:

            fig, axes = plt.subplots(3,2, figsize=(10,15))
            axes=axes.flatten()
            for i in range(6):
                force = self._force[i]
                force._data['amplitude'].plot(ax=axes[i], cmap=plt.cm.GnBu)
                axes[i].set_title(self._modes[i])
            fig.suptitle('Force RAO amplitudes')

            figs.append(fig)

        # --- RAO phass

        if phase:

            fig, axes = plt.subplots(3,2, figsize=(10,15))
            axes=axes.flatten()
            for i in range(6):
                force = self._force[i]

                force._data['phase'].plot(ax=axes[i], cmap=plt.cm.twilight_shifted)
                axes[i].set_title(self._modes[i])
            fig.suptitle('Force RAO phase [rad]')

            figs.append(fig)

        # Added mass

        if adm:
            fig, axes = plt.subplots(3, 2, figsize=(10, 15))
            axes = axes.flatten()
            for i in range(6):

                mode = self._modes[i]

                for other in range(6):
                    if i==other:
                        lw = 2
                    else:
                        lw = 1
                    self._mass.sel(radiating_dof=i, influenced_dof=other).plot(ax=axes[i], lw=lw, label = self._modes[other])

                axes[i].set_title(self._modes[i])
                if i==5:
                    axes[i].legend()
            fig.suptitle('Added mass\nDiagonal terms shown with thicker line')

            figs.append(fig)

        # Damping

        if damp:
            fig, axes = plt.subplots(3, 2, figsize=(10, 15))
            axes = axes.flatten()
            for i in range(6):

                mode = self._modes[i]

                for other in range(6):
                    if i == other:
                        lw = 2
                    else:
                        lw = 1
                    self._damping.sel(radiating_dof=i, influenced_dof=other).plot(ax=axes[i], lw=lw, label = self._modes[other])

                axes[i].set_title(self._modes[i])
                if i == 5:
                    axes[i].legend()
            fig.suptitle('Damping \nDiagonal terms shown with thicker line')

            figs.append(fig)

        if do_show:
            plt.show()