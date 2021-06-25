Conventions
============

Added mass and damping
-----------------------
Added-mass and damping are 6x6 matrices and are a function of frequency.

 - added mass : mT, mT*m
 - damping: kN/(m/s) or kN*m/(rad/s)
 - frequency : rad/s

The matrices are constructed as follows:
 M(i,j) is the force or moment that results in dof_i as result of a unit acceleration of dof_j
  i = influenced dof
  j = exited dof

Wave-Forces
------------
Wave-Forces are a function of heading, mode and frequency. Wave-forces are stored as amplitude and a phase. The phase
is stored as a complex number with unit length "complex_unit". The reason for this complexity is that is makes it easy to interpolate.

 - wave-direction : degrees going to
 - phase-angles : radians
 - frequency : rad/s


Frequency-Grid
---------------

The frequency grid for added-mass, damping and wave-forces may be the same, but this is not enforced.

 - frequency : rad/s

Definitions:
------------

 **Headings**
 Wave headings are in ```[degrees]``` counter-clockwise from positive X-axis; ```going to```.
 For a ship with stern-to-bow aligned with the x-axis this means:
 - 0 deg = stern waves,
 - 180 deg = bow waves,
 - 90 deg = waves from SB 
   
**Phase angles**

Phase angles for wave-forces are in ```[radians]``` and are ```Lagging relative to wave crest```.

pseudo code:
```python
        A = -(omega ** 2) * M_total  \
                + 1j * omega * B  \
                + K )              # mass, damping, stiffness
        re = np.real(Fwave)
        im = np.imag(Fwave)
        excitation = re - 1j * im

        rao = np.linalg.solve(A, excitation )

        # calculate in time-domain:

        time_phasor = np.exp(1j * t)
        motions = np.outer(time_phasor, rao)
        response = np.real(motions)
        wave_elevation = wave_amplitude * np.real(time_phasor)
```