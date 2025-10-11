**Ma** rine **Fre** quency **Do** main

A set of tools for working with frequency domain data for marine applications.

The purpose of this package is not to provide yet another format for storing hydrodynamic data.

The purpose is to provide an easier way to store, visualize, exchange, compare and modify this data.

This is done by creating classes for the typical data-types:

- Rao (Raos of any kind)
- Hyddb1 (First order hydrodynamic database with added mass, damping and wave-forces)

These classes supply static methods to create them from supported data-types. For example:

>>> my_vessel = Hyddb1.create_from_capytaine(filename = 'titanic.nc')

modification functions

>>> my_vessel.regrid_omega(new_omega)
>>> my_vessel.add_heading(new_heading)

and export/plot function

>>> my_vessel.plot()
>>> my_vessel.save_as_hyd('titanic.hyd')

Inspired by and build to work with:

- capytaine (BEM) [https://github.com/mancellin/capytaine]
- wavespectra ([https://github.com/wavespectra/wavespectra])
- DAVE (General marine modeller) [https://open-ocean.org/DAVE]

Install
========

Any of the following:

- Conda: `conda install mafredo -c conda-forge`

- Mamba: `mamba install mafredo -c conda-forge`

- pip: `pip install mafredo`

Contributions, compliments and complaints
================================================
https://github.com/RubendeBruin/mafredo

Docs
========
https://mafredo.readthedocs.io/en/latest/



