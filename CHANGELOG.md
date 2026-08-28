# Changelog

## Unreleased

Plotting API cleaned up so figures can be embedded in a GUI (breaking):

- `Hyddb1.plot()` now returns the list of figures and no longer shows them;
  the `do_show` argument is gone. Call `matplotlib.pyplot.show()` yourself
  when running interactively.
- `Rao.plot`, `plot_amplitude`, `plot_phase` and `plot_surface` create their
  own figure when no `ax` is given instead of drawing into pyplot's "current
  axes" (which could be a figure owned by someone else), and return the Axes.
- All figures are now created with matplotlib's object-oriented API and are
  unknown to pyplot: embedding them in a GUI has no side-effects and nothing
  accumulates in pyplot's global registry. To display them interactively use
  the new `mafredo.show(figs)`, which adopts them into pyplot;
  `matplotlib.pyplot.show()` alone will not show them.
- Fixed: `Rao.plot_surface` with a single wave-direction drew a line plot and
  then an image on top of it; it now falls back to the line plot only.
- Fixed: `Rao.plot` labeled the y-axis "Amplitude" even when plotting phase,
  and always drew a legend even for a single direction.

## Version 2025.2.0

- Removed NetCDF4 in favour of hypy5, this solves fatal crashes that occured when reading a db twice.

## Version 2025.1.0

- Converted setup.cfg to pyproject.toml
- Fixed all unit tests

## Version 0.1

- Made creation methods static and renamed from `load` to `create`
