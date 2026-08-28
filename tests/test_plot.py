import matplotlib.pyplot as plt

from mafredo import show
from mafredo.hyddb1 import FrequencyUnit, Hyddb1


def test_plot_figures_are_not_registered_with_pyplot(data_path):
    hyd = Hyddb1.create_from(data_path / "grid_t20.dhyd")
    figs = hyd.plot(unit=FrequencyUnit.seconds)

    assert len(figs) == 4
    assert plt.get_fignums() == []  # object-oriented figures: pyplot-free


def test_show_adopts_figures_into_pyplot(data_path):
    hyd = Hyddb1.create_from(data_path / "grid_t20.dhyd")
    figs = hyd.plot(unit=FrequencyUnit.seconds)

    show(figs, block=False)  # Agg backend: registers, does not open windows

    assert len(plt.get_fignums()) == len(figs)
    assert plt.gcf() is figs[-1]

    # adopting twice must not duplicate
    show(figs, block=False)
    assert len(plt.get_fignums()) == len(figs)

    # a single Axes is accepted as well
    ax = hyd._force[0].plot_amplitude()
    show(ax, block=False)
    assert len(plt.get_fignums()) == len(figs) + 1


if __name__ == "__main__":
    h = Hyddb1.create_from(r"files/grid_t20.dhyd")
    show(h.plot(unit=FrequencyUnit.seconds))
