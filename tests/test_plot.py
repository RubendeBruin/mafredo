import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba

from mafredo import show
from mafredo.hyddb1 import FrequencyUnit, Hyddb1


def test_plot_figures_are_not_registered_with_pyplot(data_path):
    hyd = Hyddb1.create_from(data_path / "grid_t20.dhyd")
    figs = hyd.plot(unit=FrequencyUnit.seconds)

    assert len(figs) == 4
    assert plt.get_fignums() == []  # object-oriented figures: pyplot-free


def test_force_amplitudes_have_one_legend_that_fits(data_path):
    """One legend, beside the bottom-right subplot, with every heading in it and nothing cut off."""
    hyd = Hyddb1.create_from(data_path / "grid_t20.dhyd")
    hyd.expand360_using_symmetry()
    fig = hyd.plot(adm=False, damp=False, phase=False, unit=FrequencyUnit.seconds)[0]

    axes = fig.get_axes()
    assert [ax.get_legend() is not None for ax in axes] == [False] * 5 + [True]

    legend = axes[5].get_legend()
    assert len(legend.get_texts()) == len(hyd.wave_directions)

    fig.canvas.draw()
    box = legend.get_window_extent()
    figure = fig.bbox
    assert figure.x0 <= box.x0 and box.x1 <= figure.x1
    assert figure.y0 <= box.y0 and box.y1 <= figure.y1


def test_line_colors_follow_heading(data_path):
    """Colors come from the heading, not from matplotlib's 10-color cycle.

    With the cycle, a database with 24 headings drew 45 and 195 deg in the same color, and 135 and
    285 deg likewise, so a legend entry could not be told from another heading's line.
    """
    hyd = Hyddb1.create_from(data_path / "grid_t20.dhyd")
    hyd.expand360_using_symmetry()
    ax = hyd._force[3].plot_amplitude()

    lines = {float(line.get_label()): line for line in ax.get_lines()}
    headings = sorted(lines)
    assert len(headings) > 10  # more headings than the default color cycle has colors

    def rgba(heading):
        return to_rgba(lines[heading].get_color())

    up_to_180 = [h for h in headings if h <= 180]
    assert len({rgba(h) for h in up_to_180}) == len(up_to_180)  # every heading up to 180 its own color

    for h in headings:
        if h > 180:
            assert lines[h].get_linestyle() == "--"
            assert rgba(h) == rgba(360 - h)  # the mirror image, same color
        else:
            assert lines[h].get_linestyle() == "-"


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
