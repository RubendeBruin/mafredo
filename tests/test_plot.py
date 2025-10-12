from mafredo.hyddb1 import Hyddb1, FrequencyUnit


def gimme():
    hyd = Hyddb1.create_from(r"files/grid_t20.dhyd")
    return hyd


if __name__ == "__main__":
    h = gimme()
    h.plot(unit=FrequencyUnit.seconds)

    import matplotlib.pyplot as plt

    plt.show()
