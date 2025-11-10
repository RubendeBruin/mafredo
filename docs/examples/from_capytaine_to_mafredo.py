import matplotlib.pyplot as plt

from mafredo import hyddb1

filename = "cheetah"

hyd = hyddb1.Hyddb1()

print("reading data")

hyd.create_from_capytaine(rf"{filename}.nc")

plt.subplots(2, 3)
for i in range(6):
    plt.subplot(2, 3, i + 1)
    rao = hyd.force_rao(i)
    rao["amplitude"].plot()

plt.show()


hyd.save_as(f"{filename}.dhyd")
