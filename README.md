# MaFreDo — Marine Frequency Domain

A set of tools for working with frequency domain data for marine applications.

The purpose of this package is **not** to provide yet another format for storing hydrodynamic data.
Instead, it aims to offer an easier way to **store, visualize, exchange, compare, and modify** this data.

---

## Overview

MaFreDo provides Python classes for common hydrodynamic data types:

- **Rao** — Response Amplitude Operators of any kind
- **Hyddb1** — First-order hydrodynamic database with added mass, damping, and wave forces

Each class includes static constructors to easily load from supported data sources. For example:

```python
from mafredo import Hyddb1

my_vessel = Hyddb1.create_from_capytaine(filename="titanic.nc")
```

You can then modify, visualize, or export your data:

```python
my_vessel.regrid_omega(new_omega)
my_vessel.add_heading(new_heading)

my_vessel.plot()
my_vessel.save_as_hyd("titanic.hyd")
```

---

## Built to integrate with

- [Capytaine (BEM)](https://github.com/mancellin/capytaine)
- [Wavespectra](https://github.com/wavespectra/wavespectra)
- [DAVE – General Marine Modeller](https://open-ocean.org/DAVE)

---

## Installation

You can install MaFreDo using **Conda**, **Mamba**, or **pip**:

```bash
# Conda
conda install mafredo -c conda-forge

# Mamba
mamba install mafredo -c conda-forge

# pip
pip install mafredo
```

---

## Documentation

📘 **Docs:** [https://mafredo.readthedocs.io/en/latest/](https://mafredo.readthedocs.io/en/latest/)

---

## Contributing

Contributions, compliments, and complaints are welcome!
👉 [https://github.com/RubendeBruin/mafredo](https://github.com/RubendeBruin/mafredo)
