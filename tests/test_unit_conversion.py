import numpy as np
from mafredo import FrequencyUnit


def test_FrequencyUnit():
    omega = 4 * np.pi

    unit = FrequencyUnit.Hz
    unit_name, x = unit.to_unit(omega)

    assert unit_name == "Hz"
    assert x == 2.0

    unit = FrequencyUnit.rad_s
    unit_name, x = unit.to_unit(omega)

    assert unit_name == "rad/s"
    assert x == omega

    unit = FrequencyUnit.seconds
    unit_name, x = unit.to_unit(omega)

    assert unit_name == "s"
    assert x == 0.5
