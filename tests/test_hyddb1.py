from mafredo.hyddb1 import Hyddb1


def test_load_dhyd():
    hyd = Hyddb1()
    hyd.load_from(r'files/barge_100_30_4.dhyd')

def test_load_nc():
    hyd = Hyddb1()
    hyd.load_from_capytaine(r"files/capytaine.nc")

    omega = 0.01

    mass = hyd.amass(omega=omega)
    damping = hyd.damping(omega=omega)
    force = hyd.force(omega=omega, wave_direction=90)
