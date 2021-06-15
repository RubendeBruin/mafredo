from mafredo import hyddb1
import numpy as np
import matplotlib.pyplot as plt

hyd = hyddb1.Hyddb1()

print('reading data')

# hyd.load_from_capytaine(r"capytaine.nc")
hyd.create_from(r"c:\data\temp.nc")

# Barge is 100m x 30m x 5m draft
Awl = 100*30
Disp = 100*30*5
KB = 2.5
KG = 5
BMt = (1/12)*100*30**3 / Disp
BMl = (1/12)*30*100**3 / Disp
GMt = KB + BMt - KG
GMl = KB + BMl - KG

mass = 1.025 * Disp
rxx = 10
ryy = 40
rzz = 40

# Intertia matrix
I = np.zeros((6,6))
I[0,0] = mass
I[1,1] = mass
I[2,2] = mass
I[3,3] = mass * rxx**2
I[4,4] = mass * ryy**2
I[5,5] = mass * rzz**2

print(I.diagonal())

# Stiffness matrix
K = np.zeros((6,6))
K[0,0] = 100
K[1,1] = 100
K[2,2] = Awl * 1.025*9.81
K[3,3] = GMt * Disp * 1.025 * 9.81
K[4,4] = GMl * Disp * 1.025 * 9.81
K[5,5] = 100

# Calculate the RAO

omegas = np.linspace(0.01,4,100)
# omegas = hyd.frequencies
heading = 90

print('regridding')

hyd.add_frequencies(omegas)

print('calculating rao')
rao = []


for omega in omegas:

    # get the hydrodynamic components
    added_mass = hyd.amass(omega)
    B_hyd = hyd.damping(omega)
    F_wave = hyd.force(omega, heading)

    # construct the system
    # -omega^2 M + omega B + K = F

    # inertia terms
    A = np.zeros((6,6),dtype=complex)

    A += -omega**2 * (I + added_mass)

    # damping terms
    A += 1j * omega * B_hyd

    # stiffness terms
    A += K

    # solve
    x = np.linalg.solve(A, F_wave)
    rao.append(x)

print('plotting')


plt.figure()
for i in range(6):

    a = np.array([r[i] for r in rao])

    ax1 = plt.subplot(3, 2, i + 1)

    amplitude = abs(a)
    if i>2:
        amplitude = np.rad2deg(amplitude)


    ax1.plot(omegas, amplitude, label = "amplitude", color = 'black', linewidth = 1)
    plt.title(hyd._modes[i])
    ax1.set_xlabel('omega [rad/s]')

    yy = plt.ylim()
    if yy[1] < 1e-5:
        plt.ylim((0,0.1))
        continue
    else:
        plt.ylim((0, yy[1]))


    if i==4:
        plt.legend()

    ax2 = ax1.twinx()
    ax2.plot(omegas, np.angle(a),label = "phase", color = 'black', linestyle = ':', linewidth = 1)


    if i==5:
        plt.legend()

plt.suptitle('Incoming wave direction = {}'.format(heading))
plt.tight_layout()

plt.show()










