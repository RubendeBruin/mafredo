import pickle
from mafredo import hyddb1
import numpy as np
import matplotlib.pyplot as plt


filename = './files/OrcaWave_FPSO.pickle'

infile = open(filename,'rb')
vesselType = pickle.load(infile)
infile.close()

filename = r'C:\data\python\rao\docs\examples\FPSO'
hyd = hyddb1.Hyddb1()
print('reading data')
hyd.create_from_capytaine(r"{}.nc".format(filename))

print('test')

d0 = vesselType['Draughts'][0]
adb = d0['FrequencyDependentAddedMassAndDamping']

iMode = 0

period = [d['AMDPeriodOrFrequency'] for d in adb]
axx =  [d['AddedMassMatrixX, AddedMassMatrixY, AddedMassMatrixZ, AddedMassMatrixRx, AddedMassMatrixRy, AddedMassMatrixRz'][0][0] for d in adb]

plt.plot(period, axx)
axx_cap = [hyd.amass(2*np.pi/p)[0][0] for p in period]

plt.plot(period, axx_cap)
plt.title('Added mass for surge, wave-direction 0')
plt.show()

# ------

period = [d['AMDPeriodOrFrequency'] for d in adb]
ayy =  [d['AddedMassMatrixX, AddedMassMatrixY, AddedMassMatrixZ, AddedMassMatrixRx, AddedMassMatrixRy, AddedMassMatrixRz'][1][0] for d in adb]

plt.plot(period, ayy)
ayy_cap = [hyd.amass(2*np.pi/p)[1][0] for p in period]

plt.plot(period, ayy_cap)
plt.title('Added mass for sway, wave-direction [1,0] (should be 0)]')
plt.show()

# ------

period = [d['AMDPeriodOrFrequency'] for d in adb]
ayy =  [float(d['AddedMassMatrixX, AddedMassMatrixY, AddedMassMatrixZ, AddedMassMatrixRx, AddedMassMatrixRy, AddedMassMatrixRz'][1][1]) for d in adb]

plt.plot(period, ayy)

ayy_cap = [hyd.amass(2*np.pi/p)[1][1] for p in period]

plt.plot(period, ayy_cap)
plt.title('Added mass [1,1]')
plt.show()
