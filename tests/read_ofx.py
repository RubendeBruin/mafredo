import pickle
from mafredo import Hyddb1
import numpy as np
import matplotlib.pyplot as plt
import yaml

filename = r'C:\Data\DAVE\mafredo\tests\files\barge_100x30x4.yml'
vessel_type_name = 'Vessel type1'
iDraught = 0

data = Hyddb1.create_from_orcaflex_yml(filename, vessel_type_name, iDraught)
data.expand360_using_symmetry()
data.plot(adm=False, damp=False, amp=False)

# data2 = Hyddb1.create_from_orcaflex_yml(filename, 'Full_directions', iDraught)
# data2.plot()


#
#
# # Orcaflex exported .yml files contain two "documents", separated by '---' , we want the second one.
# with open(filename,'r') as stream:
#     documents = yaml.safe_load_all(stream)
#     model = list(documents)[1]
#
# vessel_types = model['VesselTypes']
#
# vessel_type = None
# for vt in vessel_types:
#     if vt['Name'] == vessel_type_name:
#         vessel_type = vt
#         break
# if vessel_type is None:
#     vessel_types = [vt['Name'] for vt in vessel_types]
#     raise ValueError(f'No vessel type found with name "{vessel_type_name}", available names are {str(vessel_types)}')
#
# # conventions and conversions
#
# PeriodOrFrequency = vessel_type['WavesReferredToBy']
# def to_rads(values):
#     if PeriodOrFrequency == 'period (s)':
#         return 2*np.pi / np.array(values,dtype=float)
#     elif PeriodOrFrequency == 'frequency (rad/s)':
#         return np.array(values, dtype=float)
#     elif PeriodOrFrequency == 'frequency (Hz)':
#         return 2*np.pi * np.array(values, dtype=float)
#     else:
#         raise ValueError(f'Unknown setting for "WavesReferredToBy" : {PeriodOrFrequency}')
#
# RAOPhaseUnitsConvention = vessel_type['RAOPhaseUnitsConvention']
#
# def to_phase_rad(values):
#     if RAOPhaseUnitsConvention == 'degrees':
#         values = np.deg2rad(values)
#     elif RAOPhaseUnitsConvention == 'radians':
#         pass
#     else:
#         raise ValueError(f'Unknown setting for "RAOPhaseUnitsConvention" : {RAOPhaseUnitsConvention}')
#
#     if vessel_type['RAOPhaseConvention'] == 'leads':
#         values =-values
#
#     if vessel_type['RAOPhaseRelativeToConvention'] == 'crest':
#         pass
#     elif vessel_type['RAOPhaseRelativeToConvention'] == 'trough':
#         values += np.pi
#     elif vessel_type['RAOPhaseRelativeToConvention'] == 'zero down-crossing':
#         values += (np.pi/2)
#     elif vessel_type['RAOPhaseRelativeToConvention'] == 'zero up-crossing':
#         values += 1.5*np.pi
#     else:
#         raise ValueError(f'Unknown setting for "RAOPhaseRelativeToConvention" : {vessel_type["RAOPhaseRelativeToConvention"]}')
#
# d0 = vessel_type['Draughts'][iDraught]
#
# # Added mass and damping
# adb = d0['FrequencyDependentAddedMassAndDamping']
#
# period = [d['AMDPeriodOrFrequency'] for d in adb]
#
# freqs = to_rads(period)
# amass =  [d['AddedMassMatrixX, AddedMassMatrixY, AddedMassMatrixZ, AddedMassMatrixRx, AddedMassMatrixRy, AddedMassMatrixRz'] for d in adb]
# damp =  [d['DampingX, DampingY, DampingZ, DampingRx, DampingRy, DampingRz'] for d in adb]
#
# # Wave forces [Load RAO]
# RAOs = d0['LoadRAOs']['RAOs']
#
# for heading in RAOs:
#     direction = heading['RAODirection']
#     data = heading['RAOPeriodOrFrequency, RAOSurgeAmp, RAOSurgePhase, RAOSwayAmp, RAOSwayPhase, RAOHeaveAmp, RAOHeavePhase, RAORollAmp, RAORollPhase, RAOPitchAmp, RAOPitchPhase, RAOYawAmp, RAOYawPhase']
#     # split data into amplitude and phase
#     np_data = np.array(data, dtype=float)
#
#     period = to_rads(np_data[:,0])
#
#     for iMode in range(6):
#         iAmp = 2*iMode + 1
#         iPhase = 2*iMode + 2
#
#         amp = np_data[:,iAmp]
#         phase = np_data[:,iPhase]
#
#         phase = to_phase_rad(phase)
#
#         if iMode > 2:
#             if vessel_type['RAOResponseUnits'] == 'radians':
#                 amp = np.rad2deg(amp)
#
# # Symmetry – possible values: None; xz plane; yz plane; xz and yz planes; Circular
# if vessel_type['Symmetry'] == 'None':
#     pass
#
