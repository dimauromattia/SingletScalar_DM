r'''
Tests for the calculation of the relic density
==========================

This script performes a test to verify that the relic density calculation in the script is correct.
'''
import matplotlib.pyplot as plt
import numpy as np
from singletscalar_dm import *


# %%
# The part below performes a check for the calculation of the relic density done with MICROMEGAs
# 

print('TESTING INTERPOLATION WITH MICROMEGAS RESULTS')

table_Micro = np.loadtxt('/Users/mattiadimauro/Dropbox/MadDM/MG5_aMC_v2_9_9/SHP_paper_data_scripts/code/SingletScalar_DM/src/singletscalar_dm/data/test_omega_MicroOMEGAs_interpolation.txt')
mass_Micro = table_Micro[:,0]
lambda_Micro = table_Micro[:,1]
omega_Micro = table_Micro[:,2]

Omega_vec_int = np.zeros(len(mass_Micro))
lambdahs_vec_int = np.zeros(len(mass_Micro))
#print('TESTING INTERPOLATION WITH MICROMEGAs RESULTS')
#print('Mass   lambda     Omega_code  Omega_micromegas ratio_omega ratio_lambda')
print('Checking ', end="", flush=True)
for t in range(len(mass_Micro)):
    if t % 10 ==0 and t!=0:
        print('.', end="", flush=True)
    if lambdahs_vec_int[t]<1e0 and mass_Micro[t]>1e1:
        Omega_vec_int[t] = interpolate_Omega_MicrOMEGAs(mass_Micro[t],lambda_Micro[t])
        lambdahs_vec_int[t] = interpolate_lambda_MicrOMEGAs(mass_Micro[t],omega_Micro[t],False)
        if abs((Omega_vec_int[t]-omega_Micro[t])/omega_Micro[t])>0.05 or abs((lambda_Micro[t]-lambdahs_vec_int[t])/lambda_Micro[t])>0.05:
            print(r'WARNING IMPRECISE CALCULATION FOR')
            print(r'm_S=%.3f and lambdahs=%.3e'%(mass_Micro[t],lambda_Micro[t]))
            print(r'omegah2 real is %.3e and interpolated is %.3e'%(omega_Micro[t],Omega_vec_int[t]))
            print(r'Error is %.3f for omegah2 and %.3f for lambda'%((Omega_vec_int[t]-omega_Micro[t])/omega_Micro[t],(lambdahs_vec_int[t]-lambda_Micro[t])/lambda_Micro[t]))
print('')

# %%
# The part below performes a check for the calculation of the relic density done with DRAKE
# 

r'''
Tests for the calculation of the relic density
==========================

This script performes a test to verify that the relic density calculation in the script is correct.
'''

print('TESTING INTERPOLATION WITH DRAKE RESULTS (QCDA)')

table_Micro = np.loadtxt('/Users/mattiadimauro/Dropbox/MadDM/MG5_aMC_v2_9_9/SHP_paper_data_scripts/code/SingletScalar_DM/src/singletscalar_dm/data/test_DRAKE_omega_fBE_QCDA_testinterp.txt')
mass_Micro = table_Micro[:,0]
lambda_Micro = table_Micro[:,1]
omega_Micro = table_Micro[:,2]

Omega_vec_int = np.zeros(len(mass_Micro))
lambdahs_vec_int = np.zeros(len(mass_Micro))

print('Checking ', end="", flush=True)
for t in range(len(mass_Micro)):
    if t % 10 ==0 and t!=0:
        print('.', end="", flush=True)
    Omega_vec_int[t] = interpolate_Omega(mass_Micro[t],lambda_Micro[t],'QCDA',True)
    lambdahs_vec_int[t] = interpolate_lambda(mass_Micro[t],omega_Micro[t],'QCDA',True)
    if abs((Omega_vec_int[t]-omega_Micro[t])/omega_Micro[t])>0.05 or abs((lambda_Micro[t]-lambdahs_vec_int[t])/lambda_Micro[t])>0.05:
        print(r'WARNING IMPRECISE CALCULATION FOR')
        print(r'm_S=%.3f and lambdahs=%.3e'%(mass_Micro[t],lambda_Micro[t]))
        print(r'omegah2 real is %.3e and interpolated is %.3e'%(omega_Micro[t],Omega_vec_int[t]))
        print(r'omegah2 real is %.3e and interpolated is %.3e'%(lambda_Micro[t],lambdahs_vec_int[t]))
        print(r'Error is %.3f for omegah2 and %.3f for lambda'%((Omega_vec_int[t]-omega_Micro[t])/omega_Micro[t],(lambdahs_vec_int[t]-lambda_Micro[t])/lambda_Micro[t]))
print('')

'''
TESTING INTERPOLATION WITH DRAKE RESULTS (QCDA)
Checking ......WARNING IMPRECISE CALCULATION FOR
m_S=55.789 and lambdahs=1.117e-03
omegah2 real is 1.532e+00 and interpolated is 1.552e+00
omegah2 real is 1.117e-03 and interpolated is 0.000e+00
Error is 0.013 for omegah2 and -1.000 for lambda
.WARNING IMPRECISE CALCULATION FOR
m_S=53.684 and lambdahs=1.651e-02
omegah2 real is 8.318e-02 and interpolated is 8.370e-02
omegah2 real is 1.651e-02 and interpolated is 1.792e-02
Error is 0.006 for omegah2 and 0.085 for lambda
'''

print('TESTING INTERPOLATION WITH DRAKE RESULTS (QCDB)')

table_Micro = np.loadtxt('/Users/mattiadimauro/Dropbox/MadDM/MG5_aMC_v2_9_9/SHP_paper_data_scripts/code/SingletScalar_DM/src/singletscalar_dm/data/test_DRAKE_omega_fBE_QCDB_testinterp.txt')
mass_Micro = table_Micro[:,0]
lambda_Micro = table_Micro[:,1]
omega_Micro = table_Micro[:,2]

Omega_vec_int = np.zeros(len(mass_Micro))
lambdahs_vec_int = np.zeros(len(mass_Micro))
print('Checking ', end="", flush=True)
for t in range(len(mass_Micro)):
    if t % 10 ==0 and t!=0:
        print('.', end="", flush=True)
    Omega_vec_int[t] = interpolate_Omega(mass_Micro[t],lambda_Micro[t],'QCDB',False)
    lambdahs_vec_int[t] = interpolate_lambda(mass_Micro[t],omega_Micro[t],'QCDB',False)
    if abs((Omega_vec_int[t]-omega_Micro[t])/omega_Micro[t])>0.05 or abs((lambda_Micro[t]-lambdahs_vec_int[t])/lambda_Micro[t])>0.05:
        print(r'WARNING IMPRECISE CALCULATION FOR')
        print(r'm_S=%.3f and lambdahs=%.3e'%(mass_Micro[t],lambda_Micro[t]))
        print(r'omegah2 real is %.3e and interpolated is %.3e'%(omega_Micro[t],Omega_vec_int[t]))
        print(r'Error is %.3f for omegah2 and %.3f for lambda'%((Omega_vec_int[t]-omega_Micro[t])/omega_Micro[t],(lambdahs_vec_int[t]-lambda_Micro[t])/lambda_Micro[t]))
print('')