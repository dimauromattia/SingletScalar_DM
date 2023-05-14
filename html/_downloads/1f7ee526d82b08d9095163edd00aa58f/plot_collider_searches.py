'''
Direct detection
================

Calculate various direct detection observables.
'''

import matplotlib.pyplot as plt
import numpy as np
from singletscalar_dm import *

# %%
# The constraints from collider searches are calculated using the upper limits on the contribution of invisible Higgs decay. 
# The package calculates the branching ratio in invibile decay with the function `Gamma_inv`.

gamma_inv = Gamma_inv(50,0.01)
print(gamma_inv)

# %%
# In the following, we show how to compute the upper limits using collider data.

Gamma_inv_measured = 0.170
Gamma_inv_future = 0.025
mh = 125
Gamma_H_SM = 3.8e-3 # GeV

DMmass_vec = np.linspace(10.,mh/2.,1000)
Lambda_current_vec = np.zeros(len(DMmass_vec))
Lambda_future_vec = np.zeros(len(DMmass_vec))
Br_value = Gamma_inv_measured
for t in range(len(DMmass_vec)):
    Lambda_current_vec[t] = minimize_br_inv(DMmass_vec[t],Gamma_H_SM,Gamma_inv_measured)
    Lambda_future_vec[t] = minimize_br_inv(DMmass_vec[t],Gamma_H_SM,Gamma_inv_future)

# %%

table_newB = np.loadtxt(import_data_file('bounds_collider_SHP_current.dat'))
table_newB_p = np.loadtxt(import_data_file('bounds_collider_SHP_projection.dat'))

# %%

fig = plt.figure(figsize=(8,6))
plt.plot(DMmass_vec,Lambda_current_vec,lw=2.0,ls='--',color='blue',label=r'ATLAS+2020, Analytic')
plt.plot(table_newB[:,0],table_newB[:,1],lw=2.5,ls='--',color='black',label=r'ATLAS+2020')
plt.plot(table_newB_p[:,0],table_newB_p[:,4],lw=2.0,ls='-.',color='red',label=r'14 TeV HL-LHC')
plt.plot(table_newB_p[:,0],table_newB_p[:,5],lw=2.0,ls=':',color='orange',label=r'27 TeV HL-LHC')
plt.xlabel(r'$m_{S}$ [GeV]', fontsize=20)
plt.ylabel(r'$\lambda_{HS}$', fontsize=20)
plt.axis([10,100,1e-4,1e1])
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params('both', length=8, width=3, which='major')
plt.tick_params('both', length=6, width=3, which='minor')
plt.grid(True)
plt.yscale('log')
plt.xscale('linear')
plt.legend(loc=4,prop={'size':16},numpoints=1, scatterpoints=1, ncol=2)
fig.tight_layout(pad=0.5)
plt.show()
