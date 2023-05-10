'''
Direct detection
================

Calculate various direct detection observables.
'''

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from importlib.resources import files
from singletscalar_dm import *
_data_dir = files('singletscalar_dm.data')

# %%
# The direct detection part is calculated using the function `SI_noomega`.
# This function performs the calculation of the spin independent cross section analytically assuming that :math:`f_N=0.297`.
# In the following, we sho the calculation of the spin independent cross section for a specific value of DM mass and coupling. 

mS = 70 # GeV
lambda_hs = 0.01
SI_noomega(mS,lambda_hs)

# %%
# Instead, in order to find the upper limits for :math:`\lambda_{HS}`, you can use the function `GetUL_DD_Mine_nomega`.
# The example below calculates the upper limits for :math:`\lambda_{HS}` for the LZ and future DARWIN experiments.

DMmassDD_vec = np.logspace(np.log10(2.),4.,100)
LambdaDD_LZ_HS_vec = np.zeros(len(DMmassDD_vec))
LambdaDD_Darwin_HS_vec = np.zeros(len(DMmassDD_vec))
for t in range(len(DMmassDD_vec)):
    LambdaDD_LZ_HS_vec[t] = GetUL_DD_nomega(DMmassDD_vec[t],'LZ')
    LambdaDD_Darwin_HS_vec[t] = GetUL_DD_nomega(DMmassDD_vec[t],'DARWIN')

# %%

fig = plt.figure(figsize=(8,6))
plt.plot(DMmassDD_vec[18:],LambdaDD_LZ_HS_vec[18:],lw=3.0,ls='--',color='black',label=r'DD, LZ+2022')
plt.fill_between(DMmassDD_vec[18:],(0.3/0.2)*LambdaDD_LZ_HS_vec[18:],(0.3/0.6)*LambdaDD_LZ_HS_vec[18:],alpha=0.5,color='grey')
plt.plot(DMmassDD_vec[18:],1e-10*(0.3/0.2)*LambdaDD_LZ_HS_vec[18:],lw=5.0,alpha=0.5,color='grey',label=r'LZ, $\rho_{\odot}=[0.2,0.6]$ GeV/cm$^3$')
plt.plot(DMmassDD_vec[11:],LambdaDD_Darwin_HS_vec[11:],lw=3.0,ls='--',color='red',label=r'DD, DARWIN proj.')
plt.xlabel(r'$m_S$ [GeV]', fontsize=20)
plt.ylabel(r'$\lambda_{\rm{HS}}$', fontsize=20)
plt.axis([5,1e4,1e-4,1e1])
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params('both', length=8, width=3, which='major')
plt.tick_params('both', length=6, width=3, which='minor')
plt.grid(True)
plt.yscale('log')
plt.xscale('log')
plt.legend(loc=2,prop={'size':16},numpoints=1, scatterpoints=1, ncol=1)
fig.tight_layout(pad=0.5)
plt.show()

# %%
# It is also possible to calculate the direct detection constraints including the rescaling due to the density of dark matter of the :math:`S` particles.
# First, we calculate the :math:`\xi` term which is the relative ratio of :math:`S` particle density over the total dark matter one.

Omegah2_best = 0.120

table = np.loadtxt(_data_dir.joinpath('Omega_MICROMEGAs_zoom_paper.dat'))

Lambda_vec = np.logspace(-5,1,500)
csi_vec = np.zeros(len(MassDD_vec)*len(Lambda_vec))

cont = 0
for t in range(len(MassDD_vec)):
    for u in range(len(Lambda_vec)):
        if table[cont,2]<0.:
            table[cont,2] = 1e6
        csi_vec[cont] = table[cont,2]/Omegah2_best
        cont = cont + 1

print(csi_vec)

# %%
# The we use the function `GetUL_DD` to get the upper limits including also the relic density of :math:`S`.

DMmassDD_vec = np.logspace(np.log10(5),4.,100)
LambdaDD_LZ_HS_vec = np.zeros(len(MassDD_vec))
LambdaDD_Darwin_HS_vec = np.zeros(len(MassDD_vec))
for t in range(len(MassDD_vec)):
    LambdaDD_LZ_HS_vec[t] = GetUL_DD_withomega(MassDD_vec[t],Lambda_vec,MassDD_vec,csi_vec,'LZ')
    LambdaDD_Darwin_HS_vec[t] = GetUL_DD_withomega(MassDD_vec[t],Lambda_vec,MassDD_vec,csi_vec,'Darwin')

# %%

table_RD_FB = np.loadtxt(_data_dir.joinpath('DRAKE_omega_fBE_QCDB_paper.dat'))
mass_RD_FBQCDB = table_RD_FB[:,0]
lambda_RD_FBQCDB = table_RD_FB[:,1]

funcint_RD = interp1d(mass_RD_FBQCDB,lambda_RD_FBQCDB)
for t in range(len(MassDD_vec)):
    if funcint_RD(MassDD_vec[t])>LambdaDD_LZ_HS_vec[t]:
        LambdaDD_LZ_HS_vec[t]=funcint_RD(MassDD_vec[t])
    if funcint_RD(MassDD_vec[t])>LambdaDD_Darwin_HS_vec[t]:
        LambdaDD_Darwin_HS_vec[t]=funcint_RD(MassDD_vec[t])

# %%

fig = plt.figure(figsize=(8,6))
plt.fill_between(mass_RD_FBQCDB,lambda_RD_FBQCDB*1e-6,lambda_RD_FBQCDB*1.05,color='grey')
plt.fill_between(MassDD_vec,funcint_RD(MassDD_vec),LambdaDD_LZ_HS_vec,lw=2.0,alpha=0.5,color='green')
plt.fill_between(MassDD_vec,funcint_RD(MassDD_vec),LambdaDD_Darwin_HS_vec,lw=2.0,alpha=0.2,color='orange')
plt.plot(MassDD_vec,LambdaDD_LZ_HS_vec,lw=3.0,ls='-',color='darkgreen',label=r'LZ+2022')
plt.plot(MassDD_vec,LambdaDD_Darwin_HS_vec,lw=2.0,ls='-.',color='red',label=r'DARWIN PROJ.')
plt.plot(mass_RD_FBQCDB,lambda_RD_FBQCDB,lw=2.0,ls=':',color='blue', label=r'$\Omega_{S}h^2=0.12}$')
plt.text(54.2,0.05,r'$\Omega_{S}h^2<0.12}$',color='blue', fontsize=16)
plt.text(54.2,0.003,r'$\Omega_{S}h^2>0.12}$',color='blue', fontsize=16)
plt.xlabel(r'$m_{\rm{S}}$ [GeV]', fontsize=20)
plt.ylabel(r'$\lambda_{\rm{HS}}$', fontsize=20)
plt.axis([54,64,1e-4,1e1])
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params('both', length=8, width=3, which='major')
plt.tick_params('both', length=6, width=3, which='minor')
plt.grid(True)
plt.yscale('log')
plt.xscale('linear')
plt.legend(loc=2,prop={'size':16},numpoints=1, scatterpoints=1, ncol=1)
fig.tight_layout(pad=0.5)
plt.show()
