'''
Indirect detection calculation
==========================

This script calculates the gamma-ray flux for the Galactic center excess.
'''

import matplotlib.pyplot as plt
import numpy as np
from singletscalar_dm import *
from scipy.interpolate import interp1d

# %%
# In order to calculate the relic density it is possible to use the function `interpolate_Omega`, which uses the code DRAKE near the Higgs resonance and MicrOMEGAs elsewhere.
# For examplte, let's compute the relic density, expressed as :math:`\Omega h^2` for the QCDA model.

table_best  = np.loadtxt('/Users/mattiadimauro/Dropbox/FERMI-LAT/GC_PSRDM/ZENODO_FILES/Figures_12_and_14_GCE_Spectra/GCE_BestFitModel_flux_Inner40x40_masked_disk.dat')
table_up  = np.loadtxt('/Users/mattiadimauro/Dropbox/FERMI-LAT/GC_PSRDM/ZENODO_FILES/GCE_band_up.txt')
table_down  = np.loadtxt('/Users/mattiadimauro/Dropbox/FERMI-LAT/GC_PSRDM/ZENODO_FILES/GCE_band_down.txt')
#func_best = interp1d(table_best[:,0],table_best[:,1])
func_min = interp1d(table_down[:,0],table_down[:,1])
func_max = interp1d(table_up[:,0],table_up[:,1])
energy_C_vec = np.logspace(np.log10(0.315),np.log10(34.0),100)
energy_Cerr_vec = table_best[:,0]
#flux_C_best = func_best(energy_C_vec)
flux_C_min = func_min(energy_C_vec)
flux_C_max = func_max(energy_C_vec)
flux_Cerr_min = func_min(energy_Cerr_vec)
flux_Cerr_max = func_max(energy_Cerr_vec)
flux_av = table_best[:,1]

DMmass=62.485
lambda_hs=2e-4
energy_vec = np.logspace(-1.,2.,100)
flux_DM = np.zeros(len(energy_vec))
for i in range(len(energy_vec)):
    flux_DM[i] = flux_DM_prompt(energy_vec[i],DMmass,lambda_hs)

fig = plt.figure(figsize=(8,6))
plt.plot(energy_vec,flux_DM*np.power(energy_vec,2.), lw=5.0, ls='--', color='blue', label=r'SHP model')
plt.fill_between(energy_C_vec,flux_C_min,flux_C_max,alpha=0.5,color="grey")
plt.errorbar(table_best[:,0], table_best[:,1], yerr=[table_best[:,1]-table_best[:,2],table_best[:,3]-table_best[:,1]], fmt='*', color='black',label=r'Cholis+2022')
#plt.errorbar(energy_av, flux_ul, yerr=flux_ul*0.2, xerr=[energy_min,energy_max], fmt='*', color='black', upltims=True)
plt.ylabel(r'$E^2 \frac{dN}{dE}$ [GeV/cm$^2$/s/sr]', fontsize=18)
plt.xlabel(r'$E$ [GeV]', fontsize=18)
plt.axis([0.1,100,8.e-8,3e-6])
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.tick_params('both', length=7, width=2, which='major')
plt.tick_params('both', length=5, width=2, which='minor')
plt.grid(True)
plt.yscale('log')
plt.xscale('log')
plt.legend(loc=2,prop={'size':16},numpoints=1, scatterpoints=1, ncol=1)
fig.tight_layout(pad=0.5)
plt.show()