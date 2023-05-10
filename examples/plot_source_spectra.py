'''
Source spectra
==============

Calculate the source spectra.
'''

import matplotlib.pyplot as plt
import numpy as np
from singletscalar_dm import *

# %%
# You can use the function `DMspectra_inttable` to calculate the source spectrum.
# The example below calculates the spectrum for gamma rays with smoothing.

mS = 200 # GeV
lambda_hs = 0.01
log10x,dNdlog10x = DMspectra_inttable(mS,lambda_hs,'gammas',True)
print(log10x,dNdlog10x)

# %%
# The following examples, instead, show the source spectra for gamma rays, positrons, antiprotons and neutrinos.

mass_vec = np.array([3.,30.0,62.,100,400.,2000,10000])
colors_vec = np.array(['black','blue','green','brown','orange','purple','pink'])
ls_vec = np.array(['-','--','-.',':','-','--','-.',':','-'])

# %%

fig = plt.figure(figsize=(8,6))

for t in range(len(mass_vec)):
    x_vec,dNdlogx = DMspectra_inttable(mass_vec[t],0.1,'gammas',True)
    plt.plot(x_vec,dNdlogx,lw=2.5,ls=ls_vec[t],color=colors_vec[t],label=r'$m_s=%d$ GeV'%mass_vec[t])

plt.xlabel(r'$\log_{10}(x)$', fontsize=20)
plt.ylabel(r'$dN/d\log_{10}(x)$', fontsize=20)
plt.text(-7.6,15,r'$\gamma$ rays', fontsize=20)
plt.axis([-8,0.1,1e-4,6e1])
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params('both', length=8, width=3, which='major')
plt.tick_params('both', length=6, width=3, which='minor')
plt.grid(True)
plt.yscale('log')
plt.xscale('linear')
plt.legend(loc=4,prop={'size':16},numpoints=1, scatterpoints=1, ncol=2)
fig.tight_layout(pad=0.5)
plt.savefig('plots/DMspectra_gamma.pdf')

# %%

fig = plt.figure(figsize=(8,6))

for t in range(len(mass_vec)):
    x_vec,dNdlogx = DMspectra_inttable(mass_vec[t],0.1,'antiprotons',True)
    plt.plot(x_vec,dNdlogx,lw=2.5,ls=ls_vec[t],color=colors_vec[t],label=r'$m_s=%d$ GeV'%mass_vec[t])

plt.xlabel(r'$\log_{10}(x)$', fontsize=20)
plt.ylabel(r'$dN/d\log_{10}(x)$', fontsize=20)
plt.text(-1.9,3e-4,r'$\bar{p}$', fontsize=20)
plt.axis([-6,0.1,1e-4,5e0])
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params('both', length=8, width=3, which='major')
plt.tick_params('both', length=6, width=3, which='minor')
plt.grid(True)
plt.yscale('log')
plt.xscale('linear')
plt.legend(loc=2,prop={'size':16},numpoints=1, scatterpoints=1, ncol=1)
fig.tight_layout(pad=0.5)
plt.savefig('plots/DMspectra_antiprotons.pdf')

# %%

fig = plt.figure(figsize=(8,6))

for t in range(len(mass_vec)):
    x_vec,dNdlogx = DMspectra_inttable(mass_vec[t],0.1,'positrons',True)
    plt.plot(x_vec,dNdlogx,lw=2.5,ls=ls_vec[t],color=colors_vec[t],label=r'$m_s=%d$ GeV'%mass_vec[t])

plt.xlabel(r'$\log_{10}(x)$', fontsize=20)
plt.ylabel(r'$dN/d\log_{10}(x)$', fontsize=20)
plt.text(-7.6,3,r'$\gamma$ rays', fontsize=20)
plt.axis([-8,0.1,1e-4,2e1])
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params('both', length=8, width=3, which='major')
plt.tick_params('both', length=6, width=3, which='minor')
plt.grid(True)
plt.yscale('log')
plt.xscale('linear')
plt.legend(loc=4,prop={'size':16},numpoints=1, scatterpoints=1, ncol=1)
fig.tight_layout(pad=0.5)
plt.show()

# %%

fig = plt.figure(figsize=(8,6))

for t in range(len(mass_vec)):
    x_vec,dNdlogx_ne = DMspectra_inttable(mass_vec[t],0.1,'neutrinos_e',True)
    x_vec,dNdlogx_nmu = DMspectra_inttable(mass_vec[t],0.1,'neutrinos_mu',True)
    x_vec,dNdlogx_ntau = DMspectra_inttable(mass_vec[t],0.1,'neutrinos_tau',True)
    plt.plot(x_vec,dNdlogx_ne+dNdlogx_nmu+dNdlogx_ntau,lw=2.5,ls=ls_vec[t],color=colors_vec[t],label=r'$m_s=%d$'%mass_vec[t])

plt.xlabel(r'$\log_{10}(x)$', fontsize=20)
plt.ylabel(r'$dN/d\log_{10}(x)$', fontsize=20)
plt.text(-7.6,15,r'$\gamma$ rays', fontsize=20)
plt.axis([-8,0.1,1e-3,1e2])
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params('both', length=8, width=3, which='major')
plt.tick_params('both', length=6, width=3, which='minor')
plt.grid(True)
plt.yscale('log')
plt.xscale('linear')
plt.legend(loc=4,prop={'size':16},numpoints=1, scatterpoints=1, ncol=1)
fig.tight_layout(pad=0.5)
plt.show()
