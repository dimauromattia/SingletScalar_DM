'''
Relic density calculations
==========================

Calculate the relic density of the model.
'''

import matplotlib.pyplot as plt
import numpy as np
from singletscalar_dm import *

table_int = np.loadtxt(import_data_file('SHP_sigmav_table.dat'))
        sigmav = _lambda2sigmav(DMmass,lambda_hs,table_int)

# %%
# Here we create the likelihood profile for dwarf.

table_dwarf = np.loadtxt(import_data_file('likelihood/LogLike_stacked_paper_SHP.txt'))
LogLike_table_dwarf = np.zeros(shape=(len(sigmav_vec),len(mass_vec)))
for t in range(len(mass_vec)):
    for u in range(len(sigmav_vec)):
        LogLike_table_dwarf[u,t] = -2.*table_dwarf[t*len(sigmav_vec)+u,3]
funcint_dwarf = interp2d(mass_vec,sigmav_vec,LogLike_table_dwarf)

LogLikel_table_dwarf = np.zeros(shape=(len(lambdap_vec),len(mass_vec)))
for t in range(len(mass_vec)):
    for u in range(len(lambdap_vec)):
        sigmav_val = func_lambda2sigmav(mass_vec[t],lambdap_vec[u])
        LogLikel_table_dwarf[u,t] = funcint_dwarf(mass_vec[t],sigmav_val)


#####pbar
table_pbar = np.loadtxt(import_data_file('likelihood/LogLike_Manconi2021_pbar_paper.txt'))
LogLike_table_pbar = np.zeros(shape=(len(sigmav_vec),len(mass_vec)))
for t in range(len(mass_vec)):
    for u in range(len(sigmav_vec)):
        LogLike_table_pbar[u,t] = table_pbar[t*len(sigmav_vec)+u,2]
funcint_pbar = interp2d(mass_vec,sigmav_vec,LogLike_table_pbar)

LogLikel_table_pbar = np.zeros(shape=(len(lambdap_vec),len(mass_vec)))
for t in range(len(mass_vec)):
    for u in range(len(lambdap_vec)):
        sigmav_val = func_lambda2sigmav(mass_vec[t],lambdap_vec[u])
        LogLikel_table_pbar[u,t] = funcint_pbar(mass_vec[t],sigmav_val)
        

#####GCE
table_gce = np.loadtxt(import_data_file('likelihood/Chi_table_Cholis_GCE_MED_paper.txt'))
LogLike_table_gce = np.zeros(shape=(len(sigmav_vec),len(mass_vec)))
for t in range(len(mass_vec)):
    for u in range(len(sigmav_vec)):
        LogLike_table_gce[u,t] = -table_gce[t*len(sigmav_vec)+u,2]
funcint_gce = interp2d(mass_vec,sigmav_vec,LogLike_table_gce)

LogLikel_table_gce = np.zeros(shape=(len(lambdap_vec),len(mass_vec)))
for t in range(len(mass_vec)):
    for u in range(len(lambdap_vec)):
        sigmav_val = func_lambda2sigmav(mass_vec[t],lambdap_vec[u])
        LogLikel_table_gce[u,t] = funcint_gce(mass_vec[t],sigmav_val)
        

#####Combined
LogLike_table_combined = LogLike_table_gce+LogLike_table_pbar+LogLike_table_dwarf
LogLikel_table_combined = LogLikel_table_gce+LogLikel_table_pbar+LogLikel_table_dwarf


fig, ax = pl.subplots(figsize=(8,6))

dlin = ( 80+2 )/30.
scale_vec = np.arange( -80,2, dlin )
scale_cb = np.arange( -80,2, dlin*10.)

cf = ax.contourf(mass_vec, lambdap_vec, LogLikel_table_combined-LogLikel_table_combined.max(), 30, levels=list(scale_vec), cmap='seismic')
#ax.set_xlabel(r'$\lambda_{HS}$', fontsize=18)
#ax.set_ylabel(r'$m_{S}$ [GeV]', fontsize=18)
cbar = fig.colorbar(cf, ax = ax)
cbar.ax.set_ylabel(r'$-\Delta\chi^2$ Combined AstroP', fontsize="large")

#plt.colorbar()
pl.ylabel(r'$\lambda_{HS}$', fontsize=18)
pl.xlabel(r'$m_{S}$ [GeV]', fontsize=18)
pl.axis([35.,90,1e-4,1.0])
pl.xticks(fontsize=16)
pl.yticks(fontsize=16)
pl.grid(True)
pl.yscale('log')
pl.xscale('linear') 
pl.legend(loc=1,prop={'size':15},numpoints=1, scatterpoints=1, ncol=2)
fig.tight_layout(pad=0.5)
pl.savefig(folder+"contour_lambda_combined_test_MED_zoom_paper.pdf")



fig, ax = pl.subplots(figsize=(8,6))

dlin = ( 80+2 )/30.
scale_vec = np.arange( -80,2, dlin )
scale_cb = np.arange( -80,2, dlin*10.)

cf = ax.contourf(mass_vec-62.50, lambdap_vec, LogLikel_table_combined-LogLikel_table_combined.max(), 30, levels=list(scale_vec), cmap='seismic')
#ax.set_xlabel(r'$\lambda_{HS}$', fontsize=18)
#ax.set_ylabel(r'$m_{S}$ [GeV]', fontsize=18)
cbar = fig.colorbar(cf, ax = ax)
cbar.ax.set_ylabel(r'$-\Delta\chi^2$ Combined AstroP', fontsize="large")

pl.text(0.11,2e-5,'MED', fontsize=18)

#plt.colorbar()
pl.ylabel(r'$\lambda_{HS}$', fontsize=18)
pl.xlabel(r'$m_{S}-m_h/2$ [GeV]', fontsize=18)
pl.axis([-0.2,0.2,1e-5,0.005])
pl.xticks(fontsize=16)
pl.yticks(fontsize=16)
pl.grid(True)
pl.yscale('log')
pl.xscale('linear') 
pl.legend(loc=2,prop={'size':16},numpoints=1, scatterpoints=1, ncol=1)
fig.tight_layout(pad=0.5)
pl.savefig(folder+"contour_lambda_combined_test_MED_superzoom_paper.pdf")

