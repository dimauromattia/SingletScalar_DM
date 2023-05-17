r'''
Tests for the calculation of annihilation cross section
==========================

This script performes a test to verify that the annihilation cross section calculation in the script is correct.
'''
import matplotlib.pyplot as plt
import numpy as np
from singletscalar_dm import *


# %%
# The part below performes a check for the calculation of the annihilation cross section
# 

def returnsigmav(string):
    '''
    This function returns the sigmav value in the files MadDM_results.txt written in the output folder by MaDDM.
    '''
    sep = ","
    x = 'Total_xsec                    = ['
    f = open(string)
    sigmav_value=0
    for line in f:
      if x in line:
        matches = [l[l.find(x)+len(x):] for l in line.split(sep) if l[l.find(x)+len(x):]]
        sigmav_value = float(matches[0])
    return sigmav_value

makeplots = True
resolution = 0.05

####TEST SPECTRA
massestest_vec = np.array([7.32130641e+00, 2.68320630e+00, 6.50879081e+01, 2.79608900e+01,
       2.21004329e+00, 3.49476587e+01, 2.57498330e+00, 7.81947595e+01,
       8.20403008e+01, 1.53773661e+01, 4.22060542e+01, 3.84737274e+01,
       1.28310945e+01, 4.40556075e+00, 7.92990603e+00, 5.00431664e+01,
       6.91445083e+01, 3.68017099e+01, 2.72113158e+01, 8.23640988e+01,
       1.96941678e+01, 7.72496273e+01, 5.08068593e+01, 4.04675780e+01,
       1.52268685e+01, 6.48589445e+00, 1.71542928e+01, 2.06535496e+00,
       6.39330134e+00, 6.10213982e+00, 7.04721789e+01, 1.45082919e+01,
       9.56094277e+00, 4.64504767e+00, 3.22418969e+00, 8.84617858e+00,
       3.48176161e+00, 7.32271282e+00, 5.82477619e+01, 3.94873004e+00,
       5.43920019e+01, 3.13299477e+00, 1.51600838e+01, 7.41409970e+00,
       3.26848160e+00, 1.65308136e+01, 6.73403493e+00, 4.37227613e+00,
       4.73748809e+01, 1.52631205e+01, 4.34944985e+00, 3.85736311e+01,
       1.05633210e+01, 2.62138163e+00, 6.87298741e+00, 7.86563109e+00,
       9.37946502e+01, 6.52585996e+01, 4.19632310e+01, 8.50707166e+00,
       9.29973129e+00, 4.91579754e+00, 2.93697343e+01, 2.47273574e+01,
       2.09456463e+01, 4.18853682e+00, 3.00787528e+01, 8.50848229e+01,
       6.23985285e+00, 5.87503995e+00, 1.02283565e+01, 4.32994792e+00,
       4.11662886e+00, 2.74497912e+01, 7.72678436e+00, 3.69009771e+00,
       5.31741302e+00, 6.84304475e+01, 1.31330527e+01, 2.18302833e+01,
       4.22357787e+00, 1.12941165e+01, 4.07966129e+01, 4.68489009e+01,
       5.17724494e+01, 2.81206541e+00, 3.64798938e+01, 3.63788060e+00,
       2.37024728e+01, 6.09269107e+00, 2.65955548e+00, 2.20677863e+01,
       4.51690703e+01, 3.60098921e+01, 9.08743988e+01, 4.95564929e+00,
       3.11275198e+00, 8.75004407e+00, 1.93922098e+01, 3.97516804e+01,
       6.21337141e+01, 6.24162110e+01, 6.21746762e+01, 6.24219219e+01,
       6.23194073e+01, 6.29543202e+01, 6.23521342e+01, 6.21516508e+01,
       6.26416792e+01, 6.27716722e+01, 6.25720828e+01, 6.27213531e+01,
       6.21188968e+01, 6.25426327e+01, 6.22998865e+01, 6.22743363e+01,
       6.21245060e+01, 6.22511340e+01, 6.22068575e+01, 6.24115040e+01,
       6.24780007e+01, 6.25283143e+01, 6.25366005e+01, 6.24586209e+01,
       6.25384763e+01, 6.25078733e+01, 6.24995373e+01, 6.25496712e+01,
       6.25141107e+01, 6.25108798e+01, 6.25456293e+01, 6.25190905e+01,
       6.25279615e+01, 6.24684866e+01, 6.24528223e+01, 6.24886801e+01,
       6.25496037e+01, 6.25120192e+01, 6.25417894e+01, 6.25188399e+01,
       2.91326637e+03, 2.05005156e+03, 7.70656881e+03, 1.15602423e+04,
       1.70370207e+02, 8.89332256e+02, 1.09396294e+04, 8.06604474e+02,
       2.14141575e+03, 1.04811190e+03, 1.03360110e+03, 9.37286782e+02,
       9.29585793e+02, 1.72098327e+02, 2.12326955e+03, 6.40344208e+02,
       1.55552659e+04, 1.01816848e+02, 2.38425891e+02, 5.64297822e+03,
       3.64051255e+02, 7.64775052e+02, 2.24742156e+03, 9.81694749e+02,
       2.36202792e+02, 6.04815473e+03, 2.94259302e+03, 1.58214219e+04,
       1.26036602e+02, 1.14800260e+02, 1.77581379e+02, 1.30488425e+03,
       1.33916014e+02, 1.91483013e+03, 1.85642300e+02, 1.02863620e+04,
       3.48081726e+02, 4.96790690e+03, 3.34745678e+02, 1.63621039e+03,
       3.61557438e+02, 2.78302202e+02, 1.14258234e+04, 4.19677599e+03,
       1.20027049e+04, 1.63549464e+04, 1.16367164e+03, 8.73587677e+03,
       1.02107376e+04, 2.31115125e+03, 8.17474331e+03, 6.00862956e+02,
       6.56858857e+03, 6.65504706e+03, 3.19521575e+02, 7.53504539e+02,
       2.70524568e+02, 3.09146559e+02, 3.55894333e+02, 2.41187639e+02])
lambdahstest_vec = np.array([1.11205470e+00, 4.10261475e+01, 1.04662960e-01, 5.75400810e+00,
       5.87133660e+01, 1.63687942e-01, 1.22889178e-02, 2.31843920e-01,
       3.41433749e-06, 8.88630287e-02, 2.09498083e-06, 7.80213619e-02,
       2.57522175e-02, 1.82598373e+01, 6.11541648e-01, 1.42734124e-06,
       1.31769439e-04, 5.22930491e+00, 5.56593994e-01, 4.26006630e-06,
       1.10408611e-01, 1.73655844e-01, 3.65124929e-05, 8.06176265e-04,
       1.33252417e-03, 5.89176357e-01, 3.81587623e-05, 2.44684697e-06,
       5.21416879e-03, 3.09976187e+01, 2.23443295e-01, 2.07770852e-04,
       2.29243660e+00, 6.58016525e-01, 1.89679639e-06, 2.23221199e+00,
       7.87774525e-04, 1.84412869e-05, 1.09825089e+00, 1.40445543e-02,
       3.84182918e-02, 9.50623237e+00, 1.37161258e-06, 1.25630321e-01,
       4.23308695e-01, 7.90385031e-01, 1.15950851e-04, 6.75748795e-03,
       1.42578456e-01, 4.88277235e-04, 1.82561460e-06, 3.48500819e-03,
       2.29913191e+00, 1.57716247e-06, 9.29785970e-06, 1.43694622e-01,
       2.86305754e-04, 2.21191708e-04, 6.15274958e+00, 8.48379538e-05,
       8.26410071e-03, 1.70569818e-06, 1.25523911e-04, 2.30933995e-02,
       4.87275918e-01, 2.43808554e-04, 2.19095178e-04, 1.15948780e-04,
       8.16374199e-03, 1.82897846e-04, 5.35314470e+01, 3.44755598e-01,
       1.89335277e-02, 4.35519352e-05, 3.99032019e-01, 7.62789103e+01,
       2.53731172e-01, 4.38908975e-03, 5.52604776e-01, 2.06001818e-04,
       7.80321367e-06, 1.23557321e-06, 6.21011983e-04, 7.89562981e-02,
       5.99832067e-02, 1.64732342e-01, 8.12277816e-02, 1.30130794e-04,
       2.32516105e-02, 4.39919133e+00, 8.80779966e-03, 3.95430608e-05,
       1.89890606e-02, 3.38136019e+00, 1.11146677e+00, 4.97572376e-04,
       5.65301174e-06, 1.42481801e-04, 7.43518241e-04, 4.58678278e-04,
       2.87958554e+00, 8.95982702e+01, 1.22706971e-01, 1.76508630e+00,
       1.34769305e+00, 7.16086608e+00, 6.14053149e-03, 1.75414034e-06,
       9.99886068e-02, 1.69458522e-02, 2.11901163e-05, 1.00817883e-03,
       4.06208954e-03, 4.47517024e-02, 7.95305658e-02, 6.25085187e-02,
       4.77643953e-03, 1.36738730e-06, 9.77507312e-06, 1.79935209e-04,
       2.73180175e-03, 2.72813023e-04, 3.22517141e-04, 1.11856944e-06,
       4.78707682e+00, 1.09126397e+00, 3.72162329e+00, 1.74992132e+00,
       9.06643693e+01, 3.43358273e-02, 3.73179333e-06, 3.56176202e+00,
       5.53822255e-05, 6.24544228e-05, 6.97589970e+01, 1.18998103e-01,
       1.19610243e-04, 9.70496363e+00, 9.85314909e-06, 4.78813614e-02,
       5.43250746e+00, 1.36095112e+00, 3.99910315e+01, 7.47912817e-03,
       1.90720801e-03, 4.93449031e-03, 2.55312044e-04, 1.34961874e+00,
       4.08399636e-01, 2.24997284e+01, 3.97538380e-02, 1.51509054e-05,
       1.61651325e-05, 2.69228938e-02, 1.26981423e-06, 9.35539902e-04,
       1.11374071e-02, 7.44235103e-01, 1.27246020e-03, 1.30468875e-02,
       4.68400354e+00, 3.97382050e+01, 1.98384941e-01, 4.27975673e-01,
       2.34294136e+01, 2.17756874e+01, 2.77299383e-01, 7.69860558e-02,
       1.96334118e+01, 1.65455984e-02, 1.66915196e+01, 5.29812616e-01,
       6.14742575e-02, 7.33338459e+01, 1.15442891e-01, 1.89953412e+00,
       8.50892312e+00, 2.02521423e-03, 3.61459130e-03, 1.20085068e-04,
       5.70672845e-03, 5.73801064e-04, 1.62558854e-06, 4.01735697e+01,
       3.28114901e-06, 7.46154610e+01, 4.69580168e+01, 1.16195391e-04,
       2.00975438e-01, 5.85581809e+01, 3.66161901e-03, 2.35433441e-03,
       2.76452650e+00, 1.00730341e-01, 6.45351207e-03, 3.46067940e-04,
       2.32647936e+01, 4.71834759e-03, 4.47101924e-06, 7.59388119e+01])
logenergyx_bins = np.arange(-8.955,-0.000,0.09)
x_vec = pow(10.,logenergyx_bins)
x_bin = np.zeros(len(x_vec))
for t in range(len(x_vec)):
    x_bin[t] = pow(10.,logenergyx_bins[t]+0.045)-pow(10.,logenergyx_bins[t]-0.045)

print('Testing the Dark matter spectra for gamma rays')

print('Checking ', end="", flush=True)
for t in range(len(massestest_vec)):
    if lambdahstest_vec[t]<1e0 and massestest_vec[t]>1e1:
        if t % 10 ==0 and t!=0:
            print('.', end="", flush=True)
        fluxresult_vec = np.zeros(len(logenergyx_bins))
        if t+1<10:
            table = np.loadtxt('/Users/mattiadimauro/Dropbox/MadDM/MG5_aMC_v2_9_9_iMac/test_SHP_spectra_sigmav/output/run_0%d/gammas_spectrum_pythia8.dat'%(t+1))
        else:
            table = np.loadtxt('/Users/mattiadimauro/Dropbox/MadDM/MG5_aMC_v2_9_9_iMac/test_SHP_spectra_sigmav/output/run_%d/gammas_spectrum_pythia8.dat'%(t+1))
    
        N_table = table[:,1]*np.log10(np.e)/pow(10.,logenergyx_bins)*x_bin
        N_int = DMspectra_inttable(massestest_vec[t],lambdahstest_vec[t],'gammas')[1]*np.log10(np.e)/pow(10.,logenergyx_bins)*x_bin
        if N_table.sum()/N_int.sum()<1.-resolution or N_table.sum()/N_int.sum()>1.+resolution:
            print(r'WARNING IMPRECISE CALCULATION FOR')
            print(r'm_S=%.3f and lambdahs=%.3e'%(massestest_vec[t],lambdahstest_vec[t]))
            print(r'Error is %.3f'%(N_table.sum()/N_int.sum()))
            if makeplots==True:
                fig = plt.figure(figsize=(8,6))
                plt.plot(logenergyx_bins,table[:,1], lw=1.5, ls=':', color='red', label=r'MadDM')
                plt.plot(logenergyx_bins,DMspectra_inttable(massestest_vec[t],lambdahstest_vec[t],'gammas')[1], lw=1.5, ls='-.', color='blue', label=r'Interp, table')
                plt.ylabel(r'$dN/dlog(x)$', fontsize=18)
                plt.xlabel(r'$log(x)$', fontsize=18)
                plt.axis([-6,0.2,1e-3,1e2])
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=18)
                plt.tick_params('both', length=7, width=2, which='major')
                plt.tick_params('both', length=5, width=2, which='minor')
                plt.grid(True)
                plt.yscale('log')
                plt.xscale('linear')
                plt.legend(loc=2,prop={'size':16},numpoints=1, scatterpoints=1, ncol=1)
                fig.tight_layout(pad=0.5)
                plt.savefig("/Users/mattiadimauro/Dropbox/MadDM/MG5_aMC_v2_9_9/SHP_paper_data_scripts/code/SingletScalar_DM/tests/test_spectra/test_gammas_spectra_%d.pdf"%t)
    
                fig = plt.figure(figsize=(8,6))
                plt.plot(logenergyx_bins,table[:,1]/DMspectra_inttable(massestest_vec[t],lambdahstest_vec[t],'gammas')[1], lw=1.5, ls='-.', color='blue', label=r'MadDM/Interp table')
                plt.ylabel(r'ratio', fontsize=18)
                plt.xlabel(r'$log(x)$', fontsize=18)
                plt.axis([-6,0.2,0.5,2e0])
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=18)
                plt.tick_params('both', length=7, width=2, which='major')
                plt.tick_params('both', length=5, width=2, which='minor')
                plt.grid(True)
                plt.yscale('log')
                plt.xscale('linear')
                plt.legend(loc=2,prop={'size':16},numpoints=1, scatterpoints=1, ncol=1)
                fig.tight_layout(pad=0.5)
                plt.savefig("/Users/mattiadimauro/Dropbox/MadDM/MG5_aMC_v2_9_9/SHP_paper_data_scripts/code/SingletScalar_DM/tests/test_spectra/testratio_gammas_spectra_%d.pdf"%t)
print('')


print('Testing the Dark matter spectra for positrons')

print('Checking ', end="", flush=True)
for t in range(len(massestest_vec)):
    if lambdahstest_vec[t]<1e0 and massestest_vec[t]>1e1:
        if t % 10 ==0 and t!=0:
            print('.', end="", flush=True)
        fluxresult_vec = np.zeros(len(logenergyx_bins))
        if t+1<10:
            table = np.loadtxt('/Users/mattiadimauro/Dropbox/MadDM/MG5_aMC_v2_9_9_iMac/test_SHP_spectra_sigmav/output/run_0%d/positrons_spectrum_pythia8.dat'%(t+1))
        else:
            table = np.loadtxt('/Users/mattiadimauro/Dropbox/MadDM/MG5_aMC_v2_9_9_iMac/test_SHP_spectra_sigmav/output/run_%d/positrons_spectrum_pythia8.dat'%(t+1))
    
        N_table = table[:,1]*np.log10(np.e)/pow(10.,logenergyx_bins)*x_bin
        N_int = DMspectra_inttable(massestest_vec[t],lambdahstest_vec[t],'positrons')[1]*np.log10(np.e)/pow(10.,logenergyx_bins)*x_bin
        if N_table.sum()/N_int.sum()<1.-resolution or N_table.sum()/N_int.sum()>1.+resolution:
            print(r'WARNING IMPRECISE CALCULATION FOR')
            print(r'm_S=%.3f and lambdahs=%.3e'%(massestest_vec[t],lambdahstest_vec[t]))
            print(r'Error is %.3f'%(N_table.sum()/N_int.sum()))
            if makeplots==True:
                fig = plt.figure(figsize=(8,6))
                plt.plot(logenergyx_bins,table[:,1], lw=1.5, ls=':', color='red', label=r'MadDM')
                plt.plot(logenergyx_bins,DMspectra_inttable(massestest_vec[t],lambdahstest_vec[t],'positrons')[1], lw=1.5, ls='-.', color='blue', label=r'Interp, table')
                plt.ylabel(r'$dN/dlog(x)$', fontsize=18)
                plt.xlabel(r'$log(x)$', fontsize=18)
                plt.axis([-6,0.2,1e-3,5e1])
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=18)
                plt.tick_params('both', length=7, width=2, which='major')
                plt.tick_params('both', length=5, width=2, which='minor')
                plt.grid(True)
                plt.yscale('log')
                plt.xscale('linear')
                plt.legend(loc=2,prop={'size':16},numpoints=1, scatterpoints=1, ncol=1)
                fig.tight_layout(pad=0.5)
                plt.savefig("/Users/mattiadimauro/Dropbox/MadDM/MG5_aMC_v2_9_9/SHP_paper_data_scripts/code/SingletScalar_DM/tests/test_spectra/test_positrons_spectra_%d.pdf"%t)
    
                fig = plt.figure(figsize=(8,6))
                plt.plot(logenergyx_bins,table[:,1]/DMspectra_inttable(massestest_vec[t],lambdahstest_vec[t],'positrons')[1], lw=1.5, ls='-.', color='blue', label=r'MadDM/Interp table')
                plt.ylabel(r'ratio', fontsize=18)
                plt.xlabel(r'$log(x)$', fontsize=18)
                plt.axis([-6,0.2,0.5,2e0])
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=18)
                plt.tick_params('both', length=7, width=2, which='major')
                plt.tick_params('both', length=5, width=2, which='minor')
                plt.grid(True)
                plt.yscale('log')
                plt.xscale('linear')
                plt.legend(loc=2,prop={'size':16},numpoints=1, scatterpoints=1, ncol=1)
                fig.tight_layout(pad=0.5)
                plt.savefig("/Users/mattiadimauro/Dropbox/MadDM/MG5_aMC_v2_9_9/SHP_paper_data_scripts/code/SingletScalar_DM/tests/test_spectra/testratio_positrons_spectra_%d.pdf"%t)
print('')


print('Testing the Dark matter spectra for antiprotons')

print('Checking ', end="", flush=True)
for t in range(len(massestest_vec)):
    if lambdahstest_vec[t]<1e0 and massestest_vec[t]>1e1:
        if t % 10 ==0 and t!=0:
            print('.', end="", flush=True)
        fluxresult_vec = np.zeros(len(logenergyx_bins))
        if t+1<10:
            table = np.loadtxt('/Users/mattiadimauro/Dropbox/MadDM/MG5_aMC_v2_9_9_iMac/test_SHP_spectra_sigmav/output/run_0%d/antiprotons_spectrum_pythia8.dat'%(t+1))
        else:
            table = np.loadtxt('/Users/mattiadimauro/Dropbox/MadDM/MG5_aMC_v2_9_9_iMac/test_SHP_spectra_sigmav/output/run_%d/antiprotons_spectrum_pythia8.dat'%(t+1))

        N_table = table[:,1]*np.log10(np.e)/pow(10.,logenergyx_bins)*x_bin
        N_int = DMspectra_inttable(massestest_vec[t],lambdahstest_vec[t],'antiprotons')[1]*np.log10(np.e)/pow(10.,logenergyx_bins)*x_bin
        if N_table.sum()/N_int.sum()<1.-resolution or N_table.sum()/N_int.sum()>1.+resolution:
            print(r'WARNING IMPRECISE CALCULATION FOR')
            print(r'm_S=%.3f and lambdahs=%.3e'%(massestest_vec[t],lambdahstest_vec[t]))
            print(r'Error is %.3f'%(N_table.sum()/N_int.sum()))
            if makeplots==True:
                fig = plt.figure(figsize=(8,6))
                plt.plot(logenergyx_bins,table[:,1], lw=1.5, ls=':', color='red', label=r'MadDM')
                plt.plot(logenergyx_bins,DMspectra_inttable(massestest_vec[t],lambdahstest_vec[t],'antiprotons')[1], lw=1.5, ls='-.', color='blue', label=r'Interp, table')
                plt.ylabel(r'$dN/dlog(x)$', fontsize=18)
                plt.xlabel(r'$log(x)$', fontsize=18)
                plt.axis([-6,0.2,1e-3,3e1])
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=18)
                plt.tick_params('both', length=7, width=2, which='major')
                plt.tick_params('both', length=5, width=2, which='minor')
                plt.grid(True)
                plt.yscale('log')
                plt.xscale('linear')
                plt.legend(loc=2,prop={'size':16},numpoints=1, scatterpoints=1, ncol=1)
                fig.tight_layout(pad=0.5)
                plt.savefig("/Users/mattiadimauro/Dropbox/MadDM/MG5_aMC_v2_9_9/SHP_paper_data_scripts/code/SingletScalar_DM/tests/test_spectra/test_antiprotons_spectra_%d.pdf"%t)
    
                fig = plt.figure(figsize=(8,6))
                plt.plot(logenergyx_bins,table[:,1]/DMspectra_inttable(massestest_vec[t],lambdahstest_vec[t],'antiprotons')[1], lw=1.5, ls='-.', color='blue', label=r'MadDM/Interp table')
                plt.ylabel(r'ratio', fontsize=18)
                plt.xlabel(r'$log(x)$', fontsize=18)
                plt.axis([-6,0.2,0.5,2e0])
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=18)
                plt.tick_params('both', length=7, width=2, which='major')
                plt.tick_params('both', length=5, width=2, which='minor')
                plt.grid(True)
                plt.yscale('log')
                plt.xscale('linear')
                plt.legend(loc=2,prop={'size':16},numpoints=1, scatterpoints=1, ncol=1)
                fig.tight_layout(pad=0.5)
                plt.savefig("/Users/mattiadimauro/Dropbox/MadDM/MG5_aMC_v2_9_9/SHP_paper_data_scripts/code/SingletScalar_DM/tests/test_spectra/testratio_antiprotons_spectra_%d.pdf"%t)