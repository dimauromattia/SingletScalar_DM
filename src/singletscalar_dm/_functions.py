'''Main functions

Sub-module containing the main package functions.
Functions defined here are also present in the ``singletscalar_dm`` namespace.

'''
import numpy as np
from scipy.interpolate import (interp1d, interp2d, CubicSpline)
from scipy.optimize import minimize
from ._globals import *

__all__ = [
    'interpolate_Omega',
    'interpolate_Omega_MicrOMEGAs',
    'interpolate_lambda',
    'interpolate_lambda_MicrOMEGAs',
    '_lambda2sigmav',
    'sigmav_channels',
    'DMspectra_inttable',
    'provide_ULEXP',
    'SI_noomega',
    'SI_withomega',
    'GetUL_DD_nomega',
    'GetUL_DD_withomega',
    'minimize_br_inv',
    'Gamma_inv',
    'Br_inv',
    'Br_inv_UL',
    'func_interpolate',
    'flux_DM_prompt'
]


def interpolate_relicdensity(mass_val,QCDmodel):
    '''Calculates the lambda_hs value for which we obtain the observed DM relic density (Omegah^2=0.120).

    Parameters
    ----------
    mass_val : np.ndarray
        Dark matter mass values in GeV.
    QCDmodel : {'MICROMEGAs', 'QCDA', 'QCDB'}
        Model for the QCD phase transition.

    Return
    ------
    lambda_val : np.ndarray
        Values of the lambda_HS parameter.

    Notes
    -----
    The computation of the relic density has been obtained via the code DRAKE and MICROMEGAs.
    '''
        
    if QCDmodel=='QCDA':
        column = 1
    elif QCDmodel=='QCDB':
        column = 2
    table_RD_FB = np.loadtxt(folder+'/Omega_MicroOMEGAs_DRAKE_QCDB_QCDA.dat')
    mass_vec = table_RD_FB[:,0]
    lambda_vec = table_RD_FB[:,column]
    func = interp1d(mass_vec,lambda_vec)
    return func(mass_val)


def interpolate_Omega(mass_val,lambda_val,QCDmodel,warningprint):
    '''Calculates the relic density as Omega h^2 given the mass and lambda.

    Parameters
    ----------
    mass_val : np.ndarray
        Dark matter mass values in GeV.
    lambda_val : np.ndarray
        Values of the lambda_HS parameter.
    QCDmodel : {'QCDA', 'QCDB'}
        Model for the QCD phase transition.

    Return
    ------
    Omega_val : np.ndarray
        The values of the relic density interpolated on the 2d grid of mass
        and coupling.

    Notes
    -----
    The computation of the relic density has been obtained via the code DRAKE and MICROMEGAs.
    '''

    if mass_val>mass_vector_drake.min() and mass_val<mass_vector_drake.max():

        bins=100
        table = np.loadtxt(import_data_file('tableOmega_DRAKE_fBE_%s.dat'%QCDmodel))
        idx = np.searchsorted(mass_vector_drake, mass_val, side='right')-1

        cont1 = idx*bins
        cont2 = cont1+100
        Omega_mass1 = table[cont1:cont2,2]
        Omega_mass2 = table[(cont1+100):(cont2+100),2]
        lambda_mass1 = table[cont1:cont2,1]
        lambda_mass2 = table[(cont1+100):(cont2+100),1]

        check = 0
        if lambda_val<lambda_mass1.max() and lambda_val>lambda_mass1.min():
            func1 = interp1d(lambda_mass1,np.sqrt(Omega_mass1))
            omega1 = np.power(float(func1(lambda_val)),2.)
            #Below there are other two ways to interpolate.
            #func1 = interp1d(np.sqrt(lambda_mass1),Omega_mass1)
            #omega1 = float(func1(sqrt(lambda_val)))
            #func1 = interp1d(lambda_mass1,np.sqrt(Omega_mass1))
            #omega1 = np.power(float(func1(lambda_val)),2.)
        else:
            if warningprint==True:
                print('Warning, extrapolating.....')
                print('The relic density with DRAKE has been calculated for lambda_HS between ',lambda_mass1.min(),lambda_mass1.max())
            func1 = CubicSpline(lambda_mass1,np.sqrt(Omega_mass1), bc_type='not-a-knot')
            omega1 = np.power(float(func1(lambda_val)),2.)
            check = 1

        if lambda_val<lambda_mass2.max() and lambda_val>lambda_mass2.min():
            func2 = interp1d(lambda_mass2,np.sqrt(Omega_mass2))
            omega2 = np.power(float(func2(lambda_val)),2)
            #Below there are other two ways to interpolate.
            #func2 = interp1d(np.sqrt(lambda_mass2),Omega_mass2)
            #omega2 = float(func2(sqrt(lambda_val)))
            #func2 = interp1d(lambda_mass2,np.sqrt(Omega_mass2))
            #omega2 = np.power(float(func2(lambda_val)),2)
        else:
            func2 = CubicSpline(lambda_mass2,np.sqrt(Omega_mass2), bc_type='not-a-knot')
            omega2 = np.power(float(func2(lambda_val)),2)
            if check==0:
                if warningprint==True:
                    print('Warning, extrapolating.')
                    print('The relic density with DRAKE has been calculated for lambda_HS between ',lambda_mass1.min(),lambda_mass1.max())
        
        func_final = interp1d(np.array([mass_vector_drake[idx],mass_vector_drake[idx+1]]),np.array([np.power(omega1,0.1),np.power(omega2,0.1)]),kind='linear')
        #interpolate with CubicSpline
        #func_final = CubicSpline(np.array([mass_vector_drake[idx],mass_vector_drake[idx+1]]),np.array([np.power(omega1,0.1),np.power(omega2,0.1)]), bc_type='not-a-knot')
        return float(np.power(func_final(mass_val),1/0.1))

    else:
        if mass_val>mass_vector_micro.min() and mass_val<mass_vector_micro.max():
            table = np.loadtxt(import_data_file('tableOmega_MicroOMEGAs.dat'))
            Omega_table = table[:,2]
            lambdahs_vec = np.logspace(-5,2,500)
            func  = interp2d(mass_vector_micro,lambdahs_vec,Omega_table,kind='cubic')
            return float(func(mass_val,lambda_val))
        else:
            print('Warning, mass should be between',mass_vector_micro.min(),mass_vector_micro.max())

def interpolate_Omega_MicrOMEGAs(mass_val,lambda_val):
    '''Calculates the relic density as Omega h^2 given the mass and lambda.

    Parameters
    ----------
    mass_val : np.ndarray
        Dark matter mass values in GeV.
    lambda_val : np.ndarray
        Values of the lambda_HS parameter.

    Return
    ------
    Omega_val : np.ndarray
        The value of the relic density interpolated on the 2d grid of mass
        and coupling.

    Notes
    -----
    The computation of the relic density has been obtained via the code MicrOMEGAs.
    '''

    if mass_val>mass_vector_micro.min() and mass_val<mass_vector_micro.max():
        table = np.loadtxt(import_data_file('tableOmega_MicroOMEGAs.dat'))
        Omega_table = table[:,2]
        lambdahs_vec = np.logspace(-5,2,500)
        func  = interp2d(mass_vector_micro,lambdahs_vec,Omega_table,kind='cubic')
        return float(func(mass_val,lambda_val))
    else:
        print('Warning, mass should be between',mass_vector_micro.min(),mass_vector_micro.max())

def interpolate_lambda(mass_val,Omega_val,QCDmodel,warningprint):
    '''Calculates the lambda_HS parameter for given relic density and mass values.

    Parameters
    ----------
    mass_val : np.ndarray
        Dark matter mass values in GeV.
    Omega_val : np.ndarray
        The values of the relic density interpolated on the 2d grid of mass
    QCDmodel : {'QCDA', 'QCDB'}
        Model for the QCD phase transition.

    Return
    ------
    lambda_val : np.ndarray
        Values of the lambda_HS parameter, obtained through interpolation.

    Notes
    -----
    The computation of the parameter uses the computation of the relic density
    obtained via the code DRAKE and MICROMEGAs.
    '''

    if mass_val>mass_vector_drake.min() and mass_val<mass_vector_drake.max():

        bins=100
        table = np.loadtxt(import_data_file('tableOmega_DRAKE_fBE_%s.dat'%QCDmodel))
        idx = np.searchsorted(mass_vector_drake, mass_val, side='right')-1

        cont1 = idx*bins
        cont2 = cont1+100
        Omega_mass1 = table[cont1:cont2,2]
        Omega_mass2 = table[(cont1+100):(cont2+100),2]
        lambda_mass1 = table[cont1:cont2,1]
        lambda_mass2 = table[(cont1+100):(cont2+100),1]

        check = 0
        if Omega_val<Omega_mass1.max() and Omega_val>Omega_mass1.min():
            func1 = interp1d(Omega_mass1,lambda_mass1)
            lambda1 = float(func1(Omega_val))
        else:
            if warningprint==True:
                print('Warning, extrapolating.',mass_val)
                print('For this mass pick a range of Omegah^2 between ',Omega_mass1.min(),Omega_mass1.max())
            check = 1
            #func1 = CubicSpline(Omega_mass1,np.log10(lambda_mass1), bc_type='not-a-knot')
            #lambda1 = float(func1(Omega_val))

        if Omega_val<Omega_mass2.max() and Omega_val>Omega_mass2.min():
            func2 = interp1d(Omega_mass2,lambda_mass2)
            lambda2 = float(func2(Omega_val))
        else:
            if check==0:
                if warningprint==True:
                    print('Warning, extrapolating.',mass_val)
                    print('For this mass pick a range of Omegah^2 between ',Omega_mass2.min(),Omega_mass2.max())
                check = 1
                #func2 = CubicSpline(Omega_mass2,lambda_mass2, bc_type='not-a-knot')
                #lambda2 = float(func2(Omega_val))
        if check==0:
            func_final = interp1d(np.array([mass_vector_drake[idx],mass_vector_drake[idx+1]]),np.array([lambda1,lambda2]))
            return func_final(mass_val)
        else:
            return 0
            #func_final = interp1d(np.array([mass_vector_drake[idx],mass_vector_drake[idx+1]]),np.array([lambda1,lambda2]))
            #return func_final(mass_val)

    else:
        if mass_val>mass_vector_micro.min() and mass_val<mass_vector_micro.max():
            table = np.loadtxt(import_data_file('tableOmega_MicroOMEGAs.dat'))
            bins = 500
            lambdahs_vec = np.logspace(-5,2,bins)
            bins=len(lambdahs_vec)
            idx = np.searchsorted(mass_vector_micro, mass_val, side='right')-1

            cont1 = idx*bins
            cont2 = cont1+bins
            Omega_mass1 = table[cont1:cont2,2]
            Omega_mass2 = table[(cont1+bins):(cont2+bins),2]
            lambda_mass1 = table[cont1:cont2,1]
            lambda_mass2 = table[(cont1+bins):(cont2+bins),1]

            check = 0
            if Omega_val<Omega_mass1.max() and Omega_val>Omega_mass1.min():
                check_bin = 0
                for t in range(len(Omega_mass1)):
                    if Omega_mass1[t]<Omega_val and check_bin==0:
                        check_bin = t
                if Omega_mass1[check_bin-1]<Omega_mass1[check_bin+1]:
                    func1 = interp1d(Omega_mass1[check_bin-1:check_bin+1],lambda_mass1[check_bin-1:check_bin+1])
                else:
                    func1 = interp1d(np.array([Omega_mass1[check_bin],Omega_mass1[check_bin-1]]),np.array([lambda_mass1[check_bin],lambda_mass1[check_bin-1]]))

                lambda1 = float(func1(Omega_val))
            else:
                if warningprint==True:
                    print('Warning, extrapolating.')
                    print('For this mass pick a range of Omegah^2 between ',Omega_mass1.min(),Omega_mass1.max())
                check = 1
                #func1 = CubicSpline(Omega_mass1,np.log10(lambda_mass1), bc_type='not-a-knot')
                #lambda1 = float(func1(Omega_val))

            if Omega_val<Omega_mass2.max() and Omega_val>Omega_mass2.min():
                check_bin = 0
                for t in range(len(Omega_mass2)):
                    if Omega_mass2[t]<Omega_val and check_bin==0:
                        check_bin = t
                if Omega_mass2[check_bin-1]<Omega_mass2[check_bin+1]:
                    #print(Omega_val,Omega_mass2[check_bin-1:check_bin+1],lambda_mass2[check_bin-1:check_bin+1])
                    func2 = interp1d(Omega_mass2[check_bin-1:check_bin+1],lambda_mass2[check_bin-1:check_bin+1])
                else:
                    #print(Omega_val,np.array([Omega_mass2[check_bin],Omega_mass2[check_bin-1]]),np.array([lambda_mass2[check_bin],lambda_mass2[check_bin-1]]))
                    func2 = interp1d(np.array([Omega_mass2[check_bin],Omega_mass2[check_bin-1]]),np.array([lambda_mass2[check_bin],lambda_mass2[check_bin-1]]))
                lambda2 = float(func2(Omega_val))
            else:
                if check==0:
                    if warningprint==True:
                        print('Warning, extrapolating.')
                        print('For this mass pick a range of Omegah^2 between ',Omega_mass2.min(),Omega_mass2.max())
                    check = 1
                    #func2 = CubicSpline(Omega_mass2,lambda_mass2, bc_type='not-a-knot')
                    #lambda2 = float(func2(Omega_val))
            if check==0:
                func_final = interp1d(np.array([mass_vector_micro[idx],mass_vector_micro[idx+1]]),np.array([lambda1,lambda2]))
                lambda_val = float(func_final(mass_val))
                if mass_val<70 and lambda_val>1.0 and warningprint==True:
                    print('Warning the problem f(lambda)=Omega h^2 could have two solutions for lambda')
                return lambda_val
            else:
                return 0
                #func_final = interp1d(np.array([mass_vector_micro[idx],mass_vector_micro[idx+1]]),np.array([lambda1,lambda2]))
                #lambda_val = float(func_final(mass_val))
                #return lambda_val
        else:
            print('Warning, mass should be between',mass_vector_micro.min(),mass_vector_micro.max())


def interpolate_lambda_MicrOMEGAs(mass_val,Omega_val,warningprint):
    '''Calculates the lambda_HS parameter for given relic density and mass values.

    Parameters
    ----------
    mass_val : np.ndarray
        Dark matter mass values in GeV.
    Omega_val : np.ndarray
        The value of the relic density interpolated on the 2d grid of mass

    Return
    ------
    lambda_val : np.ndarray
        Values of the lambda_HS parameter, obtained through interpolation.

    Notes
    -----
    The computation of the parameter uses the computation of the relic density
    obtained via the code MicrOMEGAs.
    '''

    if mass_val>mass_vector_micro.min() and mass_val<mass_vector_micro.max():
        table = np.loadtxt(import_data_file('tableOmega_MicroOMEGAs.dat'))
        bins = 500
        lambdahs_vec = np.logspace(-5,2,bins)
        bins=len(lambdahs_vec)
        idx = np.searchsorted(mass_vector_micro, mass_val, side='right')-1

        cont1 = idx*bins
        cont2 = cont1+bins
        Omega_mass1 = table[cont1:cont2,2]
        Omega_mass2 = table[(cont1+bins):(cont2+bins),2]
        lambda_mass1 = table[cont1:cont2,1]
        lambda_mass2 = table[(cont1+bins):(cont2+bins),1]

        check = 0
        if Omega_val<Omega_mass1.max() and Omega_val>Omega_mass1.min():
            check_bin = 0
            for t in range(len(Omega_mass1)):
                if Omega_mass1[t]<Omega_val and check_bin==0:
                    check_bin = t
            if Omega_mass1[check_bin-1]<Omega_mass1[check_bin+1]:
                func1 = interp1d(Omega_mass1[check_bin-1:check_bin+1],lambda_mass1[check_bin-1:check_bin+1])
            else:
                func1 = interp1d(np.array([Omega_mass1[check_bin],Omega_mass1[check_bin-1]]),np.array([lambda_mass1[check_bin],lambda_mass1[check_bin-1]]))

            lambda1 = float(func1(Omega_val))
        else:
            if warningprint==True:
                print('Warning, extrapolating.')
                print('For this mass pick a range of Omegah^2 between ',Omega_mass1.min(),Omega_mass1.max())
            check = 1
            #func1 = CubicSpline(Omega_mass1,np.log10(lambda_mass1), bc_type='not-a-knot')
            #lambda1 = float(func1(Omega_val))

        if Omega_val<Omega_mass2.max() and Omega_val>Omega_mass2.min():
            check_bin = 0
            for t in range(len(Omega_mass2)):
                if Omega_mass2[t]<Omega_val and check_bin==0:
                    check_bin = t
            if Omega_mass2[check_bin-1]<Omega_mass2[check_bin+1]:
                func2 = interp1d(Omega_mass2[check_bin-1:check_bin+1],lambda_mass2[check_bin-1:check_bin+1])
            else:
                func2 = interp1d(np.array([Omega_mass2[check_bin],Omega_mass2[check_bin-1]]),np.array([lambda_mass2[check_bin],lambda_mass2[check_bin-1]]))
            lambda2 = float(func2(Omega_val))
        else:
            if check==0:
                if warningprint==True:
                    print('Warning, extrapolating.')
                    print('For this mass pick a range of Omegah^2 between ',Omega_mass2.min(),Omega_mass2.max())
                check = 1
                #func2 = CubicSpline(Omega_mass2,lambda_mass2, bc_type='not-a-knot')
                #lambda2 = float(func2(Omega_val))
        if check==0:
            func_final = interp1d(np.array([mass_vector_micro[idx],mass_vector_micro[idx+1]]),np.array([lambda1,lambda2]))
            lambda_val = float(func_final(mass_val))
            if mass_val<70 and lambda_val>1.0 and warningprint==True:
                print('Warning the problem f(lambda)=Omega h^2 could have two solutions for lambda')
            return lambda_val
        else:
            return 0
            #func_final = interp1d(np.array([mass_vector_micro[idx],mass_vector_micro[idx+1]]),np.array([lambda1,lambda2]))
            #lambda_val = float(func_final(mass_val))
            #return lambda_val
    else:
        print('Warning, mass should be between',mass_vector_micro.min(),mass_vector_micro.max())



def _lambda2sigmav(DMmass_val,lambdahs_val,table_int):
    '''Returns the interpolated :math:`\sigma v` value on mass and coupling.

    The value is given in units of cm^3/s.

    Parameters
    ----------
    DMmass_val : np.ndarray
        Dark matter mass values in GeV.
    lambdahs_val : np.ndarray
        Values of the lambda_HS parameter.
    table_int : np.ndarray
        Array containing all the computed values of the :math:`\sigma v`.

    Return
    ------
    sigmav : np.ndarray
        The values of the :math:`\sigma v`.
    '''
    if lambdahs_val < 1e-4:
        func_int = interp1d(massz_vec,table_int[:,1],kind='cubic')
        rescaling_lambda = pow(lambdahs_val/1e-4,2.)
        return float(rescaling_lambda*func_int(DMmass_val))

    elif lambdahs_val > 1e1:
        func_int_1 = interp1d(massz_vec,table_int[:,len(lambdahs_vec)-1])
        func_int_2 = interp1d(massz_vec,table_int[:,len(lambdahs_vec)])
        val1 = func_int_1(DMmass_val)
        val2 = func_int_2(DMmass_val)
        if val1==0:
            val1=1e-50
        if val2==0:
            val2=1e-50
        logres_1 = np.log10(val1)
        logres_2 = np.log10(val2)
        loglambda_1 = np.log10(lambdahs_vec[len(lambdahs_vec)-2])
        loglambda_2  = np.log10(lambdahs_vec[len(lambdahs_vec)-1])
        #rescaling_lambda = pow(lambdahs_val/1e1,2.)
        result = logres_2 + (logres_2-logres_1)*(np.log10(lambdahs_val)-loglambda_2)/(loglambda_2-loglambda_1)
        #print(pow(10.,logres_1),pow(10.,logres_2),pow(10.,loglambda_1),pow(10.,loglambda_2),pow(10.,result))
        return pow(10.,result)

    else:
        func_int = interp2d( massz_vec, lambdahs_vec , np.ndarray.flatten( table_int[:,1:(len(lambdahs_vec)+1)]),kind='cubic')
        return float(func_int( DMmass_val, lambdahs_val) )


def sigmav_channels(DMmass_val,lambdahs_val,channel):
    '''Calculate the :math:`\sigma v` for a given value of mass and lambda_HS for a certain channel.

    The value is given in units of cm^3/s.

    Parameters
    ----------
    DMmass_val : np.ndarray
        Dark matter mass values in GeV.
    lambdahs_val : np.ndarray
        Values of the lambda_HS parameter.
    channel : {'tot', 'cc', 'bb', 'tt', 'tautau', 'gg', 'ww', 'zz', 'hh', 'aa', 'za'}
        The annihilation channel of dark matter.

    Return
    ------
    sigmav : np.ndarray
        The values of the :math:`\sigma v` for the selected channel.

    See Also
    --------
    _lambda2sigma

    Notes
    -----
    This function internally calls `_lambda2sigmav` passing the ``table_int`` as
    the content of the data file related to the ``channel`` specified.
    '''
    #check = np.shape(np.where((channel_vec[:]==channel))[0])
    #print(check,np.where((channel_vec[:]==channel)),np.where((channel_vec[:]==channel)))
    #if check==1:
    if channel=='tot':
        table_int = np.loadtxt(import_data_file('SHP_sigmav_table.dat'))
        sigmav_val = _lambda2sigmav(DMmass_val,lambdahs_val,table_int)

    else:
        table_int = np.loadtxt(import_data_file('SHP_sigmav_%s.dat'%channel))
        sigmav_val = _lambda2sigmav(DMmass_val,lambdahs_val,table_int)

    return sigmav_val
    #else:
    #    print('WARNING, possible channels are:',channel_vec)


def DMspectra_inttable(DMmass_val,lambdahs_val,particle,smooth=False):
    '''Calculates the dark matter source spectra in units of dN/dlog10(x).

    The computation is done for given values of the mass and lambda_HS.

    Parameters
    ----------
    DMmass_val : np.ndarray
        Dark matter mass values in GeV.
    lambdahs_val : np.ndarray
        Values of the lambda_HS parameter.
    particle : {'gammas', 'positrons', 'antiprotons', 'neutrinos_e', 'neutrinos_mu', 'neutrinos_tau'}
        The annihilation channel of dark matter.
    smooth : bool, default=False
        Whether to consider the smoothed spectra.

    Return
    ------
    bins : ndarray
        The log-spaced bins in variable :math:`x = E/m_{DM}`.
    spectra : ndarray
        The values of the spectra.
    '''
    add = ''
    if smooth:
        add = '_smooth'
    tablespectra_int = np.loadtxt(import_data_file('SHP_spectra%s_table_%s.dat'%(add,particle)))

    if DMmass_val<=massz_vec[164]:
        func_int = interp2d(massz_vec,logenergyx_bins,tablespectra_int[:,5])
        dNdx = np.ndarray.flatten(func_int(DMmass_val,logenergyx_bins))
        return logenergyx_bins,dNdx

    elif DMmass_val>massz_vec[164]:

        if lambdahs_val <= 1e-2:
            lambda_index = 2
            func_int = interp2d(massz_vec,logenergyx_bins,tablespectra_int[:,lambda_index])
            result_int = func_int(DMmass_val,logenergyx_bins)
            dNdx = np.ndarray.flatten(result_int)
            return logenergyx_bins,dNdx
        elif lambdahs_val >= 1e1:
            lambda_index = 2 + len(lambda_vec)-1
            func_int = interp2d(massz_vec,logenergyx_bins,tablespectra_int[:,lambda_index])
            result_int = func_int(DMmass_val,logenergyx_bins)
            dNdx = np.ndarray.flatten(result_int)
            return logenergyx_bins,dNdx
        else:
            check=0
            for t in range(len(lambda_vec)):
                if lambdahs_val<lambda_vec[t] and check==0:
                    contl = t-1
                    check=1
            lambda_index = int(contl) + 2

            func_int_1 = interp2d(massz_vec,logenergyx_bins,tablespectra_int[:,lambda_index])
            func_int_2 = interp2d(massz_vec,logenergyx_bins,tablespectra_int[:,lambda_index+1])

            #print(lambda_index,table_int[len(table_int)-100:len(table_int)-1,lambda_index],lambda_index+1,table_int[len(table_int)-100:len(table_int)-1,lambda_index+1])
            flux_a = np.array(func_int_1(DMmass_val,logenergyx_bins))
            flux_b = np.array(func_int_2(DMmass_val,logenergyx_bins))
            #for t in range(len(flux_a)):
            #    print(logenergyx_bins[t],flux_a[t],flux_b[t])
            #print('a,b',lambda_index,flux_a,flux_b)
            #print(lambda_vec[lambda_index-2],lambda_vec[lambda_index-1],flux_a,flux_b)
            #func = interp1d(np.array([pow(lambda_vec[lambda_index-2],2.),pow(lambda_vec[lambda_index-1],2.)]),np.array([flux_a,flux_b]))
            result_int = flux_b + (flux_a - flux_b) * ( pow(lambdahs_val,2.) - pow(lambda_vec[lambda_index-1],2.) ) / ( pow(lambda_vec[lambda_index-2],2.) - pow(lambda_vec[lambda_index-1],2.) )

            dNdx = np.ndarray.flatten(result_int)
            return logenergyx_bins,dNdx

    else:
        print('WARNING MASSES SHOULD BE BETWEEN 2 AND 10000 GeV')
        return 0.



def _LZUL(DMmass):
    ''' Interpolates LZ data on an array of masses. '''
    table = np.loadtxt(import_data_file('LZ_SI_2022_data.dat'))
    mass_vec = table[:,0]
    sigma_SI = table[:,1]*1e-48
    if DMmass<mass_vec[0] or DMmass>mass_vec.max():
        return 0
    else:
        func_int_csi = interp1d(mass_vec,sigma_SI,kind='quadratic')
        return func_int_csi(DMmass)

def _DARWINUL(DMmass):
    ''' Interpolates DARWIN data on an array of masses. '''
    table = np.loadtxt(import_data_file('DARWIN_SI_proj.dat'))
    mass_vec = table[:,0]
    sigma_SI = table[:,1]*1e-48
    if DMmass<mass_vec[0] or DMmass>mass_vec.max():
        return 0
    else:
        func_int_csi = interp1d(mass_vec,sigma_SI,kind='quadratic')
        return func_int_csi(DMmass)

def provide_ULEXP(DMmass,Exp):
    '''Provides the direct detection sigma SI upper limits for specified experiment.

    It interpolates the limits on the array of dark matter masses provided.
    Limits are returned in units cm^2.

    Parameters
    ----------
    DMmass : np.ndarray
        Dark matter mass values in GeV.
    Exp : {'LZ', 'Darwin'}
        The experiment to consider.

    Return
    ------
    ul : ndarray
        The upper limit values sampled at the provided masses.

    Ref
    ---
    '''
    if Exp=='LZ':
        SI_exp = _LZUL(DMmass)
    elif Exp== 'Darwin':
        SI_exp = _DARWINUL(DMmass)
    return SI_exp

def SI_noomega(DMmass_val,lambdahs_val):
    '''Provides the spin-independent cross-section for direct detection in cm^2.

    The computation is done for given values of the mass and lambda_HS.

    Parameters
    ----------
    DMmass_val : np.ndarray
        Dark matter mass values in GeV.
    lambdahs_val : np.ndarray
        Values of the lambda_HS parameter.

    Return
    ------
    sigma_si : ndarray
        The values of the spin-independent cross-section.
    '''
    mu = (mN*DMmass_val)/(mN+DMmass_val)
    val = np.power(lambdahs_val*fN*mu*mN,2.)/(4.*np.pi*pow(mh,4.)*np.power(DMmass_val,2.))
    return GeVm2tocm2*val

def SI_withomega(DMmass,lambda_hs,Lambda_vec,Mass_vec,csi_vec):
    '''Provides the spin-independent cross-section for direct detection in cm^2.

    This function takes into account the relic density for the values of the parameters 
    lambda_hs and DMmass.

    Parameters
    ----------
    DMmass_val : np.ndarray
        Dark matter mass values in GeV.
    lambdahs_val : np.ndarray
        Values of the lambda_HS parameter.
    Lambda_vec : np.ndarray
        Vector of the values of lambda_hs used for the calculation of the relic density
    Mass_vec : np.ndarray
        Vector of the values of mass used for the calculation of the relic density
    csi_vec : np.ndarray
        Vector of the values of csi used for the calculation of the relic density

    Return
    ------
    sigma_si : ndarray
        The values of the spin-independent cross-section.
    '''
    mu = (mN*DMmass)/(mN+DMmass)
    val = np.power(lambda_hs*fN*mu*mN,2.)/(4.*np.pi*pow(mh,4.)*np.power(DMmass,2.))
    func_int_csi = interp2d(Mass_vec,Lambda_vec,csi_vec)
    csi_val = func_int_csi(DMmass,lambda_hs)[0]
    return GeVm2tocm2*val*csi_val

def GetUL_DD_nomega(DMmass_val,Exp):
    '''Provides the direct detection sigma SI upper limits for specified experiment.

    It interpolates the limits on the array of dark matter masses provided.
    Limits are returned in units cm^2.
    The relic relic density and local dark matter density is the one of DM.

    Parameters
    ----------
    DMmass : np.ndarray
        Dark matter mass values in GeV.
    Exp : {'LZ', 'Darwin'}
        The experiment to consider.

    Return
    ------
    ul : ndarray
        The upper limit values sampled at the provided masses.
    '''
    value = 0
    if Exp=='LZ':
        value = np.sqrt(_LZUL(DMmass_val)/SI_noomega(DMmass_val,1.0))
        return value
    elif Exp== 'DARWIN':
        value = np.sqrt(_DARWINUL(DMmass_val)/SI_noomega(DMmass_val,1.0))
        return value
    else:
        print('WARNING: Wrong parameter for the experiment.')
        print('Choose among LZ and DARWIN')

def _tominimize_DD(DMmass,Lambda_val,Lambda_vec,Mass_vec,csi_vec,Exp):
    ''' Function to minimize the difference between the theoretically calculated and measured cross section
    
    Parameters
    ----------
    DMmass_val : np.ndarray
        Dark matter mass values in GeV.
    lambdahs_val : np.ndarray
        Values of the lambda_HS parameter.
    Lambda_vec : np.ndarray
        Vector of the values of lambda_hs used for the calculation of the relic density
    Mass_vec : np.ndarray
        Vector of the values of mass used for the calculation of the relic density
    csi_vec : np.ndarray
        Vector of the values of csi used for the calculation of the relic density

    Return
    ------
    ul : ndarray
        The absolute difference between the theoretical and experimenthal cross section. 
    '''
    #print(Lambda_val)
    if Lambda_val<=Lambda_vec[0]:
        Lambda_val=Lambda_vec[0]
    elif Lambda_val>=Lambda_vec[len(Lambda_vec)-1]:
        Lambda_val=Lambda_vec[len(Lambda_vec)-1]
    SI_theory = SI_withomega(DMmass,Lambda_val,Lambda_vec,Mass_vec,csi_vec)
    SI_exp = provide_ULEXP(DMmass,Exp)

    #print(Lambda_val,SI_exp,SI_theory,abs(SI_exp-SI_theory)/SI_theory)
    #return abs(func_LuxUL(DMmass)-SI_theory)/SI_theory
    return np.abs(SI_exp-SI_theory)

def GetUL_DD_withomega(DMmass,Lambda_vec,Mass_vec,csi_vec,Exp):
    '''Provides the direct detection sigma SI upper limits for specified experiment.

    It interpolates the limits on the array of dark matter masses provided.
    Limits are returned in units cm^2.
    <missing>

    Parameters
    ----------
    DMmass : np.ndarray
        Dark matter mass values in GeV.
    Lambda_vec : np.ndarray
        Vector of the values of lambda_hs used for the calculation of the relic density
    Mass_vec : np.ndarray
        Vector of the values of mass used for the calculation of the relic density
    csi_vec : np.ndarray
        Vector of the values of csi used for the calculation of the relic density
    Exp : {'LZ', 'Darwin'}
        The experiment to consider.

    Return
    ------
    ul : ndarray
        The upper limit values sampled at the provided masses.
    '''
    SI_theory = SI_withomega(DMmass,Lambda_vec,Lambda_vec,Mass_vec,csi_vec)
    SI_exp = provide_ULEXP(DMmass,Exp)
    #print(SI_theory,SI_theory.max(),SI_theory.min(),SI_exp)
    if SI_theory.max()<SI_exp:
        return Lambda_vec.max()

    elif SI_theory.min()>SI_exp:
        return Lambda_vec.min()

    else:
        lambda_0 = 1e-2
        #print(DMmass,Lambda_vec,Mass_vec,csi_vec,Exp)
        #print(minimize(tominimize_DD, lambda_0, method='nelder-mead', bounds=(1e-5,1e1), options={'xatol': 1e-8, 'disp': True}, args=(DMmass,Lambda_vec,Mass_vec,csi_vec,Exp)))
        value = minimize(_tominimize_DD, lambda_0, method='nelder-mead', options={'xatol': 1e-8, 'disp': True}, args=(DMmass,Lambda_vec,Mass_vec,csi_vec,Exp)).x[0]
        SI_theory = SI_withomega(value,DMmass,Lambda_vec,Mass_vec,csi_vec)
        SI_exp = provide_ULEXP(DMmass,Exp)

        if abs(SI_exp-SI_theory) > 1e-47:
            lambda_0 = 1e0
            value = minimize(_tominimize_DD, lambda_0, method='nelder-mead', options={'xatol': 1e-8, 'disp': True}, args=(DMmass,Lambda_vec,Mass_vec,csi_vec,Exp)).x[0]
            SI_theory = SI_withomega(value,DMmass,Lambda_vec,Mass_vec,csi_vec)
            SI_exp = provide_ULEXP(DMmass,Exp)

            if abs(SI_exp-SI_theory) > 1e-47:
                lambda_0 = 1e-4
                value = minimize(_tominimize_DD, lambda_0, method='nelder-mead', options={'xatol': 1e-8, 'disp': True}, args=(DMmass,Lambda_vec,Mass_vec,csi_vec,Exp)).x[0]
                SI_theory = SI_withomega(value,DMmass,Lambda_vec,Mass_vec,csi_vec)
                SI_exp = provide_ULEXP(DMmass,Exp)

                if abs(SI_exp-SI_theory) > 1e-47:
                    lambda_0 = 1e1
                    value = minimize(_tominimize_DD, lambda_0, method='nelder-mead', options={'xatol': 1e-8, 'disp': True}, args=(DMmass,Lambda_vec,Mass_vec,csi_vec,Exp)).x[0]
                    SI_theory = SI_withomega(value,DMmass,Lambda_vec,Mass_vec,csi_vec)
                    SI_exp = provide_ULEXP(DMmass,Exp)

                    if abs(SI_exp-SI_theory) > 1e-47:
                        value = 0.
        return value

def Br_inv_UL(lambda_hs,DMmass,Gamma_H,Gamma_inv_measured):
    '''Calculates the upper limit on the Branching ratio of Higgs boson to invisible particles.

    Parameters
    ----------
    lambda_hs : np.ndarray
        The coupling between Higgs boson and the dark matter.
    DMmass : np.ndarray
        Dark matter mass values in GeV.
    Gamma_H : np.ndarray
        The Higgs width in GeV.
    Gamma_inv_measured : np.ndarray
        The Higgs with to invisible particles as measured experimentally.

    Return
    ------
    Br_inv_ul : ndarray
        The upper limit on Branching ratio of Higgs to invisible particles.
    '''
    return np.abs(Br_inv(lambda_hs,DMmass,Gamma_H)-Gamma_inv_measured)/Gamma_inv_measured

def minimize_br_inv(DMmass,Gamma_H,Gamma_inv_measured):
    '''Returns the upper limits for lambda_hs given the upper limits for the branching ratio of the Higgs boson into invisible particles

    This function utilizes a minimazer to find the upper limit for lambda_hs.
    In particular, it finds the value of lambda_hs for which the theoretical value of Gamma^H_invisible is equal to the observed value.

    Parameters
    ----------
    DMmass : np.ndarray
        Dark matter mass values in GeV.
    Gamma_H : np.ndarray
        The Higgs width in GeV.
    Gamma_inv_measured : np.ndarray
        The Higgs with to invisible particles as measured experimentally.

    Return
    ------
    lambda_ul : ndarray
        Upper limit for lambda_hs found with collider constraints.
    '''
    lambda_0 = 0.01
    lambda_val = minimize(Br_inv_UL, lambda_0, method='nelder-mead', options={'xatol': 1e-8, 'disp': True}, args=(DMmass,Gamma_H,Gamma_inv_measured)).x[0]
    return lambda_val 

def Gamma_inv(DMmass,lambda_hs):
    '''Calculates the invisible width of the Higgs boson in GeV.

    Parameters
    ----------
    DMmass : np.ndarray
        Dark matter mass values in GeV.
    lambda_hs : np.ndarray
        The coupling between Higgs boson and the dark matter.

    Return
    ------
    Gamma_inv : ndarray
        The invisible width of the Higgs boson.
    '''
    Gamma_inv = np.power(lambda_hs*v,2.)/(32.*np.pi*mh)*np.sqrt(1.-np.power(2.*DMmass/mh,2.))
    return Gamma_inv

def Br_inv(lambda_hs,DMmass,Gamma_H):
    '''Calculates the Branching ratio of the Standard Model Higgs boson to invisible particles.

    Parameters
    ----------
    lambda_hs : np.ndarray
        The coupling between Higgs boson and the dark matter.
    DMmass : np.ndarray
        Dark matter mass values in GeV.
    Gamma_H : np.ndarray
        The Higgs width in GeV.

    Return
    ------
    Br_inv : ndarray
        The Branching ratio of Higgs to invisible.
    '''
    Gamma_inv_value = Gamma_inv(DMmass,lambda_hs)
    Br_inv = (Gamma_inv_value)/(Gamma_H+Gamma_inv_value)
    return Br_inv

def func_interpolate(varval,variablevec,funcvec):
    '''This is a generic function to perform an interpolation

    Parameters
    ----------
    varval : np.ndarray
        Value of the variable for which we want the interpolated value
    variablevec : np.ndarray
        Vector for the independent variable
    funcvec : np.ndarray
        Vector for the dependent variable.

    Return
    ------
    value : ndarray
        Interpolated value
    '''
    result = 0.
    if varval<variablevec[0]:
        result = 0.
    elif varval>variablevec[len(variablevec)-1]:
        result = 0.
    else:
        Log10E_bin = np.log10(variablevec[1])-np.log10(variablevec[0])
        nbin = ( np.log10(varval)-np.log10(variablevec[0]) )/Log10E_bin
        binval = int(nbin)
        result = pow(10.,  np.log10(funcvec[binval]) + (np.log10(funcvec[binval+1])-np.log10(funcvec[binval]))*(np.log10(varval)-np.log10(variablevec[binval]))/(np.log10(variablevec[binval+1])-np.log10(variablevec[binval]))  )
    return result

def flux_DM_prompt(Energy,DMmass,lambda_hs):
    '''Calculates the gamma-ray flux for the Galactic center excess.

    Parameters
    ----------
    Energy : np.ndarray
        Gamma-ray energy in GeV
    DMmass : np.ndarray
        Dark matter mass values in GeV.
    lambda_hs : np.ndarray
        The coupling between Higgs boson and the dark matter.

    Return
    ------
    flux : ndarray
        gamma-ray flux for the Galactic center excess.
    '''
    
    kpc = pow(10.,3.)*3.0856775*pow(10.,18.) #cm
    rodot = 8.12*kpc #cm
    rhodot = 0.300 #GeV/cm^3
    Jfact = 117.0 #MED for Cholis
    solidangle = 0.4288
    
    if Energy>DMmass:
        return 0.
    else:
        table_int = np.loadtxt(import_data_file('SHP_sigmav_table.dat'))
        sigmav = _lambda2sigmav(DMmass,lambda_hs,table_int)
        x_vec,fluxDM_x = DMspectra_inttable(DMmass,lambda_hs,'gammas',smooth=False)
        EnergyDM_vec = pow(10.,x_vec)*DMmass
        fluxDM_vec = fluxDM_x/(np.log(10.)*np.power(10.,x_vec)*DMmass)
        for t in range(len(EnergyDM_vec)):
            if fluxDM_vec[t]>0.:
                fluxDM_vec[t]=fluxDM_vec[t]
            else:
                fluxDM_vec[t]=1e-30
        
        dNdE = func_interpolate(Energy,EnergyDM_vec,fluxDM_vec)
        
        return 0.5*(rodot/(4.*np.pi))*pow(rhodot/DMmass,2.)*Jfact*sigmav*dNdE