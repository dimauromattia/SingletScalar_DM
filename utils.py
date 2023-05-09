from math import *
import matplotlib.pyplot as pl
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.ndimage.filters import gaussian_filter1d
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize
from include import *

def func_interpolate_Omega(mass_val,lambda_val,QCDmodel):
    '''
    This script calculates the relic density as Omega h^2 given the mass and lambda
    mass_val : dark matter mass value in GeV
    lambda_val: lambda_HS value
    QCDmodel: Model for the QCD phase transition. Possible choices are 'QCDA' and 'QCDB'
    '''
    
    if mass_val>mass_vector_drake.min() and mass_val<mass_vector_drake.max():
        
        bins=100
        table = np.loadtxt('DATA/tableOmega_DRAKE_fBE_%s.txt'%QCDmodel)
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
            #func1 = interp1d(np.sqrt(lambda_mass1),Omega_mass1)
            #omega1 = float(func1(sqrt(lambda_val)))
            #func1 = interp1d(lambda_mass1,np.sqrt(Omega_mass1))
            #omega1 = np.power(float(func1(lambda_val)),2.)
        else:
            print('Warning, extrapolating.....')
            print('The relic density with DRAKE has been calculated for lambda_HS between ',lambda_mass1.min(),lambda_mass1.max())
            func1 = CubicSpline(lambda_mass1,np.sqrt(Omega_mass1), bc_type='not-a-knot')
            omega1 = np.power(float(func1(lambda_val)),2.)
            check = 1
            
        if lambda_val<lambda_mass2.max() and lambda_val>lambda_mass2.min():
            func2 = interp1d(lambda_mass2,np.sqrt(Omega_mass2))
            omega2 = np.power(float(func2(lambda_val)),2)
            #func2 = interp1d(np.sqrt(lambda_mass2),Omega_mass2)
            #omega2 = float(func2(sqrt(lambda_val)))
            #func2 = interp1d(lambda_mass2,np.sqrt(Omega_mass2))
            #omega2 = np.power(float(func2(lambda_val)),2)
        else:
            func2 = CubicSpline(lambda_mass2,np.sqrt(Omega_mass2), bc_type='not-a-knot')
            omega2 = np.power(float(func2(lambda_val)),2)
            if check==0:
                print('Warning, extrapolating.')
                print('For this mass pick a range of lambda between ',lambda_mass1.min(),lambda_mass1.max())
        
        func_final = interp1d(np.array([mass_vector_drake[idx],mass_vector_drake[idx+1]]),np.array([np.power(omega1,0.1),np.power(omega2,0.1)]))
        return float(np.power(func_final(mass_val),1/0.1))

    else:
        if mass_val>mass_vector_micro.min() and mass_val<mass_vector_micro.max():
            model = 'DATA/'
            table = np.loadtxt(model+'/tableOmega_MicroOMEGAs.txt')
            Omega_table = table[:,2]
            lambdahs_vec = np.logspace(-5,2,500)
            func  = interp2d(mass_vector_micro,lambdahs_vec,Omega_table)
            return float(func(mass_val,lambda_val))
        else:
            print('Warning, mass should be between',mass_vector_micro.min(),mass_vector_micro.max())
            

def func_interpolate_Omega_MicrOMEGAs(mass_val,lambda_val):
    '''
    This script calculates the relic density as Omega h^2 using MicrOMEGAs results and given the mass and lambda.
    mass_val : dark matter mass value in GeV
    lambda_val: lambda_HS value
    QCDmodel: Model for the QCD phase transition. Possible choices are 'QCDA' and 'QCDB'
    '''
    
    if mass_val>mass_vector_micro.min() and mass_val<mass_vector_micro.max():
        model = 'DATA/'
        table = np.loadtxt(model+'/tableOmega_MicroOMEGAs.txt')
        Omega_table = table[:,2]
        lambdahs_vec = np.logspace(-5,2,500)
        func  = interp2d(mass_vector_micro,lambdahs_vec,Omega_table)
        return float(func(mass_val,lambda_val))
    else:
        print('Warning, mass should be between',mass_vector_micro.min(),mass_vector_micro.max())


def func_interpolate_lambda(mass_val,Omega_val,QCDmodel):
    '''
    This script calculates the lambda_HS parameter for a given relic density (Omega h^2) and dark matter mass
    mass_val : dark matter mass value in GeV
    Omega_val: relic density (Omega h^2)
    QCDmodel: Model for the QCD phase transition. Possible choices are 'QCDA' and 'QCDB'
    '''
    
    if mass_val>mass_vector_drake.min() and mass_val<mass_vector_drake.max():
        
        bins=100
        table = np.loadtxt('DATA/tableOmega_DRAKE_fBE_%s.txt'%QCDmodel)
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
            print('Warning, extrapolating.')
            print('For this mass pick a range of lambda between ',Omega_mass1.min(),Omega_mass1.max())
            check = 1
            #func1 = CubicSpline(Omega_mass1,np.log10(lambda_mass1), bc_type='not-a-knot')
            #lambda1 = float(func1(Omega_val))
            
        if Omega_val<Omega_mass2.max() and Omega_val>Omega_mass2.min():
            func2 = interp1d(Omega_mass2,lambda_mass2)
            lambda2 = float(func2(Omega_val))
        else:
            #func2 = CubicSpline(Omega_mass2,lambda_mass2, bc_type='not-a-knot')
            #lambda2 = float(func2(Omega_val))
            if check==0:
                print('Warning, extrapolating.')
                print('For this mass pick a range of lambda between ',Omega_mass2.min(),Omega_mass2.max())
                check = 1
        
        #print(omega1,omega2,mass_vector[idx],mass_vector[idx+1])
        #print(np.array([mass_vector[idx],mass_vector[idx+1]]),np.array([omega1,omega2)])
        if check==0:
            func_final = interp1d(np.array([mass_vector_drake[idx],mass_vector_drake[idx+1]]),np.array([lambda1,lambda2]))
            return func_final(mass_val)
        #return idx,mass_val,mass_vector[idx],mass_vector[idx+1],Omega_mass1,Omega_mass2,omega1,omega2,func_final(mass_val)
        else:
            return 0

    else:
        if mass_val>mass_vector_micro.min() and mass_val<mass_vector_micro.max():
            model = 'DATA/'
            table = np.loadtxt(model+'/tableOmega_MicroOMEGAs.txt')
            Omega_table = table[:,2]
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
                    #print(Omega_val,Omega_mass1[check_bin-1:check_bin+1],lambda_mass1[check_bin-1:check_bin+1])
                    func1 = interp1d(Omega_mass1[check_bin-1:check_bin+1],lambda_mass1[check_bin-1:check_bin+1])
                else:
                    #print(Omega_val,np.array([Omega_mass1[check_bin],Omega_mass1[check_bin-1]]),np.array([lambda_mass1[check_bin],lambda_mass1[check_bin-1]]))
                    func1 = interp1d(np.array([Omega_mass1[check_bin],Omega_mass1[check_bin-1]]),np.array([lambda_mass1[check_bin],lambda_mass1[check_bin-1]]))
            
                lambda1 = float(func1(Omega_val))
            else:
                print('Warning, extrapolating.')
                print('For this mass pick a range of lambda between ',Omega_mass1.min(),Omega_mass1.max())
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
                #func2 = CubicSpline(Omega_mass2,lambda_mass2, bc_type='not-a-knot')
                #lambda2 = float(func2(Omega_val))
                if check==0:
                    print('Warning, extrapolating.')
                    print('For this mass pick a range of lambda between ',Omega_mass2.min(),Omega_mass2.max())
                    check = 1
        
            #print(lambda1,lambda2)
            #print(Omega_mass1,Omega_mass2,len(Omega_mass1),len(Omega_mass2))
            #print(mass_vector_micro[idx],mass_vector_micro[idx+1],mass_val,Omega_val,lambda1,lambda2)
        
            if check==0:
                func_final = interp1d(np.array([mass_vector_micro[idx],mass_vector_micro[idx+1]]),np.array([lambda1,lambda2]))
                lambda_val = float(func_final(mass_val))
                if mass_val<70 and lambda_val>1.0:
                    print('Warning the problem f(lambda)=Omega h^2 could have two solutions for lambda')
                return lambda_val
            else:
                return 0
        else:
            print('Warning, mass should be between',mass_vector_micro.min(),mass_vector_micro.max())


def func_interpolate_lambda_MicrOMEGAs(mass_val,Omega_val,QCDmodel):
    '''
    This script calculates the lambda_HS parameter for a given relic density (Omega h^2) and dark matter mass
    mass_val : dark matter mass value in GeV
    Omega_val: relic density (Omega h^2)
    QCDmodel: Model for the QCD phase transition. Possible choices are 'QCDA' and 'QCDB'
    '''
    
    if mass_val>mass_vector_micro.min() and mass_val<mass_vector_micro.max():
        model = 'DATA/'
        table = np.loadtxt(model+'/tableOmega_MicroOMEGAs.txt')
        Omega_table = table[:,2]
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
                #print(Omega_val,Omega_mass1[check_bin-1:check_bin+1],lambda_mass1[check_bin-1:check_bin+1])
                func1 = interp1d(Omega_mass1[check_bin-1:check_bin+1],lambda_mass1[check_bin-1:check_bin+1])
            else:
                #print(Omega_val,np.array([Omega_mass1[check_bin],Omega_mass1[check_bin-1]]),np.array([lambda_mass1[check_bin],lambda_mass1[check_bin-1]]))
                func1 = interp1d(np.array([Omega_mass1[check_bin],Omega_mass1[check_bin-1]]),np.array([lambda_mass1[check_bin],lambda_mass1[check_bin-1]]))
        
            lambda1 = float(func1(Omega_val))
        else:
            print('Warning, extrapolating.')
            print('For this mass pick a range of lambda between ',Omega_mass1.min(),Omega_mass1.max())
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
            #func2 = CubicSpline(Omega_mass2,lambda_mass2, bc_type='not-a-knot')
            #lambda2 = float(func2(Omega_val))
            if check==0:
                print('Warning, extrapolating.')
                print('For this mass pick a range of lambda between ',Omega_mass2.min(),Omega_mass2.max())
                check = 1
    
        #print(lambda1,lambda2)
        #print(Omega_mass1,Omega_mass2,len(Omega_mass1),len(Omega_mass2))
        #print(mass_vector_micro[idx],mass_vector_micro[idx+1],mass_val,Omega_val,lambda1,lambda2)
    
        if check==0:
            func_final = interp1d(np.array([mass_vector_micro[idx],mass_vector_micro[idx+1]]),np.array([lambda1,lambda2]))
            lambda_val = float(func_final(mass_val))
            if mass_val<70 and lambda_val>1.0:
                print('Warning the problem f(lambda)=Omega h^2 could have two solutions for lambda')
            return lambda_val
        else:
            return 0
    else:
        print('Warning, mass should be between',mass_vector_micro.min(),mass_vector_micro.max())
        


def func_lambda2sigmav(DMmass_val,lambdahs_val,table_int):
    '''
    This function calculates the annihilation cross section in units of cm^3/s for a given value of the mass and lambda_HS
    mass_val : dark matter mass value in GeV
    lambda_val: lambda_HS value
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
        logres_1 = log10(val1)
        logres_2 = log10(val2)
        loglambda_1 = log10(lambdahs_vec[len(lambdahs_vec)-2])
        loglambda_2  = log10(lambdahs_vec[len(lambdahs_vec)-1])
        #rescaling_lambda = pow(lambdahs_val/1e1,2.)
        result = logres_2 + (logres_2-logres_1)*(log10(lambdahs_val)-loglambda_2)/(loglambda_2-loglambda_1)
        #print(pow(10.,logres_1),pow(10.,logres_2),pow(10.,loglambda_1),pow(10.,loglambda_2),pow(10.,result))
        return pow(10.,result)
        
    else:
        func_int = interp2d( massz_vec, lambdahs_vec , np.ndarray.flatten( table_int[:,1:(len(lambdahs_vec)+1)]),kind='cubic')
        return float(func_int( DMmass_val, lambdahs_val) )


def func_sigmav_channels(DMmass_val,lambdahs_val,channel):
    '''
    This function calculates the annihilation cross section in units of cm^3/s for a given value of the mass and lambda_HS and DM channel
    mass_val : dark matter mass value in GeV
    lambda_val: lambda_HS value
    channel np.array(['cc','bb','tt','tautau','gg','ww','zz','hh','aa','za'])
    '''
    
    channel_vec = np.array(['cc','bb','tt','tautau','gg','ww','zz','hh','aa','za'])
    #check = np.shape(np.where((channel_vec[:]==channel))[0])
    #print(check,np.where((channel_vec[:]==channel)),np.where((channel_vec[:]==channel)))
    #if check==1:
    if channel=='tot':
        table_int = np.loadtxt('DATA/SHP_sigmav_table.txt')
        sigmav_val = func_lambda2sigmav(DMmass_val,lambdahs_val,table_int)

    else:
        table_int = np.loadtxt('DATA/SHP_sigmav_%s.txt'%channel)
        sigmav_val = func_lambda2sigmav(DMmass_val,lambdahs_val,table_int)
    
    return sigmav_val
    #else:
    #    print('WARNING, possible channels are:',channel_vec)
    

def func_DMspectra_inttable(DMmass_val,lambdahs_val,particle,smooth):
    '''
    This function calculates the DM source spectra in units of dN/dlog10(x) for a given value of the mass and lambda_HS
    mass_val : dark matter mass value in GeV
    lambda_val: lambda_HS value
    particle: gammas, positrons, antiprotons, neutrinos_e, neutrinos_mu, neutrinos_tau
    smooth: Yes or No
    '''
    add = ''
    if smooth=='Yes':
        add = '_smooth'
    else:
        add = add
    tablespectra_int = np.loadtxt('DATA/SHP_spectra%s_table_%s.txt'%(add,particle))
    
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



def func_LZUL(DMmass):
    table = np.loadtxt('DATA/LZ_SI_2022_data.dat')
    mass_vec = table[:,0]
    sigma_SI = table[:,1]*1e-48
    if DMmass<mass_vec[0] or DMmass>mass_vec.max():
        return 0
    else:
        func_int_csi = interp1d(mass_vec,sigma_SI,kind='quadratic')
        return func_int_csi(DMmass)

def func_DARWINUL(DMmass):
    table = np.loadtxt('DATA/DARWIN_SI_proj.dat')
    mass_vec = table[:,0]
    sigma_SI = table[:,1]*1e-48
    if DMmass<mass_vec[0] or DMmass>mass_vec.max():
        return 0
    else:
        func_int_csi = interp1d(mass_vec,sigma_SI,kind='quadratic')
        return func_int_csi(DMmass)

def func_provide_ULEXP(DMmass,Exp):
    if Exp=='LZ':
        SI_exp = func_LZUL(DMmass)
    elif Exp== 'Darwin':
        SI_exp = func_DARWINUL(DMmass)
    return SI_exp
    
def func_SI_noomega(DMmass_val,lambdahs_val):
    '''
    This function gives the spin independent cross section for direct detection in cm^2
    mass_val : dark matter mass value in GeV
    lambda_val: lambda_HS value
    '''
    mu = (mN*DMmass_val)/(mN+DMmass_val)
    val = np.power(lambdahs_val*fN*mu*mN,2.)/(4.*np.pi*pow(mh,4.)*np.power(DMmass_val,2.))
    return GeVm2tocm2*val

def func_SI_withomega(lambda_hs,DMmass,Lambda_vec,Mass_vec,csi_vec):
    
    mu = (mN*DMmass)/(mN+DMmass)
    val = np.power(lambda_hs*fN*mu*mN,2.)/(4.*np.pi*pow(mh,4.)*np.power(DMmass,2.))
    func_int_csi = interp2d(Mass_vec,Lambda_vec,csi_vec)
    csi_val = func_int_csi(DMmass,lambda_hs)[0]
    return GeVm2tocm2*val*csi_val

def func_GetUL_DD_nomega(DMmass_val,Exp):
    value = 0
    if Exp=='LZ':
        value = sqrt(func_LZUL(DMmass_val)/func_SI_noomega(DMmass_val,1.0))
        return value
    elif Exp== 'DARWIN':
        value = sqrt(func_DARWINUL(DMmass_val)/func_SI_noomega(DMmass_val,1.0))
        return value
    else:
        print('WARNING: Wrong parameter for the experiment.')
        print('Choose among LZ and DARWIN')

def func_tominimize_DD(Lambda_val,DMmass,Lambda_vec,Mass_vec,csi_vec,Exp):
    #print(Lambda_val)
    if Lambda_val<=Lambda_vec[0]:
        Lambda_val=Lambda_vec[0]
    elif Lambda_val>=Lambda_vec[len(Lambda_vec)-1]:
        Lambda_val=Lambda_vec[len(Lambda_vec)-1]
    SI_theory = func_SI_withomega(Lambda_val,DMmass,Lambda_vec,Mass_vec,csi_vec)
    SI_exp = func_provide_ULEXP(DMmass,Exp)
        
    #print(Lambda_val,SI_exp,SI_theory,abs(SI_exp-SI_theory)/SI_theory)
    #return abs(func_LuxUL(DMmass)-SI_theory)/SI_theory
    return abs(SI_exp-SI_theory)
    
def func_GetUL_DD_withomega(DMmass,Lambda_vec,Mass_vec,csi_vec,Exp):
    
    SI_theory = func_SI_withomega(Lambda_vec,DMmass,Lambda_vec,Mass_vec,csi_vec)
    SI_exp = func_provide_ULEXP(DMmass,Exp)
    #print(SI_theory,SI_theory.max(),SI_theory.min(),SI_exp)
    if SI_theory.max()<SI_exp:
        return Lambda_vec.max()
    
    elif SI_theory.min()>SI_exp:
        return Lambda_vec.min()
    
    else:
        lambda_0 = 1e-2
        #print(DMmass,Lambda_vec,Mass_vec,csi_vec,Exp)
        #print(minimize(func_tominimize_DD, lambda_0, method='nelder-mead', bounds=(1e-5,1e1), options={'xatol': 1e-8, 'disp': True}, args=(DMmass,Lambda_vec,Mass_vec,csi_vec,Exp)))
        value = minimize(func_tominimize_DD, lambda_0, method='nelder-mead', options={'xatol': 1e-8, 'disp': True}, args=(DMmass,Lambda_vec,Mass_vec,csi_vec,Exp)).x[0]
        SI_theory = func_SI_withomega(value,DMmass,Lambda_vec,Mass_vec,csi_vec)
        SI_exp = func_provide_ULEXP(DMmass,Exp)
        
        if abs(SI_exp-SI_theory) > 1e-47:
            lambda_0 = 1e0
            value = minimize(func_tominimize_DD, lambda_0, method='nelder-mead', options={'xatol': 1e-8, 'disp': True}, args=(DMmass,Lambda_vec,Mass_vec,csi_vec,Exp)).x[0]
            SI_theory = func_SI_withomega(value,DMmass,Lambda_vec,Mass_vec,csi_vec)
            SI_exp = func_provide_ULEXP(DMmass,Exp)
        
            if abs(SI_exp-SI_theory) > 1e-47:
                lambda_0 = 1e-4
                value = minimize(func_tominimize_DD, lambda_0, method='nelder-mead', options={'xatol': 1e-8, 'disp': True}, args=(DMmass,Lambda_vec,Mass_vec,csi_vec,Exp)).x[0]
                SI_theory = func_SI_withomega(value,DMmass,Lambda_vec,Mass_vec,csi_vec)
                SI_exp = func_provide_ULEXP(DMmass,Exp)
            
                if abs(SI_exp-SI_theory) > 1e-47:
                    lambda_0 = 1e1
                    value = minimize(func_tominimize_DD, lambda_0, method='nelder-mead', options={'xatol': 1e-8, 'disp': True}, args=(DMmass,Lambda_vec,Mass_vec,csi_vec,Exp)).x[0]
                    SI_theory = func_SI_withomega(value,DMmass,Lambda_vec,Mass_vec,csi_vec)
                    SI_exp = func_provide_ULEXP(DMmass,Exp)
                
                    if abs(SI_exp-SI_theory) > 1e-47:
                        value = 0.
        return value

def func_minimize(DMmass,Gamma_H,Gamma_inv_measured):
    lambda_0 = 0.01
    value = minimize(func_Br_inv_UL, lambda_0, method='nelder-mead', options={'xatol': 1e-8, 'disp': True}, args=(DMmass,Gamma_H,Gamma_inv_measured)).x[0]
    return value 

def func_Gamma_inv(DMmass,lambda_hs):
    return np.power(lambda_hs*v,2.)/(32.*np.pi*mh)*np.sqrt(1.-np.power(2.*DMmass/mh,2.))

def func_Br_inv(lambda_hs,DMmass,Gamma_H):
    Gamma_inv = func_Gamma_inv(DMmass,lambda_hs)
    return (Gamma_inv)/(Gamma_H+Gamma_inv)

def func_Br_inv_UL(lambda_hs,DMmass,Gamma_H,Gamma_inv_measured):
    return np.abs(func_Br_inv(lambda_hs,DMmass,Gamma_H)-Gamma_inv_measured)/Gamma_inv_measured
