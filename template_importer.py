
import matplotlib
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from scipy.optimize import curve_fit
import numpy as np
import scipy as sc
from scipy import stats
from scipy.stats import norm
from astropy  import units as u
from astropy import constants as const
import scipy as sp
from scipy.integrate import quad
from scipy import integrate
import math
from scipy.signal import convolve
from scipy import interpolate
from scipy.interpolate import interp1d
import time
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from scipy.stats import norm
import matplotlib.mlab as mlab
import sys

from config import *
from functions import *

'''
Defining constants and wavelength range for all templates
'''
c=const.c
Lsol=(const.L_sun.cgs).value


RFull=np.logspace(-4,5.5,10**4)
RFullnu=c*(RFull*10**-6)**-1
ltn=(RFull**2)*(3*10**8)**-1
#==============================================================================
diagnostic_optical=False

#OT=Table.read('templates/narrow_test.fits')
OT=Table.read('templates/narrow_test_no_lines_ext.fits')
OT_arr = (np.lib.recfunctions.structured_to_unstructured(OT.as_array())).T

#t_param=Table.read('templates/fsps_QSF_12_v3.param.fits')
t_param=Table.read('templates/fsps_QSF_12_v3_narrow.param.fits')


RFull=np.array(OT['lambda'])
ltn=(RFull**2)*(3*10**8)**-1

OT_arr_flux=OT_arr[1:,:]

GB=np.zeros((len(OT_arr_flux),(len(RFull))))
for i,_ in enumerate(OT_arr_flux[:,0]):
    #f=interp1d(OT_arr[0,:],OT_arr_flux[i,:],bounds_error=False,fill_value=0,kind='linear')
    #GB[i,:]=f(RFull)
    GB[i,:]=OT_arr_flux[i,:]

if qso:
    #qso_templ=Table.read('templates/shen2016_ext.fits')
    #qso_func=interp1d(qso_templ['wave']*10**-4,qso_templ['flux'],bounds_error=False,fill_value=0,kind='cubic')
    #qso=qso_func(RFull)*ltn
    GB[:,:]*=0
    QSO,qso_wave=get_qso_templates()
    for i,_ in enumerate(QSO[:,0]):
        qso_func=interp1d(qso_wave*10**-4,QSO[i,:],bounds_error=False,fill_value=0,kind='cubic')
        qso=qso_func(RFull)*ltn
        GB[i,:]=qso


if diagnostic_optical:
    for i,_ in enumerate(GB[:,0]):
        plot(RFull,GB[i,:])
        yscale('log')
        xscale('log')
    show()
#==============================================================================


AGN_l=Table.read('templates/AGN/LoLum_AGN_GM.fits')
AGN_h=Table.read('templates/AGN/HiLum_AGN_GM.fits')


agn_h=interp1d(AGN_h['wave']*10**-4,AGN_h['flux'],bounds_error=False,fill_value=0,kind='cubic')
agn_l=interp1d(AGN_l['wave']*10**-4,AGN_l['flux'],bounds_error=False,fill_value=0,kind='cubic')

agn_high=agn_h(RFull)*ltn
agn_low=agn_l(RFull)*ltn

agn=np.zeros((2,len(agn_low)))
agn[0,:]=agn_high
agn[1,:]=agn_low



#==============================================================================
'''DL07 CONFIGURATION PANEL'''

'Range of the gamma parameter, change at your discretion'
g1=np.array([0,0.001,0.0025,0.005,0.0075])
g2=np.arange(0.01,0.1,0.01)
g3=np.array([0.2,0.35,0.5])
g=np.concatenate((g1,g2,g3))

#==============================================================================
'''
LOAD IN DL07 Templates
'''
model1 = fits.open('templates/DL07/model1.fits', memmap=True)
model2 = fits.open('templates/DL07/model2.fits', memmap=True)
lam = fits.open('templates/DL07/lambda.fits', memmap=True)


model1_full = Table.read('templates/DL07/model1_full.fits', memmap=True)
model2_full = Table.read('templates/DL07/model2_full.fits', memmap=True)



if use_cold_dl:
    diffuse=np.array(model1_full['MODEL1'][0]) #DL diffuse part as 3d array
    pdr=np.array(model2_full['MODEL1'][0])  #DL pdr part as 3d array
    dat1=diffuse
    dat2=pdr
    Umin=np.array([0.100,0.200,0.300,0.400,0.500,0.600,0.700,0.800,1.000,1.200,
    1.700,2.000,3.000,4.000,5.000,7.000,8.000,10.00,12.00,15.00,20.00,25.00,30.00,
    35.00,40.00,50.00])
    q_indices=np.array([0,1,2,3,4,5,6,7,8,9,10])
    Umin=Umin[:]
    dat1=dat1[:,:,q_indices]
    dat2=dat2[:,:,q_indices]
else:
    dat1=model1[0].data
    dat2=model2[0].data
    Umin=np.array([0.800,1.200,2.000,3.000,4.000,5.000,7.000,8.000,10.00,12.00,
    15.00,20.00,25.00,30.00,35.00,40.00,50.00])

DG_ratio=0.01 #Template specific dust-to-gas ratio, see DL07 for more details.


dat3=lam[0].data


def S(i,j,k):
    sed=((1.-g[k])*dat1[:,i,j]+g[k]*dat2[:,i,j])
    ss=interp1d(np.log10(dat3), np.log10(sed), bounds_error=False,fill_value='extrapolate', kind='linear')
    sed=10**ss(np.log10(RFull))
    sed=np.array(sed)
    return sed

DL07_0=np.zeros((len(RFull),len(dat1[0,:,0]),len(dat1[0,0,:]),len(g)))

for i,obj in enumerate(dat1[0,:,0]):
    for j,obj in enumerate(dat1[0,0,:]):
        for k,obj in enumerate(g):
            DL07_0[:,i,j,k]=S(i,j,k)


#==============================================================================

print('Templates Imported')
