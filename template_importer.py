
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
OT=Table.read('templates/GB.fits')
OT_arr = (np.lib.recfunctions.structured_to_unstructured(OT.as_array())).T

OT_arr_flux=OT_arr[1:,:]

GB=np.zeros((len(OT_arr_flux),(len(RFull))))
for i,_ in enumerate(OT_arr_flux[:,0]):
    f=interp1d(OT_arr[0,:],OT_arr_flux[i,:],bounds_error=False,fill_value=0,kind='linear')
    GB[i,:]=f(RFull)


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
